#include "build.h"

#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>
#include <sys/prctl.h>

namespace {
void PrintUsage() {
  std::cout << "用法：./build RUN [选项]，或 ./build --run RUN [选项]\n\n"
            << "示例：./build 0616 --outdir ../rootfile/trial_0616\n"
            << "每轮同时输出 dataRUN_build_narrow.root 和 dataRUN_build_wide.root。\n"
            << "已有输出会拒绝覆盖；prompt/random 分别保存在 tr_event/tr_bg。\n\n"
            << "  --mode windows|cache|replay  默认 windows；cache 保存原始字段宽窗\n"
            << "  --input FILE            原始 C2 ROOT；replay 时指定本地宽窗 ROOT\n"
            << "  --cache-receipt FILE    replay 的缓存完成回执，默认与缓存同目录\n"
            << "  --outdir DIR            输出目录；默认项目 rootfile/\n"
            << "  --project DIR           项目根目录；通常由可执行文件位置确定\n"
            << "  --window-file FILE      逐 Ge、Si/S3、时期的窄/宽窗口表\n"
            << "  --profile-dir DIR       完整时间修正参数目录\n"
            << "  --source-manifest FILE  原始输入名单及校验记录\n"
            << "  --receipt FILE          完成回执；默认放在输出目录\n"
            << "  --entries N             只检查前 N 条；完整处理时省略\n"
            << "  --chunk-entries N       每批读取条数，默认 1000000\n"
            << "  --tag LABEL             输出名附加标签\n"
            << "  --inspect-parameters    只打印参数，不读取事件或写 ROOT\n"
            << "  --help                  显示本说明\n";
}

uint64_t ReadUnsigned(const std::string &text, const std::string &name) {
  if (text.empty() || text.find_first_not_of("0123456789") != std::string::npos)
    throw std::invalid_argument(name + " 必须是非负整数");
  std::size_t used = 0;
  auto value = std::stoull(text, &used);
  if (used != text.size())
    throw std::invalid_argument(name + " 不是完整整数");
  return value;
}

BuildOptions ReadOptions(int argc, char **argv) {
  BuildOptions options;
  std::set<std::string> seen;
  for (int i = 1; i < argc; ++i) {
    std::string key = argv[i];
    if (key == "--inspect-parameters") {
      options.inspect_parameters = true;
      continue;
    }
    std::string value;
    if (key.rfind("--", 0) != 0) {
      if (i != 1)
        throw std::invalid_argument(
            "当前按逐 Ge 窗口表选窗，不接受旧的 run win jump 参数；请用 --help 查看用法");
      value = key;
      key = "--run";
    } else {
      if (++i == argc)
        throw std::invalid_argument("缺少参数值：" + key);
      value = argv[i];
    }
    if (!seen.insert(key).second)
      throw std::invalid_argument("重复参数：" + key);

    if (key == "--run") {
      auto run = ReadUnsigned(value, key);
      if (run == 0 || run > uint64_t(std::numeric_limits<int>::max()))
        throw std::invalid_argument("无效轮次");
      options.run = int(run);
    } else if (key == "--project")
      options.project = value;
    else if (key == "--input")
      options.input = value;
    else if (key == "--outdir")
      options.outdir = value;
    else if (key == "--receipt")
      options.receipt = value;
    else if (key == "--window-file")
      options.window_file = value;
    else if (key == "--profile-dir")
      options.profile_dir = value;
    else if (key == "--source-manifest")
      options.source_manifest = value;
    else if (key == "--cache-receipt")
      options.cache_receipt = value;
    else if (key == "--mode")
      options.mode = value;
    else if (key == "--tag")
      options.tag = value;
    else if (key == "--entries")
      options.entry_limit = ReadUnsigned(value, key);
    else if (key == "--chunk-entries")
      options.chunk_entries = ReadUnsigned(value, key);
    else
      throw std::invalid_argument("未知参数：" + key);
  }
  if (!seen.count("--run"))
    throw std::invalid_argument("请指定轮次，例如 ./build 0616");
  return options;
}
} // namespace

int main(int argc, char **argv) {
  prctl(PR_SET_NAME, "build");
  for (int i = 1; i < argc; ++i) {
    if (std::string(argv[i]) == "--help") {
      PrintUsage();
      return 0;
    }
  }
  try {
    const auto options = ReadOptions(argc, argv);
    build processor(options);
    return processor.Process();
  } catch (const std::exception &error) {
    std::cerr << "ERROR: " << error.what() << '\n';
    return 1;
  }
}
