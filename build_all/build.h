#ifndef SM150_BUILD_H
#define SM150_BUILD_H

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>

// 一轮 build 的运行设置。空路径使用 build.cpp 中的项目默认位置。
// 能量刻度、时间修正和逐 Ge 符合窗从已有参数文件读取。
struct BuildOptions {
  int run = 0;
  std::filesystem::path project, input, outdir, receipt;
  std::filesystem::path window_file, profile_dir, source_manifest, cache_receipt;
  // windows：现行窄/宽；cache：原始字段宽窗；replay：从宽窗重建窄/宽。
  std::string mode = "windows", tag;

  // 0 处理整轮；非零仅用于小样检查，回执会标记覆盖范围。
  uint64_t entry_limit = 0;
  std::size_t chunk_entries = 1000000;
  bool inspect_parameters = false;
};

class build {
public:
  explicit build(const BuildOptions &options);

  // 读取 → 刻度/时间修正 → 分组/匹配 → 窄宽选窗 → 写出 → 完整读回。
  // 返回 0 表示本轮及读回核验完成；失败返回 1，并保留错误记录。
  int Process();

private:
  BuildOptions options_;
};

#endif
