// 150Sm：当前窄/宽两套事件build，按原三文件方式维护。
// main.cpp：命令行；build.h：运行设置与build接口；本文件：全部build实现。
// 阅读顺序：先看 build::Process() 的八个步骤，再按下面编号查具体处理。
// 参数与分组算法来自已核验的sm150_build；写出时去除跨参考分组的重复配对。
// 共用时间修正函数仍在 include/sm150_time_alignment.h，避免各程序出现不同版本。
// 旧的set.h固定窗/能量阈值不参与本入口，实际门列在classify及ROOT元数据中。
#include "build.h"
#include "sm150_time_alignment.h"
#include <TBranch.h>
#include <TBufferFile.h>
#include <TFile.h>
#include <TLeaf.h>
#include <TObjString.h>
#include <TROOT.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
#include <algorithm>
#include <array>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <ctime>
#include <fcntl.h>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <linux/fs.h>
#include <map>
#include <memory>
#include <nlohmann/json.hpp>
#include <openssl/evp.h>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/file.h>
#include <sys/prctl.h>
#include <sys/stat.h>
#include <sys/syscall.h>
#include <unistd.h>
#include <utility>
#include <vector>

namespace fs = std::filesystem;
using Clock = std::chrono::steady_clock;
static std::string
sha256(const fs::path &path,
       const std::function<void(uint64_t, uint64_t)> &progress = {});

// ==========================================================================
// 1. 探测器通道映射声明（完整映射表放在文件末尾）
// ==========================================================================
// Exact maps copied from build_all/build.cpp; no changes to physical topology.
namespace sm150_native {
struct Maps {
  std::map<int, int> map_ge_ring_id, map_ge_sector_id, map_spider_ring_id,
      map_spider_sector_id, map_s3_ring_id, map_s3_sector_id;
  std::map<int, int> map_clover_id = {{212, 2}, {213, 2}, {214, 2}, {215, 2},

                                      {400, 3}, {401, 3}, {402, 3}, {403, 3},

                                      {404, 4}, {405, 4}, {406, 4}, {407, 4},

                                      {409, 5}, {410, 5}, {411, 5}};
  bool init(int run);
};
} // namespace sm150_native

// ==========================================================================
// 2. 读取原有pars参数：能量刻度、逐轮增益、衰减和正背面归一
// ==========================================================================

namespace sm150_native {
using json = nlohmann::json;
inline void require(bool ok, const std::string &message) {
  if (!ok)
    throw std::runtime_error(message);
}
inline std::string run_name(int run) {
  std::ostringstream s;
  s << std::setw(4) << std::setfill('0') << run;
  return s.str();
}
inline std::vector<std::vector<double>> table(const std::string &path,
                                              bool header = true) {
  std::ifstream input(path);
  require(bool(input), "Missing parameter file: " + path);
  std::string line;
  if (header)
    std::getline(input, line);
  std::vector<std::vector<double>> rows;
  while (std::getline(input, line)) {
    const auto hash = line.find('#');
    if (hash != std::string::npos)
      line.resize(hash);
    std::istringstream stream(line);
    double x;
    std::vector<double> row;
    while (stream >> x) {
      require(std::isfinite(x), "Nonfinite parameter: " + path);
      row.push_back(x);
    }
    require(stream.eof(), "Non-numeric parameter row: " + path);
    if (!row.empty())
      rows.push_back(row);
  }
  require(!rows.empty(), "Empty parameter table: " + path);
  return rows;
}
struct Parameters {
  Maps maps;
  std::vector<double> p = std::vector<double>(16384 * 9, 0.);
  std::vector<int32_t> q = std::vector<int32_t>(16384 * 5, 0);
  std::vector<std::string> files;
  json file_hashes = json::object();
  json records = json::object();
  std::string group, run;
  std::vector<std::vector<double>> read(const std::string &path,
                                        bool header = true) {
    const auto hash = sha256(path);
    auto rows = table(path, header);
    require(sha256(path) == hash, "Parameter changed while loading: " + path);
    files.push_back(path);
    file_hashes[path] = hash;
    return rows;
  }
  using Keyed = std::map<int, std::vector<double>>;
  Keyed keyed(const std::string &path, int width) {
    Keyed result;
    for (const auto &row : read(path)) {
      require(row.size() >= size_t(3 + width),
              "Short calibration row: " + path);
      int key = int(10000 * row[0] + 100 * row[1] + row[2]);
      require(!result.count(key), "Duplicate channel: " + path);
      result[key] =
          std::vector<double>(row.begin() + 3, row.begin() + 3 + width);
    }
    return result;
  }
  Parameters(const std::string &project, int r, const json &calibration) {
    run = run_name(r);
    group = calibration.at("groups").at(run).get<std::string>();
    require(group == "low" || group == "high1" || group == "high2",
            "Unsupported run group");
    require(maps.init(r), "Native map does not support run");
    for (int k = 0; k < 16384; ++k)
      q[k * 5 + 3] = -1;
    const std::string pars = project + "/pars/";
    const std::string cal_name = group == "low" ? "cali_si_r457_low_current.txt"
                                 : group == "high1"
                                     ? "cali_si_r460.txt"
                                     : "cali_si_r457_high_current.txt";
    const std::string off_name = group == "low"     ? "ts_offset_r616_r674.txt"
                                 : group == "high1" ? "ts_offset_r472_r605.txt"
                                                    : "ts_offset_r680_r718.txt";
    const auto ge = keyed(pars + "cali_ge.txt", 3),
               si = keyed(pars + cal_name, 3),
               offset = keyed(pars + off_name, 1),
               gain = keyed(pars + "run_k/k_" + run + ".txt", 1);
    std::set<int> good;
    for (const auto &row : read(pars + "config/b.txt")) {
      require(row.size() >= 3, "Short good-Ge row");
      if (row[2] == 1)
        good.insert(int(100 * row[0] + row[1]));
    }
    auto add = [&](int key, int kind, int ring, int sector,
                   const std::vector<double> &c, double g, double kg,
                   double att, double fb0, double fb1) {
      require(key >= 0 && key < 16384 && !q[key * 5],
              "Duplicate or invalid physical channel");
      const auto &off = offset.at(key);
      double baseline = calibration.at("channel_offsets")
                            .at(group)
                            .at(std::to_string(key))
                            .get<double>();
      require(std::floor(off[0]) == std::floor(baseline),
              "Frozen native offset changed");
      const std::array<double, 9> values = {c[0], c[1], c[2], g,     kg,
                                            att,  fb0,  fb1,  off[0]};
      std::copy(values.begin(), values.end(), p.begin() + key * 9);
      q[key * 5] = kind;
      q[key * 5 + 1] = ring;
      q[key * 5 + 2] = sector;
      auto clover = maps.map_clover_id.find(key);
      q[key * 5 + 3] = clover == maps.map_clover_id.end() ? -1 : clover->second;
      q[key * 5 + 4] = good.count(key) ? 1 : 0;
      records[std::to_string(key)] = {{"run", run},
                                      {"key", key},
                                      {"kind", kind == 1   ? "ge"
                                               : kind == 2 ? "si"
                                               : kind == 3 ? "s3_ring"
                                                           : "s3_sector"},
                                      {"cid", key / 10000},
                                      {"sid", key % 10000 / 100},
                                      {"ch", key % 100},
                                      {"ring", ring},
                                      {"sector", sector},
                                      {"cal_p0_kev", c[0]},
                                      {"cal_p1_kev_per_adc", c[1]},
                                      {"cal_p2_kev_per_adc2", c[2]},
                                      {"gain", g},
                                      {"k_ge", kg},
                                      {"attenuation", att},
                                      {"fb_p0_kev", fb0},
                                      {"fb_p1", fb1},
                                      {"offset_ns", off[0]}};
    };
    for (const auto &pair : maps.map_ge_ring_id) {
      int key = pair.first;
      add(key, 1, pair.second, maps.map_ge_sector_id.at(key), ge.at(key), 1.,
          gain.at(key)[0], 1., 0., 1.);
    }
    for (const auto &pair : maps.map_spider_ring_id) {
      int key = pair.first, sid = key % 10000 / 100;
      double g = 1.;
      if (((475 < r && r < 616) || (680 <= r && r <= 718)) && sid == 2)
        g = 4.;
      if (616 <= r && r <= 674 && sid == 9)
        g = 4.;
      add(key, 2, pair.second, maps.map_spider_sector_id.at(key), si.at(key), g,
          1., 1., 0., 1.);
    }
    for (int face = 0; face < 2; ++face) {
      const std::string name = face ? "sector" : "ring";
      const int count = face ? 32 : 24;
      std::map<int, std::array<double, 2>> fb;
      for (const auto &row :
           read(pars + "run_fb_cor/correction_" + name + "_" + run + ".txt")) {
        require(row.size() >= 3 && !fb.count(int(row[0])),
                "Invalid front/back table");
        fb[int(row[0])] = {row[1], row[2]};
      }
      require(fb.size() == size_t(count), "Incomplete front/back table");
      const auto &ids = face ? maps.map_s3_sector_id : maps.map_s3_ring_id;
      for (const auto &pair : ids) {
        int key = pair.first, id = pair.second;
        std::ostringstream number;
        number << std::setw(2) << std::setfill('0') << id;
        auto rows =
            read(pars + "run_s3_att/" + name + number.str() + ".txt", false);
        double attenuation = 0.;
        int matched = 0;
        for (const auto &row : rows) {
          require(row.size() >= 2, "Short attenuation row");
          if (row[0] == r) {
            attenuation = row[1];
            ++matched;
          }
        }
        require(matched == 1 && attenuation > 0.,
                "Missing/invalid run attenuation");
        add(key, face ? 4 : 3, face ? 0 : id, face ? id : 0, si.at(key),
            key % 10000 / 100 == 11 ? 10. : 1., 1., attenuation, fb.at(id)[0],
            fb.at(id)[1]);
      }
    }
    require(records.size() == 189, "Incomplete native parameter coverage");
  }
};
} // namespace sm150_native

// ==========================================================================
// 3. 完整时间修正、排序和跨分块事件分组（时间单位ns）
// ==========================================================================
// Frozen streaming core retained during the C++ frontend port.
// Derived from codex/analysis/spectrum_all_runs_timing_20261015/src/extract.cpp
// Source SHA256:
// 6b76b98ee215b8cb375b35d1e41570f30e92d148a1c885bd765a01aeccf69944 Complete
// timing is applied before sorting and grouping; double dt retained.

// Single-thread streaming extraction. All ordering is (corrected ns, raw
// entry). The raw decoder has already ordered coarse timestamps; here an
// explicit worst-case signed CFD bound supplies the cross-chunk watermark.
namespace sm150_windowed {
struct Hit {
  double t;
  int64_t entry;
  double e;
  int key;
  uint8_t flags;
};
struct Engine {
  std::vector<double> par;
  std::vector<int32_t> info;
  std::vector<int> clovers;
  std::vector<Hit> tail;
  double consumed = -std::numeric_limits<double>::infinity();
  sm150::TimeAlignment timing;
  int64_t previous_raw = std::numeric_limits<int64_t>::min();
  int64_t max_shift = 0, read = 0, all_ge_groups = 0;
  std::string error;
  std::vector<int64_t> ref_entry, goff, poff, ge_entry, ge_second, p_entry;
  std::vector<double> ref_time;
  std::vector<int16_t> ge_key, p_key;
  std::vector<double> ge_dt, p_dt;
  std::vector<uint8_t> ge_flags, p_flags;
  std::vector<double> ge_energy, p_energy;
  Engine(const double *p, const int32_t *q, const char *directory, int run)
      : par(p, p + 16384 * 9), info(q, q + 16384 * 5) {
    timing.LoadProfile(directory, run);
    if (!timing.active() || !timing.total())
      throw std::runtime_error("Complete timing profile required");
    for (int k = 0; k < 16384; ++k)
      if (info[k * 5]) {
        max_shift =
            std::max(max_shift, int64_t(std::ceil(std::abs(par[k * 9 + 8]))));
        int cl = info[k * 5 + 3];
        if (cl >= 0 &&
            std::find(clovers.begin(), clovers.end(), cl) == clovers.end())
          clovers.push_back(cl);
      }
    std::sort(clovers.begin(), clovers.end());
    // cfds is signed int16, cfd is signed int16: <= 131080 ns from coarse Ge
    // time.
    max_shift =
        1000000; // computational watermark bound, not a coincidence gate
  }
  void clear_output() {
    ref_entry.clear();
    ref_time.clear();
    goff.assign(1, 0);
    poff.assign(1, 0);
    ge_entry.clear();
    ge_second.clear();
    p_entry.clear();
    ge_key.clear();
    p_key.clear();
    ge_dt.clear();
    p_dt.clear();
    ge_flags.clear();
    p_flags.clear();
    ge_energy.clear();
    p_energy.clear();
  }
  void emit(const Hit &a, const Hit *b, double ref) {
    if (!info[a.key * 5 + 4])
      return;
    ge_key.push_back(a.key);
    ge_entry.push_back(a.entry);
    ge_second.push_back(b ? b->entry : -1);
    ge_energy.push_back(a.e + (b ? b->e : 0.));
    ge_dt.push_back(a.t - ref);
    ge_flags.push_back(a.flags | (b ? (b->flags & 3) : 0));
  }
  int process(int64_t n, int64_t begin, const int16_t *cid, const int16_t *sid,
              const int16_t *ch, const int16_t *sr, const int64_t *ts,
              const uint16_t *adc, const int16_t *cfd, const int16_t *cfds,
              const uint8_t *cfdft, const uint8_t *pileup,
              const uint8_t *outofr, bool final,
              const int64_t *source_entry = nullptr) {
    clear_output();
    if (begin != read)
      throw std::runtime_error("Non-contiguous raw entry interval");
    for (int64_t i = 0; i < n; ++i) {
      if (sr[i] != 100 && sr[i] != 250)
        throw std::runtime_error("Unknown sampling rate");
      const int64_t coarse = ts[i] * (sr[i] == 250 ? 8 : 10);
      if (coarse < previous_raw)
        throw std::runtime_error("Coarse timestamp decreased");
      previous_raw = coarse;
      int key = 10000 * cid[i] + 100 * sid[i] + ch[i];
      if (key < 0 || key >= 16384 || !info[key * 5])
        throw std::runtime_error("Unmapped channel");
      const double *p = &par[key * 9];
      double x = adc[i], e = p[0] + p[1] * x + p[2] * x * x;
      e *= p[3] * p[4] * p[5];
      e = p[6] + p[7] * e;
      const int kind = info[key * 5];
      if (e <= (kind == 1 ? 20. : 2500.))
        continue;
      double t = coarse;
      if (cid[i] == 0) {
        if (sr[i] != 250)
          throw std::runtime_error("Ge sampling rate is not 250");
        t = cfdft[i] ? double(8 * ts[i])
                     : (2 * ts[i] - cfds[i] + cfd[i] / 16384.) * 4.;
      }
      const int64_t base = int64_t(t);
      const double correction =
          (kind == 1) ? timing.TotalGeOffset(key, cfdft[i] != 0, e)
                      : timing.TotalParticleOffsetFromAdc(key, adc[i]);
      const double corrected = double(base) + correction;
      if (std::abs(corrected - coarse) > max_shift)
        throw std::runtime_error("Timestamp bound violated");
      tail.push_back({corrected, source_entry ? source_entry[i] : begin + i, e, key,
                      uint8_t(pileup[i] | outofr[i] << 1 | cfdft[i] << 2)});
    }
    read += n;
    std::sort(tail.begin(), tail.end(), [](const Hit &a, const Hit &b) {
      return a.t < b.t || (a.t == b.t && a.entry < b.entry);
    });
    std::vector<const Hit *> ge, particles;
    ge.reserve(tail.size() / 8);
    particles.reserve(tail.size());
    for (const auto &h : tail)
      (info[h.key * 5] == 1 ? ge : particles).push_back(&h);
    const double watermark = final ? std::numeric_limits<double>::infinity()
                                   : double(previous_raw - max_shift - 1000);
    auto lower = [](const Hit *a, double t) { return a->t < t; };
    auto upper = [](double t, const Hit *a) { return t < a->t; };
    size_t pos =
        std::upper_bound(ge.begin(), ge.end(), consumed, upper) - ge.begin();
    std::vector<const Hit *> group;
    while (pos < ge.size() && ge[pos]->t <= watermark) {
      const Hit *ref = ge[pos];
      const double t = ref->t;
      size_t start =
          std::lower_bound(ge.begin(), ge.end(), t - 100, lower) - ge.begin();
      size_t stop =
          std::upper_bound(ge.begin(), ge.end(), t + 100, upper) - ge.begin();
      group.clear();
      group.push_back(ref);
      for (size_t j = pos; j > start;)
        group.push_back(ge[--j]);
      for (size_t j = pos + 1; j < stop; ++j)
        group.push_back(ge[j]);
      pos = stop;
      consumed = t + 100;
      ++all_ge_groups;
      auto pl =
          std::lower_bound(particles.begin(), particles.end(), t - 1000, lower);
      auto pr =
          std::upper_bound(particles.begin(), particles.end(), t + 1000, upper);
      if (pl == pr)
        continue;
      size_t old = ge_energy.size();
      // Match the reference helper ordering: singles, then Clover IDs. Photon
      // ordering has no influence on spectra or the frozen reference timestamp.
      for (auto a : group)
        if (info[a->key * 5 + 3] < 0)
          emit(*a, nullptr, t);
      for (int cl : clovers) {
        const Hit *a = nullptr, *b = nullptr;
        int count = 0;
        for (auto h : group)
          if (info[h->key * 5 + 3] == cl) {
            if (count == 0)
              a = h;
            if (count == 1)
              b = h;
            ++count;
          }
        if (count == 1)
          emit(*a, nullptr, t);
        else if (count == 2 && std::abs(a->key % 100 - b->key % 100) != 2) {
          if (a->e <= b->e)
            std::swap(a, b);
          emit(*a, b, t);
        }
      }
      if (old == ge_energy.size())
        continue;
      ref_entry.push_back(ref->entry);
      ref_time.push_back(t);
      for (auto it = pl; it != pr; ++it) {
        const Hit &h = **it;
        p_entry.push_back(h.entry);
        p_key.push_back(h.key);
        p_energy.push_back(h.e);
        p_dt.push_back(h.t - t);
        p_flags.push_back(h.flags);
      }
      goff.push_back(ge_energy.size());
      poff.push_back(p_energy.size());
    }
    if (final)
      tail.clear();
    else {
      auto keep =
          std::lower_bound(tail.begin(), tail.end(), watermark - 1000,
                           [](const Hit &a, double t) { return a.t < t; });
      tail.erase(tail.begin(), keep);
    }
    return 0;
  }
};
} // namespace sm150_windowed

// ==========================================================================
// 4. 逐Ge窗表、粒子匹配、prompt/random写出及读回检查
// ==========================================================================

namespace sm150_native {
struct Gate {
  double lower = 0., upper = 0., alpha = 0.;
  std::vector<std::pair<double, double>> random;
  json record;
  int side(double dt) const {
    for (size_t i = 0; i < random.size(); ++i)
      if (random[i].first <= dt && dt < random[i].second)
        return i == 0 ? -1 : 1;
    return 0;
  }
};
struct Gates {
  std::map<std::pair<int, int>, std::array<Gate, 2>> rows;
  json windows;
  explicit Gates(const json &input, const std::string &group) {
    for (const auto &row : input.at("windows"))
      if (row.at("group") == group) {
        int ge = row.at("ge"), view = row.at("view");
        require(ge >= 0 && ge < 1600 && (view == 0 || view == 1),
                "Invalid Ge/reference in windows");
        auto key = std::make_pair(ge, view);
        require(!rows.count(key), "Duplicate window combination");
        std::array<Gate, 2> pair;
        for (int label = 0; label < 2; ++label) {
          const auto &w = row.at(label ? "B" : "A");
          auto &g = pair[label];
          g.lower = w.at("lower_ns");
          g.upper = w.at("upper_ns");
          g.alpha = w.at("alpha");
          auto same_width = [](double a, double b) {
            return std::isfinite(a) && std::isfinite(b) &&
                   std::abs(a - b) <= 1e-9;
          };
          require(std::isfinite(g.lower) && std::isfinite(g.upper) &&
                      g.upper > g.lower && std::isfinite(g.alpha) &&
                      g.alpha > 0 &&
                      same_width(g.upper - g.lower,
                                 w.at("width_ns").get<double>()) &&
                      w.at("extra_factor") == 1,
                  "Invalid prompt width or extra factor");
          double width = 0.;
          for (const auto &r : w.at("random")) {
            require(r.is_array() && r.size() == 2, "Invalid random interval");
            const double lower = r.at(0), upper = r.at(1);
            require(std::isfinite(lower) && std::isfinite(upper) &&
                        upper > lower,
                    "Invalid random interval width");
            g.random.emplace_back(lower, upper);
            width += upper - lower;
          }
          require(
              g.random.size() == 2 &&
                  same_width(width, w.at("random_width_ns").get<double>()) &&
                  same_width(g.alpha * width, g.upper - g.lower),
              "Invalid random normalization");
          // random_side=-1/+1 必须分别代表左/右；两个区间不得重复或相交。
          require(g.random[0].second <= g.lower && g.random[1].first >= g.upper,
                  "Random windows must be ordered left/right, outside prompt");
          g.record = w;
          g.record["ge"] = ge;
          g.record["reference"] = view;
          g.record["selection_status"] = row.at("status");
        }
        require(pair[1].lower <= pair[0].lower &&
                    pair[1].upper >= pair[0].upper,
                "A/B windows not nested");
        rows.emplace(key, std::move(pair));
      }
    require(rows.size() == 40, "Missing Ge/reference windows");
    for (const auto &row : rows)
      require(rows.count({row.first.first, 1 - row.first.second}),
              "Missing companion Ge/reference window");
  }
  const Gate &at(int ge, int view, int label) const {
    return rows.at({ge, view})[label];
  }
  json manifest(int label) const {
    json a = json::array();
    for (const auto &p : rows)
      a.push_back(p.second[label].record);
    return a;
  }
};
struct Particle {
  int key = 0, view = 0;
  double e = 0., t = 0., dt = 0.;
  int64_t entry = 0;
  uint8_t flags = 0;
  int sector = -1;
};
struct Event {
  int key;
  double energy, time, reference_time;
  int64_t entry, native, reference;
  uint8_t flags;
  std::vector<Particle> particles, sectors;
};
// 每个输出文件（单轮、单候选）独立去重，身份为原始Ge entry + 粒子entry +
// Si/S3。
// reference_entry只是分组锚点，不能作为物理配对的身份；同一Ge的新粒子必须保留。
struct PairDeduplicator {
  struct Pair {
    Particle particle, sector;
    int side;
  };
  struct Photon {
    int key;
    double energy, time;
    int64_t native;
    uint8_t flags;
    std::vector<Pair> pairs;
  };
  std::map<int64_t, Photon> active;
  std::multimap<double, int64_t> expiry;
  double last_reference = -std::numeric_limits<double>::infinity();
  int64_t considered[2][3]{}, removed[2][3]{};
  size_t peak_photons = 0;

  void begin(const Event &e) {
    require(e.reference_time >= last_reference,
            "Reference time decreased during deduplication");
    require(std::abs(e.time - e.reference_time) <= 100.,
            "Ge outside reference grouping bound");
    last_reference = e.reference_time;
    // Ge分组为[reference-100, reference+100]。后续reference单调递增，
    // Ge时间严格小于reference-100时不可能再出现，因此缓存不随整轮长度增长。
    while (!expiry.empty() && expiry.begin()->first < e.reference_time - 100.) {
      active.erase(expiry.begin()->second);
      expiry.erase(expiry.begin());
    }
    auto result = active.emplace(
        e.entry, Photon{e.key, e.energy, e.time, e.native, e.flags, {}});
    if (result.second)
      expiry.emplace(e.time, e.entry);
    const auto &g = result.first->second;
    require(g.key == e.key && g.energy == e.energy && g.time == e.time &&
                g.native == e.native && g.flags == e.flags,
            "Repeated Ge source entry has inconsistent data");
    peak_photons = std::max(peak_photons, active.size());
  }
  static bool same_particle(const Particle &a, const Particle &b) {
    // sector是当前分组中sector向量的下标，不属于原始事件身份。
    return a.entry == b.entry && a.key == b.key && a.view == b.view &&
           a.e == b.e && a.t == b.t && a.dt == b.dt && a.flags == b.flags;
  }
  bool accept(const Event &e, const Particle &p, int side) {
    require(p.view == 0 || p.view == 1,
            "Invalid particle view during deduplication");
    require(side >= -1 && side <= 1,
            "Invalid prompt/random side during deduplication");
    const int region = side == 0 ? 0 : (side < 0 ? 1 : 2);
    ++considered[p.view][region];
    auto &pairs = active.at(e.entry).pairs;
    for (const auto &old : pairs) {
      if (old.particle.view != p.view || old.particle.entry != p.entry)
        continue;
      require(
          old.side == side && same_particle(old.particle, p),
          "Repeated Ge-particle pair has inconsistent data or window region");
      if (p.view == 1)
        require(same_particle(old.sector, e.sectors.at(p.sector)),
                "Repeated Ge-S3 ring pair has a different sector match");
      ++removed[p.view][region];
      return false;
    }
    pairs.push_back(
        {p, p.view == 1 ? e.sectors.at(p.sector) : Particle{}, side});
    return true;
  }
  json counts() const {
    json result = {
        {"key", "run,candidate,Ge_source_entry,Si_or_S3,particle_source_entry"},
        {"policy", "keep_first_accepted_pair; reject_inconsistent_repeats"},
        {"peak_cached_Ge", peak_photons}};
    for (int view = 0; view < 2; ++view)
      for (int region = 0; region < 3; ++region) {
        const char *name = region == 0
                               ? "prompt"
                               : (region == 1 ? "random_left" : "random_right");
        result[view ? "s3" : "si"][name] = {
            {"considered", considered[view][region]},
            {"duplicates_removed", removed[view][region]},
            {"written", considered[view][region] - removed[view][region]}};
      }
    return result;
  }
};
struct Family {
  int n = 0, capacity;
  std::vector<Short_t> sid, ch, ring, sector, side;
  std::vector<Double_t> energy, time, dt, alpha;
  std::vector<Long64_t> ts, entry, native;
  std::vector<UChar_t> flags;
  explicit Family(int cap)
      : capacity(cap), sid(cap), ch(cap), ring(cap), sector(cap), side(cap),
        energy(cap), time(cap), dt(cap), alpha(cap), ts(cap), entry(cap),
        native(cap), flags(cap) {}
  void branches(TTree &tree, const std::string &name) {
    const std::string count = "n_" + name;
    tree.Branch(count.c_str(), &n, (count + "/I").c_str());
    auto field = [&](const std::string &field, void *pointer,
                     const std::string &type) {
      const auto full = name + "_" + field;
      tree.Branch(full.c_str(), pointer,
                  (full + "[" + count + "]/" + type).c_str());
    };
    field("sid", sid.data(), "S");
    field("ch", ch.data(), "S");
    if (name == "ge" || name == "spider") {
      field("ring_id", ring.data(), "S");
      field("sector_id", sector.data(), "S");
    } else
      field("id", name == "s3_ring" ? ring.data() : sector.data(), "S");
    field("energy", energy.data(), "D");
    field("ts", ts.data(), "L");
    field("time_ns", time.data(), "D");
    field("source_entry", entry.data(), "L");
    field("flags", flags.data(), "b");
    if (name == "ge")
      field("native_time_ns", native.data(), "L");
    else {
      field("dt_ns", dt.data(), "D");
      field("alpha", alpha.data(), "D");
      field("random_side", side.data(), "S");
    }
  }
  void put(const Parameters &pars, int key, double e, double t, int64_t source,
           uint8_t f, double difference = 0., double a = 1., int s = 0) {
    require(n < capacity, "Output multiplicity exceeds capacity");
    const int i = n++;
    sid[i] = key % 10000 / 100;
    ch[i] = key % 100;
    ring[i] = pars.q[key * 5 + 1];
    sector[i] = pars.q[key * 5 + 2];
    energy[i] = e;
    time[i] = t;
    ts[i] = int64_t(std::floor(t));
    entry[i] = source;
    flags[i] = f;
    dt[i] = difference;
    alpha[i] = a;
    side[i] = s;
  }
};
struct TreeWriter {
  TTree *tree;
  int run = 0;
  Long64_t reference = 0;
  Double_t reference_time = 0.;
  Family ge{1}, si{96}, ring{24}, sector{32};
  int64_t events = 0, si_pairs = 0, s3_pairs = 0, max_si = 0;
  TreeWriter(TFile &file, const char *name, int r) : run(r) {
    file.cd();
    tree = new TTree(name, name);
    tree->Branch("run", &run, "run/I");
    tree->Branch("reference_entry", &reference, "reference_entry/L");
    tree->Branch("reference_group_time_ns", &reference_time,
                 "reference_group_time_ns/D");
    ge.branches(*tree, "ge");
    si.branches(*tree, "spider");
    ring.branches(*tree, "s3_ring");
    sector.branches(*tree, "s3_sector");
  }
  void fill(const Parameters &pars, const Gates &gates, int label,
            const Event &event,
            const std::vector<std::pair<size_t, int>> &selected) {
    ge.n = si.n = ring.n = sector.n = 0;
    reference = event.reference;
    reference_time = event.reference_time;
    ge.put(pars, event.key, event.energy, event.time, event.entry, event.flags);
    ge.native[0] = event.native;
    for (const auto &choice : selected) {
      const auto &p = event.particles[choice.first];
      const auto &gate = gates.at(event.key, p.view, label);
      (p.view == 0 ? si : ring)
          .put(pars, p.key, p.e, p.t, p.entry, p.flags, p.dt, gate.alpha,
               choice.second);
      if (p.view == 1) {
        const auto &s = event.sectors.at(p.sector);
        sector.put(pars, s.key, s.e, s.t, s.entry, s.flags, s.dt, gate.alpha,
                   choice.second);
      }
    }
    require(si.n + ring.n > 0 && tree->Fill() > 0, "ROOT event write failed");
    ++events;
    si_pairs += si.n;
    s3_pairs += ring.n;
    max_si = std::max(max_si, int64_t(si.n));
  }
  json counts() const {
    return {{"events", events},
            {"si_pairs", si_pairs},
            {"s3_pairs", s3_pairs},
            {"max_si_multiplicity", max_si}};
  }
};
struct FileWriter {
  std::filesystem::path path, partial;
  std::unique_ptr<TFile> file;
  std::unique_ptr<TreeWriter> prompt, random;
  const Parameters &pars;
  const Gates &gates;
  int label;
  PairDeduplicator unique_pairs;
  bool closed = false;
  FileWriter(const std::filesystem::path &target, const Parameters &p,
             const Gates &g, int l, const json &manifest)
      : path(target), partial(target.string() + ".partial"), pars(p), gates(g),
        label(l) {
    require(!std::filesystem::exists(path) && !std::filesystem::exists(partial),
            "Refuse to overwrite output or partial");
    std::filesystem::create_directories(path.parent_path());
    file.reset(TFile::Open(partial.c_str(), "CREATE"));
    require(file && !file->IsZombie(), "Cannot create ROOT output");
    file->SetCompressionSettings(504);
    prompt =
        std::make_unique<TreeWriter>(*file, "tr_event", std::stoi(pars.run));
    random = std::make_unique<TreeWriter>(*file, "tr_bg", std::stoi(pars.run));
    file->cd();
    require(TObjString(manifest.dump().c_str()).Write("sm150_build_manifest") >
                0,
            "ROOT provenance write failed");
  }
  void consume(const Event &e) {
    unique_pairs.begin(e);
    std::vector<std::pair<size_t, int>> a, b;
    for (size_t i = 0; i < e.particles.size(); ++i) {
      const auto &p = e.particles[i];
      const auto &g = gates.at(e.key, p.view, label);
      if (g.lower <= p.dt && p.dt < g.upper && unique_pairs.accept(e, p, 0))
        a.emplace_back(i, 0);
      int side = g.side(p.dt);
      if (side && unique_pairs.accept(e, p, side))
        b.emplace_back(i, side);
    }
    if (!a.empty())
      prompt->fill(pars, gates, label, e, a);
    if (!b.empty())
      random->fill(pars, gates, label, e, b);
  }
  json close() {
    require(!closed, "Output closed twice");
    file->cd();
    json c = {{"tr_event", prompt->counts()}, {"tr_bg", random->counts()}};
    require(TObjString(c.dump().c_str()).Write("sm150_build_counts") > 0 &&
                TObjString(unique_pairs.counts().dump().c_str())
                        .Write("sm150_build_pair_deduplication") > 0,
            "ROOT count metadata write failed");
    require(file->Write() > 0, "ROOT file write failed");
    file->Close();
    require(!file->TestBit(TFile::kWriteError), "ROOT file close/write error");
    closed = true;
    // 保留.partial；两套输出和冻结输入均通过核验后，Process才发布正式文件。
    return c;
  }
  ~FileWriter() {
    if (file && file->IsOpen())
      file->Close();
  }
};
inline int64_t classify(const sm150_windowed::Engine &cache,
                        const Parameters &pars, const Gates &gates,
                        FileWriter &a, FileWriter &b) {
  int64_t accepted = 0;
  for (size_t i = 0; i < cache.ref_entry.size(); ++i) {
    if (cache.goff[i + 1] - cache.goff[i] != 1)
      continue;
    const size_t gi = cache.goff[i];
    int key = cache.ge_key[gi];
    double energy = cache.ge_energy[gi];
    if (cache.ge_second[gi] >= 0 || (cache.ge_flags[gi] & 3) ||
        !gates.rows.count({key, 0}) || energy < 40. || energy >= 2000.)
      continue;
    Event e;
    e.key = key;
    e.energy = energy;
    e.entry = cache.ge_entry[gi];
    e.flags = cache.ge_flags[gi];
    e.reference = cache.ref_entry[i];
    e.reference_time = cache.ref_time[i];
    e.time = e.reference_time + cache.ge_dt[gi];
    const double correction =
        cache.timing.TotalGeOffset(key, (e.flags & 4) != 0, energy);
    e.native = std::llrint(e.time - correction);
    std::vector<Particle> rings;
    for (int64_t j = cache.poff[i]; j < cache.poff[i + 1]; ++j) {
      int pk = cache.p_key[j], kind = pars.q[pk * 5];
      double en = cache.p_energy[j];
      if ((cache.p_flags[j] & 3) || kind < 2 ||
          en <= (kind == 2 ? 15000. : 25000.) || en >= 100000.)
        continue;
      if (kind == 3 && (pars.q[pk * 5 + 1] < 1 ||
                        pars.q[pk * 5 + 1] > (pars.group == "low" ? 24 : 8)))
        continue;
      Particle p;
      p.key = pk;
      p.e = en;
      p.t = e.reference_time + cache.p_dt[j];
      p.dt = cache.p_dt[j] - cache.ge_dt[gi];
      p.entry = cache.p_entry[j];
      p.flags = cache.p_flags[j];
      p.view = kind - 2;
      if (kind == 2)
        e.particles.push_back(p);
      else if (kind == 3)
        rings.push_back(p);
      else
        e.sectors.push_back(p);
    }
    size_t match_ring = 0, match_sector = 0, matches = 0;
    for (size_t r = 0; r < rings.size(); ++r)
      for (size_t s = 0; s < e.sectors.size(); ++s)
        if (std::abs(rings[r].e - e.sectors[s].e) < 2500. &&
            std::abs(rings[r].dt - e.sectors[s].dt) < 80.) {
          ++matches;
          match_ring = r;
          match_sector = s;
        }
    if (matches == 1) {
      Particle r = rings[match_ring];
      r.view = 1;
      r.sector = match_sector;
      e.particles.push_back(r);
    }
    if (!e.particles.empty()) {
      a.consume(e);
      b.consume(e);
      ++accepted;
    }
  }
  return accepted;
}
struct ParticleReader {
  TTreeReaderValue<Int_t> n;
  TTreeReaderArray<Short_t> sid, ch, id, side;
  TTreeReaderArray<Double_t> energy, dt, alpha;
  TTreeReaderArray<Long64_t> entry;
  ParticleReader(TTreeReader &r, const std::string &name)
      : n(r, ("n_" + name).c_str()), sid(r, (name + "_sid").c_str()),
        ch(r, (name + "_ch").c_str()),
        id(r, (name + (name == "spider" ? "_ring_id" : "_id")).c_str()),
        side(r, (name + "_random_side").c_str()),
        energy(r, (name + "_energy").c_str()), dt(r, (name + "_dt_ns").c_str()),
        alpha(r, (name + "_alpha").c_str()),
        entry(r, (name + "_source_entry").c_str()) {}
};
inline json verify_root(const std::filesystem::path &path,
                        const Parameters &pars, const Gates &gates, int label,
                        const json &expected) {
  std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "READ"));
  require(file && !file->IsZombie(), "Cannot reopen output");
  json result;
  auto *metadata =
      dynamic_cast<TObjString *>(file->Get("sm150_build_manifest"));
  require(metadata, "Missing embedded provenance");
  auto manifest = json::parse(metadata->GetString().Data());
  require(manifest.at("run") == pars.run &&
              manifest.at("candidate") == std::string(label ? "B" : "A"),
          "Wrong output provenance");
  // 与写出缓存独立核验：完整读回P/R两棵树后排序原始配对身份。
  // 身份不含prompt/random类别，避免同一配对误入两类而漏检。
  std::array<std::vector<std::pair<int64_t, int64_t>>, 2> pair_ids;
  for (const std::string name : {"tr_event", "tr_bg"}) {
    auto *tree = dynamic_cast<TTree *>(file->Get(name.c_str()));
    require(tree, "Missing output tree");
    TTreeReader reader(tree);
    TTreeReaderValue<Int_t> ng(reader, "n_ge"), ns(reader, "n_s3_sector"),
        run(reader, "run");
    TTreeReaderArray<Short_t> sid(reader, "ge_sid"), ch(reader, "ge_ch");
    TTreeReaderArray<Double_t> energy(reader, "ge_energy");
    TTreeReaderArray<Long64_t> ge_entry(reader, "ge_source_entry");
    ParticleReader si(reader, "spider"), ring(reader, "s3_ring");
    TTreeReaderArray<Double_t> se(reader, "s3_sector_energy"),
        st(reader, "s3_sector_dt_ns");
    int64_t events = 0, nsi = 0, nring = 0, max_si = 0;
    while (reader.Next()) {
      require(*run == std::stoi(pars.run), "Wrong run in output row");
      require(*ng == 1 && sid.GetSize() == 1 && ch.GetSize() == 1 &&
                  energy.GetSize() == 1 && ge_entry.GetSize() == 1,
              "Ge multiplicity mismatch");
      require(ge_entry[0] >= 0, "Invalid original Ge entry");
      int key = 100 * sid[0] + ch[0];
      require(energy[0] >= 40. && energy[0] < 2000., "Ge energy gate mismatch");
      for (int view = 0; view < 2; ++view) {
        auto &p = view ? ring : si;
        const auto count = p.energy.GetSize();
        require(*p.n >= 0 && size_t(*p.n) == count && p.dt.GetSize() == count &&
                    p.entry.GetSize() == count && p.sid.GetSize() == count &&
                    p.ch.GetSize() == count && p.id.GetSize() == count &&
                    p.side.GetSize() == count && p.alpha.GetSize() == count,
                "Particle array count mismatch");
        const auto &gate = gates.at(key, view, label);
        for (int i = 0; i < *p.n; ++i) {
          require(p.entry[i] >= 0, "Invalid original particle entry");
          pair_ids[view].emplace_back(ge_entry[0], p.entry[i]);
          require(p.energy[i] > (view ? 25000. : 15000.) &&
                      p.energy[i] < 100000.,
                  "Particle energy gate mismatch");
          require(p.sid[i] >= 0 && p.sid[i] < 16 && p.ch[i] >= 0 &&
                      p.ch[i] < 16 && p.alpha[i] == gate.alpha,
                  "Particle channel/alpha mismatch");
          if (name == "tr_event")
            require(p.side[i] == 0 && p.dt[i] >= gate.lower &&
                        p.dt[i] < gate.upper,
                    "Prompt gate mismatch");
          else
            require(p.side[i] != 0 && gate.side(p.dt[i]) == p.side[i],
                    "Random gate mismatch");
          if (view)
            require(p.id[i] >= 1 && p.id[i] <= (pars.group == "low" ? 24 : 8),
                    "S3 ring selection mismatch");
        }
      }
      require(*ring.n <= 1 && *si.n + *ring.n > 0,
              "Invalid output particle multiplicity");
      require(*ring.n == *ns && se.GetSize() == size_t(*ns) &&
                  st.GetSize() == size_t(*ns),
              "S3 face multiplicity mismatch");
      for (int i = 0; i < *ns; ++i)
        require(std::abs(ring.energy[i] - se[i]) < 2500. &&
                    std::abs(ring.dt[i] - st[i]) < 80.,
                "S3 face agreement mismatch");
      ++events;
      nsi += *si.n;
      nring += *ring.n;
      max_si = std::max(max_si, int64_t(*si.n));
    }
    require(reader.GetEntryStatus() == TTreeReader::kEntryBeyondEnd,
            "ROOT readback did not finish");
    require(events == tree->GetEntries() &&
                events == expected.at(name).at("events").get<int64_t>() &&
                nsi == expected.at(name).at("si_pairs").get<int64_t>() &&
                nring == expected.at(name).at("s3_pairs").get<int64_t>() &&
                max_si ==
                    expected.at(name).at("max_si_multiplicity").get<int64_t>(),
            "Output count mismatch");
    result[name] = {{"events", events},
                    {"si_pairs", nsi},
                    {"s3_pairs", nring},
                    {"pairs", nsi + nring},
                    {"max_si_multiplicity", max_si}};
  }
  for (int view = 0; view < 2; ++view) {
    auto &ids = pair_ids[view];
    std::sort(ids.begin(), ids.end());
    auto duplicate = std::adjacent_find(ids.begin(), ids.end());
    if (duplicate != ids.end())
      throw std::runtime_error(
          "Duplicate Ge-" + std::string(view ? "S3" : "Si") +
          " pair in output: Ge source entry " +
          std::to_string(duplicate->first) + ", particle source entry " +
          std::to_string(duplicate->second));
  }
  return result;
}
} // namespace sm150_native

// ==========================================================================
// 5. 原始ROOT批量读取（保留已验证的高效basket读取方式）
// ==========================================================================
using namespace sm150_native;

static double seconds(Clock::time_point start) {
  return std::chrono::duration<double>(Clock::now() - start).count();
}
static json read_json(const fs::path &path) {
  std::ifstream s(path);
  require(bool(s), "Missing JSON: " + path.string());
  json j;
  s >> j;
  return j;
}
static void publish_file(const fs::path &partial, const fs::path &target);
static void save_json(const fs::path &path, const json &value,
                      bool no_replace = false) {
  fs::create_directories(path.parent_path());
  const auto temp = fs::path(path.string() + ".tmp");
  if (no_replace)
    require(!fs::exists(temp), "Existing partial JSON receipt");
  {
    std::ofstream s(temp);
    require(bool(s), "Cannot write JSON");
    s << value.dump(2) << '\n';
    s.flush();
    require(bool(s), "JSON write failed");
  }
  if (no_replace)
    publish_file(temp, path);
  else
    fs::rename(temp, path);
}
static std::string
sha256(const fs::path &path,
       const std::function<void(uint64_t, uint64_t)> &progress) {
  std::ifstream s(path, std::ios::binary);
  require(bool(s), "Cannot hash file: " + path.string());
  std::unique_ptr<EVP_MD_CTX, decltype(&EVP_MD_CTX_free)> context(
      EVP_MD_CTX_new(), &EVP_MD_CTX_free);
  require(bool(context) &&
              EVP_DigestInit_ex(context.get(), EVP_sha256(), nullptr) == 1,
          "SHA initialization failed");
  std::vector<char> buffer(4 * 1024 * 1024);
  uint64_t done = 0, total = fs::file_size(path);
  while (s) {
    s.read(buffer.data(), buffer.size());
    auto n = s.gcount();
    if (n) {
      require(EVP_DigestUpdate(context.get(), buffer.data(), n) == 1,
              "SHA update failed");
      done += n;
      if (progress)
        progress(done, total);
    }
  }
  require(s.eof() && done == total, "Incomplete hash read");
  unsigned char digest[EVP_MAX_MD_SIZE];
  unsigned int n = 0;
  require(EVP_DigestFinal_ex(context.get(), digest, &n) == 1 && n == 32,
          "SHA completion failed");
  std::ostringstream out;
  for (unsigned int i = 0; i < n; ++i)
    out << std::hex << std::setw(2) << std::setfill('0') << int(digest[i]);
  return out.str();
}

// 输入在加载时冻结，避免长时间原始文件SHA核验期间的参数变动被错误记为基线。
static json frozen_json(const fs::path &path, json &hashes) {
  const auto hash = sha256(path);
  auto value = read_json(path);
  require(sha256(path) == hash, "JSON changed while loading: " + path.string());
  hashes[path.string()] = hash;
  return value;
}
static void verify_frozen(const json &hashes) {
  for (auto it = hashes.begin(); it != hashes.end(); ++it)
    require(sha256(it.key()) == it.value(),
            "Frozen source/parameter changed: " + it.key());
}
// flock随进程退出自动释放；空锁文件保留，不能用删除锁文件的方式抢占活动任务。
struct OutputLock {
  int fd = -1;
  explicit OutputLock(const fs::path &path) {
    fs::create_directories(path.parent_path());
    fd = ::open(path.c_str(), O_RDWR | O_CREAT | O_CLOEXEC, 0644);
    require(fd >= 0, "Cannot open output lock: " + path.string());
    if (::flock(fd, LOCK_EX | LOCK_NB) != 0) {
      ::close(fd);
      fd = -1;
      throw std::runtime_error("Output is owned by another build: " +
                               path.string());
    }
  }
  OutputLock(const OutputLock &) = delete;
  OutputLock &operator=(const OutputLock &) = delete;
  ~OutputLock() {
    if (fd >= 0)
      ::close(fd);
  }
};
static void distinct_paths(const std::vector<fs::path> &outputs,
                           const std::vector<fs::path> &inputs) {
  auto same = [](const fs::path &a, const fs::path &b) {
    return fs::weakly_canonical(fs::absolute(a)) ==
               fs::weakly_canonical(fs::absolute(b)) ||
           (fs::exists(a) && fs::exists(b) && fs::equivalent(a, b));
  };
  for (size_t i = 0; i < outputs.size(); ++i) {
    for (size_t j = 0; j < i; ++j)
      require(!same(outputs[i], outputs[j]),
              "Output/receipt path collision: " + outputs[i].string());
    for (const auto &input : inputs)
      require(!same(outputs[i], input),
              "Output would overwrite input: " + outputs[i].string());
  }
}
static void publish_file(const fs::path &partial, const fs::path &target) {
  // 不使用可覆盖目标的普通rename；无论是否有协作锁，已有正式文件都不能覆盖。
  if (::syscall(SYS_renameat2, AT_FDCWD, partial.c_str(), AT_FDCWD,
                target.c_str(), RENAME_NOREPLACE) != 0)
    throw std::runtime_error("Cannot publish " + target.string() + ": " +
                             std::strerror(errno));
}
struct RawBlock {
  std::vector<int16_t> cid, sid, ch, sr, cfd, cfds;
  std::vector<int64_t> ts, source_entry;
  std::vector<uint16_t> evte;
  std::vector<uint8_t> cfdft, pileup, outofr;
  explicit RawBlock(size_t n)
      : cid(n), sid(n), ch(n), sr(n), cfd(n), cfds(n), ts(n), source_entry(n), evte(n), cfdft(n),
        pileup(n), outofr(n) {}
};
// Bulk scalar ROOT reader, introduced after preserving the first comparison.
// ROOT bulk API documentation: https://root.cern.ch/doc/v632/classTBranch.html
// Local ROOT 6.34.10 TBranch.h/TLeaf*.h provide the same checked interfaces.

template <class T> struct BulkColumn {
  TBranch *branch = nullptr;
  TBufferFile buffer{TBuffer::kRead};
  int64_t read_entry = 0;
  size_t available = 0, position = 0;
  const char *data = nullptr;
  void bind(TTree &tree, const char *name, const char *type) {
    branch = tree.GetBranch(name);
    require(branch, std::string("Missing branch: ") + name);
    require(branch->GetListOfLeaves()->GetEntries() == 1,
            "Bulk input requires one scalar leaf per branch");
    auto *leaf = dynamic_cast<TLeaf *>(branch->GetListOfLeaves()->At(0));
    require(leaf && std::string(leaf->GetTypeName()) == type &&
                leaf->GetLeafCount() == nullptr && leaf->GetLenStatic() == 1,
            "Wrong scalar raw branch type");
    require(branch->GetEntries() == tree.GetEntries(),
            "Raw branch entry count mismatch");
    require(branch->GetBulkRead().SupportsBulkRead(),
            "ROOT does not support bulk input for this branch");
    tree.SetBranchStatus(name, 1);
    tree.AddBranchToCache(name, true);
  }
  void read(int64_t begin, size_t count, T *output) {
    require(begin == read_entry, "Non-contiguous bulk input interval");
    size_t done = 0;
    while (done < count) {
      if (position == available) {
        const int loaded =
            branch->GetBulkRead().GetBulkEntries(read_entry, buffer);
        require(loaded > 0 && read_entry + loaded <= branch->GetEntries(),
                "Invalid bulk ROOT basket read");
        available = size_t(loaded);
        position = 0;
        data = buffer.GetCurrent();
        require(data != nullptr && buffer.Length() + available * sizeof(T) <=
                                       size_t(buffer.BufferSize()),
                "Short bulk ROOT buffer");
      }
      const size_t take = std::min(count - done, available - position);
      std::memcpy(output + done, data + position * sizeof(T), take * sizeof(T));
      done += take;
      position += take;
      read_entry += take;
    }
  }
};

struct RawReader {
  TTree &tree;
  BulkColumn<int16_t> cid, sid, ch, sr, cfd, cfds;
  BulkColumn<int64_t> ts, source_entry;
  BulkColumn<uint16_t> evte;
  BulkColumn<uint8_t> cfdft, pileup, outofr;
  bool original_ids;
  RawReader(TTree &tr, int64_t end, bool preserve_ids = false)
      : tree(tr), original_ids(preserve_ids) {
    tree.SetBranchStatus("*", 0);
    tree.SetCacheSize(64 * 1024 * 1024);
    tree.SetCacheEntryRange(0, end);
    cid.bind(tree, "cid", "Short_t");
    sid.bind(tree, "sid", "Short_t");
    ch.bind(tree, "ch", "Short_t");
    sr.bind(tree, "sr", "Short_t");
    cfd.bind(tree, "cfd", "Short_t");
    cfds.bind(tree, "cfds", "Short_t");
    ts.bind(tree, "ts", "Long64_t");
    evte.bind(tree, "evte", "UShort_t");
    cfdft.bind(tree, "cfdft", "Bool_t");
    pileup.bind(tree, "pileup", "Bool_t");
    outofr.bind(tree, "outofr", "Bool_t");
    if (original_ids)
      source_entry.bind(tree, "source_entry", "Long64_t");
    tree.StopCacheLearningPhase();
  }
  void read(int64_t begin, size_t n, RawBlock &b) {
    cid.read(begin, n, b.cid.data());
    sid.read(begin, n, b.sid.data());
    ch.read(begin, n, b.ch.data());
    sr.read(begin, n, b.sr.data());
    ts.read(begin, n, b.ts.data());
    evte.read(begin, n, b.evte.data());
    cfd.read(begin, n, b.cfd.data());
    cfds.read(begin, n, b.cfds.data());
    cfdft.read(begin, n, b.cfdft.data());
    pileup.read(begin, n, b.pileup.data());
    outofr.read(begin, n, b.outofr.data());
    if (original_ids)
      source_entry.read(begin, n, b.source_entry.data());
  }
};

// ==========================================================================
// 6. 可重放的宽窗缓存。每个原始命中只存一次，保留source_entry及全部原始字段。
// 核心窗为任意Ge附近±1000ns；另存100ns边界余量。全部Ge均保留，
// 不在此处施加能量门、质量位、单Ge、addback或S3配对筛选。
// ==========================================================================
namespace sm150_cache {
constexpr double core_ns = 1000., guard_ns = 100., bound_ns = 1000000.;
struct Record {
  int16_t cid = 0, sid = 0, ch = 0, sr = 0, cfd = 0, cfds = 0;
  int64_t ts = 0, source_entry = 0;
  uint16_t evte = 0;
  uint8_t cfdft = 0, pileup = 0, outofr = 0, kind = 0, in_core = 0;
  double energy_keV = 0., corrected_time_ns = 0.;
  int8_t decision = 0; // 0等待后续Ge；1保留；-1丢弃
};
static Record calibrate(const RawBlock &b, size_t i, int64_t entry,
                        const Parameters &pars,
                        const sm150::TimeAlignment &timing) {
  Record r;
  r.cid = b.cid[i]; r.sid = b.sid[i]; r.ch = b.ch[i]; r.sr = b.sr[i];
  r.ts = b.ts[i]; r.evte = b.evte[i]; r.cfd = b.cfd[i]; r.cfds = b.cfds[i];
  r.cfdft = b.cfdft[i]; r.pileup = b.pileup[i]; r.outofr = b.outofr[i];
  r.source_entry = entry;
  require(r.sr == 100 || r.sr == 250, "Unknown cache sampling rate");
  require(r.ts >= 0 && r.ts <= std::numeric_limits<int64_t>::max() / 10,
          "Cache timestamp outside supported range");
  const int key = 10000 * r.cid + 100 * r.sid + r.ch;
  require(key >= 0 && key < 16384 && pars.q[key * 5],
          "Unmapped cache channel");
  r.kind = pars.q[key * 5];
  const double *p = &pars.p[key * 9];
  double x = r.evte, e = p[0] + p[1] * x + p[2] * x * x;
  e *= p[3] * p[4] * p[5];
  r.energy_keV = p[6] + p[7] * e;
  const int64_t coarse = r.ts * (r.sr == 250 ? 8 : 10);
  double t = coarse;
  if (r.cid == 0) {
    require(r.sr == 250, "Ge cache sampling rate is not 250");
    t = r.cfdft ? double(8 * r.ts)
                : (2 * r.ts - r.cfds + r.cfd / 16384.) * 4.;
  }
  // 与现行build完全一致，先截断原生CFD时间，再加完整时间项。
  // 未截断的CFD、cfds、ts仍单独保存，可供后续方法研究。
  r.corrected_time_ns =
      double(int64_t(t)) +
      (r.kind == 1 ? timing.TotalGeOffset(key, r.cfdft != 0, r.energy_keV)
                   : timing.TotalParticleOffsetFromAdc(key, r.evte));
  require(std::isfinite(r.energy_keV) && std::isfinite(r.corrected_time_ns) &&
              std::abs(r.corrected_time_ns - coarse) <= bound_ns,
          "Cache calibration/timestamp bound violated");
  return r;
}
// 仅将字段本身按固定顺序送入SHA，不包含C++结构体的padding。
struct Digest {
  std::unique_ptr<EVP_MD_CTX, decltype(&EVP_MD_CTX_free)> ctx{
      EVP_MD_CTX_new(), &EVP_MD_CTX_free};
  Digest() {
    require(ctx && EVP_DigestInit_ex(ctx.get(), EVP_sha256(), nullptr) == 1,
            "Cache digest initialization failed");
  }
  void add(const Record &r) {
    std::array<char, 51> bytes{};
    size_t offset = 0;
    auto pack = [&](const auto &value) {
      std::memcpy(bytes.data() + offset, &value, sizeof(value));
      offset += sizeof(value);
    };
    pack(r.cid); pack(r.sid); pack(r.ch); pack(r.sr); pack(r.cfd); pack(r.cfds);
    pack(r.ts); pack(r.evte); pack(r.cfdft); pack(r.pileup); pack(r.outofr);
    pack(r.source_entry); pack(r.energy_keV); pack(r.corrected_time_ns);
    pack(r.kind); pack(r.in_core);
    require(offset == bytes.size() &&
                EVP_DigestUpdate(ctx.get(), bytes.data(), bytes.size()) == 1,
            "Cache payload digest failed");
  }
  std::string finish() {
    unsigned char bytes[EVP_MAX_MD_SIZE];
    unsigned int n = 0;
    require(EVP_DigestFinal_ex(ctx.get(), bytes, &n) == 1 && n == 32,
            "Cache digest completion failed");
    std::ostringstream s;
    for (unsigned i = 0; i < n; ++i)
      s << std::hex << std::setw(2) << std::setfill('0') << int(bytes[i]);
    return s.str();
  }
};
struct Selector {
  double core_half_width, guard_width;
  std::vector<Record> pending;
  std::vector<double> anchors;
  std::array<uint64_t, 5> input_counts{}, kept_counts{};
  uint64_t guard_particles = 0;
  int64_t last_coarse = std::numeric_limits<int64_t>::min(), last_entry = -1;
  // 默认行为保持1000+100ns；较窄范围供独立体积试验复用同一选择器。
  explicit Selector(double core = core_ns, double guard = guard_ns)
      : core_half_width(core), guard_width(guard) {
    require(std::isfinite(core) && std::isfinite(guard) && core > 0. &&
                guard >= 0. && core + guard <= core_ns + guard_ns,
            "Invalid cache selector core/guard width");
  }
  void add(Record r) {
    const int64_t coarse = r.ts * (r.sr == 250 ? 8 : 10);
    require(coarse >= last_coarse && r.source_entry > last_entry,
            "Cache input timestamp/order decreased or duplicate source entry");
    last_coarse = coarse; last_entry = r.source_entry;
    require(r.kind >= 1 && r.kind <= 4, "Invalid cache detector kind");
    ++input_counts[r.kind];
    if (r.kind == 1) {
      r.decision = 1;
      r.in_core = 1;
      anchors.push_back(r.corrected_time_ns);
    }
    pending.push_back(r);
  }
  template <class Sink> void flush(bool final, Sink sink) {
    std::sort(anchors.begin(), anchors.end());
    const double radius = core_half_width + guard_width;
    const double cutoff = final ? std::numeric_limits<double>::infinity()
                                : double(last_coarse) - bound_ns - radius;
    for (auto &r : pending) {
      if (r.decision || r.corrected_time_ns >= cutoff)
        continue;
      auto next = std::lower_bound(anchors.begin(), anchors.end(),
                                   r.corrected_time_ns);
      double distance = std::numeric_limits<double>::infinity();
      if (next != anchors.end())
        distance = std::abs(*next - r.corrected_time_ns);
      if (next != anchors.begin())
        distance = std::min(distance,
                            std::abs(*std::prev(next) - r.corrected_time_ns));
      r.decision = distance <= radius ? 1 : -1;
      r.in_core = distance <= core_half_width;
    }
    size_t consumed = 0;
    while (consumed < pending.size() && pending[consumed].decision) {
      const auto &r = pending[consumed++];
      if (r.decision == 1) {
        sink(r);
        ++kept_counts[r.kind];
        if (r.kind != 1 && !r.in_core)
          ++guard_particles;
      }
    }
    pending.erase(pending.begin(), pending.begin() + consumed);
    anchors.erase(anchors.begin(),
                  std::lower_bound(anchors.begin(), anchors.end(),
                                   cutoff - radius));
    if (final)
      require(pending.empty() && anchors.empty(), "Incomplete cache flush");
  }
  json counts() const {
    return {{"kind_order", {"unused", "Ge", "Si", "S3_ring", "S3_sector"}},
            {"input", input_counts}, {"retained", kept_counts},
            {"guard_only_particles", guard_particles}};
  }
};
struct Writer {
  fs::path path, partial;
  std::unique_ptr<TFile> file;
  TTree *tree = nullptr;
  Record row;
  Digest digest;
  int64_t count = 0;
  Writer(const fs::path &p) : path(p), partial(p.string() + ".partial") {
    file.reset(TFile::Open(partial.c_str(), "CREATE"));
    require(file && !file->IsZombie(), "Cannot create cache partial");
    file->SetCompressionSettings(504);
    file->cd();
    tree = new TTree("tree", "Unique raw hits: all Ge and nearby particles");
    tree->SetAutoFlush(-32 * 1024 * 1024);
#define CACHE_BRANCH(name, type) tree->Branch(#name, &row.name, #name "/" type)
    CACHE_BRANCH(cid, "S"); CACHE_BRANCH(sid, "S"); CACHE_BRANCH(ch, "S");
    CACHE_BRANCH(sr, "S"); CACHE_BRANCH(cfd, "S"); CACHE_BRANCH(cfds, "S");
    CACHE_BRANCH(ts, "L"); CACHE_BRANCH(evte, "s");
    CACHE_BRANCH(cfdft, "O"); CACHE_BRANCH(pileup, "O"); CACHE_BRANCH(outofr, "O");
    CACHE_BRANCH(source_entry, "L"); CACHE_BRANCH(energy_keV, "D");
    CACHE_BRANCH(corrected_time_ns, "D"); CACHE_BRANCH(kind, "b");
    CACHE_BRANCH(in_core, "O");
#undef CACHE_BRANCH
  }
  void write(const Record &r) {
    row = r;
    require(tree->Fill() >= 0 && !file->TestBit(TFile::kWriteError),
            "Cache ROOT write failed");
    digest.add(r); ++count;
  }
  json close(json metadata, const Selector &selector) {
    metadata["retained_entries"] = count;
    metadata["counts"] = selector.counts();
    metadata["payload_sha256"] = digest.finish();
    metadata["payload_encoding"] =
        "fixed scalar field order in source; little-endian, no padding";
    file->cd();
    require(tree->Write("", TObject::kOverwrite) > 0, "Cache tree write failed");
    TObjString text(metadata.dump().c_str());
    require(text.Write("sm150_cache_manifest") > 0, "Cache manifest write failed");
    file->Flush();
    require(!file->TestBit(TFile::kWriteError), "Cache flush failed");
    file->Close(); tree = nullptr;
    return metadata;
  }
  ~Writer() { if (file && file->IsOpen()) file->Close(); }
};
// 完整读回所有字段；核对保存前/后SHA、原始ID唯一、刻度及核心/余量归属。
static json verify(const fs::path &path, const json &metadata,
                   const Parameters &pars, const sm150::TimeAlignment &timing,
                   size_t chunk) {
  std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "READ"));
  require(file && !file->IsZombie(), "Cannot read back cache");
  auto *tree = dynamic_cast<TTree *>(file->Get("tree"));
  auto *object = dynamic_cast<TObjString *>(file->Get("sm150_cache_manifest"));
  require(tree && object && json::parse(object->GetString().Data()) == metadata,
          "Cache metadata readback mismatch");
  int64_t entries = tree->GetEntries();
  require(entries == metadata.at("retained_entries"), "Cache count mismatch");
  RawReader reader(*tree, entries, true);
  BulkColumn<double> energy, time;
  BulkColumn<uint8_t> kind, core;
  energy.bind(*tree, "energy_keV", "Double_t");
  time.bind(*tree, "corrected_time_ns", "Double_t");
  kind.bind(*tree, "kind", "UChar_t"); core.bind(*tree, "in_core", "Bool_t");
  std::vector<double> energies(chunk), times(chunk);
  std::vector<uint8_t> kinds(chunk), cores(chunk);
  RawBlock raw(chunk);
  Digest disk_digest, recomputed_digest;
  Selector selector(metadata.at("core_half_width_ns").get<double>(),
                    metadata.at("guard_ns").get<double>());
  int64_t kept = 0;
  const int64_t source_stop = metadata.at("source").at("coverage").at("stop");
  auto sink = [&](const Record &r) { recomputed_digest.add(r); ++kept; };
  for (int64_t begin = 0; begin < entries;) {
    const size_t n = std::min<int64_t>(chunk, entries - begin);
    reader.read(begin, n, raw);
    energy.read(begin, n, energies.data()); time.read(begin, n, times.data());
    kind.read(begin, n, kinds.data()); core.read(begin, n, cores.data());
    for (size_t i = 0; i < n; ++i) {
      require(raw.source_entry[i] >= 0 && raw.source_entry[i] < source_stop,
              "Cache source identity out of range");
      Record r = calibrate(raw, i, raw.source_entry[i], pars, timing);
      require(r.energy_keV == energies[i] && r.corrected_time_ns == times[i] &&
                  r.kind == kinds[i] && cores[i] <= 1,
              "Cache calibrated fields disagree with raw fields");
      r.in_core = cores[i];
      disk_digest.add(r); selector.add(r);
    }
    begin += n;
    selector.flush(begin == entries, sink);
  }
  auto disk_sha = disk_digest.finish(), selected_sha = recomputed_digest.finish();
  require(disk_sha == metadata.at("payload_sha256") &&
              selected_sha == disk_sha && kept == entries &&
              selector.kept_counts == metadata.at("counts").at("retained")
                                          .get<std::array<uint64_t, 5>>() &&
              selector.guard_particles ==
                  metadata.at("counts").at("guard_only_particles"),
          "Cache raw payload, core/guard or coverage readback mismatch");
  return {{"all_fields_readback", true}, {"payload_sha256", disk_sha},
          {"source_ids_unique_ordered", true}, {"calibration_exact", true},
          {"retained_hits_pass_window", true}, {"entries", entries}};
}
} // namespace sm150_cache

// ==========================================================================
// 7. 只读参数检查与运行路径辅助
// ==========================================================================
static json inspect_parameters(const Parameters &pars,
                               const sm150_windowed::Engine &engine,
                               const json &calibration) {
  json rows = json::object();
  for (auto it = pars.records.begin(); it != pars.records.end(); ++it) {
    int key = std::stoi(it.key());
    rows[it.key()] = {
        {"p", std::vector<double>(pars.p.begin() + 9 * key,
                                  pars.p.begin() + 9 * key + 9)},
        {"q", std::vector<int32_t>(pars.q.begin() + 5 * key,
                                   pars.q.begin() + 5 * key + 5)}};
  }
  json timing = json::array();
  std::vector<int> ge = {204, 205, 206, 209, 211, 212, 213, 214, 215, 304,
                         305, 306, 308, 310, 311, 506, 507, 508, 509, 510};
  for (int key : ge)
    for (int failed : {0, 1})
      for (double energy :
           {39., 40., 50., 100., 334., 439., 737., 1000., 1999.9, 2000.})
        timing.push_back(
            {key, failed, energy,
             engine.timing.TotalGeOffset(key, bool(failed), energy)});
  json particles = json::array();
  for (auto it = calibration.at("particle_tables").at(pars.run).begin();
       it != calibration.at("particle_tables").at(pars.run).end(); ++it)
    for (double adc : {0., 1000., 4000., 10000., 65535.})
      particles.push_back(
          {std::stoi(it.key()), adc,
           engine.timing.TotalParticleOffsetFromAdc(std::stoi(it.key()), adc)});
  return {{"run", pars.run},
          {"group", pars.group},
          {"parameters", rows},
          {"ge_timing", timing},
          {"particle_timing", particles}};
}

static fs::path default_path(const fs::path &requested,
                             const fs::path &fallback) {
  return requested.empty() ? fallback : requested;
}
build::build(const BuildOptions &options) : options_(options) {}

// ==========================================================================
// 7. 一轮处理主流程：参数 → 读取 → 构建 → 两套选窗 → 校验
// ==========================================================================
int build::Process() {
  const auto &options = options_;
  const bool caching = options.mode == "cache", replay = options.mode == "replay";
  auto start = Clock::now();
  const auto cpu_start = std::clock();
  fs::path receipt_path, status_path;
  // 锁保持到异常状态写完，防止失败者覆盖另一任务的进度。
  std::vector<std::unique_ptr<OutputLock>> output_locks;
  try {
    require(options.mode == "windows" || caching || replay, "Unknown build mode");
    require(!replay || (!options.input.empty() && !options.entry_limit),
            "Replay requires --input cache.root and does not accept --entries");
    require(replay || options.cache_receipt.empty(),
            "--cache-receipt is only used in replay mode");
    require(options.entry_limit <=
                uint64_t(std::numeric_limits<int64_t>::max()),
            "Entry limit exceeds ROOT signed entry range");
    fs::path executable = fs::canonical("/proc/self/exe");
    // 第1步：确定输入与输出位置，读取刻度、时间修正和两套逐Ge窗口。
    fs::path default_project =
        executable.parent_path().filename() == "build_all"
            ? executable.parent_path().parent_path()
            : fs::current_path();
    fs::path project =
        fs::canonical(default_path(options.project, default_project));
    int run = options.run;
    std::string rn = run_name(run);
    fs::path window_file = default_path(
        options.window_file,
        project / "codex/analysis/spectrum_ge_equal_height_background_20261022/"
                  "selected_windows.json");
    json input_hashes = json::object();
    json window_json = frozen_json(window_file, input_hashes);
    std::string profile = window_json.at("profile");
    fs::path profile_dir =
        default_path(options.profile_dir,
                     project / "pars/time_alignment/versions" / profile);
    profile_dir = fs::canonical(profile_dir);
    profile = profile_dir.filename().string();
    // 冻结完整配置目录（也覆盖共用加载器可能使用的兼容参数文件）。
    for (const auto &item : fs::directory_iterator(profile_dir))
      if (item.is_regular_file() && item.path().filename() != "README.md")
        input_hashes[item.path().string()] = sha256(item.path());
    json calibration =
        frozen_json(profile_dir / "calibration.json", input_hashes);
    Parameters pars(project.string(), run, calibration);
    input_hashes.update(pars.file_hashes);
    Gates gates(window_json, pars.group);
    if (replay)
      for (const auto &row : gates.rows)
        for (const auto &gate : row.second) {
          // 现行分组以参考Ge为中心取±1000ns，组内Ge最多偏离100ns。
          // 留出该差值；更外侧窗口需另做按目标Ge直接配对的分析。
          require(gate.lower >= -900. && gate.upper <= 900.,
                  "Replay prompt exceeds the current grouping coverage");
          for (const auto &side : gate.random)
            require(side.first >= -900. && side.second <= 900.,
                    "Replay random exceeds the current grouping coverage");
        }
    for (int key = 0; key < 16384; ++key)
      if (pars.q[key * 5] == 1 && pars.q[key * 5 + 4])
        require(gates.rows.count({key, 0}) && gates.rows.count({key, 1}),
                "Missing window for active Ge channel");
    sm150_windowed::Engine engine(pars.p.data(), pars.q.data(),
                                  profile_dir.c_str(), run);
    verify_frozen(input_hashes);
    if (options.inspect_parameters) {
      std::cout << inspect_parameters(pars, engine, calibration).dump() << '\n';
      return 0;
    }
    // 第2步：核对严格轮次名单；已有输出拒绝覆盖。
    fs::path manifest_file = default_path(
        options.source_manifest,
        project /
            "codex/analysis/spectrum_all_runs_timing_20261015/manifest.json");
    json source_manifest = frozen_json(manifest_file, input_hashes);
    for (const auto &p :
         {project / "include/sm150_time_alignment.h",
          project / "build_all/main.cpp", project / "build_all/build.h",
          project / "build_all/build.cpp", project / "build_all/Makefile"})
      input_hashes[p.string()] = sha256(p);
    const std::string binary_hash = sha256(executable);
    json expected;
    for (const auto &row : source_manifest.at("sources"))
      if (row.at("run") == rn) {
        require(expected.is_null(), "Duplicate run in source manifest");
        expected = row;
      }
    require(!expected.is_null(), "Run not in strict173 source list");
    fs::path input = fs::canonical(
        default_path(options.input, "/mnt/e/Sm150/rootfile/data_C2_" + rn +
                                        "_" + rn + ".root"));
    json cache_receipt, cache_metadata;
    fs::path cache_receipt_path;
    if (replay) {
      cache_receipt_path = fs::canonical(default_path(
          options.cache_receipt,
          input.parent_path() / (input.stem().string() + "_complete.json")));
      cache_receipt = frozen_json(cache_receipt_path, input_hashes);
      require(cache_receipt.at("status") == "complete" &&
                  cache_receipt.at("mode") == "cache" &&
                  cache_receipt.at("run") == rn &&
                  cache_receipt.at("verification").at("all_fields_readback") == true,
              "Cache has no matching completed readback receipt");
    }
    fs::path outdir =
        fs::absolute(default_path(options.outdir, caching
            ? project / "rootfile/wide_cache_1000ns" : project / "rootfile"));
    fs::create_directories(outdir);
    std::string tag = options.tag;
    require(tag.find_first_not_of("abcdefghijklmnopqrstuvwxyz0123456789_") ==
                std::string::npos,
            "Invalid output tag");
    if (!tag.empty())
      tag = "_" + tag;
    std::string basename = "data" + rn + "_build" +
                           (caching ? "_cache_1000ns" : "");
    fs::path narrow = outdir / (basename + "_narrow" + tag + ".root"),
             wide = outdir / (basename + "_wide" + tag + ".root");
    fs::path cache_path = outdir / (basename + tag + ".root");
    receipt_path = fs::absolute(default_path(
        options.receipt, outdir / (basename + tag +
            (caching ? "_complete.json" : "_cpp_complete.json"))));
    const fs::path prospective_status = receipt_path.string() + ".status.json";
    const fs::path run_lock = outdir / ("." + basename + tag + ".lock"),
                   receipt_lock =
                       receipt_path.parent_path() /
                       ("." + receipt_path.filename().string() + ".lock");
    std::vector<fs::path> protected_inputs = {input, executable};
    for (auto it = input_hashes.begin(); it != input_hashes.end(); ++it)
      protected_inputs.emplace_back(it.key());
    std::vector<fs::path> targets = caching ? std::vector<fs::path>{cache_path}
                                          : std::vector<fs::path>{narrow, wide};
    std::vector<fs::path> protected_outputs = targets;
    for (const auto &path : targets)
      protected_outputs.emplace_back(path.string() + ".partial");
    for (const auto &path :
         {receipt_path, fs::path(receipt_path.string() + ".tmp"),
          prospective_status, fs::path(prospective_status.string() + ".tmp"),
          run_lock, receipt_lock})
      protected_outputs.push_back(path);
    distinct_paths(protected_outputs, protected_inputs);
    output_locks.emplace_back(std::make_unique<OutputLock>(run_lock));
    output_locks.emplace_back(std::make_unique<OutputLock>(receipt_lock));
    for (const auto &path : targets)
      require(!fs::exists(path) && !fs::exists(path.string() + ".partial"),
              "Existing output or partial; refuse overwrite");
    require(!fs::exists(receipt_path) &&
                !fs::exists(receipt_path.string() + ".tmp"),
            "Existing completion receipt or partial receipt");
    status_path = prospective_status;
    // 第3步：记录进度；整轮运行先核验原始输入SHA。
    auto state = [&](json value) {
      value["run"] = rn;
      value["seconds"] = seconds(start);
      save_json(status_path, value);
    };
    struct stat before, after;
    require(::stat(input.c_str(), &before) == 0, "Cannot stat input");
    require(uint64_t(before.st_size) ==
                (replay ? cache_receipt.at("output").at("bytes") : expected.at("bytes"))
                    .get<uint64_t>(),
            "Wrong input size");
    uint64_t limit = options.entry_limit;
    size_t chunk = options.chunk_entries;
    require(chunk > 0 && chunk <= 2000000, "Block size outside checked range");
    json input_hash = nullptr;
    double last_report = -1., hash_started = seconds(start);
    if (!limit) {
      input_hash = sha256(input, [&](uint64_t done, uint64_t total) {
        double elapsed = seconds(start);
        if (elapsed - last_report >= 10 || done == total) {
          state({{"phase", "source_check"},
                 {"source_bytes_done", done},
                 {"source_bytes_total", total}});
          last_report = elapsed;
        }
      });
      require(input_hash == (replay ? cache_receipt.at("output").at("sha256")
                                   : expected.at("sha256")),
              "Input SHA256 mismatch");
    }
    double hash_seconds = seconds(start) - hash_started;
    verify_frozen(input_hashes);
    // 第4步：打开原始tree，保持原始entry身份和小样覆盖记录。
    std::unique_ptr<TFile> file(TFile::Open(input.c_str(), "READ"));
    require(file && !file->IsZombie(), "Cannot open input ROOT");
    auto *tree = dynamic_cast<TTree *>(file->Get("tree"));
    require(tree, "Missing raw tree");
    int64_t entries = tree->GetEntries();
    if (!replay)
      require(entries == expected.at("entries").get<int64_t>(),
              "Wrong raw tree entry count");
    int64_t stop = limit ? std::min(int64_t(limit), entries) : entries;
    json source = {
        {"path", input.string()},
        {"bytes", before.st_size},
        {"sha256", input_hash},
        {"expected_decoded_sha256", expected.at("sha256")},
        {"entries", entries},
        {"coverage",
         {{"start", 0}, {"stop", stop}, {"fraction", double(stop) / entries}}},
        {"timing_stage", "complete_before_sorting_and_grouping"},
        {"entry_limit", limit ? json(limit) : json(nullptr)}};
    bool full_run = stop == entries && !limit;
    json replay_input = nullptr;
    if (replay) {
      auto *object = dynamic_cast<TObjString *>(file->Get("sm150_cache_manifest"));
      require(object, "Missing cache manifest");
      cache_metadata = json::parse(object->GetString().Data());
      require(cache_metadata == cache_receipt.at("cache_manifest") &&
                  cache_metadata.at("schema_version") == 1 &&
                  cache_metadata.at("run") == rn &&
                  cache_metadata.at("group") == pars.group &&
                  cache_metadata.at("retained_entries") == entries &&
                  cache_metadata.at("all_Ge_retained") == true &&
                  cache_metadata.at("core_half_width_ns") == sm150_cache::core_ns &&
                  cache_metadata.at("guard_ns") == sm150_cache::guard_ns &&
                  cache_metadata.at("energy_parameters") == pars.records,
              "Cache identity/coverage/energy calibration mismatch");
      for (const auto &item : fs::directory_iterator(profile_dir))
        if (item.is_regular_file() && item.path().filename() != "README.md")
          require(cache_metadata.at("time_parameter_sha256")
                      .at(item.path().filename().string()) ==
                      input_hashes.at(item.path().string()),
                  "Time parameters changed since cache selection; re-audit coverage");
      source = cache_metadata.at("source");
      require(source.at("entries") == expected.at("entries") &&
                  source.at("expected_decoded_sha256") == expected.at("sha256") &&
                  source.at("bytes") == expected.at("bytes") &&
                  source.at("coverage").at("start") == 0 &&
                  source.at("coverage").at("stop").get<int64_t>() <=
                      expected.at("entries").get<int64_t>(),
              "Cached original-source provenance mismatch");
      full_run = cache_receipt.at("full_run");
      require(!full_run ||
                  (source.at("sha256") == expected.at("sha256") &&
                   source.at("coverage").at("stop") == expected.at("entries")),
              "Full cache lacks verified complete original-source coverage");
      replay_input = {{"path", input.string()}, {"bytes", before.st_size},
                      {"sha256", input_hash}, {"entries", entries},
                      {"receipt", cache_receipt_path.string()}};
    }
    if (caching) {
      json time_hashes = json::object();
      for (const auto &item : fs::directory_iterator(profile_dir))
        if (item.is_regular_file() && item.path().filename() != "README.md")
          time_hashes[item.path().filename().string()] =
              input_hashes.at(item.path().string());
      json metadata = {
          {"schema_version", 1}, {"run", rn}, {"group", pars.group},
          {"source", source}, {"core_half_width_ns", sm150_cache::core_ns},
          {"guard_ns", sm150_cache::guard_ns}, {"all_Ge_retained", true},
          {"selection", "All Ge, plus every Si/S3 hit within inclusive 1100ns "
                        "of any corrected Ge. No energy/quality/multiplicity cut."},
          {"in_core", "true for Ge; particle within inclusive 1000ns of any Ge"},
          {"ordering", "original source_entry, unique raw hits"},
          {"energy_parameters", pars.records},
          {"time_alignment_profile", profile},
          {"time_parameter_sha256", time_hashes},
          {"frozen_input_sha256", input_hashes},
          {"executable_sha256", binary_hash},
          {"compression", "ZSTD level4"}, {"ROOT_version", gROOT->GetVersion()}};
      sm150_cache::Writer writer(cache_path);
      sm150_cache::Selector selector;
      RawReader reader(*tree, stop);
      RawBlock raw(chunk);
      auto build_start = Clock::now();
      double raw_seconds = 0.;
      for (int64_t begin = 0; begin < stop;) {
        const size_t n = std::min<int64_t>(chunk, stop - begin);
        auto read_start = Clock::now();
        reader.read(begin, n, raw);
        raw_seconds += seconds(read_start);
        for (size_t i = 0; i < n; ++i)
          selector.add(sm150_cache::calibrate(raw, i, begin + i, pars, engine.timing));
        begin += n;
        selector.flush(begin == stop,
                       [&](const sm150_cache::Record &r) { writer.write(r); });
        state({{"phase", "caching"}, {"entries_done", begin},
               {"entries_total", stop}, {"retained_entries", writer.count}});
        if (begin % 10000000 == 0 || begin == stop)
          std::cout << "CACHE " << rn << ": " << begin << '/' << stop
                    << "; kept " << writer.count << "; " << seconds(start)
                    << "s\n" << std::flush;
        require(fs::space(outdir).available > 2ULL * 1024 * 1024 * 1024,
                "Output free-space floor");
        require(seconds(start) < 86400, "24hour per-run budget");
      }
      require(selector.pending.empty() &&
                  selector.input_counts[1] == selector.kept_counts[1],
              "Incomplete cache or missing Ge");
      metadata = writer.close(metadata, selector);
      const double build_seconds = seconds(build_start);
      state({{"phase", "validating_cache"}});
      auto verify_start = Clock::now();
      json verified = sm150_cache::verify(writer.partial, metadata, pars,
                                          engine.timing, chunk);
      json output = {{"path", cache_path.string()},
                     {"bytes", fs::file_size(writer.partial)},
                     {"sha256", sha256(writer.partial)}};
      require(::stat(input.c_str(), &after) == 0 &&
                  before.st_size == after.st_size && before.st_ino == after.st_ino &&
                  before.st_mtim.tv_sec == after.st_mtim.tv_sec &&
                  before.st_mtim.tv_nsec == after.st_mtim.tv_nsec,
              "Input identity changed during caching");
      verify_frozen(input_hashes);
      publish_file(writer.partial, writer.path);
      json result = {
          {"status", "complete"}, {"mode", "cache"}, {"run", rn},
          {"group", pars.group}, {"full_run", full_run},
          {"source", source}, {"cache_manifest", metadata},
          {"output", output}, {"verification", verified},
          {"timing_seconds",
           {{"source_hash", hash_seconds}, {"read_build_write", build_seconds},
            {"raw_root_read", raw_seconds},
            {"readback_and_final_checks", seconds(verify_start)},
            {"total", seconds(start)},
            {"cpu", double(std::clock() - cpu_start) / CLOCKS_PER_SEC}}},
          {"frozen_input_sha256", input_hashes},
          {"executable_sha256", binary_hash}, {"chunk_entries", chunk}};
      save_json(receipt_path, result, true);
      state({{"phase", "complete"}, {"receipt", receipt_path.string()}});
      std::cout << "CACHE COMPLETE " << rn << ": " << output.at("bytes")
                << " bytes; receipt " << receipt_path << '\n';
      return 0;
    }
    // 每份ROOT内嵌实际窗口、时间版本、筛选条件和输入来源。
    auto manifest = [&](int label) {
      return json{
          {"schema_version", 2},
          {"run", rn},
          {"intensity", pars.group == "low" ? "low" : "high"},
          {"group", pars.group},
          {"candidate", label ? "B" : "A"},
          {"window_file", window_file.string()},
          {"window_sha256", input_hashes.at(window_file.string())},
          {"time_alignment_profile", profile},
          {"time_alignment_profile_dir", profile_dir.string()},
          {"time_alignment_calibration_sha256",
           input_hashes.at((profile_dir / "calibration.json").string())},
          {"energy_parameters", pars.records},
          {"source", source},
          {"replay_input", replay_input},
          {"windows", gates.manifest(label)},
          {"selection",
           {{"Ge_lab_keV", {40, 2000}},
            {"Si_MeV", {15, 100}},
            {"S3_MeV", {25, 100}},
            {"particle_lower_energy_comparison", "strict greater than"},
            {"single_Ge", true},
            {"no_addback", true},
            {"S3_ring_range", {1, pars.group == "low" ? 24 : 8}}}},
          {"time_definition", "dt=t_particle_corrected-t_Ge_corrected; "
                              "lower<=dt<upper, exact double dt"},
          {"random", "Both sides unweighted; per-partneralpha=Wp/WR,extra1"},
          {"histograms",
           "No random subtraction or Doppler correction during build"},
          {"pair_deduplication",
           {{"key",
             {"run", "candidate", "Ge_source_entry", "Si_or_S3",
              "particle_source_entry"}},
            {"keep", "first accepted occurrence"},
            {"S3_sector",
             "kept with ring; repeated ring requires identical sector"},
            {"readback", "unique across both prompt and random trees"}}},
          {"implementation", "original_layout_cpp"},
          {"raw_reader", "root_bulk_baskets"},
          {"executable_sha256", binary_hash},
          {"frozen_input_sha256", input_hashes},
          {"compression", "ZSTD level4"},
          {"ROOT_version", gROOT->GetVersion()}};
    };
    // 第5步：一次读取同时构建窄窗A和宽窗B；prompt/random保留原始计数。
    FileWriter a(narrow, pars, gates, 0, manifest(0)),
        b(wide, pars, gates, 1, manifest(1));
    RawReader reader(*tree, stop, replay);
    RawBlock raw(chunk);
    int64_t selected = 0;
    double raw_read_seconds = 0, engine_seconds = 0, classify_seconds = 0;
    auto build_start = Clock::now();
    for (int64_t begin = 0; begin < stop;) {
      size_t n = std::min<int64_t>(chunk, stop - begin);
      auto part = Clock::now();
      reader.read(begin, n, raw);
      raw_read_seconds += seconds(part);
      part = Clock::now();
      engine.process(n, begin, raw.cid.data(), raw.sid.data(), raw.ch.data(),
                     raw.sr.data(), raw.ts.data(), raw.evte.data(),
                     raw.cfd.data(), raw.cfds.data(), raw.cfdft.data(),
                     raw.pileup.data(), raw.outofr.data(),
                     begin + int64_t(n) == stop,
                     replay ? raw.source_entry.data() : nullptr);
      engine_seconds += seconds(part);
      part = Clock::now();
      selected += classify(engine, pars, gates, a, b);
      classify_seconds += seconds(part);
      begin += n;
      state({{"phase", "building"},
             {"entries_done", begin},
             {"entries_total", stop},
             {"selected_Ge_events", selected}});
      if (begin % 10000000 == 0 || begin == stop)
        std::cout << "BUILD " << rn << ": " << begin << '/' << stop
                  << " records; " << seconds(start) << "s\n"
                  << std::flush;
      require(fs::space(outdir).available > 2ULL * 1024 * 1024 * 1024,
              "Output free-space floor");
      require(seconds(start) < 86400, "24hour per-run budget");
    }
    // 第6步：写完最后一块，确认事件覆盖完整，再关闭两套输出。
    require(engine.read == stop && engine.tail.empty(),
            "Incomplete streaming coverage or uncleared tail");
    auto ac = a.close(), bc = b.close();
    require(bc["tr_event"]["si_pairs"].get<int64_t>() >=
                    ac["tr_event"]["si_pairs"].get<int64_t>() &&
                bc["tr_event"]["s3_pairs"].get<int64_t>() >=
                    ac["tr_event"]["s3_pairs"].get<int64_t>(),
            "Nested window counts failed");
    // 第7步：完整读回两套ROOT，核验门、alpha、计数及P/R合计后的原始配对唯一性。
    double build_seconds = seconds(build_start);
    state({{"phase", "validating"}});
    auto verify_start = Clock::now();
    json outputs = json::array();
    for (int label = 0; label < 2; ++label) {
      const auto &writer = label ? b : a;
      const auto &path = writer.partial;
      auto checked = verify_root(path, pars, gates, label, label ? bc : ac);
      outputs.push_back(
          {{"path", writer.path.string()},
           {"bytes", fs::file_size(path)},
           {"sha256", sha256(path)},
           {"trees", checked},
           {"pair_deduplication", (label ? b : a).unique_pairs.counts()},
           {"unique_original_pairs_verified", true},
           {"all_cuts_and_alpha_verified", true}});
    }
    require(::stat(input.c_str(), &after) == 0 &&
                before.st_size == after.st_size &&
                before.st_ino == after.st_ino &&
                before.st_mtim.tv_sec == after.st_mtim.tv_sec &&
                before.st_mtim.tv_nsec == after.st_mtim.tv_nsec,
            "Input identity changed during build");
    verify_frozen(input_hashes);
    // 第8步：两套读回及冻结输入检查全部通过后才发布；完成回执是整轮完成标志。
    publish_file(a.partial, a.path);
    publish_file(b.partial, b.path);
    json result = {
        {"status", "complete"},
        {"run", rn},
        {"group", pars.group},
        {"implementation", "original_layout_cpp"},
        {"raw_reader", "root_bulk_baskets"},
        {"mode", options.mode},
        {"full_run", full_run},
        {"source", source},
        {"replay_input", replay_input},
        {"outputs", outputs},
        {"extractor_counters",
         {{"raw_entries", engine.read},
          {"ge_reference_groups", engine.all_ge_groups},
          {"tail_hits", engine.tail.size()}}},
        {"timing_seconds",
         {{"source_hash", hash_seconds},
          {"read_build_write", build_seconds},
          {"raw_root_read", raw_read_seconds},
          {"calibrate_sort_group", engine_seconds},
          {"classify_write", classify_seconds},
          {"readback_and_final_checks", seconds(verify_start)},
          {"total", seconds(start)},
          {"cpu", double(std::clock() - cpu_start) / CLOCKS_PER_SEC}}},
        {"frozen_input_sha256", input_hashes},
        {"executable_sha256", binary_hash},
        {"ROOT_version", gROOT->GetVersion()},
        {"chunk_entries", chunk}};
    save_json(receipt_path, result, true);
    state({{"phase", "complete"}, {"receipt", receipt_path.string()}});
    std::cout << result.dump() << '\n';
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "ERROR: " << error.what() << '\n';
    if (!status_path.empty())
      try {
        save_json(status_path, {{"phase", "failed"},
                                {"error", error.what()},
                                {"seconds", seconds(start)}});
      } catch (...) {
      }
    return 1;
  }
}

// ==========================================================================
// 8. Ge、Si array、S3的ring/sector通道映射（保持原映射）
// ==========================================================================
namespace sm150_native {
bool Maps::init(int run) {

  //
  map_ge_sector_id = {
      {204, 8}, {205, 9}, {206, 10}, {208, 2}, {209, 3}, {210, 7}, {211, 1},
      {212, 2}, {213, 2}, {214, 2},  {215, 2},

      {304, 2}, {305, 3}, {306, 5},  {307, 3}, {308, 8}, {309, 6}, {310, 4},
      {311, 6},

      {400, 3}, {401, 3}, {402, 3},  {403, 3}, {404, 4}, {405, 4}, {406, 4},
      {407, 4}, {409, 5}, {410, 5},  {411, 5},

      {504, 1}, {505, 3}, {506, 2},  {507, 4}, {508, 6}, {509, 5}, {510, 7}};

  //
  map_ge_ring_id = {{204, 3}, {205, 3}, {206, 3}, {208, 2}, {209, 5}, {210, 2},
                    {211, 4}, {212, 3}, {213, 3}, {214, 3}, {215, 3},

                    {304, 4}, {305, 4}, {306, 4}, {307, 1}, {308, 4}, {309, 1},
                    {310, 5}, {311, 5},

                    {400, 3}, {401, 3}, {402, 3}, {403, 3}, {404, 3}, {405, 3},
                    {406, 3}, {407, 3}, {409, 3}, {410, 3}, {411, 3},

                    {504, 2}, {505, 2}, {506, 5}, {507, 4}, {508, 4}, {509, 5},
                    {510, 4}};

  //
  if ((run >= 472 && run <= 605) || (run >= 680 && run <= 718)) {
    //
    map_spider_sector_id = {{10200, 2},  {10201, 2},  {10202, 2},  {10203, 2},
                            {10204, 2},  {10205, 2},  {10206, 2},  {10207, 2},
                            {10208, 1},  {10209, 1},  {10210, 1},  {10211, 1},
                            {10212, 1},  {10213, 1},  {10214, 1},  {10215, 1},

                            {10300, 4},  {10301, 4},  {10302, 4},  {10303, 4},
                            {10304, 4},  {10305, 4},  {10306, 4},  {10307, 4},
                            {10308, 3},  {10309, 3},  {10310, 3},  {10311, 3},
                            {10312, 3},  {10313, 3},  {10314, 3},  {10315, 3},

                            {10400, 6},  {10401, 6},  {10402, 6},  {10403, 6},
                            {10404, 6},  {10405, 6},  {10406, 6},  {10407, 6},
                            {10408, 5},  {10409, 5},  {10410, 5},  {10411, 5},
                            {10412, 5},  {10413, 5},  {10414, 5},  {10415, 5},

                            {10600, 7},  {10601, 7},  {10602, 7},  {10603, 7},
                            {10604, 7},  {10605, 7},  {10606, 7},  {10607, 7},
                            {10608, 8},  {10609, 8},  {10610, 8},  {10611, 8},
                            {10612, 8},  {10613, 8},  {10614, 8},  {10615, 8},

                            {10700, 9},  {10701, 9},  {10702, 9},  {10703, 9},
                            {10704, 9},  {10705, 9},  {10706, 9},  {10707, 9},
                            {10708, 10}, {10709, 10}, {10710, 10}, {10711, 10},
                            {10712, 10}, {10713, 10}, {10714, 10}, {10715, 10},

                            {10800, 11}, {10801, 11}, {10802, 11}, {10803, 11},
                            {10804, 11}, {10805, 11}, {10806, 11}, {10807, 11},
                            {10808, 12}, {10809, 12}, {10810, 12}, {10811, 12},
                            {10812, 12}, {10813, 12}, {10814, 12}, {10815, 12}};

    //
    map_spider_ring_id = {
        {10200, 8}, {10201, 7}, {10202, 6}, {10203, 5}, {10204, 4}, {10205, 3},
        {10206, 2}, {10207, 1}, {10208, 8}, {10209, 7}, {10210, 6}, {10211, 5},
        {10212, 4}, {10213, 3}, {10214, 2}, {10215, 1},

        {10300, 1}, {10301, 2}, {10302, 3}, {10303, 4}, {10304, 5}, {10305, 6},
        {10306, 7}, {10307, 8}, {10308, 1}, {10309, 2}, {10310, 3}, {10311, 4},
        {10312, 5}, {10313, 6}, {10314, 7}, {10315, 8},

        {10400, 8}, {10401, 7}, {10402, 6}, {10403, 5}, {10404, 4}, {10405, 3},
        {10406, 2}, {10407, 1}, {10408, 8}, {10409, 7}, {10410, 6}, {10411, 5},
        {10412, 4}, {10413, 3}, {10414, 2}, {10415, 1},

        {10600, 1}, {10601, 2}, {10602, 3}, {10603, 4}, {10604, 5}, {10605, 6},
        {10606, 7}, {10607, 8}, {10608, 1}, {10609, 2}, {10610, 3}, {10611, 4},
        {10612, 5}, {10613, 6}, {10614, 7}, {10615, 8},

        {10700, 1}, {10701, 2}, {10702, 3}, {10703, 4}, {10704, 5}, {10705, 6},
        {10706, 7}, {10707, 8}, {10708, 1}, {10709, 2}, {10710, 3}, {10711, 4},
        {10712, 5}, {10713, 6}, {10714, 7}, {10715, 8},

        {10800, 1}, {10801, 2}, {10802, 3}, {10803, 4}, {10804, 5}, {10805, 6},
        {10806, 7}, {10807, 8}, {10808, 1}, {10809, 2}, {10810, 3}, {10811, 4},
        {10812, 5}, {10813, 6}, {10814, 7}, {10815, 8}};

    //
    map_s3_sector_id = {{10900, 10}, {10901, 12}, {10902, 14}, {10903, 16},
                        {10904, 18}, {10905, 20}, {10906, 22}, {10907, 24},
                        {10908, 26}, {10909, 28}, {10910, 30}, {10911, 32},
                        {10912, 2},  {10913, 4},  {10914, 6},  {10915, 8},

                        {11000, 11}, {11001, 13}, {11002, 15}, {11003, 17},
                        {11004, 19}, {11005, 21}, {11006, 23}, {11007, 25},
                        {11008, 27}, {11009, 29}, {11010, 31}, {11011, 1},
                        {11012, 3},  {11013, 5},  {11014, 7},  {11015, 9}};
  } else if ((run >= 616 && run <= 674) || (run >= 721 && run <= 723)) {
    //
    map_spider_sector_id = {{10900, 2},  {10901, 2},  {10902, 2},  {10903, 2},
                            {10904, 2},  {10905, 2},  {10906, 2},  {10907, 2},
                            {10908, 1},  {10909, 1},  {10910, 1},  {10911, 1},
                            {10912, 1},  {10913, 1},  {10914, 1},  {10915, 1},

                            {11000, 4},  {11001, 4},  {11002, 4},  {11003, 4},
                            {11004, 4},  {11005, 4},  {11006, 4},  {11007, 4},
                            {11008, 3},  {11009, 3},  {11010, 3},  {11011, 3},
                            {11012, 3},  {11013, 3},  {11014, 3},  {11015, 3},

                            {10400, 6},  {10401, 6},  {10402, 6},  {10403, 6},
                            {10404, 6},  {10405, 6},  {10406, 6},  {10407, 6},
                            {10408, 5},  {10409, 5},  {10410, 5},  {10411, 5},
                            {10412, 5},  {10413, 5},  {10414, 5},  {10415, 5},

                            {10600, 7},  {10601, 7},  {10602, 7},  {10603, 7},
                            {10604, 7},  {10605, 7},  {10606, 7},  {10607, 7},
                            {10608, 8},  {10609, 8},  {10610, 8},  {10611, 8},
                            {10612, 8},  {10613, 8},  {10614, 8},  {10615, 8},

                            {10700, 9},  {10701, 9},  {10702, 9},  {10703, 9},
                            {10704, 9},  {10705, 9},  {10706, 9},  {10707, 9},
                            {10708, 10}, {10709, 10}, {10710, 10}, {10711, 10},
                            {10712, 10}, {10713, 10}, {10714, 10}, {10715, 10},

                            {10800, 11}, {10801, 11}, {10802, 11}, {10803, 11},
                            {10804, 11}, {10805, 11}, {10806, 11}, {10807, 11},
                            {10808, 12}, {10809, 12}, {10810, 12}, {10811, 12},
                            {10812, 12}, {10813, 12}, {10814, 12}, {10815, 12}};

    //
    map_spider_ring_id = {
        {10900, 8}, {10901, 7}, {10902, 6}, {10903, 5}, {10904, 4}, {10905, 3},
        {10906, 2}, {10907, 1}, {10908, 8}, {10909, 7}, {10910, 6}, {10911, 5},
        {10912, 4}, {10913, 3}, {10914, 2}, {10915, 1},

        {11000, 1}, {11001, 2}, {11002, 3}, {11003, 4}, {11004, 5}, {11005, 6},
        {11006, 7}, {11007, 8}, {11008, 1}, {11009, 2}, {11010, 3}, {11011, 4},
        {11012, 5}, {11013, 6}, {11014, 7}, {11015, 8},

        {10400, 8}, {10401, 7}, {10402, 6}, {10403, 5}, {10404, 4}, {10405, 3},
        {10406, 2}, {10407, 1}, {10408, 8}, {10409, 7}, {10410, 6}, {10411, 5},
        {10412, 4}, {10413, 3}, {10414, 2}, {10415, 1},

        {10600, 1}, {10601, 2}, {10602, 3}, {10603, 4}, {10604, 5}, {10605, 6},
        {10606, 7}, {10607, 8}, {10608, 1}, {10609, 2}, {10610, 3}, {10611, 4},
        {10612, 5}, {10613, 6}, {10614, 7}, {10615, 8},

        {10700, 1}, {10701, 2}, {10702, 3}, {10703, 4}, {10704, 5}, {10705, 6},
        {10706, 7}, {10707, 8}, {10708, 1}, {10709, 2}, {10710, 3}, {10711, 4},
        {10712, 5}, {10713, 6}, {10714, 7}, {10715, 8},

        {10800, 1}, {10801, 2}, {10802, 3}, {10803, 4}, {10804, 5}, {10805, 6},
        {10806, 7}, {10807, 8}, {10808, 1}, {10809, 2}, {10810, 3}, {10811, 4},
        {10812, 5}, {10813, 6}, {10814, 7}, {10815, 8}};

    //
    map_s3_sector_id = {{10200, 10}, {10201, 12}, {10202, 14}, {10203, 16},
                        {10204, 18}, {10205, 20}, {10206, 22}, {10207, 24},
                        {10208, 26}, {10209, 28}, {10210, 30}, {10211, 32},
                        {10212, 2},  {10213, 4},  {10214, 6},  {10215, 8},

                        {10300, 11}, {10301, 13}, {10302, 15}, {10303, 17},
                        {10304, 19}, {10305, 21}, {10306, 23}, {10307, 25},
                        {10308, 27}, {10309, 29}, {10310, 31}, {10311, 1},
                        {10312, 3},  {10313, 5},  {10314, 7},  {10315, 9}};
  } else {

    return 0;
  }

  //
  map_s3_ring_id = {{11102, 1},  {11103, 3},  {11104, 5},  {11105, 7},
                    {11106, 9},  {11107, 11}, {11108, 13}, {11109, 15},
                    {11110, 17}, {11111, 19}, {11112, 21}, {11113, 23},

                    {11200, 2},  {11201, 4},  {11202, 6},  {11203, 8},
                    {11204, 10}, {11205, 12}, {11206, 14}, {11207, 16},
                    {11208, 18}, {11209, 20}, {11210, 22}, {11211, 24}};

  return 1;
}
} // namespace sm150_native
