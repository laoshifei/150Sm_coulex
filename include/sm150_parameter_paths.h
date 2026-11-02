#ifndef SM150_PARAMETER_PATHS_H
#define SM150_PARAMETER_PATHS_H

// Resolve the adopted runtime files without embedding calibration values in C++.
// Explicit candidate/version arguments remain the responsibility of each caller.
#include <nlohmann/json.hpp>
#include <openssl/evp.h>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace sm150 {
class ParameterPaths {
 public:
  explicit ParameterPaths(const std::filesystem::path& pars)
      : root_(std::filesystem::absolute(pars).lexically_normal()) {
    const auto file = root_ / "active.json";
    if (!std::filesystem::exists(file)) return;
    active_hash_ = sha256(file);
    active_ = read_json(file);
    check(active_.at("schema") == "sm150_active_parameters_v1", "Unsupported runtime parameter manifest");
    catalog_file_ = root_ / relative(active_.at("catalog").get<std::string>());
    catalog_hash_ = sha256(catalog_file_);
    check(catalog_hash_ == active_.at("catalog_sha256").get<std::string>(), "Parameter catalog SHA256 mismatch");
    catalog_ = read_json(catalog_file_);
    check(catalog_.at("schema") == "sm150_parameter_catalog_v1", "Unsupported parameter catalog");
    check(active_.at("paths").is_object() && catalog_.at("files").is_object() &&
          catalog_.at("compatibility_paths").is_object(), "Invalid parameter catalog structure");
    verify_unchanged();
  }

  bool active() const { return !active_.is_null(); }
  std::filesystem::path manifest_path() const { return root_ / "active.json"; }
  std::vector<std::filesystem::path> selection_files() const {
    return active() ? std::vector<std::filesystem::path>{manifest_path(), catalog_file_}
                    : std::vector<std::filesystem::path>{};
  }
  nlohmann::json identity() const {
    if (!active()) return {{"selection", "legacy_paths"}};
    return {{"path", manifest_path().string()}, {"sha256", active_hash_},
            {"version", active_.at("version")}, {"catalog_sha256", catalog_hash_}};
  }

  std::filesystem::path resource(const std::string& name,
                                const std::filesystem::path& fallback = {}) const {
    if (!active()) {
      check(!fallback.empty(), "No active parameter manifest or fallback for " + name);
      return fallback;
    }
    const auto path = root_ / relative(active_.at("paths").at(name).at("path").get<std::string>());
    check(std::filesystem::exists(path), "Missing parameter resource: " + path.string());
    if (std::filesystem::is_directory(path)) {
      const auto prefix = path.lexically_relative(root_).generic_string() + "/";
      bool found = false;
      for (auto it = catalog_.at("files").begin(); it != catalog_.at("files").end(); ++it)
        if (it.key().rfind(prefix, 0) == 0) {
          verify_file(root_ / relative(it.key()), true);
          found = true;
        }
      check(found, "Parameter directory has no catalog entries: " + path.string());
    } else {
      verify_file(path, true);
    }
    return path;
  }

  // Map old pars-relative names to the same frozen numerical files. Unknown
  // names retain historical lookup, so explicit old profiles remain usable.
  std::filesystem::path legacy(const std::filesystem::path& name) const {
    auto rel = relative(name.generic_string());
    if (active()) {
      const auto input = rel.generic_string();
      std::size_t length = 0;
      for (auto it = catalog_.at("compatibility_paths").begin(); it != catalog_.at("compatibility_paths").end(); ++it) {
        const auto& old = it.key();
        if (old.size() > length && (input == old || input.rfind(old + "/", 0) == 0)) {
          rel = relative(it.value().get<std::string>() + input.substr(old.size()));
          length = old.size();
        }
      }
    }
    const auto path = root_ / rel;
    if (active() && std::filesystem::is_regular_file(path)) verify_file(path, false);
    return path;
  }

  std::filesystem::path resolve(const std::filesystem::path& path) const {
    const auto absolute = std::filesystem::absolute(path).lexically_normal();
    const auto rel = absolute.lexically_relative(root_);
    if (rel.empty() || *rel.begin() == "..") return path;
    return legacy(rel);
  }

  std::filesystem::path timing_directory() const {
    if (active()) return resource("timing");
    std::ifstream stream(root_ / "time_alignment/active.txt");
    if (!stream) return {};
    std::string profile, extra;
    check(bool(stream >> profile) && !(stream >> extra), "Invalid legacy active timing selector");
    if (profile == "legacy") return {};
    check(profile.find_first_not_of("abcdefghijklmnopqrstuvwxyz0123456789_") == std::string::npos,
          "Invalid timing profile name");
    return root_ / "time_alignment/versions" / profile;
  }

  void verify_unchanged() const {
    if (!active()) return;
    check(sha256(manifest_path()) == active_hash_ && sha256(catalog_file_) == catalog_hash_,
          "Runtime parameter selection changed while loading");
  }

  static std::string sha256(const std::filesystem::path& path) {
    std::ifstream input(path, std::ios::binary);
    check(bool(input), "Cannot read parameter file: " + path.string());
    std::unique_ptr<EVP_MD_CTX, decltype(&EVP_MD_CTX_free)> context(EVP_MD_CTX_new(), EVP_MD_CTX_free);
    check(context && EVP_DigestInit_ex(context.get(), EVP_sha256(), nullptr) == 1, "Parameter SHA256 initialization failed");
    std::vector<char> buffer(1024 * 1024);
    while (input) {
      input.read(buffer.data(), buffer.size());
      check(EVP_DigestUpdate(context.get(), buffer.data(), input.gcount()) == 1, "Parameter SHA256 update failed");
    }
    check(input.eof(), "Parameter read failed: " + path.string());
    unsigned char bytes[EVP_MAX_MD_SIZE];
    unsigned size = 0;
    check(EVP_DigestFinal_ex(context.get(), bytes, &size) == 1, "Parameter SHA256 finish failed");
    std::ostringstream out;
    for (unsigned i = 0; i < size; ++i) out << std::hex << std::setw(2) << std::setfill('0') << unsigned(bytes[i]);
    return out.str();
  }

 private:
  std::filesystem::path root_, catalog_file_;
  nlohmann::json active_, catalog_;
  std::string active_hash_, catalog_hash_;
  static void check(bool value, const std::string& message) {
    if (!value) throw std::runtime_error(message);
  }
  static std::filesystem::path relative(const std::string& value) {
    const std::filesystem::path path(value);
    check(!path.empty() && !path.is_absolute(), "Parameter path must be relative: " + value);
    for (const auto& part : path) check(part != "..", "Parameter path escapes pars: " + value);
    return path.lexically_normal();
  }
  static nlohmann::json read_json(const std::filesystem::path& path) {
    std::ifstream stream(path);
    check(bool(stream), "Cannot read parameter JSON: " + path.string());
    nlohmann::json value; stream >> value; return value;
  }
  void verify_file(const std::filesystem::path& path, bool required) const {
    const auto name = path.lexically_relative(root_).generic_string();
    const auto& files = catalog_.at("files");
    if (!files.contains(name)) {
      check(!required, "Parameter file missing from catalog: " + name);
      return;
    }
    check(sha256(path) == files.at(name).at("sha256").get<std::string>(), "Parameter SHA256 mismatch: " + name);
  }
};
}  // namespace sm150

#endif
