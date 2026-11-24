#include "atomic_file.hpp"

#include <filesystem>
#include <fstream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <utility>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif

namespace nova::io::detail {
namespace {
namespace fs = std::filesystem;

void reject_unsafe_destination(const fs::path& path) {
    const auto status = fs::symlink_status(path);
    if (fs::is_symlink(status) ||
        (fs::exists(status) && !fs::is_regular_file(status))) {
        throw std::invalid_argument(
            "NOVA file destination must not be a symlink or non-file.");
    }
}

fs::path reserve_staging_directory(const fs::path& parent) {
    std::random_device randomness;
    for (unsigned attempt = 0; attempt < 64U; ++attempt) {
        std::ostringstream name;
        name << ".nova-staging-" << std::hex
             << randomness() << '-' << randomness() << '-' << attempt;
        const auto candidate = parent / name.str();
        if (fs::create_directory(candidate)) {
            return candidate;
        }
    }
    throw std::runtime_error("Cannot reserve a unique NOVA staging directory.");
}

class StagingDirectory {
public:
    explicit StagingDirectory(fs::path directory)
        : directory_(std::move(directory)) {}

    StagingDirectory(const StagingDirectory&) = delete;
    StagingDirectory& operator=(const StagingDirectory&) = delete;

    ~StagingDirectory() {
        std::error_code ignored;
        fs::remove_all(directory_, ignored);
    }

    [[nodiscard]] fs::path record() const {
        return directory_ / "record";
    }

private:
    fs::path directory_;
};

void replace_file(const fs::path& staged, const fs::path& destination) {
    reject_unsafe_destination(destination);
#ifdef _WIN32
    // ReplaceFileW replaces an existing file without an application-level
    // remove gap. MoveFileExW handles first-time creation. Both files share
    // a parent filesystem; neither call promises power-failure durability.
    const bool existing = fs::exists(destination);
    const BOOL success = existing
        ? ::ReplaceFileW(destination.c_str(), staged.c_str(), nullptr,
                         0, nullptr, nullptr)
        : ::MoveFileExW(staged.c_str(), destination.c_str(),
                        MOVEFILE_WRITE_THROUGH);
    if (!success) {
        throw std::system_error(
            static_cast<int>(::GetLastError()), std::system_category(),
            "Failed to publish staged NOVA file");
    }
#else
    fs::rename(staged, destination);
#endif
}
} // namespace

void write_atomic_file(const fs::path& destination,
                       const std::function<void(std::ostream&)>& write,
                       const std::function<void(const fs::path&)>& verify) {
    if (destination.empty() || destination.filename().empty()) {
        throw std::invalid_argument("NOVA file destination requires a filename.");
    }
    reject_unsafe_destination(destination);
    const auto parent = destination.has_parent_path() ? destination.parent_path() : fs::path(".");
    StagingDirectory temporary(reserve_staging_directory(parent));
    const auto record = temporary.record();
    std::ofstream file(record, std::ios::binary | std::ios::trunc);
    if (!file) throw std::runtime_error("Cannot open NOVA staging file.");
    write(file);
    file.flush();
    if (!file) throw std::runtime_error("Cannot flush NOVA staging file.");
    file.close();
    if (!file) throw std::runtime_error("Cannot close NOVA staging file.");
    if (verify) verify(record);
    replace_file(record, destination);
}

} // namespace nova::io::detail
