#pragma once

#include <filesystem>
#include <functional>
#include <iosfwd>

namespace nova::io::detail {

// Stage on the destination filesystem. A writer or verifier exception leaves
// the previous destination intact. Writers targeting the same name must each
// use a separate staging directory.
void write_atomic_file(const std::filesystem::path& destination,
                       const std::function<void(std::ostream&)>& write,
                       const std::function<void(const std::filesystem::path&)>& verify = {});
} // namespace nova::io::detail
