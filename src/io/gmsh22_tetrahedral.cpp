#include "nova/io/gmsh22_tetrahedral.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <istream>
#include <limits>
#include <locale>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace nova::io {
namespace {

using FaceKey = std::array<std::size_t, 3>;

[[nodiscard]] std::string line_or_throw(std::istream& input) {
    std::string line;
    if (!std::getline(input, line)) {
        throw std::invalid_argument("Unexpected end of Gmsh MSH 2.2 file.");
    }
    if (!line.empty() && line.back() == '\r') {
        line.pop_back();
    }
    return line;
}

void expect_line(std::istream& input, const char* expected) {
    if (line_or_throw(input) != expected) {
        throw std::invalid_argument(
            std::string("Expected Gmsh section marker: ") + expected
        );
    }
}

[[nodiscard]] std::vector<std::string> words(const std::string& line) {
    std::istringstream stream(line);
    stream.imbue(std::locale::classic());
    std::vector<std::string> result;
    std::string word;
    while (stream >> word) {
        result.push_back(std::move(word));
    }
    return result;
}

[[nodiscard]] std::uint64_t positive_integer(const std::string& token) {
    std::uint64_t result{};
    const auto converted = std::from_chars(
        token.data(), token.data() + token.size(), result
    );
    if (converted.ec != std::errc{} ||
        converted.ptr != token.data() + token.size() || result == 0) {
        throw std::invalid_argument("Invalid positive Gmsh integer.");
    }
    return result;
}

[[nodiscard]] std::size_t count_integer(const std::string& token) {
    std::uint64_t result{};
    const auto converted = std::from_chars(
        token.data(), token.data() + token.size(), result
    );
    if (converted.ec != std::errc{} ||
        converted.ptr != token.data() + token.size() ||
        result > std::numeric_limits<std::size_t>::max()) {
        throw std::invalid_argument("Invalid Gmsh count or node index.");
    }
    return static_cast<std::size_t>(result);
}

[[nodiscard]] int signed_tag(const std::string& token) {
    int result{};
    const auto converted = std::from_chars(
        token.data(), token.data() + token.size(), result
    );
    if (converted.ec != std::errc{} ||
        converted.ptr != token.data() + token.size()) {
        throw std::invalid_argument("Invalid Gmsh element tag.");
    }
    return result;
}

[[nodiscard]] double finite_coordinate(const std::string& token) {
    std::istringstream stream(token);
    stream.imbue(std::locale::classic());
    double result{};
    stream >> std::noskipws >> result;
    if (!stream || !stream.eof() || !std::isfinite(result)) {
        throw std::invalid_argument("Invalid Gmsh node coordinate.");
    }
    return result;
}

[[nodiscard]] std::size_t single_count(std::istream& input) {
    const auto data = words(line_or_throw(input));
    if (data.size() != 1U) {
        throw std::invalid_argument("Invalid Gmsh section count.");
    }
    return count_integer(data[0]);
}

[[nodiscard]] FaceKey canonical_face(FaceKey face) {
    std::sort(face.begin(), face.end());
    return face;
}

} // namespace

GmshTetrahedralMesh3D read_gmsh22_tetrahedra(std::istream& input) {
    expect_line(input, "$MeshFormat");
    const auto format = words(line_or_throw(input));
    if (format.size() != 3U || format[0] != "2.2" ||
        format[1] != "0" || format[2] != "8") {
        throw std::invalid_argument("Only ASCII Gmsh MSH 2.2 with 8-byte doubles is supported.");
    }
    expect_line(input, "$EndMeshFormat");

    std::string section = line_or_throw(input);
    if (section == "$PhysicalNames") {
        const std::size_t n = single_count(input);
        for (std::size_t i = 0; i < n; ++i) {
            static_cast<void>(line_or_throw(input));
        }
        expect_line(input, "$EndPhysicalNames");
        section = line_or_throw(input);
    }
    if (section != "$Nodes") {
        throw std::invalid_argument("Expected $Nodes in Gmsh MSH 2.2 file.");
    }

    GmshTetrahedralMesh3D result;
    std::map<std::uint64_t, std::size_t> ids;
    const std::size_t node_count = single_count(input);
    for (std::size_t i = 0; i < node_count; ++i) {
        const auto data = words(line_or_throw(input));
        if (data.size() != 4U) {
            throw std::invalid_argument("Gmsh node record must contain ID and three coordinates.");
        }
        const auto id = positive_integer(data[0]);
        if (!ids.emplace(id, result.mesh.nodes.size()).second) {
            throw std::invalid_argument("Duplicate Gmsh node ID.");
        }
        result.mesh.nodes.push_back({
            finite_coordinate(data[1]),
            finite_coordinate(data[2]),
            finite_coordinate(data[3])
        });
    }
    expect_line(input, "$EndNodes");
    expect_line(input, "$Elements");

    const std::size_t element_count = single_count(input);
    std::set<std::uint64_t> element_ids;
    std::map<FaceKey, int> boundary_tags;
    for (std::size_t i = 0; i < element_count; ++i) {
        const auto data = words(line_or_throw(input));
        if (data.size() < 3U) {
            throw std::invalid_argument("Incomplete Gmsh element record.");
        }
        const auto element_id = positive_integer(data[0]);
        if (!element_ids.insert(element_id).second) {
            throw std::invalid_argument("Duplicate Gmsh element ID.");
        }
        const std::size_t type = count_integer(data[1]);
        const std::size_t tag_count = count_integer(data[2]);
        if (tag_count > data.size() - 3U) {
            throw std::invalid_argument("Incomplete Gmsh element tags.");
        }
        std::size_t nodes_per_element{};
        switch (type) {
        case 1U: nodes_per_element = 2U; break; // lower-dimensional line
        case 2U: nodes_per_element = 3U; break; // boundary triangle
        case 4U: nodes_per_element = 4U; break; // P1 tetrahedron
        case 15U: nodes_per_element = 1U; break; // geometric point
        default:
            throw std::invalid_argument("Unsupported Gmsh element type (linear tetrahedra only).");
        }
        if (data.size() - 3U - tag_count != nodes_per_element) {
            throw std::invalid_argument("Incorrect Gmsh element node count.");
        }
        int physical_tag = 0;
        for (std::size_t j = 0; j < tag_count; ++j) {
            const int parsed_tag = signed_tag(data[3U + j]);
            if (j == 0U) {
                physical_tag = parsed_tag;
            }
        }
        std::array<std::size_t, 4> vertices{};
        for (std::size_t j = 0; j < nodes_per_element; ++j) {
            const auto id = positive_integer(data[3U + tag_count + j]);
            const auto node = ids.find(id);
            if (node == ids.end()) {
                throw std::invalid_argument("Gmsh element references an undefined node.");
            }
            vertices[j] = node->second;
        }
        if (type == 4U) {
            result.mesh.elements.push_back({
                vertices[0], vertices[1], vertices[2], vertices[3]
            });
            result.volume_physical_tags.push_back(physical_tag);
        } else if (type == 2U) {
            const FaceKey face = canonical_face({
                vertices[0], vertices[1], vertices[2]
            });
            if (face[0] == face[1] || face[1] == face[2]) {
                throw std::invalid_argument("Repeated node in Gmsh boundary triangle.");
            }
            if (!boundary_tags.emplace(face, physical_tag).second) {
                throw std::invalid_argument("Duplicate Gmsh boundary triangle.");
            }
        }
    }
    expect_line(input, "$EndElements");
    std::string trailing;
    while (std::getline(input, trailing)) {
        if (trailing.find_first_not_of(" \t\r") != std::string::npos) {
            throw std::invalid_argument("Unexpected data after $EndElements.");
        }
    }

    // Extraction validates the complete tetrahedral topology in the same pass.
    result.exterior_faces = fem::extract_boundary_faces_3d(result.mesh);
    result.exterior_physical_tags.reserve(result.exterior_faces.size());
    for (const fem::BoundaryFace3D& face : result.exterior_faces) {
        const FaceKey key = canonical_face(face);
        const auto found = boundary_tags.find(key);
        if (found == boundary_tags.end()) {
            result.exterior_physical_tags.push_back(0);
        } else {
            result.exterior_physical_tags.push_back(found->second);
            boundary_tags.erase(found);
        }
    }
    if (!boundary_tags.empty()) {
        throw std::invalid_argument("Gmsh triangle is not an exterior tetrahedral face.");
    }
    return result;
}

} // namespace nova::io
