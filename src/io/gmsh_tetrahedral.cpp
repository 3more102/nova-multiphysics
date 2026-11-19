#include "nova/io/gmsh_tetrahedral.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <istream>
#include <limits>
#include <locale>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace nova::io {
namespace {

[[nodiscard]] std::string token(std::istream& stream) {
    std::string value;
    if (!(stream >> value)) {
        throw std::invalid_argument("Truncated Gmsh MSH input.");
    }
    return value;
}

void expect(std::istream& stream, const char* expected) {
    if (token(stream) != expected) {
        throw std::invalid_argument(std::string("Expected Gmsh section marker: ") + expected);
    }
}

[[nodiscard]] std::uint64_t unsigned_number(const std::string& text, bool positive = false) {
    std::uint64_t value = 0;
    const auto parsed = std::from_chars(text.data(), text.data() + text.size(), value);
    if (parsed.ec != std::errc{} || parsed.ptr != text.data() + text.size() ||
        (positive && value == 0)) {
        throw std::invalid_argument("Invalid Gmsh unsigned integer.");
    }
    return value;
}

[[nodiscard]] std::size_t count(std::istream& stream) {
    const auto value = unsigned_number(token(stream));
    if (value > std::numeric_limits<std::size_t>::max()) {
        throw std::invalid_argument("Gmsh section count exceeds the platform limit.");
    }
    return static_cast<std::size_t>(value);
}

[[nodiscard]] std::uint64_t tag(std::istream& stream) {
    return unsigned_number(token(stream), true);
}

[[nodiscard]] int signed_number(const std::string& text) {
    int value = 0;
    const auto parsed = std::from_chars(text.data(), text.data() + text.size(), value);
    if (parsed.ec != std::errc{} || parsed.ptr != text.data() + text.size()) {
        throw std::invalid_argument("Invalid Gmsh signed integer.");
    }
    return value;
}

[[nodiscard]] int dimension(std::istream& stream) {
    const int value = signed_number(token(stream));
    if (value < 0 || value > 3) {
        throw std::invalid_argument("Invalid Gmsh entity dimension.");
    }
    return value;
}

void finite_coordinate(const std::string& text) {
    std::istringstream number(text);
    number.imbue(std::locale::classic());
    double value = 0.0;
    number >> std::noskipws >> value;
    if (!number || !number.eof() || !std::isfinite(value)) {
        throw std::invalid_argument("Invalid Gmsh coordinate.");
    }
}

[[nodiscard]] std::vector<std::string> fields(const std::string& line) {
    std::istringstream stream(line);
    stream.imbue(std::locale::classic());
    std::vector<std::string> result;
    std::string item;
    while (stream >> item) result.push_back(item);
    return result;
}

void parse_entities(
    std::istream& stream, std::map<std::pair<int, std::uint64_t>, int>& physical
) {
    std::size_t blocks[4]{};
    for (auto& block : blocks) block = count(stream);
    for (int dim = 0; dim < 4; ++dim) {
        for (std::size_t row = 0; row < blocks[dim]; ++row) {
            std::string line;
            if (!std::getline(stream >> std::ws, line)) {
                throw std::invalid_argument("Truncated Gmsh entity record.");
            }
            const auto data = fields(line);
            const std::size_t coordinate_count = (dim == 0) ? 3U : 6U;
            if (data.size() < coordinate_count + 2U) {
                throw std::invalid_argument("Incomplete Gmsh entity record.");
            }
            const auto entity_tag = unsigned_number(data[0], true);
            for (std::size_t i = 1; i <= coordinate_count; ++i) {
                finite_coordinate(data[i]);
            }
            std::size_t offset = coordinate_count + 1U;
            const auto num_physical = unsigned_number(data[offset++]);
            if (num_physical > 1U) {
                throw std::invalid_argument(
                    "Gmsh entity belongs to multiple physical groups; NOVA needs a single group."
                );
            }
            const std::size_t remaining = data.size() - offset;
            if (num_physical > remaining) {
                throw std::invalid_argument("Truncated Gmsh physical entity tags.");
            }
            int physical_tag = 0;
            if (num_physical != 0U) {
                physical_tag = signed_number(data[offset++]);
                if (physical_tag <= 0) {
                    throw std::invalid_argument("Gmsh physical group tags must be positive.");
                }
            }
            if (dim == 0) {
                if (offset != data.size()) {
                    throw std::invalid_argument("Unexpected point entity fields.");
                }
            } else {
                if (offset >= data.size()) {
                    throw std::invalid_argument("Missing Gmsh entity boundary count.");
                }
                const auto boundaries = unsigned_number(data[offset++]);
                if (boundaries != data.size() - offset) {
                    throw std::invalid_argument("Invalid Gmsh entity boundary count.");
                }
                for (; offset < data.size(); ++offset) {
                    if (signed_number(data[offset]) == 0) {
                        throw std::invalid_argument("Invalid zero bounding entity tag.");
                    }
                }
            }
            if (!physical.emplace(std::make_pair(dim, entity_tag), physical_tag).second) {
                throw std::invalid_argument("Duplicate Gmsh geometric entity.");
            }
        }
    }
    expect(stream, "$EndEntities");
}

[[nodiscard]] std::string convert_msh41(std::istream& stream) {
    expect(stream, "$MeshFormat");
    if (token(stream) != "4.1" || token(stream) != "0" || token(stream) != "8") {
        throw std::invalid_argument("Only ASCII Gmsh MSH 4.1 with 8-byte doubles is supported.");
    }
    expect(stream, "$EndMeshFormat");

    bool has_entities = false, has_names = false, has_nodes = false, has_elements = false;
    std::map<std::pair<int, std::uint64_t>, int> physical;
    std::ostringstream names, nodes, elements;
    std::size_t total_nodes = 0, total_elements = 0;

    std::string section;
    while (stream >> section) {
        if (section == "$PhysicalNames" && !has_names && !has_nodes) {
            has_names = true;
            const auto n = count(stream);
            names << "$PhysicalNames\n" << n << '\n';
            for (std::size_t i = 0; i < n; ++i) {
                std::string entry;
                if (!std::getline(stream >> std::ws, entry)) {
                    throw std::invalid_argument("Truncated Gmsh physical name.");
                }
                names << entry << '\n';
            }
            expect(stream, "$EndPhysicalNames");
            names << "$EndPhysicalNames\n";
        } else if (section == "$Entities" && !has_entities && !has_nodes) {
            has_entities = true;
            parse_entities(stream, physical);
        } else if (section == "$Nodes" && !has_nodes && !has_elements) {
            has_nodes = true;
            const auto nblocks = count(stream);
            total_nodes = count(stream);
            const auto minimum_tag = unsigned_number(token(stream));
            const auto maximum_tag = unsigned_number(token(stream));
            std::size_t seen = 0;
            std::uint64_t observed_min = std::numeric_limits<std::uint64_t>::max();
            std::uint64_t observed_max = 0;
            for (std::size_t b = 0; b < nblocks; ++b) {
                const int dim = dimension(stream);
                const auto entity = tag(stream);
                const auto parametric = count(stream);
                const auto nblock = count(stream);
                if (parametric > 1U || nblock > total_nodes - seen) {
                    throw std::invalid_argument("Invalid Gmsh node block header.");
                }
                if (has_entities && !physical.contains({dim, entity})) {
                    throw std::invalid_argument("Gmsh node block references missing entity.");
                }
                std::vector<std::uint64_t> ids;
                ids.reserve(nblock);
                for (std::size_t i = 0; i < nblock; ++i) {
                    const auto id = tag(stream);
                    ids.push_back(id);
                    observed_min = std::min(observed_min, id);
                    observed_max = std::max(observed_max, id);
                }
                for (std::size_t i = 0; i < nblock; ++i) {
                    nodes << ids[i];
                    for (int coordinate = 0; coordinate < 3; ++coordinate) {
                        const auto value = token(stream);
                        finite_coordinate(value);
                        nodes << ' ' << value;
                    }
                    for (int param = 0; param < (parametric ? dim : 0); ++param) {
                        finite_coordinate(token(stream));
                    }
                    nodes << '\n';
                }
                seen += nblock;
            }
            if (seen != total_nodes ||
                (seen != 0 && (minimum_tag != observed_min || maximum_tag != observed_max)) ||
                (seen == 0 && (minimum_tag != 0 || maximum_tag != 0))) {
                throw std::invalid_argument("Incorrect Gmsh node block totals or tag range.");
            }
            expect(stream, "$EndNodes");
        } else if (section == "$Elements" && !has_elements && has_nodes) {
            has_elements = true;
            const auto nblocks = count(stream);
            total_elements = count(stream);
            const auto minimum_tag = unsigned_number(token(stream));
            const auto maximum_tag = unsigned_number(token(stream));
            std::size_t seen = 0;
            std::uint64_t observed_min = std::numeric_limits<std::uint64_t>::max();
            std::uint64_t observed_max = 0;
            for (std::size_t b = 0; b < nblocks; ++b) {
                const int dim = dimension(stream);
                const auto entity = tag(stream);
                const auto type = count(stream);
                const auto nblock = count(stream);
                if (nblock > total_elements - seen) {
                    throw std::invalid_argument("Invalid Gmsh element block count.");
                }
                if (has_entities && !physical.contains({dim, entity})) {
                    throw std::invalid_argument("Gmsh element block references missing entity.");
                }
                std::size_t nvertices = 0;
                if (dim == 0 && type == 15U) nvertices = 1;
                else if (dim == 1 && type == 1U) nvertices = 2;
                else if (dim == 2 && type == 2U) nvertices = 3;
                else if (dim == 3 && type == 4U) nvertices = 4;
                else throw std::invalid_argument(
                    "Unsupported Gmsh element type or dimension; only first-order simplices."
                );
                const int material = has_entities ? physical.at({dim, entity}) : 0;
                for (std::size_t i = 0; i < nblock; ++i) {
                    const auto id = tag(stream);
                    observed_min = std::min(observed_min, id);
                    observed_max = std::max(observed_max, id);
                    elements << id << ' ' << type << " 2 " << material << ' ' << entity;
                    for (std::size_t j = 0; j < nvertices; ++j) {
                        elements << ' ' << tag(stream);
                    }
                    elements << '\n';
                }
                seen += nblock;
            }
            if (seen != total_elements ||
                (seen != 0 && (minimum_tag != observed_min || maximum_tag != observed_max)) ||
                (seen == 0 && (minimum_tag != 0 || maximum_tag != 0))) {
                throw std::invalid_argument("Incorrect Gmsh element totals or tag range.");
            }
            expect(stream, "$EndElements");
        } else {
            throw std::invalid_argument(
                "Unsupported, repeated or out-of-order Gmsh MSH 4.1 section: " + section
            );
        }
    }
    if (!has_nodes || !has_elements) {
        throw std::invalid_argument("Gmsh MSH 4.1 requires both $Nodes and $Elements.");
    }
    std::ostringstream normalized;
    normalized << "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n" << names.str()
               << "$Nodes\n" << total_nodes << '\n' << nodes.str()
               << "$EndNodes\n$Elements\n" << total_elements << '\n'
               << elements.str() << "$EndElements\n";
    return normalized.str();
}

} // namespace

GmshTetrahedralMesh3D read_gmsh_tetrahedra(std::istream& input) {
    std::ostringstream bytes;
    bytes << input.rdbuf();
    if (input.bad()) {
        throw std::invalid_argument("Could not read Gmsh mesh stream.");
    }
    const std::string content = bytes.str();
    std::istringstream probe(content);
    if (token(probe) != "$MeshFormat") {
        throw std::invalid_argument("Expected Gmsh $MeshFormat.");
    }
    const std::string version = token(probe);
    if (version == "2.2") {
        std::istringstream original(content);
        return read_gmsh22_tetrahedra(original);
    }
    if (version == "4.1") {
        std::istringstream original(content);
        std::istringstream normalized(convert_msh41(original));
        return read_gmsh22_tetrahedra(normalized);
    }
    throw std::invalid_argument("Only Gmsh MSH 2.2 and 4.1 ASCII are supported.");
}

} // namespace nova::io
