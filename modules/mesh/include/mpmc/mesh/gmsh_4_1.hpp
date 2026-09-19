#ifndef MPMC_MESH_GMSH_4_1_HPP
#define MPMC_MESH_GMSH_4_1_HPP

#include <mpmc/mesh/csr_adjacency.hpp>
#include <mpmc/mesh/entity.hpp>
#include <mpmc/mesh/face_boundary.hpp>
#include <mpmc/mesh/geometry_2d.hpp>
#include <mpmc/mesh/topology.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

namespace mpmc::mesh {

struct GmshPhysicalName {
    int dimension;
    std::uint32_t tag;
    std::string name;
};

struct GmshCellPhysicalGroups {
    GlobalEntityId cell_global_id;
    std::vector<std::uint32_t> physical_tags;
};

struct Gmsh41ImportResult {
    Topology topology;
    Geometry2D geometry;
    FaceBoundarySnapshot face_boundary;
    std::vector<GmshPhysicalName> physical_names;
    std::vector<GmshCellPhysicalGroups> cell_physical_groups;
};

namespace gmsh41_detail {

struct EntityKey {
    int dimension;
    std::int64_t tag;

    [[nodiscard]] friend bool operator<(const EntityKey& left,
                                        const EntityKey& right) noexcept {
        return std::tie(left.dimension, left.tag) <
               std::tie(right.dimension, right.tag);
    }
};

struct EntityInfo {
    std::vector<std::uint32_t> physical_tags;
};

struct NodeRecord {
    std::uint64_t tag;
    double x;
    double y;
    double z;
};

struct ElementRecord {
    std::uint64_t tag;
    int dimension;
    std::int64_t entity_tag;
    int element_type;
    std::vector<std::uint64_t> node_tags;
};

struct CellRecord {
    std::uint64_t tag;
    std::int64_t entity_tag;
    std::vector<std::uint64_t> node_tags;
};

struct LineRecord {
    std::uint64_t tag;
    std::int64_t entity_tag;
    std::array<std::uint64_t, 2> node_tags;
};

struct EdgeKey {
    std::uint64_t first;
    std::uint64_t second;

    [[nodiscard]] friend bool operator<(const EdgeKey& left,
                                        const EdgeKey& right) noexcept {
        return std::tie(left.first, left.second) <
               std::tie(right.first, right.second);
    }
};

struct FaceBuild {
    EdgeKey key;
    std::vector<std::size_t> adjacent_cells;
    std::uint64_t global_id{0U};
    std::uint32_t physical_tag{0U};
    bool has_line_element{false};
};

struct CellBuild {
    std::uint64_t global_id;
    std::int64_t entity_tag;
    std::vector<std::uint64_t> node_tags;
    std::vector<LocalIndex> vertex_locals;
    std::vector<EdgeKey> edge_keys;
    Coordinate2D centroid;
    double area_m2;
};

template <typename T>
inline void read_value(std::istream& input, T& value, const char* message) {
    if (!(input >> value)) {
        throw std::invalid_argument(message);
    }
}

inline void require_token(std::istream& input,
                          std::string_view expected,
                          const char* message) {
    std::string actual;
    read_value(input, actual, message);
    if (actual != expected) {
        throw std::invalid_argument(message);
    }
}

[[nodiscard]] inline std::size_t checked_size(std::uint64_t value,
                                               const char* message) {
    if (value >
        static_cast<std::uint64_t>(
            std::numeric_limits<std::size_t>::max())) {
        throw std::length_error(message);
    }
    return static_cast<std::size_t>(value);
}

[[nodiscard]] inline LocalIndex checked_local(std::size_t value,
                                               const char* message) {
    if (value >
        static_cast<std::size_t>(
            std::numeric_limits<LocalIndex::value_type>::max())) {
        throw std::length_error(message);
    }
    return LocalIndex{
        static_cast<LocalIndex::value_type>(value)};
}

[[nodiscard]] inline CsrAdjacency::Offset checked_offset(
    std::size_t value,
    const char* message) {
    if (value >
        static_cast<std::size_t>(
            std::numeric_limits<CsrAdjacency::Offset>::max())) {
        throw std::length_error(message);
    }
    return static_cast<CsrAdjacency::Offset>(value);
}

[[nodiscard]] inline std::uint32_t checked_physical_tag(
    std::int64_t value,
    const char* message) {
    if (value <= 0 ||
        static_cast<std::uint64_t>(value) >
            static_cast<std::uint64_t>(
                std::numeric_limits<std::uint32_t>::max())) {
        throw std::invalid_argument(message);
    }
    return static_cast<std::uint32_t>(value);
}

inline void require_positive_tag(std::int64_t value,
                                 const char* message) {
    if (value <= 0) {
        throw std::invalid_argument(message);
    }
}

[[nodiscard]] inline EdgeKey edge_key(std::uint64_t left,
                                      std::uint64_t right) {
    if (left == right) {
        throw std::invalid_argument(
            "mpmc::mesh::import_gmsh_4_1_ascii: degenerate edge with repeated node tag");
    }
    return left < right ? EdgeKey{left, right}
                        : EdgeKey{right, left};
}

[[nodiscard]] inline double cross(Coordinate2D a,
                                  Coordinate2D b,
                                  Coordinate2D c) {
    return (b.x_m - a.x_m) * (c.y_m - a.y_m) -
           (b.y_m - a.y_m) * (c.x_m - a.x_m);
}

[[nodiscard]] inline bool segments_strictly_intersect(
    Coordinate2D a,
    Coordinate2D b,
    Coordinate2D c,
    Coordinate2D d,
    double tolerance) {
    const double c1 = cross(a, b, c);
    const double c2 = cross(a, b, d);
    const double c3 = cross(c, d, a);
    const double c4 = cross(c, d, b);
    if (std::abs(c1) <= tolerance ||
        std::abs(c2) <= tolerance ||
        std::abs(c3) <= tolerance ||
        std::abs(c4) <= tolerance) {
        return false;
    }
    return (c1 > 0.0) != (c2 > 0.0) &&
           (c3 > 0.0) != (c4 > 0.0);
}

struct PolygonMetric {
    Coordinate2D centroid;
    double area;
};

[[nodiscard]] inline PolygonMetric polygon_metric(
    const std::vector<Coordinate2D>& coordinates) {
    if (coordinates.size() != 3U &&
        coordinates.size() != 4U) {
        throw std::invalid_argument(
            "mpmc::mesh::import_gmsh_4_1_ascii: only linear triangles and quadrilaterals are supported");
    }

    double magnitude_scale = 1.0;
    for (const auto coordinate : coordinates) {
        if (!std::isfinite(coordinate.x_m) ||
            !std::isfinite(coordinate.y_m)) {
            throw std::invalid_argument(
                "mpmc::mesh::import_gmsh_4_1_ascii: non-finite cell coordinate");
        }
        magnitude_scale = std::max(
            magnitude_scale,
            std::max(
                std::abs(coordinate.x_m),
                std::abs(coordinate.y_m)));
    }
    const double tolerance =
        256.0 * std::numeric_limits<double>::epsilon() *
        magnitude_scale * magnitude_scale;

    if (coordinates.size() == 4U) {
        if (segments_strictly_intersect(
                coordinates[0], coordinates[1],
                coordinates[2], coordinates[3],
                tolerance) ||
            segments_strictly_intersect(
                coordinates[1], coordinates[2],
                coordinates[3], coordinates[0],
                tolerance)) {
            throw std::invalid_argument(
                "mpmc::mesh::import_gmsh_4_1_ascii: self-intersecting quadrilateral is unsupported");
        }

        double turn_sign = 0.0;
        for (std::size_t i = 0U; i < 4U; ++i) {
            const double turn = cross(
                coordinates[i],
                coordinates[(i + 1U) % 4U],
                coordinates[(i + 2U) % 4U]);
            if (std::abs(turn) <= tolerance) {
                throw std::invalid_argument(
                    "mpmc::mesh::import_gmsh_4_1_ascii: degenerate quadrilateral corner");
            }
            if (turn_sign == 0.0) {
                turn_sign = turn;
            } else if ((turn > 0.0) !=
                       (turn_sign > 0.0)) {
                throw std::invalid_argument(
                    "mpmc::mesh::import_gmsh_4_1_ascii: concave quadrilateral is unsupported by the minimal baseline");
            }
        }
    }

    double twice_area = 0.0;
    double centroid_numerator_x = 0.0;
    double centroid_numerator_y = 0.0;
    for (std::size_t i = 0U;
         i < coordinates.size();
         ++i) {
        const auto& a = coordinates[i];
        const auto& b =
            coordinates[
                (i + 1U) % coordinates.size()];
        const double edge_cross =
            a.x_m * b.y_m - b.x_m * a.y_m;
        twice_area += edge_cross;
        centroid_numerator_x +=
            (a.x_m + b.x_m) * edge_cross;
        centroid_numerator_y +=
            (a.y_m + b.y_m) * edge_cross;
    }
    if (!std::isfinite(twice_area) ||
        std::abs(twice_area) <= tolerance) {
        throw std::invalid_argument(
            "mpmc::mesh::import_gmsh_4_1_ascii: degenerate cell area");
    }

    const double centroid_denominator =
        3.0 * twice_area;
    const Coordinate2D centroid{
        centroid_numerator_x / centroid_denominator,
        centroid_numerator_y / centroid_denominator};
    const double area = 0.5 * std::abs(twice_area);
    if (!std::isfinite(centroid.x_m) ||
        !std::isfinite(centroid.y_m) ||
        !std::isfinite(area) ||
        area <= 0.0) {
        throw std::invalid_argument(
            "mpmc::mesh::import_gmsh_4_1_ascii: invalid polygon metric");
    }
    return PolygonMetric{centroid, area};
}

[[nodiscard]] inline const EntityInfo* find_entity(
    const std::map<EntityKey, EntityInfo>& entities,
    int dimension,
    std::int64_t tag) {
    const auto found =
        entities.find(EntityKey{dimension, tag});
    return found == entities.end()
               ? nullptr
               : &found->second;
}

inline void parse_physical_names(
    std::istream& input,
    std::map<EntityKey, std::string>& physical_names) {
    std::uint64_t raw_count = 0U;
    read_value(
        input, raw_count,
        "mpmc::mesh::import_gmsh_4_1_ascii: invalid $PhysicalNames count");
    const std::size_t count = checked_size(
        raw_count,
        "mpmc::mesh::import_gmsh_4_1_ascii: physical-name count overflow");

    for (std::size_t i = 0U; i < count; ++i) {
        int dimension = -1;
        std::int64_t raw_tag = 0;
        std::string name;
        read_value(
            input, dimension,
            "mpmc::mesh::import_gmsh_4_1_ascii: invalid physical-name dimension");
        read_value(
            input, raw_tag,
            "mpmc::mesh::import_gmsh_4_1_ascii: invalid physical-name tag");
        if (!(input >> std::quoted(name))) {
            throw std::invalid_argument(
                "mpmc::mesh::import_gmsh_4_1_ascii: invalid quoted physical-group name");
        }
        if (dimension < 0 || dimension > 2) {
            throw std::invalid_argument(
                "mpmc::mesh::import_gmsh_4_1_ascii: only 0D/1D/2D physical names are supported");
        }
        const std::uint32_t tag = checked_physical_tag(
            raw_tag,
            "mpmc::mesh::import_gmsh_4_1_ascii: physical group tag must be positive uint32");
        const EntityKey key{dimension, raw_tag};
        if (!physical_names.emplace(key, name).second) {
            throw std::invalid_argument(
                "mpmc::mesh::import_gmsh_4_1_ascii: duplicate physical-group name key");
        }
        (void)tag;
    }
    require_token(
        input, "$EndPhysicalNames",
        "mpmc::mesh::import_gmsh_4_1_ascii: missing $EndPhysicalNames");
}

inline std::vector<std::uint32_t> read_physical_tags(
    std::istream& input) {
    std::uint64_t raw_count = 0U;
    read_value(
        input, raw_count,
        "mpmc::mesh::import_gmsh_4_1_ascii: invalid entity physical-tag count");
    const std::size_t count = checked_size(
        raw_count,
        "mpmc::mesh::import_gmsh_4_1_ascii: entity physical-tag count overflow");

    std::vector<std::uint32_t> tags;
    tags.reserve(count);
    for (std::size_t i = 0U; i < count; ++i) {
        std::int64_t raw_tag = 0;
        read_value(
            input, raw_tag,
            "mpmc::mesh::import_gmsh_4_1_ascii: invalid entity physical tag");
        tags.push_back(checked_physical_tag(
            raw_tag,
            "mpmc::mesh::import_gmsh_4_1_ascii: entity physical tag must be positive uint32"));
    }
    std::sort(tags.begin(), tags.end());
    if (std::adjacent_find(
            tags.begin(), tags.end()) !=
        tags.end()) {
        throw std::invalid_argument(
            "mpmc::mesh::import_gmsh_4_1_ascii: duplicate physical tag on entity");
    }
    return tags;
}

inline void read_signed_tags(std::istream& input,
                             std::size_t count) {
    for (std::size_t i = 0U; i < count; ++i) {
        std::int64_t ignored = 0;
        read_value(
            input, ignored,
            "mpmc::mesh::import_gmsh_4_1_ascii: invalid bounding entity tag");
        if (ignored == 0) {
            throw std::invalid_argument(
                "mpmc::mesh::import_gmsh_4_1_ascii: bounding entity tag cannot be zero");
        }
    }
}

inline void parse_entities(
    std::istream& input,
    std::map<EntityKey, EntityInfo>& entities) {
    std::array<std::uint64_t, 4> raw_counts{};
    for (auto& count : raw_counts) {
        read_value(
            input, count,
            "mpmc::mesh::import_gmsh_4_1_ascii: invalid $Entities header");
    }
    if (raw_counts[3] != 0U) {
        throw std::invalid_argument(
            "mpmc::mesh::import_gmsh_4_1_ascii: volume entities are unsupported in the 2D baseline");
    }

    for (int dimension = 0; dimension <= 2; ++dimension) {
        const std::size_t count = checked_size(
            raw_counts[static_cast<std::size_t>(dimension)],
            "mpmc::mesh::import_gmsh_4_1_ascii: entity count overflow");
        for (std::size_t i = 0U; i < count; ++i) {
            std::int64_t tag = 0;
            read_value(
                input, tag,
                "mpmc::mesh::import_gmsh_4_1_ascii: invalid entity tag");
            require_positive_tag(
                tag,
                "mpmc::mesh::import_gmsh_4_1_ascii: entity tag must be positive");

            if (dimension == 0) {
                double x = 0.0;
                double y = 0.0;
                double z = 0.0;
                read_value(input, x,
                           "mpmc::mesh::import_gmsh_4_1_ascii: invalid point entity x");
                read_value(input, y,
                           "mpmc::mesh::import_gmsh_4_1_ascii: invalid point entity y");
                read_value(input, z,
                           "mpmc::mesh::import_gmsh_4_1_ascii: invalid point entity z");
                if (!std::isfinite(x) ||
                    !std::isfinite(y) ||
                    !std::isfinite(z)) {
                    throw std::invalid_argument(
                        "mpmc::mesh::import_gmsh_4_1_ascii: non-finite point entity coordinate");
                }
            } else {
                for (std::size_t coordinate = 0U;
                     coordinate < 6U;
                     ++coordinate) {
                    double value = 0.0;
                    read_value(
                        input, value,
                        "mpmc::mesh::import_gmsh_4_1_ascii: invalid entity bounding box");
                    if (!std::isfinite(value)) {
                        throw std::invalid_argument(
                            "mpmc::mesh::import_gmsh_4_1_ascii: non-finite entity bounding box");
                    }
                }
            }

            EntityInfo info;
            info.physical_tags =
                read_physical_tags(input);
            if (dimension > 0) {
                std::uint64_t raw_bounding_count = 0U;
                read_value(
                    input, raw_bounding_count,
                    "mpmc::mesh::import_gmsh_4_1_ascii: invalid entity bounding count");
                read_signed_tags(
                    input,
                    checked_size(
                        raw_bounding_count,
                        "mpmc::mesh::import_gmsh_4_1_ascii: entity bounding count overflow"));
            }

            if (!entities
                     .emplace(
                         EntityKey{dimension, tag},
                         std::move(info))
                     .second) {
                throw std::invalid_argument(
                    "mpmc::mesh::import_gmsh_4_1_ascii: duplicate model entity");
            }
        }
    }
    require_token(
        input, "$EndEntities",
        "mpmc::mesh::import_gmsh_4_1_ascii: missing $EndEntities");
}

inline void parse_nodes(
    std::istream& input,
    std::vector<NodeRecord>& nodes) {
    std::uint64_t raw_blocks = 0U;
    std::uint64_t raw_nodes = 0U;
    std::uint64_t declared_min = 0U;
    std::uint64_t declared_max = 0U;
    read_value(input, raw_blocks,
               "mpmc::mesh::import_gmsh_4_1_ascii: invalid $Nodes block count");
    read_value(input, raw_nodes,
               "mpmc::mesh::import_gmsh_4_1_ascii: invalid $Nodes node count");
    read_value(input, declared_min,
               "mpmc::mesh::import_gmsh_4_1_ascii: invalid $Nodes minimum tag");
    read_value(input, declared_max,
               "mpmc::mesh::import_gmsh_4_1_ascii: invalid $Nodes maximum tag");

    const std::size_t block_count = checked_size(
        raw_blocks,
        "mpmc::mesh::import_gmsh_4_1_ascii: node-block count overflow");
    const std::size_t node_count = checked_size(
        raw_nodes,
        "mpmc::mesh::import_gmsh_4_1_ascii: node count overflow");
    nodes.reserve(node_count);
    std::set<std::uint64_t> seen_tags;

    for (std::size_t block = 0U;
         block < block_count;
         ++block) {
        int entity_dimension = -1;
        std::int64_t entity_tag = 0;
        int parametric = -1;
        std::uint64_t raw_block_nodes = 0U;
        read_value(input, entity_dimension,
                   "mpmc::mesh::import_gmsh_4_1_ascii: invalid node-block entity dimension");
        read_value(input, entity_tag,
                   "mpmc::mesh::import_gmsh_4_1_ascii: invalid node-block entity tag");
        read_value(input, parametric,
                   "mpmc::mesh::import_gmsh_4_1_ascii: invalid node-block parametric flag");
        read_value(input, raw_block_nodes,
                   "mpmc::mesh::import_gmsh_4_1_ascii: invalid node-block size");
        if (entity_dimension < 0 ||
            entity_dimension > 2) {
            throw std::invalid_argument(
                "mpmc::mesh::import_gmsh_4_1_ascii: node blocks above dimension 2 are unsupported");
        }
        require_positive_tag(
            entity_tag,
            "mpmc::mesh::import_gmsh_4_1_ascii: node-block entity tag must be positive");
        if (parametric != 0) {
            throw std::invalid_argument(
                "mpmc::mesh::import_gmsh_4_1_ascii: parametric node coordinates are unsupported");
        }

        const std::size_t block_nodes = checked_size(
            raw_block_nodes,
            "mpmc::mesh::import_gmsh_4_1_ascii: node-block size overflow");
        std::vector<std::uint64_t> tags(block_nodes);
        for (auto& tag : tags) {
            read_value(
                input, tag,
                "mpmc::mesh::import_gmsh_4_1_ascii: invalid node tag");
            if (tag == 0U ||
                !seen_tags.insert(tag).second) {
                throw std::invalid_argument(
                    "mpmc::mesh::import_gmsh_4_1_ascii: node tags must be unique and positive");
            }
        }
        for (const auto tag : tags) {
            NodeRecord record{tag, 0.0, 0.0, 0.0};
            read_value(input, record.x,
                       "mpmc::mesh::import_gmsh_4_1_ascii: invalid node x coordinate");
            read_value(input, record.y,
                       "mpmc::mesh::import_gmsh_4_1_ascii: invalid node y coordinate");
            read_value(input, record.z,
                       "mpmc::mesh::import_gmsh_4_1_ascii: invalid node z coordinate");
            if (!std::isfinite(record.x) ||
                !std::isfinite(record.y) ||
                !std::isfinite(record.z)) {
                throw std::invalid_argument(
                    "mpmc::mesh::import_gmsh_4_1_ascii: node coordinates must be finite");
            }
            nodes.push_back(record);
        }
    }

    if (nodes.size() != node_count) {
        throw std::invalid_argument(
            "mpmc::mesh::import_gmsh_4_1_ascii: parsed node count does not match $Nodes header");
    }
    if (nodes.empty()) {
        throw std::invalid_argument(
            "mpmc::mesh::import_gmsh_4_1_ascii: mesh contains no nodes");
    }
    const auto [minimum, maximum] =
        std::minmax_element(
            nodes.begin(),
            nodes.end(),
            [](const NodeRecord& left,
               const NodeRecord& right) {
                return left.tag < right.tag;
            });
    if (minimum->tag != declared_min ||
        maximum->tag != declared_max) {
        throw std::invalid_argument(
            "mpmc::mesh::import_gmsh_4_1_ascii: node tag range does not match $Nodes header");
    }
    require_token(
        input, "$EndNodes",
        "mpmc::mesh::import_gmsh_4_1_ascii: missing $EndNodes");
}

[[nodiscard]] inline std::size_t element_node_count(
    int dimension,
    int element_type) {
    if (dimension == 0 && element_type == 15) return 1U;
    if (dimension == 1 && element_type == 1) return 2U;
    if (dimension == 2 && element_type == 2) return 3U;
    if (dimension == 2 && element_type == 3) return 4U;
    throw std::invalid_argument(
        "mpmc::mesh::import_gmsh_4_1_ascii: unsupported element type; baseline accepts point, 2-node line, 3-node triangle and 4-node quad only");
}

inline void parse_elements(
    std::istream& input,
    std::vector<ElementRecord>& elements) {
    std::uint64_t raw_blocks = 0U;
    std::uint64_t raw_elements = 0U;
    std::uint64_t declared_min = 0U;
    std::uint64_t declared_max = 0U;
    read_value(input, raw_blocks,
               "mpmc::mesh::import_gmsh_4_1_ascii: invalid $Elements block count");
    read_value(input, raw_elements,
               "mpmc::mesh::import_gmsh_4_1_ascii: invalid $Elements element count");
    read_value(input, declared_min,
               "mpmc::mesh::import_gmsh_4_1_ascii: invalid $Elements minimum tag");
    read_value(input, declared_max,
               "mpmc::mesh::import_gmsh_4_1_ascii: invalid $Elements maximum tag");

    const std::size_t block_count = checked_size(
        raw_blocks,
        "mpmc::mesh::import_gmsh_4_1_ascii: element-block count overflow");
    const std::size_t element_count = checked_size(
        raw_elements,
        "mpmc::mesh::import_gmsh_4_1_ascii: element count overflow");
    elements.reserve(element_count);
    std::set<std::uint64_t> seen_tags;

    for (std::size_t block = 0U;
         block < block_count;
         ++block) {
        int entity_dimension = -1;
        std::int64_t entity_tag = 0;
        int element_type = 0;
        std::uint64_t raw_block_elements = 0U;
        read_value(input, entity_dimension,
                   "mpmc::mesh::import_gmsh_4_1_ascii: invalid element-block entity dimension");
        read_value(input, entity_tag,
                   "mpmc::mesh::import_gmsh_4_1_ascii: invalid element-block entity tag");
        read_value(input, element_type,
                   "mpmc::mesh::import_gmsh_4_1_ascii: invalid element type");
        read_value(input, raw_block_elements,
                   "mpmc::mesh::import_gmsh_4_1_ascii: invalid element-block size");
        if (entity_dimension < 0 ||
            entity_dimension > 2) {
            throw std::invalid_argument(
                "mpmc::mesh::import_gmsh_4_1_ascii: element blocks above dimension 2 are unsupported");
        }
        require_positive_tag(
            entity_tag,
            "mpmc::mesh::import_gmsh_4_1_ascii: element-block entity tag must be positive");

        const std::size_t nodes_per_element =
            element_node_count(
                entity_dimension,
                element_type);
        const std::size_t block_elements = checked_size(
            raw_block_elements,
            "mpmc::mesh::import_gmsh_4_1_ascii: element-block size overflow");

        for (std::size_t element = 0U;
             element < block_elements;
             ++element) {
            ElementRecord record;
            record.dimension = entity_dimension;
            record.entity_tag = entity_tag;
            record.element_type = element_type;
            read_value(
                input, record.tag,
                "mpmc::mesh::import_gmsh_4_1_ascii: invalid element tag");
            if (record.tag == 0U ||
                !seen_tags.insert(record.tag).second) {
                throw std::invalid_argument(
                    "mpmc::mesh::import_gmsh_4_1_ascii: element tags must be unique and positive");
            }
            record.node_tags.resize(
                nodes_per_element);
            for (auto& node_tag : record.node_tags) {
                read_value(
                    input, node_tag,
                    "mpmc::mesh::import_gmsh_4_1_ascii: invalid element node tag");
                if (node_tag == 0U) {
                    throw std::invalid_argument(
                        "mpmc::mesh::import_gmsh_4_1_ascii: element node tag must be positive");
                }
            }
            elements.push_back(std::move(record));
        }
    }

    if (elements.size() != element_count) {
        throw std::invalid_argument(
            "mpmc::mesh::import_gmsh_4_1_ascii: parsed element count does not match $Elements header");
    }
    if (elements.empty()) {
        throw std::invalid_argument(
            "mpmc::mesh::import_gmsh_4_1_ascii: mesh contains no elements");
    }
    const auto [minimum, maximum] =
        std::minmax_element(
            elements.begin(),
            elements.end(),
            [](const ElementRecord& left,
               const ElementRecord& right) {
                return left.tag < right.tag;
            });
    if (minimum->tag != declared_min ||
        maximum->tag != declared_max) {
        throw std::invalid_argument(
            "mpmc::mesh::import_gmsh_4_1_ascii: element tag range does not match $Elements header");
    }
    require_token(
        input, "$EndElements",
        "mpmc::mesh::import_gmsh_4_1_ascii: missing $EndElements");
}

[[nodiscard]] inline std::vector<std::uint32_t>
entity_physical_tags(
    const std::map<EntityKey, EntityInfo>& entities,
    bool entities_present,
    int dimension,
    std::int64_t tag) {
    const auto* entity =
        find_entity(entities, dimension, tag);
    if (entity == nullptr) {
        if (entities_present) {
            throw std::invalid_argument(
                "mpmc::mesh::import_gmsh_4_1_ascii: element references missing model entity");
        }
        return {};
    }
    return entity->physical_tags;
}

} // namespace gmsh41_detail

[[nodiscard]] inline Gmsh41ImportResult import_gmsh_4_1_ascii(
    std::string_view content,
    double coordinate_scale_to_m) {
    using namespace gmsh41_detail;

    if (!std::isfinite(coordinate_scale_to_m) ||
        coordinate_scale_to_m <= 0.0) {
        throw std::invalid_argument(
            "mpmc::mesh::import_gmsh_4_1_ascii: coordinate_scale_to_m must be finite and positive");
    }

    std::istringstream input{std::string(content)};
    bool mesh_format_seen = false;
    bool physical_names_seen = false;
    bool entities_seen = false;
    bool nodes_seen = false;
    bool elements_seen = false;

    std::map<EntityKey, std::string> physical_name_map;
    std::map<EntityKey, EntityInfo> entities;
    std::vector<NodeRecord> nodes;
    std::vector<ElementRecord> elements;

    std::string section;
    while (input >> section) {
        if (section == "$MeshFormat") {
            if (mesh_format_seen) {
                throw std::invalid_argument(
                    "mpmc::mesh::import_gmsh_4_1_ascii: duplicate $MeshFormat section");
            }
            double version = 0.0;
            int file_type = -1;
            int data_size = 0;
            read_value(
                input, version,
                "mpmc::mesh::import_gmsh_4_1_ascii: invalid MSH version");
            read_value(
                input, file_type,
                "mpmc::mesh::import_gmsh_4_1_ascii: invalid MSH file type");
            read_value(
                input, data_size,
                "mpmc::mesh::import_gmsh_4_1_ascii: invalid MSH data size");
            if (std::abs(version - 4.1) >
                8.0 * std::numeric_limits<double>::epsilon()) {
                throw std::invalid_argument(
                    "mpmc::mesh::import_gmsh_4_1_ascii: only MSH version 4.1 is supported");
            }
            if (file_type != 0) {
                throw std::invalid_argument(
                    "mpmc::mesh::import_gmsh_4_1_ascii: binary MSH is unsupported; ASCII is required");
            }
            if (data_size != 8) {
                throw std::invalid_argument(
                    "mpmc::mesh::import_gmsh_4_1_ascii: only 8-byte MSH data-size files are supported");
            }
            require_token(
                input, "$EndMeshFormat",
                "mpmc::mesh::import_gmsh_4_1_ascii: missing $EndMeshFormat");
            mesh_format_seen = true;
        } else if (section == "$PhysicalNames") {
            if (physical_names_seen) {
                throw std::invalid_argument(
                    "mpmc::mesh::import_gmsh_4_1_ascii: duplicate $PhysicalNames section");
            }
            parse_physical_names(
                input, physical_name_map);
            physical_names_seen = true;
        } else if (section == "$Entities") {
            if (entities_seen) {
                throw std::invalid_argument(
                    "mpmc::mesh::import_gmsh_4_1_ascii: duplicate $Entities section");
            }
            parse_entities(input, entities);
            entities_seen = true;
        } else if (section == "$Nodes") {
            if (nodes_seen) {
                throw std::invalid_argument(
                    "mpmc::mesh::import_gmsh_4_1_ascii: duplicate $Nodes section");
            }
            parse_nodes(input, nodes);
            nodes_seen = true;
        } else if (section == "$Elements") {
            if (elements_seen) {
                throw std::invalid_argument(
                    "mpmc::mesh::import_gmsh_4_1_ascii: duplicate $Elements section");
            }
            parse_elements(input, elements);
            elements_seen = true;
        } else {
            throw std::invalid_argument(
                "mpmc::mesh::import_gmsh_4_1_ascii: unsupported or unexpected MSH section");
        }
    }

    if (!mesh_format_seen ||
        !nodes_seen ||
        !elements_seen) {
        throw std::invalid_argument(
            "mpmc::mesh::import_gmsh_4_1_ascii: $MeshFormat, $Nodes and $Elements are required");
    }
    if (physical_names_seen && !entities_seen) {
        throw std::invalid_argument(
            "mpmc::mesh::import_gmsh_4_1_ascii: $Entities is required when Physical Groups are declared");
    }

    std::sort(
        nodes.begin(),
        nodes.end(),
        [](const NodeRecord& left,
           const NodeRecord& right) {
            return left.tag < right.tag;
        });

    std::map<std::uint64_t, std::size_t>
        node_local_by_tag;
    std::vector<Coordinate2D> vertex_coordinates;
    vertex_coordinates.reserve(nodes.size());
    Topology::EntityIds entity_ids;
    entity_ids.vertices.reserve(nodes.size());

    for (std::size_t local = 0U;
         local < nodes.size();
         ++local) {
        const auto& node = nodes[local];
        const double x_m =
            node.x * coordinate_scale_to_m;
        const double y_m =
            node.y * coordinate_scale_to_m;
        const double z_m =
            node.z * coordinate_scale_to_m;
        if (!std::isfinite(x_m) ||
            !std::isfinite(y_m) ||
            !std::isfinite(z_m)) {
            throw std::invalid_argument(
                "mpmc::mesh::import_gmsh_4_1_ascii: scaled node coordinate is non-finite");
        }
        const double planar_scale = std::max(
            1.0,
            std::max(
                std::abs(x_m),
                std::abs(y_m)));
        const double planar_tolerance =
            256.0 *
            std::numeric_limits<double>::epsilon() *
            planar_scale;
        if (std::abs(z_m) > planar_tolerance) {
            throw std::invalid_argument(
                "mpmc::mesh::import_gmsh_4_1_ascii: Geometry2D baseline requires nodes in the XY plane");
        }
        node_local_by_tag.emplace(
            node.tag, local);
        entity_ids.vertices.emplace_back(
            node.tag);
        vertex_coordinates.push_back(
            Coordinate2D{x_m, y_m});
    }

    std::vector<CellRecord> cell_records;
    std::vector<LineRecord> line_records;
    std::uint64_t maximum_element_tag = 0U;
    for (const auto& element : elements) {
        maximum_element_tag =
            std::max(
                maximum_element_tag,
                element.tag);
        if (entities_seen &&
            find_entity(
                entities,
                element.dimension,
                element.entity_tag) == nullptr) {
            throw std::invalid_argument(
                "mpmc::mesh::import_gmsh_4_1_ascii: element block references entity absent from $Entities");
        }

        if (element.dimension == 1) {
            line_records.push_back(
                LineRecord{
                    element.tag,
                    element.entity_tag,
                    {element.node_tags[0],
                     element.node_tags[1]}});
        } else if (element.dimension == 2) {
            cell_records.push_back(
                CellRecord{
                    element.tag,
                    element.entity_tag,
                    element.node_tags});
        }
    }
    if (cell_records.empty()) {
        throw std::invalid_argument(
            "mpmc::mesh::import_gmsh_4_1_ascii: no 2D triangle/quad cells found");
    }

    std::sort(
        cell_records.begin(),
        cell_records.end(),
        [](const CellRecord& left,
           const CellRecord& right) {
            return left.tag < right.tag;
        });
    std::sort(
        line_records.begin(),
        line_records.end(),
        [](const LineRecord& left,
           const LineRecord& right) {
            return left.tag < right.tag;
        });

    std::vector<CellBuild> cells;
    cells.reserve(cell_records.size());
    std::map<EdgeKey, FaceBuild> faces_by_edge;

    for (std::size_t cell_local = 0U;
         cell_local < cell_records.size();
         ++cell_local) {
        const auto& record =
            cell_records[cell_local];
        std::set<std::uint64_t> unique_nodes(
            record.node_tags.begin(),
            record.node_tags.end());
        if (unique_nodes.size() !=
            record.node_tags.size()) {
            throw std::invalid_argument(
                "mpmc::mesh::import_gmsh_4_1_ascii: cell contains repeated node tags");
        }

        CellBuild cell;
        cell.global_id = record.tag;
        cell.entity_tag = record.entity_tag;
        cell.node_tags = record.node_tags;
        cell.vertex_locals.reserve(
            record.node_tags.size());
        cell.edge_keys.reserve(
            record.node_tags.size());

        std::vector<Coordinate2D> coordinates;
        coordinates.reserve(
            record.node_tags.size());
        for (const auto node_tag :
             record.node_tags) {
            const auto found =
                node_local_by_tag.find(node_tag);
            if (found == node_local_by_tag.end()) {
                throw std::invalid_argument(
                    "mpmc::mesh::import_gmsh_4_1_ascii: cell references missing node");
            }
            cell.vertex_locals.push_back(
                checked_local(
                    found->second,
                    "mpmc::mesh::import_gmsh_4_1_ascii: vertex local index overflow"));
            coordinates.push_back(
                vertex_coordinates[
                    found->second]);
        }
        const auto metric =
            polygon_metric(coordinates);
        cell.centroid = metric.centroid;
        cell.area_m2 = metric.area;

        for (std::size_t edge = 0U;
             edge < record.node_tags.size();
             ++edge) {
            const EdgeKey key = edge_key(
                record.node_tags[edge],
                record.node_tags[
                    (edge + 1U) %
                    record.node_tags.size()]);
            cell.edge_keys.push_back(key);
            auto [found, inserted] =
                faces_by_edge.emplace(
                    key,
                    FaceBuild{
                        key, {}, 0U, 0U, false});
            found->second.adjacent_cells.push_back(
                cell_local);
            if (found->second.adjacent_cells.size() >
                2U) {
                throw std::invalid_argument(
                    "mpmc::mesh::import_gmsh_4_1_ascii: non-manifold edge belongs to more than two cells");
            }
            (void)inserted;
        }
        cells.push_back(std::move(cell));
    }

    for (const auto& line : line_records) {
        const EdgeKey key = edge_key(
            line.node_tags[0],
            line.node_tags[1]);
        auto found = faces_by_edge.find(key);
        if (found == faces_by_edge.end()) {
            throw std::invalid_argument(
                "mpmc::mesh::import_gmsh_4_1_ascii: line element does not coincide with a triangle/quad cell edge");
        }
        if (found->second.has_line_element) {
            throw std::invalid_argument(
                "mpmc::mesh::import_gmsh_4_1_ascii: multiple line elements map to the same face");
        }

        const auto physical_tags =
            entity_physical_tags(
                entities,
                entities_seen,
                1,
                line.entity_tag);
        if (physical_tags.size() > 1U) {
            throw std::invalid_argument(
                "mpmc::mesh::import_gmsh_4_1_ascii: a boundary curve belongs to multiple Physical Groups but FaceBoundarySnapshot stores one PhysicalTag");
        }
        const std::uint32_t physical_tag =
            physical_tags.empty()
                ? 0U
                : physical_tags.front();
        if (found->second.adjacent_cells.size() ==
                2U &&
            physical_tag != 0U) {
            throw std::invalid_argument(
                "mpmc::mesh::import_gmsh_4_1_ascii: physical line element maps to an interior face");
        }

        found->second.has_line_element = true;
        found->second.global_id = line.tag;
        found->second.physical_tag =
            physical_tag;
    }

    std::uint64_t next_generated_face_id =
        maximum_element_tag;
    for (auto& [key, face] : faces_by_edge) {
        if (face.has_line_element) {
            continue;
        }
        if (next_generated_face_id ==
            std::numeric_limits<std::uint64_t>::max()) {
            throw std::length_error(
                "mpmc::mesh::import_gmsh_4_1_ascii: cannot allocate generated internal-face GlobalEntityId");
        }
        ++next_generated_face_id;
        face.global_id =
            next_generated_face_id;
        face.physical_tag = 0U;
        (void)key;
    }

    std::vector<FaceBuild> faces;
    faces.reserve(faces_by_edge.size());
    for (auto& [key, face] : faces_by_edge) {
        std::sort(
            face.adjacent_cells.begin(),
            face.adjacent_cells.end());
        faces.push_back(std::move(face));
        (void)key;
    }
    std::sort(
        faces.begin(),
        faces.end(),
        [](const FaceBuild& left,
           const FaceBuild& right) {
            return left.global_id <
                   right.global_id;
        });

    std::map<EdgeKey, std::size_t>
        face_local_by_edge;
    entity_ids.faces.reserve(faces.size());
    for (std::size_t face_local = 0U;
         face_local < faces.size();
         ++face_local) {
        if (!face_local_by_edge
                 .emplace(
                     faces[face_local].key,
                     face_local)
                 .second) {
            throw std::invalid_argument(
                "mpmc::mesh::import_gmsh_4_1_ascii: duplicate final face edge");
        }
        entity_ids.faces.emplace_back(
            faces[face_local].global_id);
    }
    entity_ids.cells.reserve(cells.size());
    for (const auto& cell : cells) {
        entity_ids.cells.emplace_back(
            cell.global_id);
    }

    std::vector<CsrAdjacency::Offset>
        cell_vertex_offsets;
    std::vector<LocalIndex>
        cell_vertices;
    std::vector<CsrAdjacency::Offset>
        cell_face_offsets;
    std::vector<LocalIndex>
        cell_faces;
    cell_vertex_offsets.reserve(cells.size() + 1U);
    cell_face_offsets.reserve(cells.size() + 1U);
    cell_vertex_offsets.push_back(
        CsrAdjacency::Offset{0U});
    cell_face_offsets.push_back(
        CsrAdjacency::Offset{0U});

    for (const auto& cell : cells) {
        cell_vertices.insert(
            cell_vertices.end(),
            cell.vertex_locals.begin(),
            cell.vertex_locals.end());
        for (const auto& key :
             cell.edge_keys) {
            const auto found =
                face_local_by_edge.find(key);
            if (found ==
                face_local_by_edge.end()) {
                throw std::logic_error(
                    "mpmc::mesh::import_gmsh_4_1_ascii: internal face lookup failed");
            }
            cell_faces.push_back(
                checked_local(
                    found->second,
                    "mpmc::mesh::import_gmsh_4_1_ascii: face local index overflow"));
        }
        cell_vertex_offsets.push_back(
            checked_offset(
                cell_vertices.size(),
                "mpmc::mesh::import_gmsh_4_1_ascii: cell->vertex CSR overflow"));
        cell_face_offsets.push_back(
            checked_offset(
                cell_faces.size(),
                "mpmc::mesh::import_gmsh_4_1_ascii: cell->face CSR overflow"));
    }

    std::vector<CsrAdjacency::Offset>
        face_vertex_offsets;
    std::vector<LocalIndex>
        face_vertices;
    std::vector<CsrAdjacency::Offset>
        face_cell_offsets;
    std::vector<LocalIndex>
        face_cells;
    face_vertex_offsets.reserve(faces.size() + 1U);
    face_cell_offsets.reserve(faces.size() + 1U);
    face_vertex_offsets.push_back(
        CsrAdjacency::Offset{0U});
    face_cell_offsets.push_back(
        CsrAdjacency::Offset{0U});

    for (const auto& face : faces) {
        for (const auto node_tag :
             {face.key.first,
              face.key.second}) {
            const auto found =
                node_local_by_tag.find(node_tag);
            if (found ==
                node_local_by_tag.end()) {
                throw std::logic_error(
                    "mpmc::mesh::import_gmsh_4_1_ascii: internal face node lookup failed");
            }
            face_vertices.push_back(
                checked_local(
                    found->second,
                    "mpmc::mesh::import_gmsh_4_1_ascii: face vertex local index overflow"));
        }
        for (const auto cell_local :
             face.adjacent_cells) {
            face_cells.push_back(
                checked_local(
                    cell_local,
                    "mpmc::mesh::import_gmsh_4_1_ascii: face cell local index overflow"));
        }
        face_vertex_offsets.push_back(
            checked_offset(
                face_vertices.size(),
                "mpmc::mesh::import_gmsh_4_1_ascii: face->vertex CSR overflow"));
        face_cell_offsets.push_back(
            checked_offset(
                face_cells.size(),
                "mpmc::mesh::import_gmsh_4_1_ascii: face->cell CSR overflow"));
    }

    std::vector<CsrAdjacency> relations;
    relations.reserve(4U);
    relations.emplace_back(
        EntityKind::cell,
        EntityKind::vertex,
        entity_ids.vertices.size(),
        std::move(cell_vertex_offsets),
        std::move(cell_vertices));
    relations.emplace_back(
        EntityKind::cell,
        EntityKind::face,
        entity_ids.faces.size(),
        std::move(cell_face_offsets),
        std::move(cell_faces));
    relations.emplace_back(
        EntityKind::face,
        EntityKind::vertex,
        entity_ids.vertices.size(),
        std::move(face_vertex_offsets),
        std::move(face_vertices));
    relations.emplace_back(
        EntityKind::face,
        EntityKind::cell,
        entity_ids.cells.size(),
        std::move(face_cell_offsets),
        std::move(face_cells));

    Topology topology{
        std::move(entity_ids),
        std::move(relations)};

    std::vector<Coordinate2D> cell_centroids;
    std::vector<double> cell_areas;
    cell_centroids.reserve(cells.size());
    cell_areas.reserve(cells.size());
    for (const auto& cell : cells) {
        cell_centroids.push_back(
            cell.centroid);
        cell_areas.push_back(
            cell.area_m2);
    }

    std::vector<Coordinate2D> face_centroids;
    std::vector<double> face_lengths;
    std::vector<LocalIndex> face_owners;
    std::vector<UnitVector2D>
        face_owner_normals;
    std::vector<PhysicalTag> physical_tags;
    face_centroids.reserve(faces.size());
    face_lengths.reserve(faces.size());
    face_owners.reserve(faces.size());
    face_owner_normals.reserve(faces.size());
    physical_tags.reserve(faces.size());

    for (std::size_t face_local = 0U;
         face_local < faces.size();
         ++face_local) {
        const auto& face = faces[face_local];
        if (face.adjacent_cells.empty() ||
            face.adjacent_cells.size() > 2U) {
            throw std::logic_error(
                "mpmc::mesh::import_gmsh_4_1_ascii: invalid finalized face adjacency");
        }
        const auto first_vertex =
            node_local_by_tag.at(
                face.key.first);
        const auto second_vertex =
            node_local_by_tag.at(
                face.key.second);
        const auto a =
            vertex_coordinates[
                first_vertex];
        const auto b =
            vertex_coordinates[
                second_vertex];
        const double dx = b.x_m - a.x_m;
        const double dy = b.y_m - a.y_m;
        const double length =
            std::hypot(dx, dy);
        if (!std::isfinite(length) ||
            length <= 0.0) {
            throw std::invalid_argument(
                "mpmc::mesh::import_gmsh_4_1_ascii: face length must be finite and positive");
        }
        const Coordinate2D centroid{
            std::midpoint(a.x_m, b.x_m),
            std::midpoint(a.y_m, b.y_m)};

        const std::size_t owner =
            face.adjacent_cells.front();
        const auto owner_centroid =
            cells[owner].centroid;
        UnitVector2D normal{
            -dy / length,
            dx / length};
        const double outward_dot =
            normal.x *
                (centroid.x_m -
                 owner_centroid.x_m) +
            normal.y *
                (centroid.y_m -
                 owner_centroid.y_m);
        const double direction_tolerance =
            256.0 *
            std::numeric_limits<double>::epsilon() *
            std::max(
                1.0,
                length);
        if (std::abs(outward_dot) <=
            direction_tolerance) {
            throw std::invalid_argument(
                "mpmc::mesh::import_gmsh_4_1_ascii: cannot orient face normal from owner centroid");
        }
        if (outward_dot < 0.0) {
            normal.x = -normal.x;
            normal.y = -normal.y;
        }

        face_centroids.push_back(centroid);
        face_lengths.push_back(length);
        face_owners.push_back(
            checked_local(
                owner,
                "mpmc::mesh::import_gmsh_4_1_ascii: face owner local index overflow"));
        face_owner_normals.push_back(normal);
        physical_tags.emplace_back(
            face.physical_tag);
    }

    Geometry2D geometry{
        std::move(vertex_coordinates),
        std::move(cell_centroids),
        std::move(cell_areas),
        std::move(face_centroids),
        std::move(face_lengths),
        std::move(face_owners),
        std::move(face_owner_normals)};

    FaceBoundarySnapshot face_boundary =
        make_face_boundary_snapshot(
            topology, physical_tags);

    std::vector<GmshPhysicalName>
        physical_names;
    physical_names.reserve(
        physical_name_map.size());
    for (const auto& [key, name] :
         physical_name_map) {
        physical_names.push_back(
            GmshPhysicalName{
                key.dimension,
                checked_physical_tag(
                    key.tag,
                    "mpmc::mesh::import_gmsh_4_1_ascii: physical-name tag overflow"),
                name});
    }

    std::vector<GmshCellPhysicalGroups>
        cell_physical_groups;
    cell_physical_groups.reserve(
        cells.size());
    for (const auto& cell : cells) {
        cell_physical_groups.push_back(
            GmshCellPhysicalGroups{
                GlobalEntityId{
                    cell.global_id},
                entity_physical_tags(
                    entities,
                    entities_seen,
                    2,
                    cell.entity_tag)});
    }

    return Gmsh41ImportResult{
        std::move(topology),
        std::move(geometry),
        std::move(face_boundary),
        std::move(physical_names),
        std::move(cell_physical_groups)};
}

} // namespace mpmc::mesh

#endif // MPMC_MESH_GMSH_4_1_HPP
