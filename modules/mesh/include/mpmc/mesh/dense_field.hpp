#ifndef MPMC_MESH_DENSE_FIELD_HPP
#define MPMC_MESH_DENSE_FIELD_HPP

#include <mpmc/mesh/entity.hpp>
#include <mpmc/mesh/topology.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mpmc::mesh {

enum class FieldSourceKind : std::uint8_t {
    unspecified = 0,
    file_import = 1,
    user_supplied = 2,
    generated = 3,
    synthetic_test = 4,
};

struct FieldSourceMetadata {
    FieldSourceKind kind = FieldSourceKind::unspecified;
    std::string reference;
    std::string revision;
    std::string locator;
};

struct DenseFieldMetadata {
    std::string id;
    std::string unit;
    FieldSourceMetadata source;
};

/// Immutable dense numeric field aligned to one Topology entity kind.
///
/// Storage is entity-major and component-interleaved:
/// values[entity * component_count + component].
///
/// This generic layer assigns no physical meaning to the field ID, components or
/// unit string and performs no unit conversion. Every stored numeric value must
/// be finite; missing-value semantics require a future explicit contract.
class DenseFieldSnapshot {
public:
    DenseFieldSnapshot(const DenseFieldSnapshot&) = default;
    DenseFieldSnapshot(DenseFieldSnapshot&&) noexcept = default;
    DenseFieldSnapshot& operator=(const DenseFieldSnapshot&) = delete;
    DenseFieldSnapshot& operator=(DenseFieldSnapshot&&) = delete;
    ~DenseFieldSnapshot() = default;

    [[nodiscard]] static DenseFieldSnapshot create(
        const Topology& topology,
        EntityKind location,
        std::size_t component_count,
        std::vector<double> values,
        DenseFieldMetadata metadata) {
        validate_location(location);
        validate_metadata(metadata);

        if (component_count == 0U) {
            throw std::invalid_argument(
                "mpmc::mesh::DenseFieldSnapshot: component_count must be positive");
        }

        const std::size_t entity_count = topology.entity_count(location);
        if (entity_count != 0U &&
            component_count > std::numeric_limits<std::size_t>::max() / entity_count) {
            throw std::length_error(
                "mpmc::mesh::DenseFieldSnapshot: dense field size overflow");
        }
        const std::size_t expected_values = entity_count * component_count;
        if (values.size() != expected_values) {
            throw std::invalid_argument(
                "mpmc::mesh::DenseFieldSnapshot: value count does not match topology and components");
        }
        if (!std::all_of(values.begin(), values.end(),
                         [](double value) { return std::isfinite(value); })) {
            throw std::invalid_argument(
                "mpmc::mesh::DenseFieldSnapshot: values must be finite");
        }

        return DenseFieldSnapshot{
            location, entity_count, component_count, std::move(values), std::move(metadata)};
    }

    [[nodiscard]] EntityKind location() const noexcept { return location_; }
    [[nodiscard]] std::size_t entity_count() const noexcept { return entity_count_; }
    [[nodiscard]] std::size_t component_count() const noexcept { return component_count_; }
    [[nodiscard]] std::size_t value_count() const noexcept { return values_.size(); }

    [[nodiscard]] const DenseFieldMetadata& metadata() const & noexcept { return metadata_; }
    const DenseFieldMetadata& metadata() const && = delete;

    [[nodiscard]] std::span<const double> values() const noexcept { return values_; }

    [[nodiscard]] std::span<const double> entity_values(LocalIndex entity) const {
        const std::size_t index = static_cast<std::size_t>(entity.value());
        if (index >= entity_count_) {
            throw std::out_of_range(
                "mpmc::mesh::DenseFieldSnapshot: entity index out of range");
        }
        return std::span<const double>{values_}.subspan(
            index * component_count_, component_count_);
    }

    [[nodiscard]] double value(LocalIndex entity, std::size_t component) const {
        if (component >= component_count_) {
            throw std::out_of_range(
                "mpmc::mesh::DenseFieldSnapshot: component index out of range");
        }
        const auto components = entity_values(entity);
        return components[component];
    }

private:
    DenseFieldSnapshot(EntityKind location,
                       std::size_t entity_count,
                       std::size_t component_count,
                       std::vector<double> values,
                       DenseFieldMetadata metadata)
        : location_(location),
          entity_count_(entity_count),
          component_count_(component_count),
          values_(std::move(values)),
          metadata_(std::move(metadata)) {}

    static void validate_location(EntityKind location) {
        switch (location) {
        case EntityKind::vertex:
        case EntityKind::face:
        case EntityKind::cell:
            return;
        case EntityKind::edge:
            throw std::invalid_argument(
                "mpmc::mesh::DenseFieldSnapshot: edge fields are not supported by this contract");
        }
        throw std::invalid_argument(
            "mpmc::mesh::DenseFieldSnapshot: invalid field location");
    }

    static void require_text(const std::string& text, const char* message) {
        const bool nonblank =
            std::any_of(text.begin(), text.end(),
                        [](unsigned char character) { return character > 32U; });
        if (!nonblank || text.find('\0') != std::string::npos) {
            throw std::invalid_argument(message);
        }
    }

    static void validate_metadata(const DenseFieldMetadata& metadata) {
        require_text(metadata.id,
                     "mpmc::mesh::DenseFieldSnapshot: nonblank field id is required");
        require_text(metadata.unit,
                     "mpmc::mesh::DenseFieldSnapshot: explicit nonblank unit is required");
        switch (metadata.source.kind) {
        case FieldSourceKind::file_import:
        case FieldSourceKind::user_supplied:
        case FieldSourceKind::generated:
        case FieldSourceKind::synthetic_test:
            break;
        default:
            throw std::invalid_argument(
                "mpmc::mesh::DenseFieldSnapshot: invalid source kind");
        }
        require_text(
            metadata.source.reference,
            "mpmc::mesh::DenseFieldSnapshot: nonblank source reference is required");
        require_text(
            metadata.source.revision,
            "mpmc::mesh::DenseFieldSnapshot: nonblank source revision is required");
        // Locator may be empty when the source has no finer-grained record address.
        if (metadata.source.locator.find('\0') != std::string::npos) {
            throw std::invalid_argument(
                "mpmc::mesh::DenseFieldSnapshot: source locator cannot contain NUL");
        }
    }

    EntityKind location_;
    std::size_t entity_count_;
    std::size_t component_count_;
    std::vector<double> values_;
    DenseFieldMetadata metadata_;
};

} // namespace mpmc::mesh

#endif // MPMC_MESH_DENSE_FIELD_HPP
