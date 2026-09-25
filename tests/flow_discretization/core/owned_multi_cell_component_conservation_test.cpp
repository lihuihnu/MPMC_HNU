#include <mpmc/flow_discretization/owned_multi_cell_component_conservation.hpp>
#include <mpmc/mesh/topology.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <source_location>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

bool owned_multi_cell_component_conservation_header();

namespace {

namespace fd = mpmc::flow_discretization;
namespace flow = mpmc::flow;
namespace mesh = mpmc::mesh;

void require(
    bool condition,
    std::string_view message,
    std::source_location where =
        std::source_location::current()) {
    if (!condition) {
        throw std::runtime_error(
            std::string{where.file_name()} +
            ":" +
            std::to_string(where.line()) +
            ": " +
            std::string{message});
    }
}

void near(
    double actual,
    double expected,
    double relative = 5.0e-13,
    double absolute = 5.0e-13,
    std::source_location where =
        std::source_location::current()) {
    if (!std::isfinite(actual) ||
        !std::isfinite(expected) ||
        std::abs(actual - expected) >
            absolute +
                relative *
                    std::max(
                        std::abs(actual),
                        std::abs(expected))) {
        std::cerr
            << "actual=" << actual
            << " expected=" << expected
            << '\n';
        require(false, "numeric mismatch", where);
    }
}

template <class Function>
void expect_invalid(
    Function&& function,
    std::string_view fragment) {
    try {
        function();
    } catch (const std::invalid_argument& error) {
        require(
            std::string_view{error.what()}.find(fragment) !=
                std::string_view::npos,
            "invalid_argument diagnostic changed");
        return;
    }
    throw std::runtime_error(
        "expected std::invalid_argument");
}

flow::NaturalVariableStateIdentity3P
identity(std::size_t cell) {
    const std::array<
        std::array<std::size_t, 3>,
        3>
        pivots{{
            {1U, 0U, 2U},
            {0U, 2U, 1U},
            {2U, 1U, 0U}}};
    const std::array<double, 3>
        pressures{1.00e6, 0.99e6, 1.02e6};
    const std::array<double, 3>
        temperatures{350.0, 351.0, 349.0};
    const std::array<
        std::array<double, 3>,
        3>
        saturations{{
            {0.20, 0.30, 0.50},
            {0.25, 0.35, 0.40},
            {0.30, 0.20, 0.50}}};
    const std::array<
        std::array<std::vector<double>, 3>,
        3>
        compositions{{
            {
                std::vector<double>{0.10, 0.70, 0.20},
                std::vector<double>{0.60, 0.20, 0.20},
                std::vector<double>{0.20, 0.30, 0.50}},
            {
                std::vector<double>{0.55, 0.25, 0.20},
                std::vector<double>{0.25, 0.25, 0.50},
                std::vector<double>{0.20, 0.60, 0.20}},
            {
                std::vector<double>{0.20, 0.20, 0.60},
                std::vector<double>{0.30, 0.50, 0.20},
                std::vector<double>{0.65, 0.15, 0.20}}
        }};

    if (cell >= pivots.size()) {
        throw std::invalid_argument(
            "cell fixture index out of range");
    }

    return {
        flow::NaturalVariableLayout3P{
            flow::NaturalVariableCompositionPivot3P::
                from_dependent_components(
                    3U,
                    pivots[cell])},
        {"A", "B", "C"},
        pressures[cell],
        temperatures[cell],
        saturations[cell],
        compositions[cell]};
}

std::vector<double> component_jacobian(
    std::size_t component_count,
    std::size_t column_count,
    double scale) {
    std::vector<double> values(
        component_count * column_count,
        0.0);
    for (std::size_t component = 0U;
         component < component_count;
         ++component) {
        for (std::size_t column = 0U;
             column < column_count;
             ++column) {
            values[
                component * column_count +
                column] =
                scale *
                static_cast<double>(
                    1U +
                    component * 10U +
                    column);
        }
    }
    return values;
}

std::vector<double> total_gradient(
    const std::vector<double>& values,
    std::size_t component_count,
    std::size_t column_count) {
    std::vector<double> total(
        column_count,
        0.0);
    for (std::size_t column = 0U;
         column < column_count;
         ++column) {
        for (std::size_t component = 0U;
             component < component_count;
             ++component) {
            total[column] +=
                values[
                    component * column_count +
                    column];
        }
    }
    return total;
}

flow::BackwardEulerComponentAccumulationResidual3P
accumulation(
    std::size_t cell,
    double scale) {
    auto state = identity(cell);
    const std::size_t q =
        state.layout.unknown_count();

    std::vector<double> residual{
        0.10 * scale,
        -0.04 * scale,
        0.06 * scale};
    double total = 0.0;
    for (double value : residual) {
        total += value;
    }

    const auto jacobian =
        component_jacobian(
            3U,
            q,
            0.001 * scale);

    return {
        state.layout,
        0.20 + 0.01 * static_cast<double>(cell),
        120.0,
        {"A", "B", "C"},
        std::move(residual),
        total,
        q,
        jacobian,
        total_gradient(
            jacobian,
            3U,
            q)};
}

fd::ConservativeComponentFaceRateScatterLinearization3D
scatter(
    std::size_t face,
    std::size_t owner_cell,
    std::size_t neighbour_cell,
    std::vector<double> owner_rate,
    double owner_scale,
    double neighbour_scale) {
    auto owner_identity =
        identity(owner_cell);
    auto neighbour_identity =
        identity(neighbour_cell);
    const std::size_t n =
        owner_rate.size();
    const std::size_t owner_q =
        owner_identity.layout.unknown_count();
    const std::size_t neighbour_q =
        neighbour_identity.layout.unknown_count();

    std::vector<double> neighbour_rate =
        owner_rate;
    for (double& value : neighbour_rate) {
        value = -value;
    }

    const auto owner_owner =
        component_jacobian(
            n,
            owner_q,
            owner_scale);
    const auto owner_neighbour =
        component_jacobian(
            n,
            neighbour_q,
            neighbour_scale);

    std::vector<double> neighbour_owner(
        owner_owner.size());
    std::vector<double> neighbour_neighbour(
        owner_neighbour.size());
    for (std::size_t index = 0U;
         index < owner_owner.size();
         ++index) {
        neighbour_owner[index] =
            -owner_owner[index];
    }
    for (std::size_t index = 0U;
         index < owner_neighbour.size();
         ++index) {
        neighbour_neighbour[index] =
            -owner_neighbour[index];
    }

    double owner_total = 0.0;
    for (double value : owner_rate) {
        owner_total += value;
    }

    const auto owner_total_owner =
        total_gradient(
            owner_owner,
            n,
            owner_q);
    const auto owner_total_neighbour =
        total_gradient(
            owner_neighbour,
            n,
            neighbour_q);

    std::vector<double> neighbour_total_owner(
        owner_q);
    std::vector<double> neighbour_total_neighbour(
        neighbour_q);
    for (std::size_t column = 0U;
         column < owner_q;
         ++column) {
        neighbour_total_owner[column] =
            -owner_total_owner[column];
    }
    for (std::size_t column = 0U;
         column < neighbour_q;
         ++column) {
        neighbour_total_neighbour[column] =
            -owner_total_neighbour[column];
    }

    return {
        mesh::LocalIndex{
            static_cast<
                mesh::LocalIndex::value_type>(
                    face)},
        {"A", "B", "C"},
        std::move(owner_identity),
        std::move(neighbour_identity),
        std::move(owner_rate),
        std::move(neighbour_rate),
        owner_owner,
        owner_neighbour,
        std::move(neighbour_owner),
        std::move(neighbour_neighbour),
        owner_total,
        -owner_total,
        owner_total_owner,
        owner_total_neighbour,
        std::move(neighbour_total_owner),
        std::move(neighbour_total_neighbour)};
}

struct SyntheticConnectionRow {
    mesh::LocalIndex face{
        mesh::LocalIndex::value_type{0}};
    mesh::GlobalEntityId face_global{
        mesh::GlobalEntityId::value_type{0}};
    mesh::LocalIndex owner_cell{
        mesh::LocalIndex::value_type{0}};
    mesh::GlobalEntityId owner_cell_global{
        mesh::GlobalEntityId::value_type{0}};
    mesh::LocalIndex neighbour_cell{
        mesh::LocalIndex::value_type{0}};
    mesh::GlobalEntityId neighbour_cell_global{
        mesh::GlobalEntityId::value_type{0}};
    double transmissibility_m3{};
};

class SyntheticSchedule {
public:
    SyntheticSchedule(
        mesh::PartitionRank local_rank,
        std::uint32_t rank_count,
        std::vector<SyntheticConnectionRow> rows)
        : local_rank_(local_rank),
          rank_count_(rank_count),
          rows_(std::move(rows)) {}

    [[nodiscard]] mesh::PartitionRank
    local_rank() const noexcept {
        return local_rank_;
    }

    [[nodiscard]] std::uint32_t
    rank_count() const noexcept {
        return rank_count_;
    }

    [[nodiscard]] std::span<
        const SyntheticConnectionRow>
    assembly_rows() const noexcept {
        return rows_;
    }

private:
    mesh::PartitionRank local_rank_;
    std::uint32_t rank_count_;
    std::vector<SyntheticConnectionRow> rows_;
};

mesh::Topology make_topology() {
    mesh::Topology::EntityIds ids;
    ids.faces = {
        mesh::GlobalEntityId{0U},
        mesh::GlobalEntityId{1U},
        mesh::GlobalEntityId{2U},
        mesh::GlobalEntityId{3U}};
    ids.cells = {
        mesh::GlobalEntityId{0U},
        mesh::GlobalEntityId{1U},
        mesh::GlobalEntityId{2U}};
    return mesh::Topology{
        std::move(ids),
        {}};
}

struct Fixture {
    mesh::Topology topology;
    mesh::PartitionSnapshot partition;
    SyntheticSchedule schedule;
    std::array<
        flow::BackwardEulerComponentAccumulationResidual3P,
        3>
        accumulations;
    std::array<
        fd::NormalizedComponentFaceContributionLinearization3D,
        2>
        faces;
    std::array<fd::OwnedCellAccumulationBinding3D, 3>
        cell_bindings;
    std::array<
        fd::AuthoritativeNormalizedFaceContributionBinding3D,
        2>
        face_bindings;

    Fixture()
        : topology(make_topology()),
          partition(
              mesh::make_serial_partition_snapshot(
                  topology)),
          schedule(
              mesh::PartitionRank{0U},
              1U,
              {
                  {
                      mesh::LocalIndex{1U},
                      mesh::GlobalEntityId{1U},
                      mesh::LocalIndex{0U},
                      mesh::GlobalEntityId{0U},
                      mesh::LocalIndex{1U},
                      mesh::GlobalEntityId{1U},
                      2.0},
                  {
                      mesh::LocalIndex{2U},
                      mesh::GlobalEntityId{2U},
                      mesh::LocalIndex{1U},
                      mesh::GlobalEntityId{1U},
                      mesh::LocalIndex{2U},
                      mesh::GlobalEntityId{2U},
                      3.0}
              }),
          accumulations{
              accumulation(0U, 1.0),
              accumulation(1U, 1.5),
              accumulation(2U, 2.0)},
          faces{
              fd::normalize_component_face_rate_by_bulk_volume(
                  scatter(
                      1U,
                      0U,
                      1U,
                      {2.0, -1.0, 3.0},
                      0.01,
                      -0.02),
                  {2.0, 3.0}),
              fd::normalize_component_face_rate_by_bulk_volume(
                  scatter(
                      2U,
                      1U,
                      2U,
                      {-4.0, 2.0, 1.0},
                      0.03,
                      -0.04),
                  {3.0, 5.0})},
          cell_bindings{
              fd::OwnedCellAccumulationBinding3D{
                  mesh::LocalIndex{0U},
                  mesh::GlobalEntityId{0U},
                  2.0,
                  identity(0U),
                  &accumulations[0]},
              fd::OwnedCellAccumulationBinding3D{
                  mesh::LocalIndex{1U},
                  mesh::GlobalEntityId{1U},
                  3.0,
                  identity(1U),
                  &accumulations[1]},
              fd::OwnedCellAccumulationBinding3D{
                  mesh::LocalIndex{2U},
                  mesh::GlobalEntityId{2U},
                  5.0,
                  identity(2U),
                  &accumulations[2]}},
          face_bindings{
              fd::AuthoritativeNormalizedFaceContributionBinding3D{
                  mesh::LocalIndex{1U},
                  mesh::GlobalEntityId{1U},
                  &faces[0]},
              fd::AuthoritativeNormalizedFaceContributionBinding3D{
                  mesh::LocalIndex{2U},
                  mesh::GlobalEntityId{2U},
                  &faces[1]}} {}
};

fd::SerialOwnedMultiCellComponentConservationSnapshot3D
build(Fixture& fixture) {
    return fd::
        build_serial_owned_multi_cell_component_conservation_snapshot(
            fixture.schedule,
            fixture.partition,
            fixture.cell_bindings,
            fixture.face_bindings);
}

void serial_patch_rows_and_conservation() {
    Fixture fixture;
    const auto result = build(fixture);

    require(
        result.local_cell_count() == 3U &&
            result.authoritative_face_count() == 2U &&
            result.component_ids().size() == 3U &&
            result.rows().size() == 3U,
        "serial patch snapshot cardinality mismatch");

    for (double value :
         result.volume_weighted_spatial_component_balance_mol_per_s()) {
        near(value, 0.0, 0.0, 2.0e-15);
    }
    near(
        result.volume_weighted_spatial_total_balance_mol_per_s(),
        0.0,
        0.0,
        2.0e-15);

    const auto& row0 =
        result.row(mesh::LocalIndex{0U});
    const auto& row1 =
        result.row(mesh::LocalIndex{1U});
    const auto& row2 =
        result.row(mesh::LocalIndex{2U});

    require(
        row0.cell_global == mesh::GlobalEntityId{0U} &&
            row1.cell_global == mesh::GlobalEntityId{1U} &&
            row2.cell_global == mesh::GlobalEntityId{2U},
        "stable cell row identity changed");

    for (std::size_t component = 0U;
         component < 3U;
         ++component) {
        near(
            row0.local_residual.residual(component),
            fixture.accumulations[0].residual(component) +
                fixture.faces[0].owner_contribution(component));
        near(
            row1.local_residual.residual(component),
            fixture.accumulations[1].residual(component) +
                fixture.faces[0].neighbour_contribution(component) +
                fixture.faces[1].owner_contribution(component));
        near(
            row2.local_residual.residual(component),
            fixture.accumulations[2].residual(component) +
                fixture.faces[1].neighbour_contribution(component));
    }

    require(
        row0.off_diagonal_cell_pair_blocks.size() == 1U &&
            row1.off_diagonal_cell_pair_blocks.size() == 2U &&
            row2.off_diagonal_cell_pair_blocks.size() == 1U,
        "cell-pair block cardinality mismatch");

    require(
        row0.off_diagonal_cell_pair_blocks[0]
                .column_cell_global ==
            mesh::GlobalEntityId{1U} &&
        row1.off_diagonal_cell_pair_blocks[0]
                .column_cell_global ==
            mesh::GlobalEntityId{0U} &&
        row1.off_diagonal_cell_pair_blocks[1]
                .column_cell_global ==
            mesh::GlobalEntityId{2U} &&
        row2.off_diagonal_cell_pair_blocks[0]
                .column_cell_global ==
            mesh::GlobalEntityId{1U},
        "off-diagonal blocks are not keyed by stable neighbour cell");

    require(
        row1.off_diagonal_cell_pair_blocks[0]
                .contributing_face_global_ids ==
            std::vector<mesh::GlobalEntityId>{
                mesh::GlobalEntityId{1U}} &&
        row1.off_diagonal_cell_pair_blocks[1]
                .contributing_face_global_ids ==
            std::vector<mesh::GlobalEntityId>{
                mesh::GlobalEntityId{2U}},
        "cell-pair block lost stable face provenance");
}

void jacobian_mapping() {
    Fixture fixture;
    const auto result = build(fixture);

    const auto& row1 =
        result.row(mesh::LocalIndex{1U});
    const std::size_t q1 =
        fixture.cell_bindings[1]
            .state_identity.layout.unknown_count();

    for (std::size_t component = 0U;
         component < 3U;
         ++component) {
        for (std::size_t column = 0U;
             column < q1;
             ++column) {
            near(
                row1.local_residual.d_local(
                    component,
                    column),
                fixture.accumulations[1]
                        .d_residual(
                            component,
                            column) +
                    fixture.faces[0]
                        .d_neighbour_contribution_wrt_neighbour(
                            component,
                            column) +
                    fixture.faces[1]
                        .d_owner_contribution_wrt_owner(
                            component,
                            column));
        }
    }

    const auto& to_cell0 =
        row1.off_diagonal_cell_pair_blocks[0];
    const std::size_t q0 =
        fixture.cell_bindings[0]
            .state_identity.layout.unknown_count();
    for (std::size_t component = 0U;
         component < 3U;
         ++component) {
        for (std::size_t column = 0U;
             column < q0;
             ++column) {
            near(
                to_cell0.d_component(
                    component,
                    column),
                fixture.faces[0]
                    .d_neighbour_contribution_wrt_owner(
                        component,
                        column));
        }
    }

    const auto& to_cell2 =
        row1.off_diagonal_cell_pair_blocks[1];
    const std::size_t q2 =
        fixture.cell_bindings[2]
            .state_identity.layout.unknown_count();
    for (std::size_t component = 0U;
         component < 3U;
         ++component) {
        for (std::size_t column = 0U;
             column < q2;
             ++column) {
            near(
                to_cell2.d_component(
                    component,
                    column),
                fixture.faces[1]
                    .d_owner_contribution_wrt_neighbour(
                        component,
                        column));
        }
    }
}

void invalid_inputs() {
    {
        Fixture fixture;
        auto rows =
            std::vector<SyntheticConnectionRow>{
                fixture.schedule.assembly_rows().begin(),
                fixture.schedule.assembly_rows().end()};
        rows.push_back(rows.front());
        SyntheticSchedule duplicate{
            mesh::PartitionRank{0U},
            1U,
            std::move(rows)};
        expect_invalid(
            [&] {
                (void)fd::
                    build_serial_owned_multi_cell_component_conservation_snapshot(
                        duplicate,
                        fixture.partition,
                        fixture.cell_bindings,
                        fixture.face_bindings);
            },
            "duplicate authoritative face");
    }

    {
        Fixture fixture;
        const std::array<
            fd::AuthoritativeNormalizedFaceContributionBinding3D,
            1>
            missing{{fixture.face_bindings[0]}};
        expect_invalid(
            [&] {
                (void)fd::
                    build_serial_owned_multi_cell_component_conservation_snapshot(
                        fixture.schedule,
                        fixture.partition,
                        fixture.cell_bindings,
                        missing);
            },
            "binding counts");
    }

    {
        Fixture fixture;
        auto wrong_cells =
            fixture.cell_bindings;
        wrong_cells[1].cell_global =
            mesh::GlobalEntityId{99U};
        expect_invalid(
            [&] {
                (void)fd::
                    build_serial_owned_multi_cell_component_conservation_snapshot(
                        fixture.schedule,
                        fixture.partition,
                        wrong_cells,
                        fixture.face_bindings);
            },
            "invalid or duplicate owned-cell");
    }

    {
        Fixture fixture;
        auto wrong_faces =
            fixture.face_bindings;
        wrong_faces[0].face_global =
            mesh::GlobalEntityId{99U};
        expect_invalid(
            [&] {
                (void)fd::
                    build_serial_owned_multi_cell_component_conservation_snapshot(
                        fixture.schedule,
                        fixture.partition,
                        fixture.cell_bindings,
                        wrong_faces);
            },
            "invalid or duplicate authoritative");
    }

    {
        Fixture fixture;
        auto wrong_cells =
            fixture.cell_bindings;
        wrong_cells[1].bulk_volume_m3 = 4.0;
        expect_invalid(
            [&] {
                (void)fd::
                    build_serial_owned_multi_cell_component_conservation_snapshot(
                        fixture.schedule,
                        fixture.partition,
                        wrong_cells,
                        fixture.face_bindings);
            },
            "bulk volume");
    }

    {
        Fixture fixture;
        mesh::EntityOwnerRanks owners;
        owners.vertices.assign(
            fixture.topology.entity_count(
                mesh::EntityKind::vertex),
            mesh::PartitionRank{0U});
        owners.edges.assign(
            fixture.topology.entity_count(
                mesh::EntityKind::edge),
            mesh::PartitionRank{0U});
        owners.faces.assign(
            fixture.topology.entity_count(
                mesh::EntityKind::face),
            mesh::PartitionRank{0U});
        owners.cells.assign(
            fixture.topology.entity_count(
                mesh::EntityKind::cell),
            mesh::PartitionRank{0U});
        auto partitioned =
            mesh::PartitionSnapshot::create(
                fixture.topology,
                mesh::PartitionRank{0U},
                2U,
                std::move(owners));
        SyntheticSchedule schedule{
            mesh::PartitionRank{0U},
            2U,
            std::vector<SyntheticConnectionRow>{
                fixture.schedule.assembly_rows().begin(),
                fixture.schedule.assembly_rows().end()}};
        expect_invalid(
            [&] {
                (void)fd::
                    build_serial_owned_multi_cell_component_conservation_snapshot(
                        schedule,
                        partitioned,
                        fixture.cell_bindings,
                        fixture.face_bindings);
            },
            "single-rank");
    }
}

void headers() {
    require(
        owned_multi_cell_component_conservation_header(),
        "owned multi-cell conservation header probe failed");
}

using Test =
    std::pair<std::string_view, void (*)()>;

constexpr Test tests[]{
    {"serial_patch_rows_and_conservation",
     serial_patch_rows_and_conservation},
    {"jacobian_mapping", jacobian_mapping},
    {"invalid_inputs", invalid_inputs},
    {"headers", headers}};

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            throw std::invalid_argument(
                "one test name required");
        }
        for (const auto& [name, run] : tests) {
            if (name == argv[1]) {
                run();
                std::cout
                    << "[PASS] "
                    << name
                    << '\n';
                return 0;
            }
        }
        throw std::invalid_argument(
            "unknown test");
    } catch (const std::exception& error) {
        std::cerr
            << "[FAIL] "
            << error.what()
            << '\n';
        return 1;
    }
}
