#include "cpa_thermopack_parameter_snapshot.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string_view>

namespace {
namespace th = mpmc::thermodynamics;
namespace snapshot = cpa_thermopack_snapshot;

void require(bool value, const char* message) {
    if (!value) { throw std::runtime_error(message); }
}

bool same_provenance(const th::Provenance& first, const th::Provenance& second) {
    return first.kind == second.kind &&
           first.reference == second.reference &&
           first.revision == second.revision &&
           first.locator == second.locator &&
           first.note == second.note &&
           first.acquisition == second.acquisition &&
           first.usage_terms == second.usage_terms;
}

void require_same_sourced_value(
    const th::CpaSourcedValue& first,
    const th::CpaSourcedValue& second,
    const char* message) {
    require(first.value == second.value &&
                same_provenance(first.source, second.source),
            message);
}

void require_same_sites(
    const std::vector<th::CpaAssociationSiteInput>& first,
    const std::vector<th::CpaAssociationSiteInput>& second) {
    require(first.size() == second.size(), "CPA snapshot site count drifted");
    for (std::size_t i = 0; i < first.size(); ++i) {
        require(first[i].id == second[i].id &&
                    first[i].multiplicity == second[i].multiplicity,
                "CPA snapshot site topology drifted");
    }
}

void audit_snapshot_identity(bool swapped) {
    const auto literature = cpa_physical_test::parameters(swapped);
    const auto parity = snapshot::parameters(swapped);

    require(literature.dataset_id() ==
                "literature-cpa-water-methanol-cr1-333.15K",
            "literature CPA dataset identity drifted");
    require(literature.revision() ==
                "Kontogeorgis-2008-pure__Folas-CR1-kij__Kurihara-1995-VLE",
            "literature CPA revision identity drifted");
    require(parity.dataset_id() == snapshot::parity_dataset_id,
            "ThermoPack parity CPA dataset identity drifted");
    require(parity.revision() == snapshot::parity_revision,
            "ThermoPack parity CPA revision identity drifted");
    require(parity.dataset_id() != literature.dataset_id(),
            "literature and ThermoPack parity snapshots must remain distinct datasets");
    require(parity.revision() != literature.revision(),
            "literature and ThermoPack parity snapshots must remain distinct revisions");

    require(parity.components().size() == literature.components().size(),
            "CPA snapshot component count drifted");
    for (std::size_t i = 0; i < literature.components().size(); ++i) {
        require(parity.components().at(i).id == literature.components().at(i).id,
                "CPA snapshot component ordering drifted");
    }

    const auto literature_pure = literature.pure_records();
    const auto parity_pure = parity.pure_records();
    require(literature_pure.size() == 2U && parity_pure.size() == 2U,
            "methanol-water CPA snapshot must contain exactly two pure records");

    for (std::size_t i = 0; i < literature_pure.size(); ++i) {
        const auto& baseline = literature_pure[i];
        const auto& aligned = parity_pure[i];
        require(baseline.component_id == aligned.component_id,
                "CPA pure-record component identity drifted");
        require_same_sourced_value(
            baseline.a0_pa_m6_per_mol2, aligned.a0_pa_m6_per_mol2,
            "CPA parity snapshot must preserve literature a0 and provenance");
        require_same_sourced_value(
            baseline.b_m3_per_mol, aligned.b_m3_per_mol,
            "CPA parity snapshot must preserve literature b and provenance");
        require_same_sourced_value(
            baseline.c1_dimensionless, aligned.c1_dimensionless,
            "CPA parity snapshot must preserve literature c1 and provenance");
        require_same_sites(baseline.sites, aligned.sites);

        require(baseline.critical_temperature_k.source.kind ==
                    th::SourceKind::literature,
                "literature CPA Tc must retain literature provenance");
        require(aligned.critical_temperature_k.source.kind ==
                    th::SourceKind::database,
                "ThermoPack parity Tc must carry database provenance");
        require(aligned.critical_temperature_k.source.reference ==
                    "thermotools/thermopack pinned component database",
                "ThermoPack parity Tc source identity drifted");
        require(aligned.critical_temperature_k.source.revision ==
                    snapshot::thermopack_commit,
                "ThermoPack parity Tc source revision drifted");
        require(!same_provenance(
                    baseline.critical_temperature_k.source,
                    aligned.critical_temperature_k.source),
                "literature and ThermoPack Tc provenance must not collapse");

        if (aligned.component_id == "METHANOL") {
            require(baseline.critical_temperature_k.value == 512.64,
                    "literature methanol Tc drifted");
            require(aligned.critical_temperature_k.value ==
                        cpa_thermopack_phase_kernel::methanol_alpha_critical_temperature_k &&
                        aligned.critical_temperature_k.value == 512.6,
                    "ThermoPack parity methanol Tc drifted");
            require(aligned.critical_temperature_k.source.locator ==
                        "fluids/Methanol.json: critical.temperature = 512.6 K",
                    "ThermoPack methanol Tc locator drifted");
        } else if (aligned.component_id == "WATER") {
            require(baseline.critical_temperature_k.value == 647.29,
                    "literature water Tc drifted");
            require(aligned.critical_temperature_k.value ==
                        cpa_thermopack_phase_kernel::water_alpha_critical_temperature_k &&
                        aligned.critical_temperature_k.value == 647.3,
                    "ThermoPack parity water Tc drifted");
            require(aligned.critical_temperature_k.source.locator ==
                        "fluids/Water.json: critical.temperature = 647.3 K",
                    "ThermoPack water Tc locator drifted");
        } else {
            throw std::runtime_error("unexpected component in CPA snapshot identity test");
        }
    }

    const auto literature_binary = literature.binary_records();
    const auto parity_binary = parity.binary_records();
    require(literature_binary.size() == 1U && parity_binary.size() == 1U,
            "methanol-water CPA snapshot must contain exactly one binary record");
    require(literature_binary[0].first_component_id ==
                parity_binary[0].first_component_id &&
                literature_binary[0].second_component_id ==
                parity_binary[0].second_component_id,
            "CPA binary component identity drifted");
    require_same_sourced_value(
        literature_binary[0].kij_dimensionless,
        parity_binary[0].kij_dimensionless,
        "ThermoPack parity snapshot must preserve Folas kij and provenance");
    require(parity_binary[0].kij_dimensionless.value == -0.055,
            "ThermoPack parity snapshot must not silently use ThermoPack default kij=-0.09");
    require(parity_binary[0].kij_dimensionless.source.kind ==
                th::SourceKind::literature,
            "ThermoPack parity kij must retain Folas literature provenance");
    require(parity_binary[0].kij_dimensionless.source.reference.find("Folas") !=
                std::string::npos,
            "ThermoPack parity kij provenance lost the Folas source identity");

    const auto literature_assoc = literature.association_records();
    const auto parity_assoc = parity.association_records();
    require(literature_assoc.size() == parity_assoc.size(),
            "CPA association-record count drifted between snapshots");
    for (std::size_t i = 0; i < literature_assoc.size(); ++i) {
        const auto& baseline = literature_assoc[i];
        const auto& aligned = parity_assoc[i];
        require(baseline.first_component_id == aligned.first_component_id &&
                    baseline.first_site_id == aligned.first_site_id &&
                    baseline.second_component_id == aligned.second_component_id &&
                    baseline.second_site_id == aligned.second_site_id,
                "CPA association site-pair identity drifted");
        require_same_sourced_value(
            baseline.epsilon_j_per_mol, aligned.epsilon_j_per_mol,
            "ThermoPack parity snapshot must preserve association epsilon provenance");
        require_same_sourced_value(
            baseline.beta_dimensionless, aligned.beta_dimensionless,
            "ThermoPack parity snapshot must preserve association beta provenance");
    }

    std::cout << "CPA_PARAMETER_SNAPSHOT_IDENTITY_OK"
              << " order=" << (swapped ? "water,methanol" : "methanol,water")
              << " literature_dataset=" << literature.dataset_id()
              << " parity_dataset=" << parity.dataset_id()
              << " Tc_MeOH_literature=512.64"
              << " Tc_MeOH_parity=512.6"
              << " Tc_H2O_literature=647.29"
              << " Tc_H2O_parity=647.3"
              << " kij_shared=-0.055"
              << " changed_scientific_fields=critical_temperature_k_only"
              << '\n';
}

} // namespace

int main() {
    try {
        audit_snapshot_identity(false);
        audit_snapshot_identity(true);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
