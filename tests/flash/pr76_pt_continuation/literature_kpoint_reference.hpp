#ifndef MPMC_TEST_PR76_LITERATURE_KPOINT_REFERENCE_HPP
#define MPMC_TEST_PR76_LITERATURE_KPOINT_REFERENCE_HPP

#include <mpmc/thermodynamics/pr76_phase.hpp>

#include <array>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace pr76_literature_reference {
namespace th = mpmc::thermodynamics;

inline constexpr const char* smits_reference =
    "J.C. Smits, Phase Behaviour in Certain CO2 + n-Alkane + 1-Alkanol Systems: Experiments and Modelling, Delft University of Technology, 1996, repository UUID ab0a04f3-3491-4314-a95c-a943aff2f955";
inline constexpr const char* mushrif_reference =
    "Samir H. Mushrif, Determining Equation of State Binary Interaction Parameters Using K- and L-Points, University of Saskatchewan, 2004, hdl:10388/etd-10212004-233350";
inline constexpr const char* redistribution_note =
    "Bibliographic and numeric facts are transcribed for scientific regression metadata only; source copyright remains with the author/publisher and no source text is redistributed.";

inline th::Provenance smits_source(std::string locator) {
    return {th::SourceKind::literature,
            smits_reference,
            "1996 master thesis",
            std::move(locator),
            "Experimental ternary critical-endpoint record; not a distinct-three-phase composition oracle.",
            "Manual transcription from the TU Delft repository thesis record/scan.",
            redistribution_note};
}

inline th::Provenance mushrif_source(std::string locator) {
    return {th::SourceKind::literature,
            mushrif_reference,
            "October 2004 thesis",
            std::move(locator),
            "Peng-Robinson model input used only with the exact cited locator; no fitted accuracy is inferred beyond the source statement.",
            "Manual transcription from the University of Saskatchewan / Library and Archives Canada thesis copy.",
            redistribution_note};
}

inline th::SourcedScalar kelvin(double value, const th::Provenance& source) {
    return {value, th::Unit::kelvin, source, "K", "identity"};
}

inline th::SourcedScalar pressure_from_bar(
    double value_bar, const th::Provenance& source) {
    return {value_bar * 1.0e5, th::Unit::pascal, source,
            "bar", "bar * 1e5 -> Pa"};
}

inline th::SourcedScalar dimensionless(
    double value, const th::Provenance& source) {
    return {value, th::Unit::dimensionless, source,
            "dimensionless", "identity"};
}

struct CriticalEndpointAnchor {
    double pressure_pa{9.22e6};
    double temperature_k{317.43};
    double co2_mole_fraction{0.9501};
    double tridecane_fraction_on_co2_free_basis{0.8603};
    bool experimental{true};
    bool critical_endpoint{true};
    bool suitable_as_distinct_three_phase_flash_oracle{false};
    th::Provenance source{smits_source(
        "Chapter 5, Table 5.1, system 4: CO2 + n-tridecane + 1-pentanol, (L=VL) UCEP")};

    [[nodiscard]] std::array<double, 3> ordered_feed() const noexcept {
        const double co2_free = 1.0 - co2_mole_fraction;
        const double tridecane =
            co2_free * tridecane_fraction_on_co2_free_basis;
        const double pentanol = co2_free - tridecane;
        // Project order: CO2, 1-pentanol, n-tridecane.
        return {co2_mole_fraction, pentanol, tridecane};
    }
};

inline CriticalEndpointAnchor experimental_kpoint() {
    return {};
}

inline th::Pr76Phase<double> illustrative_pr76_model() {
    const auto pure_source = mushrif_source(
        "Appendix B, Table B1: critical properties and acentric factors (citing Yaws, 1999)");
    const auto bip_source = mushrif_source(
        "Section 3.5.1, Figure 3.4: Gibbs energy plot at 317.43 K and 9.22 MPa; delta12=0.17, delta13=0.15, delta23=0.10");
    const std::array<std::string, 3> ids{
        "CO2", "1-pentanol", "n-tridecane"};

    std::vector<th::Component> catalog;
    catalog.reserve(ids.size());
    catalog.push_back({ids[0], "carbon dioxide", th::ComponentKind::pure,
                       pure_source, {}});
    catalog.push_back({ids[1], "1-pentanol", th::ComponentKind::pure,
                       pure_source, {}});
    catalog.push_back({ids[2], "n-tridecane", th::ComponentKind::pure,
                       pure_source, {}});

    th::PrParameterInput input;
    input.model_id = std::string(th::pr76_profile);
    input.dataset_id = "Mushrif-2004-Fig3.4-PR76-illustrative-kpoint-surface";
    input.revision = "literature-transcription-v1";
    input.applicability = {
        std::nullopt, std::nullopt,
        mushrif_source(
            "Section 3.5.1, Figure 3.4; the dataset is a traceable point/model anchor, not a validated broad applicability range")};

    input.pure.push_back({
        ids[0], kelvin(304.1, pure_source),
        pressure_from_bar(73.8, pure_source),
        dimensionless(0.2390, pure_source)});
    input.pure.push_back({
        ids[1], kelvin(588.2, pure_source),
        pressure_from_bar(39.1, pure_source),
        dimensionless(0.5784, pure_source)});
    input.pure.push_back({
        ids[2], kelvin(676.0, pure_source),
        pressure_from_bar(17.2, pure_source),
        dimensionless(0.6203, pure_source)});

    input.binary.push_back({
        ids[0], ids[1], dimensionless(0.17, bip_source)});
    input.binary.push_back({
        ids[0], ids[2], dimensionless(0.15, bip_source)});
    input.binary.push_back({
        ids[1], ids[2], dimensionless(0.10, bip_source)});

    return th::Pr76Phase<double>::from_parameters(
        th::PrParameterSet::create(catalog, ids, input, th::DataPolicy::ordinary));
}

} // namespace pr76_literature_reference

#endif // MPMC_TEST_PR76_LITERATURE_KPOINT_REFERENCE_HPP
