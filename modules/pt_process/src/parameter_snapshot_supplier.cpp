#include <mpmc/pt_process/parameter_snapshot_supplier.hpp>

#include <mpmc/flash/sw92_profile_c_pt_flash_backend.hpp>
#include <mpmc/pt_process/configured_backends.hpp>
#include <mpmc/thermodynamics/components.hpp>
#include <mpmc/thermodynamics/cpa_parameters.hpp>
#include <mpmc/thermodynamics/pr_parameters.hpp>
#include <mpmc/thermodynamics/sw92_parameters.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <optional>
#include <ostream>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mpmc::pt_process {
namespace {

namespace th = thermodynamics;
namespace fl = flash;

constexpr std::string_view pr76_configured_id =
    "pr76.methane-ethane-propane.literature-r1";
constexpr std::string_view sw92_configured_id =
    "sw92.carbon-dioxide-water.freshwater.literature-r1";
constexpr std::string_view cpa_configured_id =
    "cpa.methanol-water-333.15k.cr1.literature-r1";

bool valid_text(std::string_view value) noexcept {
    return !value.empty() && value.find('\0') == std::string_view::npos &&
           std::any_of(value.begin(), value.end(), [](unsigned char character) {
               return character > 0x20U;
           });
}

bool valid_identifier(std::string_view value) noexcept {
    return !value.empty() &&
           std::all_of(value.begin(), value.end(), [](unsigned char character) {
               return character >= 0x21U && character <= 0x7eU;
           });
}

th::Provenance literature_source(
    std::string reference, std::string revision, std::string locator,
    std::string note, std::string acquisition) {
    return {th::SourceKind::literature,
            std::move(reference),
            std::move(revision),
            std::move(locator),
            std::move(note),
            std::move(acquisition),
            "Limited attributed numerical facts and bibliographic metadata; "
            "source text is not redistributed"};
}

th::SourcedScalar sourced_scalar(
    double value, th::Unit unit, const th::Provenance& source,
    std::string original_unit = "SI or dimensionless",
    std::string conversion = "identity") {
    return {value, unit, source, std::move(original_unit),
            std::move(conversion)};
}

th::PrParameterSet pr76_snapshot() {
    const auto pure_source = literature_source(
        "https://doi.org/10.1002/aic.16730",
        "AIChE Journal 65(11), e16730; first published 26 July 2019",
        "Table 1, methane/ethane/propane Peng-Robinson parameters",
        "Tc/Pc/omega only; the fitted methane-propane interaction value in "
        "that paper is not used",
        "Manually transcribed from the open Wiley full text and Table 1");
    const auto pair_source = literature_source(
        "https://doi.org/10.1021/i160057a011",
        "Industrial & Engineering Chemistry Fundamentals 15(1), 1976",
        "Journal page 62, methane/ethane/propane discussion and Figure 4",
        "The original ternary example states that no interaction "
        "coefficients were used, represented here as three explicit zeros",
        "Manually transcribed from the project-supplied 1976 PDF");

    const std::array<th::Provenance, 3> molar_mass_sources{{
        literature_source(
            "https://webbook.nist.gov/cgi/cbook.cgi?ID=C74828",
            "NIST Chemistry WebBook SRD 69; accessed 2026-09-22",
            "Methane molecular weight 16.0425 g/mol",
            "Molar mass only; converted to kg/mol for the flow property contract",
            "Read from the public NIST Chemistry WebBook species record"),
        literature_source(
            "https://webbook.nist.gov/cgi/cbook.cgi?ID=C74840",
            "NIST Chemistry WebBook SRD 69; accessed 2026-09-22",
            "Ethane molecular weight 30.0690 g/mol",
            "Molar mass only; converted to kg/mol for the flow property contract",
            "Read from the public NIST Chemistry WebBook species record"),
        literature_source(
            "https://webbook.nist.gov/cgi/cbook.cgi?ID=C74986",
            "NIST Chemistry WebBook SRD 69; accessed 2026-09-22",
            "Propane molecular weight 44.0956 g/mol",
            "Molar mass only; converted to kg/mol for the flow property contract",
            "Read from the public NIST Chemistry WebBook species record")}};
    constexpr std::array<double, 3> molar_mass_kg_per_mol{
        0.0160425, 0.0300690, 0.0440956};

    constexpr std::array<std::string_view, 3> ids{
        "methane", "ethane", "propane"};
    constexpr std::array<std::array<double, 3>, 3> values{{
        {190.555, 4.595e6, 0.0},
        {305.4, 4.88e6, 0.099},
        {369.825, 4.248e6, 0.15308}}};

    std::vector<th::Component> catalog;
    th::PrParameterInput input;
    input.model_id = std::string(th::pr76_profile);
    input.dataset_id = "DeitersBell-aic16730-PengRobinson1976-ternary";
    input.revision =
        "DeitersBell-2019-Table1__PengRobinson-1976-zero-kij-r1";
    input.applicability = {
        std::nullopt, std::nullopt,
        literature_source(
            pair_source.reference, pair_source.revision,
            "The cited methane/ethane/propane example; no broad T/P bounds "
            "are inferred",
            "Unknown bounds are intentional and must not be read as global "
            "physical validation",
            pair_source.acquisition)};

    for (std::size_t index = 0; index < ids.size(); ++index) {
        const std::string id(ids[index]);
        catalog.push_back({
            id,
            id,
            th::ComponentKind::pure,
            pure_source,
            sourced_scalar(
                molar_mass_kg_per_mol[index],
                th::Unit::kilogram_per_mole,
                molar_mass_sources[index],
                "g/mol",
                "g/mol * 1e-3 -> kg/mol")});
        input.pure.push_back({
            id,
            sourced_scalar(values[index][0], th::Unit::kelvin, pure_source,
                            "K"),
            sourced_scalar(values[index][1], th::Unit::pascal, pure_source,
                            "MPa", "MPa * 1e6 -> Pa"),
            sourced_scalar(values[index][2], th::Unit::dimensionless,
                            pure_source, "dimensionless")});
    }
    const auto zero = sourced_scalar(
        0.0, th::Unit::dimensionless, pair_source, "dimensionless");
    input.binary.push_back({"methane", "ethane", zero});
    input.binary.push_back({"methane", "propane", zero});
    input.binary.push_back({"ethane", "propane", zero});

    const std::array<std::string, 3> order{
        "methane", "ethane", "propane"};
    return th::PrParameterSet::create(
        catalog, order, input, th::DataPolicy::ordinary);
}

th::Sw92ParameterSet sw92_snapshot() {
    const auto model_source = literature_source(
        "https://doi.org/10.1016/0378-3812(92)85105-H",
        "Fluid Phase Equilibria 77 (1992) 217-240; authors' errata; "
        "source-pdf-sha256:cb5b1d5034d78d934e887449ce0c89692d43431835d1371da5606d95b2c6bf58",
        "Corrected equations (9)-(17) and declared fresh-water slice",
        "Corrected-original SW92 CO2/water snapshot at fixed fresh-water "
        "molality; no global T/P accuracy range is inferred",
        "Manually transcribed from the project-supplied PDF and authors' errata");
    auto pure_source = model_source;
    pure_source.locator = "Table 3: CO2 and water pure properties";
    auto pair_source = model_source;
    pair_source.locator =
        "Table 5: non-aqueous CO2/water interaction parameter";

    std::vector<th::Component> catalog{
        {"carbon-dioxide", "Carbon dioxide", th::ComponentKind::pure,
         pure_source, std::nullopt},
        {"water", "Water", th::ComponentKind::pure,
         pure_source, std::nullopt}};

    th::Sw92ParameterInput input;
    input.model_id = std::string(th::sw92_corrected_profile);
    input.dataset_id = "SW92-Table3-Table5-corrected-CO2-water-freshwater";
    input.revision =
        "source-pdf-sha256-cb5b1d5034d78d934e887449ce0c89692d43431835d1371da5606d95b2c6bf58";
    input.applicability.state = {std::nullopt, std::nullopt, model_source};
    input.applicability.nacl_molality_mol_per_kg_water =
        th::ClosedInterval{0.0, 0.0};
    input.pure.push_back({
        "carbon-dioxide", th::Sw92Species::carbon_dioxide,
        sourced_scalar(304.2, th::Unit::kelvin, pure_source, "K"),
        sourced_scalar(7.38e6, th::Unit::pascal, pure_source, "bar",
                       "bar * 100000 -> Pa"),
        sourced_scalar(0.2273, th::Unit::dimensionless, pure_source,
                       "dimensionless")});
    input.pure.push_back({
        "water", th::Sw92Species::water,
        sourced_scalar(647.3, th::Unit::kelvin, pure_source, "K"),
        sourced_scalar(22.12e6, th::Unit::pascal, pure_source, "bar",
                       "bar * 100000 -> Pa"),
        sourced_scalar(0.3434, th::Unit::dimensionless, pure_source,
                       "dimensionless")});
    input.water_binary.push_back({
        "carbon-dioxide", th::Sw92NonAqueousWaterRule::constant,
        sourced_scalar(0.1896, th::Unit::dimensionless, pair_source,
                       "dimensionless")});

    const std::array<std::string, 2> order{"carbon-dioxide", "water"};
    return th::Sw92ParameterSet::create(
        catalog, order, input, th::DataPolicy::ordinary);
}

th::CpaSourcedValue cpa_value(
    double value, const th::Provenance& source) {
    return {value, source};
}

th::Provenance cpa_parameter_source(std::string locator) {
    return literature_source(
        "https://doi.org/10.2516/ogst:2008025",
        "Oil & Gas Science and Technology 63 (2008) 305-319",
        std::move(locator),
        "Methanol 2B and water 4C pure CPA parameters",
        "Manually transcribed from the cited open-access review");
}

th::Provenance cpa_cross_source(std::string locator) {
    return literature_source(
        "https://doi.org/10.2516/ogst:2008025",
        "Oil & Gas Science and Technology 63 (2008) 305-319",
        std::move(locator),
        "CR-1 cross value explicitly evaluated from cited pure parameters; "
        "no hidden runtime combining rule is used",
        "Arithmetic and geometric CR-1 rules evaluated from the cited Table "
        "1a values");
}

th::CpaParameterSet cpa_snapshot() {
    const auto kij_source = literature_source(
        "Georgios K. Folas, Modeling of Complex Mixtures Containing "
        "Hydrogen Bonding Molecules, PhD thesis, DTU",
        "Table 2.1",
        "Water/methanol at 333.15 K with CR-1: k12=-0.055",
        "The interaction parameter is used only for this fixed-temperature "
        "snapshot",
        "Manually transcribed from the thesis VLE correlation table");
    const auto applicability_source = literature_source(
        "https://doi.org/10.1021/je00019a033",
        "Journal of Chemical & Engineering Data 40 (1995) 679-684",
        "Methanol/water P-x-y at 333.15 K, interior validation points from "
        "39.223 kPa through 72.832 kPa",
        "Bounds identify the repository's validated slice and are not a claim "
        "outside that slice",
        "Manually transcribed from the published experimental dataset");
    std::vector<th::Component> catalog{
        {"METHANOL", "Methanol", th::ComponentKind::pure,
         cpa_parameter_source("Table 1a: methanol chemical identity"),
         std::nullopt},
        {"WATER", "Water", th::ComponentKind::pure,
         cpa_parameter_source("Table 1a: water chemical identity"),
         std::nullopt}};
    const std::array<std::string, 2> order{"METHANOL", "WATER"};

    th::CpaParameterInput input;
    input.dataset_id = "literature-cpa-water-methanol-cr1-333.15K";
    input.revision =
        "Kontogeorgis-2008-pure__Folas-CR1-kij__Kurihara-1995-VLE";
    input.applicability = {
        th::ClosedInterval{333.15, 333.15},
        th::ClosedInterval{3.9223e4, 7.2832e4},
        applicability_source};
    input.pure.push_back({
        "METHANOL",
        cpa_value(
            512.64,
            cpa_parameter_source("Table 1a: methanol Tc=512.64 K")),
        cpa_value(
            0.40531,
            cpa_parameter_source(
                "Table 1a: methanol a0=4.0531 bar L^2/mol^2 -> "
                "0.40531 Pa m^6/mol^2")),
        cpa_value(
            3.0978e-5,
            cpa_parameter_source(
                "Table 1a: methanol b=0.030978 L/mol -> "
                "3.0978e-5 m^3/mol")),
        cpa_value(
            0.43102,
            cpa_parameter_source(
                "Table 1a: methanol c1=0.43102")),
        {{"H", 1U}, {"e", 1U}}});
    input.pure.push_back({
        "WATER",
        cpa_value(
            647.29,
            cpa_parameter_source("Table 1a: water Tc=647.29 K")),
        cpa_value(
            0.12277,
            cpa_parameter_source(
                "Table 1a: water a0=1.2277 bar L^2/mol^2 -> "
                "0.12277 Pa m^6/mol^2")),
        cpa_value(
            1.4515e-5,
            cpa_parameter_source(
                "Table 1a: water b=0.014515 L/mol -> "
                "1.4515e-5 m^3/mol")),
        cpa_value(
            0.67359,
            cpa_parameter_source("Table 1a: water c1=0.67359")),
        {{"H", 2U}, {"e", 2U}}});
    input.binary.push_back(
        {"METHANOL", "WATER", cpa_value(-0.055, kij_source)});
    input.association_pairs.push_back({
        "METHANOL", "H", "METHANOL", "e",
        cpa_value(
            24591.0,
            cpa_parameter_source(
                "Table 1a: methanol epsilon=245.91 bar L/mol -> "
                "24591 J/mol")),
        cpa_value(
            0.0161,
            cpa_parameter_source("Table 1a: methanol beta=0.0161"))});
    input.association_pairs.push_back({
        "WATER", "H", "WATER", "e",
        cpa_value(
            16655.0,
            cpa_parameter_source(
                "Table 1a: water epsilon=166.55 bar L/mol -> "
                "16655 J/mol")),
        cpa_value(
            0.0692,
            cpa_parameter_source("Table 1a: water beta=0.0692"))});
    input.association_pairs.push_back({
        "METHANOL", "H", "WATER", "e",
        cpa_value(
            20623.0,
            cpa_cross_source(
                "CR-1 epsilon_cross arithmetic mean = "
                "206.23 bar L/mol -> 20623 J/mol")),
        cpa_value(
            0.03337843615270194,
            cpa_cross_source(
                "CR-1 beta_cross geometric mean "
                "sqrt(0.0161*0.0692)"))});
    input.association_pairs.push_back({
        "METHANOL", "e", "WATER", "H",
        cpa_value(
            20623.0,
            cpa_cross_source(
                "CR-1 epsilon_cross arithmetic mean = "
                "206.23 bar L/mol -> 20623 J/mol")),
        cpa_value(
            0.03337843615270194,
            cpa_cross_source(
                "CR-1 beta_cross geometric mean "
                "sqrt(0.0161*0.0692)"))});

    return th::CpaParameterSet::create(
        catalog, order, input, th::DataPolicy::ordinary);
}

std::vector<std::string> component_ids(
    const fl::PtFlashBackendCapability& capability) {
    return capability.component_ids;
}

PtParameterSnapshotDescriptor descriptor(
    const OwnedConfiguredPtBackend& backend, std::string location,
    std::string applicability,
    std::vector<PtParameterSnapshotSource> sources) {
    if (!backend.backend) {
        throw std::logic_error("PT snapshot supplier created a null backend");
    }
    const auto& capability = backend.backend->capability();
    return {backend.configured_backend_id,
            std::move(location),
            capability.model_profile,
            capability.dataset_id,
            capability.revision,
            component_ids(capability),
            std::move(applicability),
            std::move(sources)};
}

void write_json_string(std::ostream& output, std::string_view value) {
    output << '"';
    for (const char raw_character : value) {
        const auto character = static_cast<unsigned char>(raw_character);
        switch (character) {
        case '"':
            output << "\\\"";
            break;
        case '\\':
            output << "\\\\";
            break;
        case '\b':
            output << "\\b";
            break;
        case '\f':
            output << "\\f";
            break;
        case '\n':
            output << "\\n";
            break;
        case '\r':
            output << "\\r";
            break;
        case '\t':
            output << "\\t";
            break;
        default:
            if (character < 0x20U) {
                const auto flags = output.flags();
                const auto fill = output.fill();
                output << "\\u" << std::hex << std::setw(4)
                       << std::setfill('0')
                       << static_cast<unsigned int>(character);
                output.flags(flags);
                output.fill(fill);
            } else {
                output << raw_character;
            }
        }
    }
    output << '"';
}

bool capability_matches(
    const OwnedConfiguredPtBackend& backend,
    const PtParameterSnapshotDescriptor& snapshot) noexcept {
    if (!backend.backend ||
        backend.configured_backend_id != snapshot.configured_backend_id) {
        return false;
    }
    try {
        const auto& capability = backend.backend->capability();
        return capability.model_profile == snapshot.model_profile &&
               capability.dataset_id == snapshot.dataset_id &&
               capability.revision == snapshot.revision &&
               capability.component_ids == snapshot.component_ids;
    } catch (...) {
        return false;
    }
}

} // namespace

bool PtParameterSnapshotSource::structurally_valid() const noexcept {
    return valid_text(reference) && valid_text(revision) && valid_text(scope);
}

bool PtParameterSnapshotDescriptor::structurally_valid() const noexcept {
    if (!valid_identifier(configured_backend_id) ||
        !valid_text(snapshot_location) || !valid_identifier(model_profile) ||
        !valid_text(dataset_id) || !valid_text(revision) ||
        component_ids.empty() || !valid_text(applicability_scope) ||
        sources.empty()) {
        return false;
    }
    std::set<std::string, std::less<>> unique_components;
    for (const auto& component_id : component_ids) {
        if (!valid_identifier(component_id) ||
            !unique_components.insert(component_id).second) {
            return false;
        }
    }
    return std::all_of(
        sources.begin(), sources.end(),
        [](const PtParameterSnapshotSource& source) {
            return source.structurally_valid();
        });
}

bool PtParameterSnapshotBundle::structurally_valid() const noexcept {
    if (convention != pt_parameter_snapshot_bundle_convention ||
        !valid_identifier(bundle_id) || !valid_text(revision) ||
        !valid_text(location) || snapshots.empty() ||
        snapshots.size() != backends.size()) {
        return false;
    }
    std::set<std::string, std::less<>> unique_ids;
    for (std::size_t index = 0; index < snapshots.size(); ++index) {
        if (!snapshots[index].structurally_valid() ||
            !unique_ids.insert(snapshots[index].configured_backend_id).second ||
            !capability_matches(backends[index], snapshots[index])) {
            return false;
        }
    }
    return true;
}

RepositoryCuratedPtParameterSnapshotsV1
load_repository_curated_pt_parameters_v1() {
    return {pr76_snapshot(), sw92_snapshot(), cpa_snapshot()};
}

PtParameterSnapshotBundle
load_repository_curated_pt_parameter_snapshots_v1() {
    const auto parameters = load_repository_curated_pt_parameters_v1();

    PtParameterSnapshotBundle bundle;
    bundle.convention = std::string(pt_parameter_snapshot_bundle_convention);
    bundle.bundle_id = std::string(repository_curated_pt_bundle_id);
    bundle.revision = std::string(repository_curated_pt_bundle_revision);
    bundle.location =
        "builtin://mpmc/pt/repository-curated-literature-snapshots/v1";
    bundle.backends.push_back(make_pr76_configured_backend(
        std::string(pr76_configured_id), parameters.pr76));

    fl::Sw92ProfileCPtFlashBackendOptions sw92_options;
    sw92_options.nacl_molality_mol_per_kg_water = 0.0;
    bundle.backends.push_back(make_sw92_configured_backend(
        std::string(sw92_configured_id), parameters.sw92,
        std::move(sw92_options)));
    bundle.backends.push_back(make_cpa_configured_backend(
        std::string(cpa_configured_id), parameters.cpa));

    bundle.snapshots.push_back(descriptor(
        bundle.backends[0],
        "builtin://mpmc/pt/repository-curated-literature-snapshots/v1/pr76",
        "No broad T/P bounds are claimed; use only where the cited "
        "methane/ethane/propane PR76 parameterization is appropriate",
        {{"https://doi.org/10.1002/aic.16730",
          "AIChE Journal 65(11), e16730",
          "Table 1 pure-component Tc/Pc/omega"},
         {"https://doi.org/10.1021/i160057a011",
          "Industrial & Engineering Chemistry Fundamentals 15(1), 1976",
          "Explicit zero binary interactions for the cited ternary example"},
         {"https://webbook.nist.gov/",
          "NIST Chemistry WebBook SRD 69; methane/ethane/propane species records; accessed 2026-09-22",
          "Molecular weights used by the selected-phase mass-density/property closure"}}));
    bundle.snapshots.push_back(descriptor(
        bundle.backends[1],
        "builtin://mpmc/pt/repository-curated-literature-snapshots/v1/sw92",
        "CO2/water only; corrected-original SW92; fixed 0 mol NaCl per kg "
        "water; no broad T/P validation range is claimed",
        {{"https://doi.org/10.1016/0378-3812(92)85105-H",
          "Fluid Phase Equilibria 77 (1992) 217-240 with authors' errata; "
          "source-pdf-sha256:cb5b1d5034d78d934e887449ce0c89692d43431835d1371da5606d95b2c6bf58",
          "Corrected equations, Table 3 and Table 5"}}));
    bundle.snapshots.push_back(descriptor(
        bundle.backends[2],
        "builtin://mpmc/pt/repository-curated-literature-snapshots/v1/cpa",
        "Methanol/water CR-1 only; T=333.15 K and "
        "P in [39223,72832] Pa",
        {{"https://doi.org/10.2516/ogst:2008025",
          "Oil & Gas Science and Technology 63 (2008) 305-319",
          "Table 1a pure and association parameters plus CR-1 definition"},
         {"Georgios K. Folas, Modeling of Complex Mixtures Containing "
          "Hydrogen Bonding Molecules, PhD thesis, DTU",
          "Table 2.1",
          "Water/methanol binary interaction record at 333.15 K using CR-1"},
         {"https://doi.org/10.1021/je00019a033",
          "Journal of Chemical & Engineering Data 40 (1995) 679-684",
          "333.15 K methanol/water physical-validation pressure interval"}}));

    if (!bundle.structurally_valid()) {
        throw std::logic_error(
            "repository-curated PT snapshot bundle failed its integrity check");
    }
    return bundle;
}

void write_pt_parameter_snapshot_manifest_json(
    std::ostream& output, const PtParameterSnapshotBundle& bundle) {
    if (!bundle.structurally_valid()) {
        throw std::invalid_argument(
            "cannot serialize an invalid PT parameter snapshot bundle");
    }
    output << "{\"convention\":";
    write_json_string(output, bundle.convention);
    output << ",\"bundle_id\":";
    write_json_string(output, bundle.bundle_id);
    output << ",\"revision\":";
    write_json_string(output, bundle.revision);
    output << ",\"location\":";
    write_json_string(output, bundle.location);
    output << ",\"snapshots\":[";
    for (std::size_t index = 0; index < bundle.snapshots.size(); ++index) {
        if (index != 0U) {
            output << ',';
        }
        const auto& snapshot = bundle.snapshots[index];
        output << "{\"configured_backend_id\":";
        write_json_string(output, snapshot.configured_backend_id);
        output << ",\"snapshot_location\":";
        write_json_string(output, snapshot.snapshot_location);
        output << ",\"model_profile\":";
        write_json_string(output, snapshot.model_profile);
        output << ",\"dataset_id\":";
        write_json_string(output, snapshot.dataset_id);
        output << ",\"revision\":";
        write_json_string(output, snapshot.revision);
        output << ",\"component_ids\":[";
        for (std::size_t component = 0;
             component < snapshot.component_ids.size(); ++component) {
            if (component != 0U) {
                output << ',';
            }
            write_json_string(output, snapshot.component_ids[component]);
        }
        output << "],\"applicability_scope\":";
        write_json_string(output, snapshot.applicability_scope);
        output << ",\"sources\":[";
        for (std::size_t source = 0; source < snapshot.sources.size();
             ++source) {
            if (source != 0U) {
                output << ',';
            }
            output << "{\"reference\":";
            write_json_string(output, snapshot.sources[source].reference);
            output << ",\"revision\":";
            write_json_string(output, snapshot.sources[source].revision);
            output << ",\"scope\":";
            write_json_string(output, snapshot.sources[source].scope);
            output << '}';
        }
        output << "]}";
    }
    output << "]}\n";
}

} // namespace mpmc::pt_process
