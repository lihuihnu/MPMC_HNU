"""Validate that every selected CI owner completed successfully."""
import json
import os
from pathlib import Path

SIMPLE = {
    "flow-core": "flow_core",
    "flow-discretization": "flow_discretization",
    "flow-discretization-petsc": "flow_discretization_petsc",
    "pt-flash-backend": "pt_flash_backend",
    "pt-stability": "pt_stability",
    "pt-split": "pt_split",
    "pr76-pt-continuation": "pr76_pt_continuation",
    "pr76-max3": "pr76_max3",
    "cpa-stability": "cpa_stability",
    "cpa-split": "cpa_split",
    "cpa-max3": "cpa_max3",
    "cpa-baseline": "cpa_baseline",
    "cpa-pt-phase": "cpa_pt_phase",
    "cpa-physical-validation": "cpa_physical_validation",
    "thermodynamics": "thermo_has_work",
    "physics-closure": "physics_closure",
    "model-configuration": "model_configuration",
    "sw92-thermodynamics": "sw92_thermodynamics",
    "sw92-family-vle": "sw92_family_vle",
    "sw92-profile-c-phase-set": "sw92_profile_c_phase_set",
    "sw92-profile-c-sensitivity": "sw92_profile_c_sensitivity",
    "sw92-physics-closure": "sw92_physics_closure",
    "ad": "ad_has_work",
}

def yes(impact, key):
    return str(impact.get(key, "")).lower() == "true"

def required_jobs(impact, catalog):
    required = {job for job, key in SIMPLE.items() if yes(impact, key)}

    high_joint = any(yes(impact, key) for key in (
        "sw92_phase_assigned_h_side_witness",
        "sw92_profile_c_pt",
        "sw92_phase_assigned_three_phase",
        "sw92_phase_assigned_three_phase_closure",
        "sw92_phase_assigned_boundary",
    ))
    if yes(impact, "sw92_phase_assigned_joint") and not high_joint:
        required.add("sw92-phase-assigned-joint")

    high_h_side = any(yes(impact, key) for key in (
        "sw92_profile_c_pt",
        "sw92_phase_assigned_three_phase",
        "sw92_phase_assigned_three_phase_closure",
        "sw92_phase_assigned_boundary",
    ))
    if yes(impact, "sw92_phase_assigned_h_side_witness") and not high_h_side:
        required.add("sw92-phase-assigned-h-side-witness")

    if yes(impact, "sw92_phase_assigned_no_w") and not (
        yes(impact, "sw92_profile_c_pt") or yes(impact, "sw92_phase_assigned_boundary")
    ):
        required.add("sw92-phase-assigned-no-w")

    if any(yes(impact, key) for key in (
        "sw92_profile_c_pt",
        "sw92_profile_c_phase_set",
        "sw92_phase_assigned_three_phase",
        "sw92_phase_assigned_three_phase_closure",
        "sw92_phase_assigned_boundary",
    )):
        required.add("sw92-phase-assigned-three-phase-topology")

    for spec in catalog["workflows"].values():
        key = spec.get("route_key")
        if key and yes(impact, key):
            required.update(spec.get("central_hashes", {}).keys())
    return required

def validate(needs, impact, catalog):
    assert needs.get("impact", {}).get("result") == "success", "impact did not complete"
    assert yes(impact, "trusted"), "trusted execution required"
    failed = {
        job: record.get("result")
        for job, record in needs.items()
        if record.get("result") not in ("success", "skipped")
    }
    assert not failed, ("failed-or-cancelled", failed)

    required = required_jobs(impact, catalog)
    missing = {
        job: needs.get(job, {}).get("result", "missing")
        for job in sorted(required)
        if needs.get(job, {}).get("result") != "success"
    }
    assert not missing, ("selected-but-not-successful", missing)
    return required

def main():
    needs = json.loads(os.environ["NEEDS_JSON"])
    impact = json.loads(os.environ["IMPACT_JSON"])
    catalog = json.loads(Path(".github/ci/workflow_map.json").read_text(encoding="utf-8"))
    required = validate(needs, impact, catalog)
    print("Required selected jobs completed:", len(required))

if __name__ == "__main__":
    main()
