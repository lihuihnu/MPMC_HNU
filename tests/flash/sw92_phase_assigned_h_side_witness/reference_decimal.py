"""Independent Decimal(80) reference for Profile-C Gate C2a1.

This script imports only the already-independent stdlib Decimal C1 oracle.  It
never imports production C++.  The traceable C1 states are used to verify that
the retained H(NA) phase lies on the reconstructed common tangent.  A separate
explicit synthetic-test variant changes only the CO2/water nonaqueous BIP to
zero in order to create both (a) a negative NA point at the retained W
composition that must fail the physical-role guard and (b) a distinct,
water-poorer negative NA point usable only as an H-split witness seed.

The synthetic case is structural evidence only; it is not a physical SW92
parameter claim or an accepted three-phase state.
"""

from decimal import Decimal as D, localcontext
from pathlib import Path
import runpy

BASE = Path(__file__).parents[1] / "sw92_phase_assigned_joint" / "reference_decimal.py"
ref = runpy.run_path(str(BASE), run_name="sw92_phase_assigned_c1_reference")
minimum_phase = ref["minimum_phase"]
case = ref["case"]
NA_KIJ = ref["NA_KIJ"]


def tpd(gas_fraction: D, gas: str, family: str,
        pressure: D, temperature: D, molality: D, common):
    phase = minimum_phase(gas_fraction, gas, family, pressure, temperature, molality)
    composition = [gas_fraction, D(1) - gas_fraction]
    return sum(
        composition[i] *
        (composition[i].ln() + phase[2][i] - common[i])
        for i in range(2)
    )


def retained_h_tpd(gas: str, pressure: str, temperature: str,
                   molality: str, feed: str, x_seed: str, y_seed: str):
    values = case(gas, pressure, temperature, molality, feed, x_seed, y_seed)
    common = [values[6], values[7]]
    value = tpd(values[1], gas, "NA", D(pressure), D(temperature),
                D(molality), common)
    if abs(value) > D("2e-48"):
        raise AssertionError(f"retained H is not on C1 tangent: {value}")
    return value


if __name__ == "__main__":
    with localcontext() as context:
        context.prec = 80

        co2_h = retained_h_tpd("CO2", "3e6", "340", "0", "0.7", "0.006", "0.989")
        ch4_h = retained_h_tpd("CH4", "1e7", "350", "1", "0.5", "0.0009", "0.995")

        # Explicit synthetic-test-only NA interaction.  AQ remains the audited
        # corrected-original CO2 formula; only the NA constant is changed.
        old_na = NA_KIJ["CO2"]
        NA_KIJ["CO2"] = D("0")
        try:
            synthetic = case("CO2", "3e6", "340", "0", "0.5", "0.002", "0.988")
            x_w = synthetic[0]
            y_h = synthetic[1]
            common = [synthetic[6], synthetic[7]]
            tpd_at_w = tpd(x_w, "CO2", "NA", D("3e6"), D("340"), D("0"), common)
            trial_gas = D("0.003")
            tpd_distinct = tpd(trial_gas, "CO2", "NA", D("3e6"), D("340"), D("0"), common)

            # The retained W point is negative under NA, but it is not
            # water-poorer than W itself and therefore cannot seed H splitting.
            if not tpd_at_w < D("-1e-5"):
                raise AssertionError("synthetic NA-at-W point lost robust negative TPD")

            # The explicit trial is distinct from H and water-poorer than W.
            if not tpd_distinct < D("-1e-5"):
                raise AssertionError("synthetic distinct NA witness lost robust negative TPD")
            if not trial_gas > x_w:
                raise AssertionError("synthetic witness is not water-poorer than retained W")
            if not abs((trial_gas / y_h).ln()) > D("1"):
                raise AssertionError("synthetic witness is not compositionally distinct from H")

            print("traceable retained-H TPD CO2=", format(co2_h, ".35g"))
            print("traceable retained-H TPD CH4=", format(ch4_h, ".35g"))
            print("synthetic x_CO2^W=", format(x_w, ".35g"),
                  "y_CO2^H=", format(y_h, ".35g"))
            print("synthetic D_NA(W)=", format(tpd_at_w, ".35g"))
            print("synthetic D_NA(w=0.003)=", format(tpd_distinct, ".35g"))
        finally:
            NA_KIJ["CO2"] = old_na

        print("Independent Decimal(80) Profile-C C2a1 NA-witness reference passed")
