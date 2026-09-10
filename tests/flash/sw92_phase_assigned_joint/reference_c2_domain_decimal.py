"""Independent Decimal(80) blocker audit for SW92 Profile-C Gate C2.

This script imports only the already independent stdlib Decimal C1 reference
implementation in the same test directory.  It does not import production C++
or any Profile-A/Profile-B result.  Starting from the directly solved C1 AQ+NA
coexistence states, it evaluates the opposite family at each assigned physical
phase composition and proves the identity

    D_AQ(y_NA) = g_AQ(y_NA) - g_NA(y_NA)
    D_NA(x_AQ) = g_NA(x_AQ) - g_AQ(x_AQ)

against the C1 common tangent.  A robust negative D_AQ(y_NA) means an
*unrestricted* two-family final TPD review would reintroduce Xu lower-envelope
family competition and is therefore not an authorized Profile-C C2 finalizer.
"""

from decimal import Decimal as D, localcontext
from pathlib import Path
import runpy

BASE = Path(__file__).with_name("reference_decimal.py")
ref = runpy.run_path(str(BASE), run_name="sw92_phase_assigned_c1_reference")
minimum_phase = ref["minimum_phase"]
case = ref["case"]


def tpd(gas_fraction: D, gas: str, family: str,
        pressure: D, temperature: D, molality: D,
        common_log_activity):
    phase = minimum_phase(gas_fraction, gas, family,
                          pressure, temperature, molality)
    composition = [gas_fraction, D(1) - gas_fraction]
    return sum(
        composition[i] *
        (composition[i].ln() + phase[2][i] - common_log_activity[i])
        for i in range(2)
    )


def family_gap(gas_fraction: D, gas: str,
               pressure: D, temperature: D, molality: D):
    aq = minimum_phase(gas_fraction, gas, "AQ",
                       pressure, temperature, molality)
    na = minimum_phase(gas_fraction, gas, "NA",
                       pressure, temperature, molality)
    composition = [gas_fraction, D(1) - gas_fraction]
    return sum(
        composition[i] * (aq[2][i] - na[2][i])
        for i in range(2)
    )


def audit_case(gas: str, pressure_value: str, temperature_value: str,
               molality_value: str, feed_value: str,
               x_seed: str, y_seed: str):
    pressure = D(pressure_value)
    temperature = D(temperature_value)
    molality = D(molality_value)
    values = case(gas, pressure_value, temperature_value,
                  molality_value, feed_value, x_seed, y_seed)
    x_aq = values[0]
    y_na = values[1]
    common = [values[6], values[7]]

    d_aq_at_na = tpd(y_na, gas, "AQ", pressure, temperature, molality, common)
    d_na_at_aq = tpd(x_aq, gas, "NA", pressure, temperature, molality, common)
    gap_at_na = family_gap(y_na, gas, pressure, temperature, molality)
    gap_at_aq = family_gap(x_aq, gas, pressure, temperature, molality)

    # At exact C1 coexistence, the candidate tangent equals the assigned-family
    # chemical potential at each retained phase.  Therefore the opposite-family
    # TPD is exactly the family Gibbs gap at that same composition.
    if abs(d_aq_at_na - gap_at_na) > D("2e-48"):
        raise AssertionError(
            f"AQ-at-NA TPD/gap identity failed: {d_aq_at_na} vs {gap_at_na}")
    if abs(d_na_at_aq + gap_at_aq) > D("2e-48"):
        raise AssertionError(
            f"NA-at-AQ TPD/gap identity failed: {d_na_at_aq} vs {-gap_at_aq}")

    return d_aq_at_na, d_na_at_aq, gap_at_na, gap_at_aq


if __name__ == "__main__":
    with localcontext() as context:
        context.prec = 80

        co2 = audit_case("CO2", "3e6", "340", "0", "0.7", "0.006", "0.989")
        ch4 = audit_case("CH4", "1e7", "350", "1", "0.5", "0.0009", "0.995")

        # These literals are human-readable regression summaries, not new model
        # parameters. The decisive high-precision check is the TPD/family-gap
        # identity above; this looser tail check only catches accidental drift.
        expected = [
            (D("-0.0013073695742952295175108491551836972166774496602498"),
             D("0.017565890359096352643045761530903617546075009971271")),
            (D("-0.003109292386691609213101625542793603606717073421159"),
             D("0.004526216421619292162543436283399494459601967755158")),
        ]
        for label, values, anchors in (("CO2", co2, expected[0]),
                                       ("CH4", ch4, expected[1])):
            for actual, target in zip(values[:2], anchors):
                if abs(actual - target) > D("1e-30"):
                    raise AssertionError(
                        f"{label} C2 domain anchor drift: {actual} vs {target}")
            if not values[0] < D("-1e-4"):
                raise AssertionError(
                    f"{label} unrestricted AQ-at-NA TPD is no longer robustly negative")
            if not values[1] > D("1e-4"):
                raise AssertionError(
                    f"{label} NA-at-AQ control TPD lost positive separation")
            print(label,
                  "D_AQ(y_NA)=", format(values[0], ".40g"),
                  "D_NA(x_AQ)=", format(values[1], ".40g"))

        print("Independent Decimal(80) Profile-C C2 role-domain blocker audit passed")