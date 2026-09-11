"""Independent Decimal(80) structural checks for Profile-C no-W routing.

Stdlib only; no production imports. This combines the already-independent SW92
one-family VLE transcription with the independent Profile-C C1 minimum-root
transcription. It verifies that a wet fixed-NA candidate has a negative AQ
phase-addition direction toward a composition that is water-richer than the
water-poor retained NA branch. A fixed-NA pair may also contain a W-like
water-rich mathematical phase; that must not be promoted to physical H role.

A dry feed beyond the W+H coexistence H endpoint is checked independently
against the same W trial. These are model numerical checks, not experimental
validation or global stability proofs.
"""

from decimal import Decimal as D, localcontext
from pathlib import Path
import runpy

HERE = Path(__file__).resolve()
VLE = HERE.parents[1] / "sw92_family_vle" / "reference_decimal.py"
C1 = HERE.parents[1] / "sw92_phase_assigned_joint" / "reference_decimal.py"
vle = runpy.run_path(str(VLE), run_name="sw92_no_w_vle_reference")
c1 = runpy.run_path(str(C1), run_name="sw92_no_w_c1_reference")

phase = vle["phase"]
solve_case = vle["solve_case"]
minimum_phase = c1["minimum_phase"]
case = c1["case"]


def log_activity_binary(gas_fraction: D, selected):
    composition = [gas_fraction, D(1) - gas_fraction]
    return [composition[i].ln() + selected[1][i] for i in range(2)]


def tpd_aq(gas_fraction: D, common, pressure: D, temperature: D):
    selected = minimum_phase(gas_fraction, "CO2", "AQ", pressure, temperature, D(0))
    composition = [gas_fraction, D(1) - gas_fraction]
    return sum(composition[i] *
               (composition[i].ln() + selected[2][i] - common[i])
               for i in range(2))


if __name__ == "__main__":
    with localcontext() as context:
        context.prec = 80
        pressure = D("3e6")
        temperature = D("340")

        # Wet feed: independent fixed-NA two-phase equilibrium. The low-Z NA
        # mathematical branch is itself extremely water-rich; only the
        # water-poor NA branch is a required H contrast for W candidate creation.
        x, y, beta, _, _ = solve_case(
            "CO2", "NA", "3e6", "340", "0", "0.7", "0.00028", "0.9887")
        low = phase(x, "CO2", "NA", pressure, temperature, D(0), "liquid")
        high = phase(y, "CO2", "NA", pressure, temperature, D(0), "vapor")
        mu_low = log_activity_binary(x, low)
        mu_high = log_activity_binary(y, high)
        common_pair = [(mu_low[i] + mu_high[i]) / D(2) for i in range(2)]

        c1_values = case("CO2", "3e6", "340", "0", "0.7", "0.006", "0.989")
        w_gas = c1_values[0]
        wet_tpd = tpd_aq(w_gas, common_pair, pressure, temperature)
        if not wet_tpd < D("-1e-6"):
            raise AssertionError(f"wet fixed-NA state lost AQ phase-addition direction: {wet_tpd}")
        w_water = D(1) - w_gas
        min_na_water = min(D(1) - x, D(1) - y)
        max_na_water = max(D(1) - x, D(1) - y)
        if not w_water > min_na_water:
            raise AssertionError("independent AQ witness lost W/H water-richness contrast")
        if not max_na_water > w_water:
            raise AssertionError("fixture no longer exposes the W-like fixed-NA mathematical branch")

        # Dry feed lies beyond the H endpoint of the independent W+H tie line.
        # Check the same independently solved W composition against the minimum-
        # Gibbs NA feed tangent. This is a targeted structural sign check only.
        dry_feed = D("0.995")
        dry_na = minimum_phase(dry_feed, "CO2", "NA", pressure, temperature, D(0))
        dry_common = [dry_feed.ln() + dry_na[2][0],
                      (D(1) - dry_feed).ln() + dry_na[2][1]]
        dry_tpd = tpd_aq(w_gas, dry_common, pressure, temperature)
        if not dry_tpd > D("1e-6"):
            raise AssertionError(f"dry H-only tangent unexpectedly favors the W trial: {dry_tpd}")

        print("wet fixed-NA x_CO2=", format(x, ".35g"),
              "y_CO2=", format(y, ".35g"), "beta=", format(beta, ".35g"))
        print("wet NA water range=", format(min_na_water, ".35g"),
              format(max_na_water, ".35g"), "AQ-W water=", format(w_water, ".35g"))
        print("wet targeted AQ TPD=", format(wet_tpd, ".35g"))
        print("dry targeted AQ TPD=", format(dry_tpd, ".35g"))
        print("Independent Decimal(80) Profile-C no-W structural checks passed")
