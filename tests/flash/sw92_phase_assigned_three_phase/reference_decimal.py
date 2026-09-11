"""Independent Decimal(80) structural reference for Profile-C Gate C2b.1.

Stdlib only; no production imports. The SW92 pure-component and water-pair
formulas are reused from the already-independent ternary Decimal oracle. The
only synthetic datum is CH4/CO2 kij=0 in both AQ and NA families. SW92 does not
supply a general CH4/CO2 non-water BIP, so this fixture is deliberately tagged
as structural evidence only: it validates the unordered W(AQ)+H0(NA)+H1(NA)
three-phase equations and phase-addition path, not a physical three-phase
prediction or experimental datum.

State: p=3 MPa, T=260 K, fresh water. A metastable C1 W+H tie-line is used to
construct a feed. The C1 common tangent has a robustly negative NA TPD toward
the second hydrocarbon branch. The same feed lies strictly inside the
independently solved W+H0+H1 triangle.
"""

from decimal import Decimal as D, localcontext
from pathlib import Path
import re
import runpy

BASE = (Path(__file__).parents[1] /
        "sw92_asymmetric_max2_physical_ternary" / "reference_decimal.py")
ref = runpy.run_path(str(BASE), run_name="sw92_c2b1_base_reference")
minimum_phase = ref["minimum_phase"]
solve_linear = ref["solve_linear"]

# Mutate only the isolated oracle's globals. These names are read by its
# functions at call time; no production code or repository parameter object is
# imported.
g = minimum_phase.__globals__
g["PRESSURE"] = D("3e6")
g["TEMPERATURE"] = D("260")
g["MOLALITY"] = D("0")
g["GAS_GAS_KIJ"] = D("0")
g["FEED"] = [D("0.1"), D("0.6"), D("0.3")]  # placeholder; replaced below

N = 3


def composition(a: D, b: D):
    return [a, b, D(1) - a - b]


def mu(values, family: str):
    phase = minimum_phase(values, family)
    return [values[i].ln() + phase[2][i] for i in range(N)], phase


def valid_composition(values) -> bool:
    return all(D(0) < value < D(1) for value in values) and sum(values) == D(1)


def residual_three(unknown):
    w = composition(unknown[0], unknown[1])
    h0 = composition(unknown[2], unknown[3])
    h1 = composition(unknown[4], unknown[5])
    if not (valid_composition(w) and valid_composition(h0) and valid_composition(h1)):
        raise ArithmeticError("invalid three-phase Decimal iterate")
    mu_w, _ = mu(w, "AQ")
    mu_h0, _ = mu(h0, "NA")
    mu_h1, _ = mu(h1, "NA")
    return [mu_w[i] - mu_h0[i] for i in range(N)] + [
        mu_w[i] - mu_h1[i] for i in range(N)]


def solve_three():
    # Deliberately rough values near the independently discovered branches; the
    # Newton solve, not these seeds, defines the reference.
    unknown = [D("0.0005"), D("0.056"),
               D("0.04"), D("0.959"),
               D("0.18"), D("0.819")]
    base_step = D("1e-21")
    for _ in range(60):
        residual = residual_three(unknown)
        norm = max(abs(value) for value in residual)
        if norm < D("1e-52"):
            break
        jacobian = [[D(0) for _ in range(6)] for _ in range(6)]
        for column in range(6):
            step = min(base_step,
                       max(D("1e-29"), abs(unknown[column]) * D("1e-18")))
            plus = list(unknown)
            minus = list(unknown)
            plus[column] += step
            minus[column] -= step
            rp = residual_three(plus)
            rm = residual_three(minus)
            for row in range(6):
                jacobian[row][column] = (rp[row] - rm[row]) / (D(2) * step)
        delta = solve_linear(jacobian, [-value for value in residual])
        damping = D(1)
        accepted = False
        while damping > D("1e-18"):
            candidate = [unknown[i] + damping * delta[i] for i in range(6)]
            try:
                next_norm = max(abs(value) for value in residual_three(candidate))
            except (ArithmeticError, ValueError):
                next_norm = D("Infinity")
            if next_norm < norm:
                unknown = candidate
                accepted = True
                break
            damping /= D(2)
        if not accepted:
            raise ArithmeticError("three-phase Decimal Newton line search failed")
    else:
        raise ArithmeticError("three-phase Decimal Newton iteration limit")
    if max(abs(value) for value in residual_three(unknown)) > D("1e-48"):
        raise ArithmeticError("three-phase Decimal residual did not converge")
    return (composition(unknown[0], unknown[1]),
            composition(unknown[2], unknown[3]),
            composition(unknown[4], unknown[5]))


def residual_c1(unknown):
    # Fix x_CH4^W=0.001 to choose one metastable W(AQ)+H(NA) tie-line.
    w = [D("0.001"), unknown[0], D(1) - D("0.001") - unknown[0]]
    h = composition(unknown[1], unknown[2])
    if not (valid_composition(w) and valid_composition(h)):
        raise ArithmeticError("invalid C1 Decimal iterate")
    mu_w, _ = mu(w, "AQ")
    mu_h, _ = mu(h, "NA")
    return [mu_w[i] - mu_h[i] for i in range(N)]


def solve_c1():
    unknown = [D("0.052"), D("0.09"), D("0.909")]
    base_step = D("1e-21")
    for _ in range(60):
        residual = residual_c1(unknown)
        norm = max(abs(value) for value in residual)
        if norm < D("1e-52"):
            break
        jacobian = [[D(0) for _ in range(3)] for _ in range(3)]
        for column in range(3):
            step = min(base_step,
                       max(D("1e-29"), abs(unknown[column]) * D("1e-18")))
            plus = list(unknown)
            minus = list(unknown)
            plus[column] += step
            minus[column] -= step
            rp = residual_c1(plus)
            rm = residual_c1(minus)
            for row in range(3):
                jacobian[row][column] = (rp[row] - rm[row]) / (D(2) * step)
        delta = solve_linear(jacobian, [-value for value in residual])
        damping = D(1)
        accepted = False
        while damping > D("1e-18"):
            candidate = [unknown[i] + damping * delta[i] for i in range(3)]
            try:
                next_norm = max(abs(value) for value in residual_c1(candidate))
            except (ArithmeticError, ValueError):
                next_norm = D("Infinity")
            if next_norm < norm:
                unknown = candidate
                accepted = True
                break
            damping /= D(2)
        if not accepted:
            raise ArithmeticError("C1 Decimal Newton line search failed")
    else:
        raise ArithmeticError("C1 Decimal Newton iteration limit")
    if max(abs(value) for value in residual_c1(unknown)) > D("1e-48"):
        raise ArithmeticError("C1 Decimal residual did not converge")
    w = [D("0.001"), unknown[0], D(1) - D("0.001") - unknown[0]]
    h = composition(unknown[1], unknown[2])
    return w, h


def tpd(point, common, family="NA"):
    phase = minimum_phase(point, family)
    return sum(point[i] * (point[i].ln() + phase[2][i] - common[i])
               for i in range(N))


def fractions_for_feed(w, h0, h1, feed):
    matrix = [[w[i], h0[i], h1[i]] for i in range(N)]
    fractions = solve_linear(matrix, feed)
    if any(value <= 0 for value in fractions):
        raise AssertionError("reference feed is not inside the three-phase triangle")
    if abs(sum(fractions) - D(1)) > D("1e-48"):
        raise AssertionError("three-phase fraction sum changed")
    return fractions


def anchors():
    w, h0, h1 = solve_three()
    mu_w, phase_w = mu(w, "AQ")
    mu_h0, phase_h0 = mu(h0, "NA")
    mu_h1, phase_h1 = mu(h1, "NA")
    if max(abs(mu_w[i] - mu_h0[i]) for i in range(N)) > D("1e-48"):
        raise AssertionError("W/H0 chemical potentials differ")
    if max(abs(mu_w[i] - mu_h1[i]) for i in range(N)) > D("1e-48"):
        raise AssertionError("W/H1 chemical potentials differ")

    c1_w, c1_h = solve_c1()
    beta_c1 = D("0.73")
    feed = [(D(1) - beta_c1) * c1_w[i] + beta_c1 * c1_h[i]
            for i in range(N)]
    g["FEED"] = feed
    fractions = fractions_for_feed(w, h0, h1, feed)
    for i in range(N):
        recovered = (fractions[0] * w[i] + fractions[1] * h0[i] +
                     fractions[2] * h1[i])
        if abs(recovered - feed[i]) > D("2e-48"):
            raise AssertionError("three-phase reference material balance changed")

    c1_mu_w, _ = mu(c1_w, "AQ")
    c1_mu_h, _ = mu(c1_h, "NA")
    common_c1 = [(c1_mu_w[i] + c1_mu_h[i]) / D(2) for i in range(N)]
    retained_tpd = tpd(c1_h, common_c1)
    witness_tpd = tpd(h1, common_c1)
    if abs(retained_tpd) > D("2e-48"):
        raise AssertionError("retained C1 H is not on its common tangent")
    if not witness_tpd < D("-0.05"):
        raise AssertionError("synthetic C1 state lost robust additional-NA witness")

    return {
        "w": w, "h0": h0, "h1": h1,
        "fractions": fractions,
        "z": [phase_w[1], phase_h0[1], phase_h1[1]],
        "common": mu_w,
        "feed": feed,
        "c1_w": c1_w, "c1_h": c1_h,
        "beta_c1": beta_c1,
        "witness_tpd": witness_tpd,
    }


if __name__ == "__main__":
    with localcontext() as context:
        context.prec = 80
        values = anchors()

        cpp = Path(__file__).with_name(
            "sw92_phase_assigned_three_phase_test.cpp").read_text(encoding="utf-8")
        block = cpp.split("constexpr ThreePhaseGolden golden{", 1)[1].split("};", 1)[0]
        literals = [D(value) for value in re.findall(
            r"([-+]?[0-9]+\.[0-9]*(?:[eE][-+]?[0-9]+)?)L", block)]
        expected = (values["w"] + values["h0"] + values["h1"] +
                    values["fractions"] + values["z"] + values["common"] +
                    values["feed"] + values["c1_w"] + values["c1_h"] +
                    [values["beta_c1"], values["witness_tpd"]])
        if len(literals) != len(expected):
            raise AssertionError(
                f"C2b.1 golden count changed: {len(literals)} != {len(expected)}")
        for literal, computed in zip(literals, expected):
            if abs(literal - computed) > D("3e-30") * max(D(1), abs(computed)):
                raise AssertionError(
                    f"C2b.1 anchor mismatch: {literal} vs {computed}")

        print("synthetic structural W=", *[format(v, ".35g") for v in values["w"]])
        print("synthetic structural H0=", *[format(v, ".35g") for v in values["h0"]])
        print("synthetic structural H1=", *[format(v, ".35g") for v in values["h1"]])
        print("three-phase fractions=", *[format(v, ".35g") for v in values["fractions"]])
        print("selected Z=", *[format(v, ".35g") for v in values["z"]])
        print("C1->H1 TPD=", format(values["witness_tpd"], ".35g"))
        print("Independent Decimal(80) Profile-C C2b.1 structural reference passed")
