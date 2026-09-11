"""Independent Decimal(80) physical SW92 three-phase oracle.

Primary topology/data source: Mortezazadeh & Rasaei, Fluid Phase Equilibria 450
(2017) 160-174, DOI 10.1016/j.fluid.2017.07.007. The selected gas-condensate
Sample 6 state is P=10 MPa, T=350 K and fresh water. Figure 7 places that state
inside the water+oil+gas region; Tables 4-5 provide pure properties and feed and
state that non-water BIPs are zero. Appendix A supplies the NA water-hydrocarbon
constants used for this physical reproduction, including 0.5 for normal
hydrocarbons heavier than C4.

The numerics intentionally follow MPMC_HNU's exact
SW92/corrected-original/PR76-base profile: corrected original SW92 water alpha
and AQ BIP correlation, plus the original PR76 quadratic kappa for every
non-water component. Therefore the 2017 publication is a physical
source/topology oracle, not a source of copied numerical phase-composition
golden values.

Stdlib only; no production imports.
"""

from decimal import Decimal as D, localcontext
import math
from pathlib import Path
import re

R = D("8.31446261815324")
PRESSURE = D("1e7")
TEMPERATURE = D("350")
MOLALITY = D(0)
NAMES = ("H2O", "C1", "C2", "C3", "C4", "C5", "C6", "C7+")
N = len(NAMES)

SPEC = {
    "H2O": (D("647.30"), D("22.048e6"), D("0.344")),
    "C1": (D("190.60"), D("4.6042e6"), D("0.013")),
    "C2": (D("305.43"), D("4.8839e6"), D("0.0986")),
    "C3": (D("369.80"), D("4.2455e6"), D("0.1524")),
    "C4": (D("419.50"), D("3.747e6"), D("0.1956")),
    "C5": (D("465.90"), D("3.3589e6"), D("0.2413")),
    "C6": (D("507.50"), D("3.0104e6"), D("0.2990")),
    "C7+": (D("655.04"), D("2.2305e6"), D("0.50879")),
}
NA_WATER_KIJ = {
    "C1": D("0.4850"), "C2": D("0.4920"), "C3": D("0.5525"),
    "C4": D("0.5091"), "C5": D("0.5"), "C6": D("0.5"),
    "C7+": D("0.5"),
}

# Table-5 values are rounded and sum to 0.9999. Normalization here is an
# explicit source-data conversion for the regression fixture, not a production
# flash-input repair rule.
PRINTED_FEED = [D("0.5085"), D("0.4249"), D("0.0214"), D("0.0112"),
                D("0.0085"), D("0.0040"), D("0.0029"), D("0.0185")]
PRINTED_SUM = sum(PRINTED_FEED)
FEED = [value / PRINTED_SUM for value in PRINTED_FEED]


def powd(value: D, exponent: D) -> D:
    if value == 0:
        if exponent > 0:
            return D(0)
        raise ValueError("zero base outside reference power domain")
    return (exponent * value.ln()).exp()


def water_alpha() -> D:
    tr = TEMPERATURE / SPEC["H2O"][0]
    # Corrected-original SW92 Eq.(9): the 0.0034 term is OUTSIDE the
    # 0.4530 multiplier. Keep this parenthesization explicit because a prior
    # regression draft caught exactly this transcription hazard.
    q = (D(1) + D("0.4530") *
         (D(1) - tr * (D(1) - D("0.0103") * powd(MOLALITY, D("1.1")))) +
         D("0.0034") * (D(1) / (tr ** 3) - D(1)))
    return q * q


def standard_alpha(name: str) -> D:
    tc, _, omega = SPEC[name]
    # Exact MPMC_HNU PR76-base convention, deliberately not the later
    # piecewise kappa printed in Mortezazadeh-Rasaei Appendix A.
    kappa = D("0.37464") + D("1.54226") * omega - D("0.26992") * omega * omega
    return (D(1) + kappa * (D(1) - (TEMPERATURE / tc).sqrt())) ** 2


def pure_ab(name: str):
    tc, pc, _ = SPEC[name]
    alpha = water_alpha() if name == "H2O" else standard_alpha(name)
    return (D("0.45724") * R * R * tc * tc / pc * alpha,
            D("0.07780") * R * tc / pc)


def aqueous_water_kij(name: str) -> D:
    tc, _, omega = SPEC[name]
    tr = TEMPERATURE / tc
    a0 = D("1.1120") - D("1.7369") * powd(omega, D("-0.1"))
    a1 = D("1.1001") + D("0.8360") * omega
    a2 = -D("0.15742") - D("1.0988") * omega
    return (a0 * (D(1) + D("0.017407") * MOLALITY) +
            a1 * tr * (D(1) + D("0.033516") * MOLALITY) +
            a2 * tr * tr * (D(1) + D("0.011478") * MOLALITY))


def cubic_roots(a: D, b: D):
    c2 = b - D(1)
    c1 = a - D(3) * b * b - D(2) * b
    c0 = b ** 3 + b * b - a * b
    q = (D(3) * c1 - c2 * c2) / D(9)
    r = (D(9) * c2 * c1 - D(27) * c0 - D(2) * c2 ** 3) / D(54)
    discriminant = q ** 3 + r * r

    def cbrt(value: float) -> float:
        return math.copysign(abs(value) ** (1.0 / 3.0), value)

    if discriminant >= 0:
        sd = math.sqrt(float(discriminant))
        seeds = [cbrt(float(r) + sd) + cbrt(float(r) - sd) - float(c2) / 3.0]
    else:
        argument = float(r / (-q ** 3).sqrt())
        theta = math.acos(max(-1.0, min(1.0, argument)))
        seeds = [2.0 * math.sqrt(float(-q)) *
                 math.cos((theta + 2.0 * k * math.pi) / 3.0) - float(c2) / 3.0
                 for k in range(3)]

    roots = []
    for seed in seeds:
        z = D(str(seed))
        for _ in range(160):
            residual = ((z + c2) * z + c1) * z + c0
            derivative = (D(3) * z + D(2) * c2) * z + c1
            next_z = z - residual / derivative
            if abs(next_z - z) < D("1e-68"):
                z = next_z
                break
            z = next_z
        residual = ((z + c2) * z + c1) * z + c0
        scale = D(1) + abs(z) ** 3 + abs(c2 * z * z) + abs(c1 * z) + abs(c0)
        if abs(residual) > D("1e-66") * scale:
            raise ArithmeticError("physical Sample-6 cubic root did not converge")
        if z > b and all(abs(z - old) > D("1e-42") for old in roots):
            roots.append(z)
    return sorted(roots)


PURE = [pure_ab(name) for name in NAMES]


def minimum_phase(composition, family: str):
    ai = [entry[0] for entry in PURE]
    bi = [entry[1] for entry in PURE]
    kij = [[D(0) for _ in range(N)] for _ in range(N)]
    for j, name in enumerate(NAMES[1:], 1):
        value = aqueous_water_kij(name) if family == "AQ" else NA_WATER_KIJ[name]
        kij[0][j] = kij[j][0] = value
    # All non-water BIPs are zero by the 2017 Sample-6 source contract.
    aij = [[(ai[i] * ai[j]).sqrt() * (D(1) - kij[i][j])
            for j in range(N)] for i in range(N)]
    rows = [sum(composition[j] * aij[i][j] for j in range(N))
            for i in range(N)]
    attraction = sum(composition[i] * rows[i] for i in range(N))
    covolume = sum(composition[i] * bi[i] for i in range(N))
    aa = attraction * PRESSURE / (R * TEMPERATURE) ** 2
    bb = covolume * PRESSURE / (R * TEMPERATURE)
    c2 = bb - D(1)
    c1 = aa - D(3) * bb * bb - D(2) * bb
    sqrt2 = D(2).sqrt()

    candidates = []
    for z in cubic_roots(aa, bb):
        derivative = (D(3) * z + D(2) * c2) * z + c1
        if derivative <= 0:
            continue
        log_ratio = ((z + (D(1) + sqrt2) * bb) /
                     (z + (D(1) - sqrt2) * bb)).ln()
        ln_phi = []
        for i in range(N):
            ratio = bi[i] / covolume
            expression = D(2) * rows[i] / attraction - ratio
            ln_phi.append(ratio * (z - D(1)) - (z - bb).ln() -
                          aa / (D(2) * sqrt2 * bb) * expression * log_ratio)
        root_gibbs = sum(composition[i] * ln_phi[i] for i in range(N))
        candidates.append((root_gibbs, z, ln_phi))
    if not candidates:
        raise ArithmeticError("physical Sample-6 has no mechanically admissible root")
    return min(candidates, key=lambda item: item[0])


def rr_compositions(log_k_h0, log_k_h1, beta_h0: D, beta_h1: D):
    if beta_h0 <= 0 or beta_h1 <= 0 or beta_h0 + beta_h1 >= 1:
        raise ArithmeticError("phase fractions outside three-phase interior")
    k0 = [value.exp() for value in log_k_h0]
    k1 = [value.exp() for value in log_k_h1]
    beta_w = D(1) - beta_h0 - beta_h1
    w, h0, h1 = [], [], []
    f0 = D(0)
    f1 = D(0)
    for zi, k0i, k1i in zip(FEED, k0, k1):
        denominator = beta_w + beta_h0 * k0i + beta_h1 * k1i
        if denominator <= 0:
            raise ArithmeticError("generalized RR denominator outside domain")
        wi = zi / denominator
        h0i = k0i * wi
        h1i = k1i * wi
        if min(wi, h0i, h1i) <= 0:
            raise ArithmeticError("three-phase composition outside positive interior")
        w.append(wi)
        h0.append(h0i)
        h1.append(h1i)
        f0 += h0i - wi
        f1 += h1i - wi
    return w, h0, h1, f0, f1


def residual(unknown):
    w, h0, h1, f0, f1 = rr_compositions(
        unknown[:N], unknown[N:2 * N], unknown[-2], unknown[-1])
    pw = minimum_phase(w, "AQ")
    p0 = minimum_phase(h0, "NA")
    p1 = minimum_phase(h1, "NA")
    mu_w = [w[i].ln() + pw[2][i] for i in range(N)]
    mu_0 = [h0[i].ln() + p0[2][i] for i in range(N)]
    mu_1 = [h1[i].ln() + p1[2][i] for i in range(N)]
    return ([mu_w[i] - mu_0[i] for i in range(N)] +
            [mu_w[i] - mu_1[i] for i in range(N)] + [f0, f1])


def solve_linear(matrix, rhs):
    n = len(rhs)
    work = [list(matrix[i]) + [rhs[i]] for i in range(n)]
    for column in range(n):
        pivot = max(range(column, n), key=lambda row: abs(work[row][column]))
        if work[pivot][column] == 0:
            raise ArithmeticError("singular physical Sample-6 Newton matrix")
        work[column], work[pivot] = work[pivot], work[column]
        scale = work[column][column]
        for j in range(column, n + 1):
            work[column][j] /= scale
        for row in range(n):
            if row == column:
                continue
            factor = work[row][column]
            if factor == 0:
                continue
            for j in range(column, n + 1):
                work[row][j] -= factor * work[column][j]
    return [work[i][n] for i in range(n)]


def normalized(values):
    total = sum(values)
    return [value / total for value in values]


def log_ratio(numerator, denominator):
    return [(numerator[i] / denominator[i]).ln() for i in range(N)]


def solve_flash():
    # Deliberately rounded independent seeds, not the C++ golden literals.
    w0 = normalized([D(".999143"), D(".000801499"), D(".0000408836"),
                     D(".0000114941"), D(".00000311120"), D("3.75863e-7"),
                     D("2.79473e-8"), D("7.45412e-12")])
    h00 = normalized([D(".00317609"), D(".305047"), D(".0397125"),
                       D(".0390992"), D(".0497509"), D(".0375450"),
                       D(".0402421"), D(".485427")])
    h10 = normalized([D(".00552981"), D(".904284"), D(".0435996"),
                       D(".0213672"), D(".0146194"), D(".00575261"),
                       D(".00313126"), D(".00171569")])
    unknown = (log_ratio(h00, w0) + log_ratio(h10, w0) +
               [D(".036499"), D(".457160")])
    step = D("1e-20")

    for _ in range(12):
        current = residual(unknown)
        norm = max(abs(value) for value in current)
        if norm < D("1e-52"):
            break
        n = len(unknown)
        jacobian = [[D(0) for _ in range(n)] for _ in range(n)]
        for column in range(n):
            plus = list(unknown)
            minus = list(unknown)
            plus[column] += step
            minus[column] -= step
            rp = residual(plus)
            rm = residual(minus)
            for row in range(n):
                jacobian[row][column] = (rp[row] - rm[row]) / (D(2) * step)
        delta = solve_linear(jacobian, [-value for value in current])
        damping = D(1)
        accepted = False
        while damping > D("1e-16"):
            candidate = [unknown[i] + damping * delta[i] for i in range(n)]
            try:
                next_norm = max(abs(value) for value in residual(candidate))
            except (ArithmeticError, ValueError, OverflowError):
                next_norm = D("Infinity")
            if next_norm < norm:
                unknown = candidate
                accepted = True
                break
            damping /= D(2)
        if not accepted:
            raise ArithmeticError("physical Sample-6 Newton line search failed")
    else:
        raise ArithmeticError("physical Sample-6 Newton iteration limit")

    final_residual = residual(unknown)
    if max(abs(value) for value in final_residual) > D("1e-48"):
        raise ArithmeticError("physical Sample-6 residual did not converge")

    w, h0, h1, f0, f1 = rr_compositions(
        unknown[:N], unknown[N:2 * N], unknown[-2], unknown[-1])
    pw = minimum_phase(w, "AQ")
    p0 = minimum_phase(h0, "NA")
    p1 = minimum_phase(h1, "NA")
    fractions = [D(1) - unknown[-2] - unknown[-1], unknown[-2], unknown[-1]]
    common = [w[i].ln() + pw[2][i] for i in range(N)]

    if not all(value > 0 for value in fractions):
        raise AssertionError("physical Sample-6 lost an interior phase")
    if not w[0] > max(h0[0], h1[0]):
        raise AssertionError("physical Sample-6 water role ordering changed")
    if not p0[1] < p1[1]:
        raise AssertionError("physical Sample-6 H slots are not in canonical lower-Z order")
    recovered = [fractions[0] * w[i] + fractions[1] * h0[i] +
                 fractions[2] * h1[i] for i in range(N)]
    if max(abs(recovered[i] - FEED[i]) for i in range(N)) > D("1e-48"):
        raise AssertionError("physical Sample-6 material balance changed")
    if max(abs(value) for value in (f0, f1)) > D("1e-48"):
        raise AssertionError("physical Sample-6 generalized RR residual changed")

    return {
        "feed": FEED, "w": w, "h0": h0, "h1": h1,
        "fractions": fractions, "z": [pw[1], p0[1], p1[1]],
        "common": common,
    }


if __name__ == "__main__":
    with localcontext() as context:
        context.prec = 80
        values = solve_flash()

        header = Path(__file__).with_name("physical_sample6.hpp").read_text(encoding="utf-8")
        block = header.split("inline constexpr PhysicalGolden golden{", 1)[1].split("};", 1)[0]
        literals = [D(value) for value in re.findall(
            r"([-+]?[0-9]+(?:\.[0-9]*)?(?:[eE][-+]?[0-9]+)?)L", block)]
        expected = (values["feed"] + values["w"] + values["h0"] + values["h1"] +
                    values["fractions"] + values["z"] + values["common"])
        if len(literals) != len(expected):
            raise AssertionError(
                f"physical Sample-6 golden count changed: {len(literals)} != {len(expected)}")
        for literal, computed in zip(literals, expected):
            if abs(literal - computed) > D("3e-27") * max(D(1), abs(computed)):
                raise AssertionError(
                    f"physical Sample-6 golden mismatch: {literal} vs {computed}")

        print("MR2017 Sample-6 explicit printed-feed sum=", PRINTED_SUM)
        print("physical W fraction=", format(values["fractions"][0], ".35g"))
        print("physical H0 fraction=", format(values["fractions"][1], ".35g"))
        print("physical H1 fraction=", format(values["fractions"][2], ".35g"))
        print("physical Z=", *[format(value, ".35g") for value in values["z"]])
        print("Independent Decimal(80) physical SW92 Sample-6 three-phase oracle passed")
