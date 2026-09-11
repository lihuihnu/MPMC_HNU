"""Independent Decimal(80) SW92 n-butane/H2O three-phase-line oracle.

Stdlib only; no production imports.

Thermodynamic model data:
- Soreide & Whitson 1992, Fluid Phase Equilibria 77, 217-240,
  DOI 10.1016/0378-3812(92)85105-H.
- n-butane SW92 Table 3: Tc=425.2 K, Pc=38.0 bar, omega=0.1931.
- n-butane/water SW92 Table 5 non-aqueous BIP: kij=0.5091.
- aqueous BIP uses the corrected-original SW92 hydrocarbon/water correlation.

Experimental comparison:
H. H. Reamer, R. H. Olds, B. H. Sage and W. N. Lacey,
"Phase Equilibria in Hydrocarbon Systems. n-Butane-Water System in Three-Phase
Region", Ind. Eng. Chem. 36 (1944) 381-383,
DOI 10.1021/ie50412a024.

The binary system has one degree of freedom on the three-phase line. At each
reported temperature this oracle independently solves pressure plus the AQ and
two NA coexistence compositions from four chemical-potential equalities. Binary
three-phase phase fractions are intentionally not validated because they are not
unique at fixed T/P/feed on a three-phase tie-line.

The experimental envelope below is a model-regression criterion, not a solver
convergence tolerance and not a fitted parameter target. No production formula,
BIP or numerical tolerance is adjusted to satisfy it.
"""

from decimal import Decimal as D, localcontext
from pathlib import Path
import math
import re

R = D("8.31446261815324")
TC_B = D("425.2")
PC_B = D("38.0e5")
OMEGA_B = D("0.1931")
TC_W = D("647.3")
PC_W = D("221.2e5")
OMEGA_W = D("0.3434")
NA_KIJ = D("0.5091")

EXPERIMENT = [
    (D("310.93"), D("3.62"), D("0.9995"), D("0.9833")),
    (D("327.59"), D("5.71"), D("0.9990"), D("0.9755")),
    (D("344.26"), D("8.65"), D("0.9979"), D("0.9662")),
    (D("360.93"), D("12.60"), D("0.9957"), D("0.9561")),
    (D("377.59"), D("17.88"), D("0.9915"), D("0.9459")),
    (D("394.26"), D("24.82"), D("0.9843"), D("0.9361")),
    (D("410.93"), D("33.85"), D("0.9732"), D("0.9292")),
    (D("416.48"), D("37.40"), D("0.9683"), D("0.9302")),
]

SEEDS = [
    (D("3.66115e5"), D("0.00006145"), D("0.999289"), D("0.983471")),
    (D("5.77024e5"), D("0.00006854"), D("0.998486"), D("0.975183")),
    (D("8.72813e5"), D("0.00008005"), D("0.997002"), D("0.965178")),
    (D("1.27636e6"), D("0.00009729"), D("0.994377"), D("0.953957")),
    (D("1.81518e6"), D("0.00012232"), D("0.989856"), D("0.942406")),
    (D("2.52408e6"), D("0.00015821"), D("0.982020"), D("0.932037")),
    (D("3.44540e6"), D("0.00020927"), D("0.967373"), D("0.926363")),
    (D("3.80820e6"), D("0.00023061"), D("0.958961"), D("0.927508")),
]


def powd(value: D, exponent: D) -> D:
    if value == 0:
        if exponent > 0:
            return D(0)
        raise ValueError("zero base outside Decimal power domain")
    return (exponent * value.ln()).exp()


def water_alpha(temperature: D) -> D:
    tr = temperature / TC_W
    root = (D(1) + D("0.4530") * (D(1) - tr) +
            D("0.0034") * (tr ** D(-3) - D(1)))
    return root * root


def standard_alpha(temperature: D) -> D:
    kappa = (D("0.37464") + D("1.54226") * OMEGA_B -
             D("0.26992") * OMEGA_B * OMEGA_B)
    return (D(1) + kappa * (D(1) - (temperature / TC_B).sqrt())) ** 2


def pure_ab(temperature: D):
    butane = (D("0.45724") * R * R * TC_B * TC_B / PC_B *
              standard_alpha(temperature),
              D("0.07780") * R * TC_B / PC_B)
    water = (D("0.45724") * R * R * TC_W * TC_W / PC_W *
             water_alpha(temperature),
             D("0.07780") * R * TC_W / PC_W)
    return butane, water


def aqueous_kij(temperature: D) -> D:
    tr = temperature / TC_B
    a0 = D("1.1120") - D("1.7369") * powd(OMEGA_B, D("-0.1"))
    a1 = D("1.1001") + D("0.8360") * OMEGA_B
    a2 = -D("0.15742") - D("1.0988") * OMEGA_B
    return a0 + a1 * tr + a2 * tr * tr


def cubic_roots(a: D, b: D):
    c2 = b - D(1)
    c1 = a - D(3) * b * b - D(2) * b
    c0 = b ** 3 + b * b - a * b
    q = (D(3) * c1 - c2 * c2) / D(9)
    r = (D(9) * c2 * c1 - D(27) * c0 - D(2) * c2 ** 3) / D(54)
    discriminant = q ** 3 + r * r
    if discriminant >= 0:
        sd = math.sqrt(float(discriminant))

        def cbrt(value: float) -> float:
            return math.copysign(abs(value) ** (1.0 / 3.0), value)

        seeds = [cbrt(float(r) + sd) + cbrt(float(r) - sd) - float(c2) / 3.0]
    else:
        theta = math.acos(float(r / (-q ** 3).sqrt()))
        seeds = [2.0 * math.sqrt(float(-q)) *
                 math.cos((theta + 2.0 * k * math.pi) / 3.0) - float(c2) / 3.0
                 for k in range(3)]

    roots = []
    for seed in seeds:
        z = D(str(seed))
        for _ in range(160):
            residual = ((z + c2) * z + c1) * z + c0
            derivative = (D(3) * z + D(2) * c2) * z + c1
            if derivative == 0:
                raise ArithmeticError("reference cubic root derivative vanished")
            next_z = z - residual / derivative
            if abs(next_z - z) < D("1e-68"):
                z = next_z
                break
            z = next_z
        residual = ((z + c2) * z + c1) * z + c0
        scale = D(1) + abs(z) ** 3 + abs(c2 * z * z) + abs(c1 * z) + abs(c0)
        if abs(residual) > D("1e-64") * scale:
            raise ArithmeticError("reference cubic root did not converge")
        if z > b and all(abs(z - old) > D("1e-45") for old in roots):
            roots.append(z)
    return sorted(roots)


def minimum_phase(temperature: D, pressure: D, x_butane: D, family: str):
    if not D(0) < x_butane < D(1):
        raise ArithmeticError("reference composition outside open simplex")
    x = [x_butane, D(1) - x_butane]
    pure = pure_ab(temperature)
    ai = [pure[0][0], pure[1][0]]
    bi = [pure[0][1], pure[1][1]]
    kij = aqueous_kij(temperature) if family == "AQ" else NA_KIJ
    aij = [[D(0), D(0)], [D(0), D(0)]]
    for i in range(2):
        for j in range(2):
            pair = D(0) if i == j else kij
            aij[i][j] = (ai[i] * ai[j]).sqrt() * (D(1) - pair)
    rows = [sum(x[j] * aij[i][j] for j in range(2)) for i in range(2)]
    attraction = sum(x[i] * rows[i] for i in range(2))
    covolume = sum(x[i] * bi[i] for i in range(2))
    aa = attraction * pressure / (R * temperature) ** 2
    bb = covolume * pressure / (R * temperature)
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
        for i in range(2):
            ratio = bi[i] / covolume
            expression = D(2) * rows[i] / attraction - ratio
            ln_phi.append(ratio * (z - D(1)) - (z - bb).ln() -
                          aa / (D(2) * sqrt2 * bb) * expression * log_ratio)
        root_gibbs = sum(x[i] * ln_phi[i] for i in range(2))
        candidates.append((root_gibbs, z, ln_phi))
    if not candidates:
        raise ArithmeticError("no mechanically admissible minimum-Gibbs root")
    return min(candidates, key=lambda item: item[0])


def chemical_potentials(temperature: D, pressure: D,
                        x_butane: D, family: str):
    phase = minimum_phase(temperature, pressure, x_butane, family)
    return ([x_butane.ln() + phase[2][0],
             (D(1) - x_butane).ln() + phase[2][1]], phase)


def residuals(temperature: D, unknown):
    pressure, x_w, x_h0, x_h1 = unknown
    if pressure <= 0 or not all(D(0) < value < D(1)
                                for value in (x_w, x_h0, x_h1)):
        raise ArithmeticError("invalid three-phase-line iterate")
    mu_w, _ = chemical_potentials(temperature, pressure, x_w, "AQ")
    mu_h0, _ = chemical_potentials(temperature, pressure, x_h0, "NA")
    mu_h1, _ = chemical_potentials(temperature, pressure, x_h1, "NA")
    return [mu_w[0] - mu_h0[0], mu_w[1] - mu_h0[1],
            mu_w[0] - mu_h1[0], mu_w[1] - mu_h1[1]]


def solve_linear(matrix, rhs):
    n = len(rhs)
    work = [list(matrix[i]) + [rhs[i]] for i in range(n)]
    for column in range(n):
        pivot = max(range(column, n), key=lambda row: abs(work[row][column]))
        if work[pivot][column] == 0:
            raise ArithmeticError("singular Decimal Newton matrix")
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


def solve_case(temperature: D, seed):
    unknown = list(seed)
    for _ in range(60):
        residual = residuals(temperature, unknown)
        norm = max(abs(value) for value in residual)
        if norm < D("1e-50"):
            break
        jacobian = [[D(0) for _ in range(4)] for _ in range(4)]
        for column in range(4):
            relative = abs(unknown[column]) * D("1e-18")
            floor = D("1e-8") if column == 0 else D("1e-28")
            step = max(floor, relative)
            plus = list(unknown)
            minus = list(unknown)
            plus[column] += step
            minus[column] -= step
            rp = residuals(temperature, plus)
            rm = residuals(temperature, minus)
            for row in range(4):
                jacobian[row][column] = (rp[row] - rm[row]) / (D(2) * step)
        delta = solve_linear(jacobian, [-value for value in residual])
        damping = D(1)
        accepted = False
        while damping > D("1e-18"):
            candidate = [unknown[i] + damping * delta[i] for i in range(4)]
            try:
                next_norm = max(abs(value) for value in residuals(temperature, candidate))
            except (ArithmeticError, ValueError):
                next_norm = D("Infinity")
            if next_norm < norm:
                unknown = candidate
                accepted = True
                break
            damping /= D(2)
        if not accepted:
            raise ArithmeticError("triple-line Decimal Newton line search failed")
    else:
        raise ArithmeticError("triple-line Decimal Newton iteration limit")
    if max(abs(value) for value in residuals(temperature, unknown)) > D("1e-46"):
        raise ArithmeticError("triple-line Decimal residual did not converge")

    pressure, x_w, x_h0, x_h1 = unknown
    _, phase_w = chemical_potentials(temperature, pressure, x_w, "AQ")
    _, phase_h0 = chemical_potentials(temperature, pressure, x_h0, "NA")
    _, phase_h1 = chemical_potentials(temperature, pressure, x_h1, "NA")
    return (temperature, pressure, x_w, x_h0, x_h1,
            phase_w[1], phase_h0[1], phase_h1[1])


def anchors():
    values = []
    max_pressure_relative = D(0)
    max_h0_absolute = D(0)
    max_h1_absolute = D(0)
    for experiment, seed in zip(EXPERIMENT, SEEDS):
        temperature, p_bar_exp, x_h0_exp, x_h1_exp = experiment
        model = solve_case(temperature, seed)
        _, pressure, _, x_h0, x_h1, _, _, _ = model
        p_bar = pressure / D("1e5")
        max_pressure_relative = max(max_pressure_relative,
                                    abs(p_bar - p_bar_exp) / p_bar_exp)
        max_h0_absolute = max(max_h0_absolute, abs(x_h0 - x_h0_exp))
        max_h1_absolute = max(max_h1_absolute, abs(x_h1 - x_h1_exp))
        values.append(model + (p_bar_exp, x_h0_exp, x_h1_exp))

    if max_pressure_relative > D("0.02"):
        raise AssertionError(f"SW92 pressure regression left 2% experimental envelope: {max_pressure_relative}")
    if max_h0_absolute > D("0.01"):
        raise AssertionError(f"SW92 dense-NA composition left 0.01 envelope: {max_h0_absolute}")
    if max_h1_absolute > D("0.005"):
        raise AssertionError(f"SW92 second-NA composition left 0.005 envelope: {max_h1_absolute}")
    return values, max_pressure_relative, max_h0_absolute, max_h1_absolute


if __name__ == "__main__":
    with localcontext() as context:
        context.prec = 80
        values, p_error, h0_error, h1_error = anchors()

        cpp = Path(__file__).with_name(
            "sw92_three_phase_experimental_line_test.cpp").read_text(encoding="utf-8")
        block = cpp.split("constexpr TripleLineGolden golden[] = {", 1)[1].split("};", 1)[0]
        literals = [D(value) for value in re.findall(
            r"([-+]?[0-9]+(?:\.[0-9]*)?(?:[eE][-+]?[0-9]+)?)L", block)]
        expected = [value for row in values for value in row]
        if len(literals) != len(expected):
            raise AssertionError(
                f"n-butane triple-line golden count changed: {len(literals)} != {len(expected)}")
        for literal, computed in zip(literals, expected):
            if abs(literal - computed) > D("5e-26") * max(D(1), abs(computed)):
                raise AssertionError(
                    f"n-butane triple-line anchor mismatch: {literal} vs {computed}")

        for row in values:
            print("T=", format(row[0], ".8g"),
                  "Pbar=", format(row[1] / D("1e5"), ".16g"),
                  "xW=", format(row[2], ".16g"),
                  "xH0=", format(row[3], ".16g"),
                  "xH1=", format(row[4], ".16g"))
        print("max relative pressure error=", format(p_error, ".12g"))
        print("max H0 absolute composition error=", format(h0_error, ".12g"))
        print("max H1 absolute composition error=", format(h1_error, ".12g"))
        print("Independent Decimal(80) SW92 n-butane/H2O experimental triple-line regression passed")
