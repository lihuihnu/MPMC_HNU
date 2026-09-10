"""Independent PR76 ternary VLE references using Decimal(80), stdlib only.

Scientific data:
* methane/ethane/propane Tc, Pc and omega: Deiters & Bell, AIChE J.,
  DOI 10.1002/aic.16730, Table 1.
* Peng & Robinson, Ind. Eng. Chem. Fundam. 15 (1976),
  DOI 10.1021/i160057a011, states that no interaction coefficients were
  used for its methane/ethane/propane ternary example. This regression
  therefore records all three off-diagonal kij explicitly as zero.

The reference fixes T=220 K, p=4 MPa and x_CH4=0.5, then solves the three
component fugacity-equality equations directly in phase compositions.
It does not call production TPD, logK-SSI, Rachford-Rice, or C++ output.
Three feed states are formed only after convergence by the lever rule.
This is a model/software regression, not experimental validation.
"""

from decimal import Decimal as D, localcontext
import math
from pathlib import Path
import re

R = D("8.31446261815324")
P = D("4000000")
T = D("220")
SPECS = (
    (D("190.555"), D("4595000"), D("0.0")),
    (D("305.4"), D("4880000"), D("0.099")),
    (D("369.825"), D("4248000"), D("0.15308")),
)
KIJ = (
    (D(0), D(0), D(0)),
    (D(0), D(0), D(0)),
    (D(0), D(0), D(0)),
)
BETAS = (D(".25"), D(".50"), D(".75"))


def roots(a, b):
    c2, c1, c0 = b - 1, a - 3 * b * b - 2 * b, b * b * b + b * b - a * b
    q = (3 * c1 - c2 * c2) / 9
    r = (9 * c2 * c1 - 27 * c0 - 2 * c2**3) / 54
    disc = q**3 + r * r
    if disc >= 0:
        sd = math.sqrt(float(disc))

        def cbrt(x):
            return math.copysign(abs(x) ** (1 / 3), x)

        seeds = [cbrt(float(r) + sd) + cbrt(float(r) - sd) - float(c2) / 3]
    else:
        theta = math.acos(float(r / (-q**3).sqrt()))
        seeds = [
            2 * math.sqrt(float(-q)) * math.cos((theta + 2 * k * math.pi) / 3)
            - float(c2) / 3
            for k in range(3)
        ]
    found = []
    for seed in seeds:
        z = D(str(seed))
        for _ in range(40):
            residual = ((z + c2) * z + c1) * z + c0
            slope = (3 * z + 2 * c2) * z + c1
            z -= residual / slope
        residual = ((z + c2) * z + c1) * z + c0
        scale = 1 + abs(z) ** 3 + abs(c2 * z * z) + abs(c1 * z) + abs(c0)
        if abs(residual) > D("1e-70") * scale:
            raise ArithmeticError("reference root did not converge")
        if z > b and (3 * z + 2 * c2) * z + c1 > 0:
            if not any(abs(z - old) < D("1e-50") for old in found):
                found.append(z)
    if not found:
        raise ArithmeticError("no admissible reference root")
    return sorted(found)


def phase(composition, vapor):
    ai, bi = [], []
    for tc, pc, omega in SPECS:
        kappa = D(".37464") + D("1.54226") * omega - D(".26992") * omega**2
        ai.append(
            D(".45724")
            * R
            * R
            * tc
            * tc
            / pc
            * (1 + kappa * (1 - (T / tc).sqrt())) ** 2
        )
        bi.append(D(".07780") * R * tc / pc)

    count = len(composition)
    aij = [
        [
            (ai[i] * ai[j]).sqrt() * (1 - KIJ[i][j])
            for j in range(count)
        ]
        for i in range(count)
    ]
    row_sums = [
        sum((composition[j] * aij[i][j] for j in range(count)), D(0))
        for i in range(count)
    ]
    a_mix = sum((composition[i] * row_sums[i] for i in range(count)), D(0))
    b_mix = sum((composition[i] * bi[i] for i in range(count)), D(0))
    aa = a_mix * P / (R * T) ** 2
    bb = b_mix * P / (R * T)
    root = roots(aa, bb)[-1 if vapor else 0]
    sqrt2 = D(2).sqrt()
    log_ratio = (
        (root + (1 + sqrt2) * bb) / (root + (1 - sqrt2) * bb)
    ).ln()
    ln_phi = [
        bi[i] / b_mix * (root - 1)
        - (root - bb).ln()
        - aa
        / (2 * sqrt2 * bb)
        * (2 * row_sums[i] / a_mix - bi[i] / b_mix)
        * log_ratio
        for i in range(count)
    ]
    chemical = [composition[i].ln() + ln_phi[i] for i in range(count)]
    return root, chemical


def residual(state):
    x_ethane, y_methane, y_ethane = state
    liquid = [D(".5"), x_ethane, D(".5") - x_ethane]
    vapor = [y_methane, y_ethane, D(1) - y_methane - y_ethane]
    _, liquid_mu = phase(liquid, False)
    _, vapor_mu = phase(vapor, True)
    return [liquid_mu[i] - vapor_mu[i] for i in range(3)]


def solve_linear(matrix, rhs):
    count = len(rhs)
    augmented = [row[:] + [rhs[i]] for i, row in enumerate(matrix)]
    for column in range(count):
        pivot = max(range(column, count), key=lambda row: abs(augmented[row][column]))
        if pivot != column:
            augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        scale = augmented[column][column]
        if scale == 0:
            raise ArithmeticError("singular reference Jacobian")
        for j in range(column, count + 1):
            augmented[column][j] /= scale
        for row in range(count):
            if row == column:
                continue
            factor = augmented[row][column]
            if factor == 0:
                continue
            for j in range(column, count + 1):
                augmented[row][j] -= factor * augmented[column][j]
    return [augmented[i][count] for i in range(count)]


def valid_state(state):
    x_ethane, y_methane, y_ethane = state
    return (
        0 < x_ethane < D(".5")
        and 0 < y_methane < 1
        and 0 < y_ethane < 1
        and y_methane + y_ethane < 1
    )


def tie_line():
    state = [D(".33"), D(".895"), D(".093")]
    step = D("1e-24")
    for _ in range(40):
        values = residual(state)
        norm = max(abs(value) for value in values)
        if norm < D("1e-60"):
            break
        columns = []
        for column in range(3):
            plus, minus = state[:], state[:]
            plus[column] += step
            minus[column] -= step
            upper, lower = residual(plus), residual(minus)
            columns.append(
                [(upper[row] - lower[row]) / (2 * step) for row in range(3)]
            )
        jacobian = [[columns[column][row] for column in range(3)] for row in range(3)]
        direction = solve_linear(jacobian, [-value for value in values])
        alpha = D(1)
        for _ in range(50):
            candidate = [
                state[i] + alpha * direction[i]
                for i in range(3)
            ]
            if valid_state(candidate):
                if max(abs(value) for value in residual(candidate)) < norm:
                    state = candidate
                    break
            alpha /= 2
        else:
            raise ArithmeticError("reference Newton backtracking failed")
    else:
        raise ArithmeticError("reference equilibrium did not converge")

    x_ethane, y_methane, y_ethane = state
    liquid = [D(".5"), x_ethane, D(".5") - x_ethane]
    vapor = [y_methane, y_ethane, D(1) - y_methane - y_ethane]
    liquid_z, liquid_mu = phase(liquid, False)
    vapor_z, vapor_mu = phase(vapor, True)
    chemical_residual = max(
        abs(liquid_mu[i] - vapor_mu[i]) for i in range(3)
    )
    if chemical_residual >= D("1e-60"):
        raise ArithmeticError("reference chemical equilibrium residual too large")
    return liquid, vapor, liquid_z, vapor_z


def anchors():
    liquid, vapor, liquid_z, vapor_z = tie_line()
    rows = []
    for beta in BETAS:
        feed = [
            (1 - beta) * liquid[i] + beta * vapor[i]
            for i in range(3)
        ]
        mass_residual = max(
            abs((1 - beta) * liquid[i] + beta * vapor[i] - feed[i])
            for i in range(3)
        )
        if mass_residual != 0:
            raise ArithmeticError("reference lever rule lost exact Decimal balance")
        rows.append((feed, liquid, vapor, beta, liquid_z, vapor_z))
    return rows


if __name__ == "__main__":
    with localcontext() as context:
        context.prec = 80
        computed = anchors()
        for feed, liquid, vapor, beta, liquid_z, vapor_z in computed:
            values = [*feed, *liquid, *vapor, beta, liquid_z, vapor_z]
            print(", ".join(format(value, ".40g") + "L" for value in values))

        header = Path(__file__).with_name("ternary_references.hpp").read_text(
            encoding="utf-8"
        )
        block = header.split(
            "inline constexpr std::array<TernaryReference, 3> ternary_references{{", 1
        )[1].split("}};", 1)[0]
        literals = re.findall(
            r"([-+]?[0-9]+(?:\.[0-9]*)?(?:[eE][-+]?[0-9]+)?)L", block
        )
        flat = [
            value
            for feed, liquid, vapor, beta, liquid_z, vapor_z in computed
            for value in [*feed, *liquid, *vapor, beta, liquid_z, vapor_z]
        ]
        if len(literals) != len(flat):
            raise AssertionError((len(literals), len(flat)))
        for literal, value in zip(literals, flat):
            tolerance = D("1e-35") * max(D(1), abs(value))
            if abs(D(literal) - value) >= tolerance:
                raise AssertionError((literal, value))
        print("Independent Decimal(80) ternary split anchors: 36/36 checked")
