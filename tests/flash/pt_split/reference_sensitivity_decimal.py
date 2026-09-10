"""Independent Decimal(80) derivatives for accepted interior PR76 PT VLE states.

The reference does NOT differentiate production logK-SSI, Rachford-Rice,
line search, TPD, or C++ output. It solves a different physical coordinate
system directly:

  s = (x_0..x_{N-2}, y_0..y_{N-2}, beta)
  q = (p_Pa, T_K, z_0..z_{N-2})

with N fugacity-equality residuals and N-1 component material balances.
A high-precision central-difference Jacobian of those independent residuals is
then used only in this reference program to solve G_s ds/dq = -G_q. A second,
finite perturb-and-resolve calculation cross-checks the implicit derivatives.
"""

from decimal import Decimal as D, localcontext
import math
from pathlib import Path
import re

R = D("8.31446261815324")

BINARY = {
    "specs": (
        (D("126.2"), D("3390000"), D("0.04")),
        (D("305.4"), D("4880000"), D("0.098")),
    ),
    "kij": ((D(0), D("0.08")), (D("0.08"), D(0))),
    "q": (D("7600000"), D("270"), D("0.30")),
    "initial": (D("0.17"), D("0.48"), D("0.43")),
}

TERNARY = {
    "specs": (
        (D("190.555"), D("4595000"), D("0.0")),
        (D("305.4"), D("4880000"), D("0.099")),
        (D("369.825"), D("4248000"), D("0.15308")),
    ),
    "kij": (
        (D(0), D(0), D(0)),
        (D(0), D(0), D(0)),
        (D(0), D(0), D(0)),
    ),
    "q": (
        D("4000000"),
        D("220"),
        D("0.6971797413795296714538190394514197508909"),
        D("0.2136022423731941810042697872110338951037"),
    ),
    "initial": (D("0.50"), D("0.333"), D("0.894"), D("0.094"), D("0.50")),
}


def roots(a, b):
    c2 = b - 1
    c1 = a - 3 * b * b - 2 * b
    c0 = b * b * b + b * b - a * b
    q = (3 * c1 - c2 * c2) / 9
    r = (9 * c2 * c1 - 27 * c0 - 2 * c2**3) / 54
    disc = q**3 + r * r
    if disc >= 0:
        sd = math.sqrt(float(disc))

        def cbrt(value):
            return math.copysign(abs(value) ** (1 / 3), value)

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
        for _ in range(48):
            residual = ((z + c2) * z + c1) * z + c0
            slope = (3 * z + 2 * c2) * z + c1
            z -= residual / slope
        residual = ((z + c2) * z + c1) * z + c0
        scale = 1 + abs(z) ** 3 + abs(c2 * z * z) + abs(c1 * z) + abs(c0)
        if abs(residual) > D("1e-70") * scale:
            raise ArithmeticError("reference PR root did not converge")
        if z > b and (3 * z + 2 * c2) * z + c1 > 0:
            if not any(abs(z - old) < D("1e-50") for old in found):
                found.append(z)
    if not found:
        raise ArithmeticError("no admissible reference PR root")
    return sorted(found)


def phase(composition, specs, kij, p, t, vapor):
    count = len(composition)
    ai = []
    bi = []
    for tc, pc, omega in specs:
        kappa = D("0.37464") + D("1.54226") * omega - D("0.26992") * omega**2
        ai.append(
            D("0.45724") * R * R * tc * tc / pc
            * (1 + kappa * (1 - (t / tc).sqrt())) ** 2
        )
        bi.append(D("0.07780") * R * tc / pc)
    aij = [
        [
            (ai[i] * ai[j]).sqrt() * (1 - kij[i][j])
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
    aa = a_mix * p / (R * t) ** 2
    bb = b_mix * p / (R * t)
    z = roots(aa, bb)[-1 if vapor else 0]
    sqrt2 = D(2).sqrt()
    log_ratio = ((z + (1 + sqrt2) * bb) / (z + (1 - sqrt2) * bb)).ln()
    ln_phi = [
        bi[i] / b_mix * (z - 1)
        - (z - bb).ln()
        - aa / (2 * sqrt2 * bb)
        * (2 * row_sums[i] / a_mix - bi[i] / b_mix)
        * log_ratio
        for i in range(count)
    ]
    return [composition[i].ln() + ln_phi[i] for i in range(count)]


def compositions(state, count):
    independent = count - 1
    liquid = list(state[:independent])
    vapor = list(state[independent:2 * independent])
    liquid.append(D(1) - sum(liquid, D(0)))
    vapor.append(D(1) - sum(vapor, D(0)))
    beta = state[-1]
    return liquid, vapor, beta


def feed_from_q(q, count):
    feed = list(q[2:])
    if len(feed) != count - 1:
        raise ValueError("reference feed coordinate mismatch")
    feed.append(D(1) - sum(feed, D(0)))
    return feed


def valid_state(state, count):
    liquid, vapor, beta = compositions(state, count)
    return (
        D(0) < beta < D(1)
        and all(D(0) < value < D(1) for value in liquid)
        and all(D(0) < value < D(1) for value in vapor)
    )


def residual(state, q, case):
    specs = case["specs"]
    kij = case["kij"]
    count = len(specs)
    liquid, vapor, beta = compositions(state, count)
    feed = feed_from_q(q, count)
    p, t = q[0], q[1]
    mu_l = phase(liquid, specs, kij, p, t, False)
    mu_v = phase(vapor, specs, kij, p, t, True)
    values = [mu_l[i] - mu_v[i] for i in range(count)]
    for i in range(count - 1):
        values.append((1 - beta) * liquid[i] + beta * vapor[i] - feed[i])
    return values


def solve_linear(matrix, rhs):
    count = len(rhs)
    augmented = [matrix[i][:] + [rhs[i]] for i in range(count)]
    for column in range(count):
        pivot = max(range(column, count), key=lambda row: abs(augmented[row][column]))
        if augmented[pivot][column] == 0:
            raise ArithmeticError("singular independent reference Jacobian")
        if pivot != column:
            augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        scale = augmented[column][column]
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


def numerical_jacobian(function, point, relative_step):
    rows = len(function(point))
    columns = len(point)
    matrix = [[D(0) for _ in range(columns)] for _ in range(rows)]
    for column in range(columns):
        h = relative_step * max(D(1), abs(point[column]))
        plus = point[:]
        minus = point[:]
        plus[column] += h
        minus[column] -= h
        upper = function(plus)
        lower = function(minus)
        for row in range(rows):
            matrix[row][column] = (upper[row] - lower[row]) / (2 * h)
    return matrix


def solve_equilibrium(case, q=None, seed=None):
    count = len(case["specs"])
    q = list(case["q"] if q is None else q)
    state = list(case["initial"] if seed is None else seed)
    expected_size = 2 * count - 1
    if len(state) != expected_size or len(q) != count + 1:
        raise ValueError("reference state dimension mismatch")
    for _ in range(48):
        values = residual(state, q, case)
        norm = max(abs(value) for value in values)
        if norm < D("1e-58"):
            break
        jacobian = numerical_jacobian(
            lambda trial: residual(trial, q, case), state, D("1e-24")
        )
        direction = solve_linear(jacobian, [-value for value in values])
        alpha = D(1)
        for _ in range(60):
            candidate = [state[i] + alpha * direction[i] for i in range(expected_size)]
            if valid_state(candidate, count):
                candidate_norm = max(abs(value) for value in residual(candidate, q, case))
                if candidate_norm < norm:
                    state = candidate
                    break
            alpha /= 2
        else:
            raise ArithmeticError("independent equilibrium Newton backtracking failed")
    else:
        raise ArithmeticError("independent equilibrium did not converge")
    if max(abs(value) for value in residual(state, q, case)) >= D("1e-58"):
        raise ArithmeticError("independent equilibrium residual too large")
    return state


def implicit_derivatives(case):
    count = len(case["specs"])
    q = list(case["q"])
    state = solve_equilibrium(case, q)
    gs = numerical_jacobian(
        lambda trial: residual(trial, q, case), state, D("1e-24")
    )
    gq = numerical_jacobian(
        lambda trial_q: residual(state, trial_q, case), q, D("1e-20")
    )
    q_count = len(q)
    state_count = len(state)
    ds = [[D(0) for _ in range(q_count)] for _ in range(state_count)]
    for column in range(q_count):
        rhs = [-gq[row][column] for row in range(state_count)]
        solution = solve_linear(gs, rhs)
        for row in range(state_count):
            ds[row][column] = solution[row]

    independent = count - 1
    liquid, vapor, _ = compositions(state, count)
    beta_gradient = ds[-1][:]
    liquid_jacobian = [row[:] for row in ds[:independent]]
    liquid_jacobian.append([
        -sum((liquid_jacobian[i][column] for i in range(independent)), D(0))
        for column in range(q_count)
    ])
    vapor_jacobian = [row[:] for row in ds[independent:2 * independent]]
    vapor_jacobian.append([
        -sum((vapor_jacobian[i][column] for i in range(independent)), D(0))
        for column in range(q_count)
    ])
    log_k_jacobian = [
        [
            vapor_jacobian[i][column] / vapor[i]
            - liquid_jacobian[i][column] / liquid[i]
            for column in range(q_count)
        ]
        for i in range(count)
    ]
    return state, beta_gradient, log_k_jacobian, liquid_jacobian, vapor_jacobian


def perturbation_crosscheck(case, reference):
    state, beta_gradient, log_k_jacobian, liquid_jacobian, vapor_jacobian = reference
    count = len(case["specs"])
    q = list(case["q"])
    steps = [D("1"), D("1e-4")] + [D("1e-6")] * (count - 1)
    max_scaled_error = D(0)
    for column, h in enumerate(steps):
        plus_q = q[:]
        minus_q = q[:]
        plus_q[column] += h
        minus_q[column] -= h
        plus = solve_equilibrium(case, plus_q, state)
        minus = solve_equilibrium(case, minus_q, state)
        lp, vp, bp = compositions(plus, count)
        lm, vm, bm = compositions(minus, count)
        targets = [((bp - bm) / (2 * h), beta_gradient[column])]
        for i in range(count):
            targets.extend([
                ((lp[i] - lm[i]) / (2 * h), liquid_jacobian[i][column]),
                ((vp[i] - vm[i]) / (2 * h), vapor_jacobian[i][column]),
                (((vp[i].ln() - lp[i].ln()) -
                  (vm[i].ln() - lm[i].ln())) / (2 * h),
                 log_k_jacobian[i][column]),
            ])
        for finite, implicit in targets:
            scale = max(D("1e-10"), abs(implicit))
            max_scaled_error = max(max_scaled_error, abs(finite - implicit) / scale)
    if max_scaled_error >= D("5e-8"):
        raise AssertionError(("perturbation cross-check", max_scaled_error))
    return max_scaled_error


def flatten(reference):
    _, beta_gradient, log_k_jacobian, liquid_jacobian, vapor_jacobian = reference
    return [
        *beta_gradient,
        *(value for row in log_k_jacobian for value in row),
        *(value for row in liquid_jacobian for value in row),
        *(value for row in vapor_jacobian for value in row),
    ]


if __name__ == "__main__":
    with localcontext() as context:
        context.prec = 80
        binary = implicit_derivatives(BINARY)
        ternary = implicit_derivatives(TERNARY)
        binary_error = perturbation_crosscheck(BINARY, binary)
        ternary_error = perturbation_crosscheck(TERNARY, ternary)
        print("binary perturb-and-resolve max scaled error:", binary_error)
        print("ternary perturb-and-resolve max scaled error:", ternary_error)

        header = Path(__file__).with_name("sensitivity_references.hpp").read_text(
            encoding="utf-8"
        )
        literals = re.findall(
            r"([-+]?[0-9]+(?:\.[0-9]*)?(?:[eE][-+]?[0-9]+)?)L", header
        )
        expected = flatten(binary) + flatten(ternary)
        if len(literals) != len(expected):
            raise AssertionError((len(literals), len(expected)))
        for literal, value in zip(literals, expected):
            tolerance = D("1e-30") * max(D(1), abs(value))
            if abs(D(literal) - value) >= tolerance:
                raise AssertionError((literal, value))
        print("Independent Decimal(80) PT sensitivity anchors: 61/61 checked")
