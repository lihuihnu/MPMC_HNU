"""Independent Decimal(80) PR76 closure references.

This extends the already independent PT-equilibrium reference used by the flash
sensitivity regression. It never calls production C++, AD, RR/SSI, TPD, or the
physics adapter. Phase Z/c derivatives are assembled in independent equilibrium
coordinates and then cross-checked by finite perturb-and-resolve calculations.
"""

from decimal import Decimal as D, localcontext
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[3]
SOURCE = ROOT / "tests/flash/pt_split/reference_sensitivity_decimal.py"
SPEC = spec_from_file_location("mpmc_independent_sensitivity_reference", SOURCE)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("cannot load independent PT sensitivity reference")
REF = module_from_spec(SPEC)
SPEC.loader.exec_module(REF)


def phase_z(composition, case, pressure, temperature, vapor):
    specs = case["specs"]
    kij = case["kij"]
    count = len(composition)
    ai = []
    bi = []
    for tc, pc, omega in specs:
        kappa = D("0.37464") + D("1.54226") * omega - D("0.26992") * omega**2
        ai.append(
            D("0.45724") * REF.R * REF.R * tc * tc / pc
            * (1 + kappa * (1 - (temperature / tc).sqrt())) ** 2
        )
        bi.append(D("0.07780") * REF.R * tc / pc)
    aij = [
        [
            (ai[i] * ai[j]).sqrt() * (1 - kij[i][j])
            for j in range(count)
        ]
        for i in range(count)
    ]
    rows = [
        sum((composition[j] * aij[i][j] for j in range(count)), D(0))
        for i in range(count)
    ]
    a_mix = sum((composition[i] * rows[i] for i in range(count)), D(0))
    b_mix = sum((composition[i] * bi[i] for i in range(count)), D(0))
    aa = a_mix * pressure / (REF.R * temperature) ** 2
    bb = b_mix * pressure / (REF.R * temperature)
    return REF.roots(aa, bb)[-1 if vapor else 0]


def scalar_gradient(function, point, relative_step):
    values = []
    for column in range(len(point)):
        h = relative_step * max(D(1), abs(point[column]))
        plus = point[:]
        minus = point[:]
        plus[column] += h
        minus[column] -= h
        values.append((function(plus) - function(minus)) / (2 * h))
    return values


def phase_z_partials(composition, case, pressure, temperature, vapor):
    point = [pressure, temperature, *composition[:-1]]

    def evaluate(values):
        phase = list(values[2:])
        phase.append(D(1) - sum(phase, D(0)))
        return phase_z(phase, case, values[0], values[1], vapor)

    return evaluate(point), scalar_gradient(evaluate, point, D("1e-20"))


def closure_reference(case):
    state, _, _, liquid_jacobian, vapor_jacobian = REF.implicit_derivatives(case)
    count = len(case["specs"])
    q = list(case["q"])
    liquid, vapor, _ = REF.compositions(state, count)
    zl, partial_l = phase_z_partials(liquid, case, q[0], q[1], False)
    zv, partial_v = phase_z_partials(vapor, case, q[0], q[1], True)
    q_count = len(q)

    dzl = []
    dzv = []
    for column in range(q_count):
        direct_l = partial_l[column] if column < 2 else D(0)
        direct_v = partial_v[column] if column < 2 else D(0)
        for component in range(count - 1):
            direct_l += partial_l[2 + component] * liquid_jacobian[component][column]
            direct_v += partial_v[2 + component] * vapor_jacobian[component][column]
        dzl.append(direct_l)
        dzv.append(direct_v)

    cl = q[0] / (zl * REF.R * q[1])
    cv = q[0] / (zv * REF.R * q[1])
    dcl = [
        cl * ((D(1) / q[0] if column == 0 else D(0))
              - dzl[column] / zl
              - (D(1) / q[1] if column == 1 else D(0)))
        for column in range(q_count)
    ]
    dcv = [
        cv * ((D(1) / q[0] if column == 0 else D(0))
              - dzv[column] / zv
              - (D(1) / q[1] if column == 1 else D(0)))
        for column in range(q_count)
    ]
    return state, [zl, zv, cl, cv], dzl, dzv, dcl, dcv


def perturbation_crosscheck(case, reference):
    state, _, dzl, dzv, dcl, dcv = reference
    count = len(case["specs"])
    q = list(case["q"])
    steps = [D("1"), D("1e-4")] + [D("1e-6")] * (count - 1)
    max_scaled_error = D(0)
    for column, h in enumerate(steps):
        plus_q = q[:]
        minus_q = q[:]
        plus_q[column] += h
        minus_q[column] -= h
        plus_state = REF.solve_equilibrium(case, plus_q, state)
        minus_state = REF.solve_equilibrium(case, minus_q, state)
        lp, vp, _ = REF.compositions(plus_state, count)
        lm, vm, _ = REF.compositions(minus_state, count)
        zlp = phase_z(lp, case, plus_q[0], plus_q[1], False)
        zlm = phase_z(lm, case, minus_q[0], minus_q[1], False)
        zvp = phase_z(vp, case, plus_q[0], plus_q[1], True)
        zvm = phase_z(vm, case, minus_q[0], minus_q[1], True)
        clp = plus_q[0] / (zlp * REF.R * plus_q[1])
        clm = minus_q[0] / (zlm * REF.R * minus_q[1])
        cvp = plus_q[0] / (zvp * REF.R * plus_q[1])
        cvm = minus_q[0] / (zvm * REF.R * minus_q[1])
        finite = [
            (zlp - zlm) / (2 * h),
            (zvp - zvm) / (2 * h),
            (clp - clm) / (2 * h),
            (cvp - cvm) / (2 * h),
        ]
        implicit = [dzl[column], dzv[column], dcl[column], dcv[column]]
        for fd, exact in zip(finite, implicit):
            scale = max(D("1e-10"), abs(exact))
            max_scaled_error = max(max_scaled_error, abs(fd - exact) / scale)
    if max_scaled_error >= D("5e-8"):
        raise AssertionError(("closure perturbation cross-check", max_scaled_error))
    return max_scaled_error


def flatten(reference):
    _, primal, dzl, dzv, dcl, dcv = reference
    return [*primal, *dzl, *dzv, *dcl, *dcv]


if __name__ == "__main__":
    with localcontext() as context:
        context.prec = 80
        binary = closure_reference(REF.BINARY)
        ternary = closure_reference(REF.TERNARY)
        binary_error = perturbation_crosscheck(REF.BINARY, binary)
        ternary_error = perturbation_crosscheck(REF.TERNARY, ternary)
        print("binary closure perturb-and-resolve max scaled error:", binary_error)
        print("ternary closure perturb-and-resolve max scaled error:", ternary_error)

        header = Path(__file__).with_name("closure_references.hpp").read_text(
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
        print("Independent Decimal(80) thermodynamic-closure anchors: 36/36 checked")
