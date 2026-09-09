"""Independent PR76 binary boundary paths; no production imports.

Reuses the existing independent Decimal phase FORMULA, not C++ split output.
Newton solves two compositions at fixed p/T, or incipient composition and log(p)
with x=z (bubble) / y=z (dew). Negative lever coordinates are reference extensions,
not physical phase fractions. See boundary_regression.md for acceptance regions.
"""
from decimal import Decimal as D, localcontext
from pathlib import Path
from functools import lru_cache

from reference_split_decimal import phase, binary_equilibrium


def specification(carbon):
    raw = [('304.2','7380000','.225'), ('190.6','4600000','.008')] if carbon else [
           ('126.2','3390000','.04'), ('305.4','4880000','.098')]
    return ([[D(v) for v in row] for row in raw],
            D('.095' if carbon else '.08'), D('220' if carbon else '270'),
            D('6080000' if carbon else '7600000'))


def newton2(residual, initial, admissible):
    """Test-only high precision Newton; central derivatives never enter production."""
    x, y = initial
    h = D('1e-25')
    for _ in range(48):
        f = residual(x, y)
        norm = max(abs(v) for v in f)
        if norm < D('1e-60'):
            return x, y
        xp, xm = residual(x+h, y), residual(x-h, y)
        yp, ym = residual(x, y+h), residual(x, y-h)
        a, c = [(xp[i]-xm[i])/(2*h) for i in range(2)]
        b, d = [(yp[i]-ym[i])/(2*h) for i in range(2)]
        det = a*d-b*c
        if det == 0:
            raise ArithmeticError('singular reference Jacobian')
        dx, dy = (-d*f[0]+b*f[1])/det, (c*f[0]-a*f[1])/det
        alpha = D(1)
        for _ in range(48):
            xn, yn = x+alpha*dx, y+alpha*dy
            if admissible(xn, yn):
                if max(abs(v) for v in residual(xn, yn)) < norm:
                    x, y = xn, yn
                    break
            alpha /= 2
        else:
            raise ArithmeticError('reference boundary backtracking failed')
    raise ArithmeticError('reference boundary Newton did not converge')


@lru_cache(maxsize=16384)
def phase_at(carbon, p, w, vapor):
    specs, kij, t, _ = specification(carbon)
    return phase(w, specs, kij, p, t, vapor)


def chemical_residual(carbon, p, x, y):
    zl, ml = phase_at(carbon, p, x, False)
    zv, mv = phase_at(carbon, p, y, True)
    return [ml[i]-mv[i] for i in range(2)], zl, zv


def equilibrium(carbon, p, guess):
    def res(x, y):
        return chemical_residual(carbon, p, x, y)[0]
    x, y = newton2(res, guess, lambda x, y: 0 < x < 1 and 0 < y < 1 and abs(x-y) > D('.03'))
    f, zl, zv = chemical_residual(carbon, p, x, y)
    if not (zl < zv and max(abs(v) for v in f) < D('1e-60')):
        raise ArithmeticError('invalid equilibrium reference')
    return x, y, zl, zv


def saturation(carbon, bubble, feed=D('.30')):
    _, _, _, p0 = specification(carbon)
    x0, y0, _, _ = binary_equilibrium(carbon)
    def res(w, q):
        p = p0*q.exp()
        return chemical_residual(carbon, p, feed if bubble else w, w if bubble else feed)[0]
    # Start on the already checked coexistence branch, then bracket the requested
    # feed composition in pressure. This avoids Newton's trivial same-root basin.
    p = p0
    xy = (x0, y0)
    target = 0 if bubble else 1
    initial_sign = xy[target]-feed
    step = D('.03') if bubble else D('-.03')
    for _ in range(160):
        pn = p*(1+step)
        try:
            trial = equilibrium(carbon, pn, xy)
        except ArithmeticError:
            step /= 2
            if abs(step) < D('1e-8'):
                raise ArithmeticError('reference saturation bracket stalled')
            continue
        p, xy = pn, trial[:2]
        if (xy[target]-feed)*initial_sign < 0:
            break
    else:
        raise ArithmeticError('reference saturation bracket not found')
    w, q = newton2(res, (xy[1 if bubble else 0], (p/p0).ln()),
        lambda w, q: 0 < w < 1 and abs(w-feed) > D('.03') and abs(q) < 2)
    p = p0*q.exp()
    x, y = (feed, w) if bubble else (w, feed)
    f, zl, zv = chemical_residual(carbon, p, x, y)
    if not (zl < zv and max(abs(v) for v in f) < D('1e-60')):
        raise ArithmeticError('invalid saturation reference')
    return p, x, y, zl, zv


def minimum_phase(carbon, p, w):
    choices = [phase_at(carbon, p, w, vapor) for vapor in (False, True)]
    return min(choices, key=lambda result: w*result[1][0]+(1-w)*result[1][1])


def tpd(w, mu, ref):
    return w*(mu[0]-ref[0])+(1-w)*(mu[1]-ref[1])


def verify_tie(carbon, p, x, y):
    # Finite independent branch/grid sanity check, NOT a global certificate.
    _, common = phase_at(carbon, p, x, False)
    for w in [D(j)/32 for j in range(1, 32)]+[x, y]:
        _, mu = minimum_phase(carbon, p, w)
        if tpd(w, mu, common) < D('-1e-45'):
            raise ArithmeticError('coexistence reference has a lower sampled Gibbs branch')


def row(carbon, p, feed, tie, coordinate):
    x, y, zl, zv = tie
    _, _, t, _ = specification(carbon)
    beta = (feed-x)/(y-x)  # May lie outside [0,1]: algebraic extension only.
    if abs((1-beta)*x+beta*y-feed) > D('1e-65'):
        raise ArithmeticError('reference material balance')
    zg, mug = minimum_phase(carbon, p, feed)
    witnesses = []
    for w, vapor in [(x, False), (y, True)]:
        _, mu = phase_at(carbon, p, w, vapor)
        witnesses.append(tpd(w, mu, mug))
    distance = min(witnesses)
    # Do not publish the sign of reference Newton noise at an exact boundary.
    # This reporting convention applies ONLY to the diagnostic TPD column;
    # compositions and the signed lever coordinate are never clipped.
    if abs(distance) < D('1e-50'):
        distance = D(0)
    return [int(carbon), coordinate, p, t, feed, x, y, beta, zl, zv, zg, distance]


def check_conditioning(carbon, p, x, y, beta):
    # Local linearized sensitivity of independent composition equations.
    # These diagnostic ceilings justify test comparison budgets, not a certified
    # nonlinear error enclosure and not a changed production stopping tolerance.
    h = D('1e-25')
    xp = chemical_residual(carbon, p, x+h, y)[0]
    xm = chemical_residual(carbon, p, x-h, y)[0]
    yp = chemical_residual(carbon, p, x, y+h)[0]
    ym = chemical_residual(carbon, p, x, y-h)[0]
    a, c = [(xp[i]-xm[i])/(2*h) for i in range(2)]
    b, d = [(yp[i]-ym[i])/(2*h) for i in range(2)]
    determinant = a*d-b*c
    inverse = ((d/determinant, -b/determinant), (-c/determinant, a/determinant))
    composition_scale = max(sum(abs(v) for v in row) for row in inverse)
    dbdx, dbdy = (beta-1)/(y-x), -beta/(y-x)
    beta_scale = sum(abs(dbdx*inverse[0][i]+dbdy*inverse[1][i]) for i in range(2))
    if composition_scale >= 50 or beta_scale >= 500:
        raise ArithmeticError('reference path conditioning exceeds declared comparison budget')


def generate():
    # Cache is test-only, bounded, and never survives a precision change.
    phase_at.cache_clear()
    ties, boundaries, pressure_rows = [], [], []
    coordinates = [D(s) for s in (
        '-.001', '-.000001', '-.0000000001', '0', '.0000000001', '.000001',
        '.0001', '.01', '.1', '.25', '.5', '.75', '.9', '.99', '.9999',
        '.999999', '.9999999999', '1', '1.0000000001', '1.000001', '1.001')]
    for carbon in (False, True):
        _, _, t, p0 = specification(carbon)
        tie = binary_equilibrium(carbon)
        verify_tie(carbon, p0, *tie[:2])
        ties.append([int(carbon), p0, t, *tie])
        bubble = saturation(carbon, True)
        dew = saturation(carbon, False)
        if not (dew[0] < p0 < bubble[0]):
            raise ArithmeticError('reference isotherm ordering changed')
        for flag, data in [(1, bubble), (0, dew)]:
            boundaries.append([int(carbon), flag, data[0], t, *data[1:]])
        guess = dew[1:3]
        for q in coordinates:
            p = dew[0]+q*(bubble[0]-dew[0])
            # Repeated smaller reference-only steps prevent a large jump from
            # finding the trivial same-root solution. Not production continuation.
            if pressure_rows and pressure_rows[-1][0] == int(carbon):
                oldp = pressure_rows[-1][2]
                for j in range(1, 3):
                    intermediate = oldp+(p-oldp)*D(j)/2
                    guess = equilibrium(carbon, intermediate, guess)[:2]
            tie = dew[1:] if q == 0 else (bubble[1:] if q == 1 else equilibrium(carbon, p, guess))
            guess = tie[:2]
            verify_tie(carbon, p, *tie[:2])
            entry = row(carbon, p, D('.30'), tie, q)
            if q < 0 and not entry[7] > 1:
                raise ArithmeticError('dew exterior is not vapor-side reference')
            if q > 1 and not entry[7] < 0:
                raise ArithmeticError('bubble exterior is not liquid-side reference')
            check_conditioning(carbon, p, tie[0], tie[1], entry[7])
            pressure_rows.append(entry)
    return ties, boundaries, pressure_rows


def literal(value):
    if isinstance(value, int):
        return str(value)
    if value == 0:
        return '0.0L'  # Canonical exact zero, independent of Decimal exponent/sign.
    text = format(value, '.40g')
    if '.' not in text and 'e' not in text.lower():
        text += '.0'
    return text+'L'


def render(ties, boundaries, rows):
    text = ("// Generated ONLY by reference_boundary_decimal.py; PR76 model regression.\n"
            "// Hua 1997 revised author manuscript parameters; not experimental data.\n"
            "#ifndef MPMC_TEST_BOUNDARY_REFERENCES_HPP\n#define MPMC_TEST_BOUNDARY_REFERENCES_HPP\n"
            "namespace boundary_reference {\n"
            "struct Tie { int carbon; long double p, t, x, y, zl, zv; };\n"
            "struct Saturation { int carbon, bubble; long double p, t, x, y, zl, zv; };\n"
            "struct State { int carbon; long double coordinate, p, t, feed, x, y, beta, zl, zv, feed_z, witness_tpd; };\n")
    for typename, name, data in [('Tie','ties',ties), ('Saturation','saturation',boundaries),
                                 ('State','pressure_path',rows)]:
        text += 'inline constexpr '+typename+' '+name+'[] = {\n'
        for values in data:
            text += '    {'+', '.join(literal(v) for v in values)+'},\n'
        text += '};\n'
    return text+'} // namespace boundary_reference\n#endif\n'


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true', help='Explicitly regenerate test-only fixture')
    parser.add_argument('--precision', type=int, default=80, choices=(80, 96))
    args = parser.parse_args()
    with localcontext() as ctx:
        ctx.prec = args.precision
        ties, boundaries, rows = generate()
        text = render(ties, boundaries, rows)
        path = Path(__file__).with_name('boundary_references.hpp')
        if args.write:
            path.write_text(text, encoding='utf-8')
        elif path.read_text(encoding='utf-8') != text:
            raise AssertionError('boundary fixture mismatch; investigate, do not fit production output')
        print(f'Independent Decimal({args.precision}): {len(ties)} tie lines, '
              f'{len(boundaries)} saturation states, {len(rows)} pressure states checked')
        print('Coexistence chemical residual < 1e-60; 33-point minimum-branch checks are NOT a global proof')
        for data in boundaries:
            print(('CO2/methane' if data[0] else 'N2/ethane'),
                  ('bubble' if data[1] else 'dew'), 'Pa=', format(data[2], '.20g'))
