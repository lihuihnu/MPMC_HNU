"""Independent audit of the PR11 CO2/methane near-dew input, not a new model.

Only the frozen PR11 header and its independent Decimal composition equations are
read. No C++/RR/SSI output is imported. Exact binary64 pressure and normalized
binary64 feed are distinguished from the nominal decimal input. Existing oracles
are never rewritten. The precision option only controls test-side arithmetic.
"""
from decimal import Decimal as D, localcontext
from pathlib import Path
import math
import re
from reference_boundary_decimal import equilibrium, phase_at, chemical_residual


def table(name):
    text = Path(__file__).with_name('boundary_references.hpp').read_text(encoding='utf-8')
    block = text.split(name+'[] = {', 1)[1].split('};', 1)[0]
    return [[D(v.strip().removesuffix('L')) for v in row.split(',')]
            for row in re.findall(r'\{([^{}]+)\}', block)]


def compute():
    phase_at.cache_clear()
    sat = [row for row in table('saturation') if row[0] == 1]
    dew = next(row for row in sat if row[1] == 0)
    bubble = next(row for row in sat if row[1] == 1)
    original = next(row for row in table('pressure_path')
                    if row[0] == 1 and row[1] == D('1e-10'))
    p0 = float(original[2])
    z = D.from_float(.3)/(D.from_float(.3)+D.from_float(.7))
    nominal = equilibrium(True, original[2], (original[5], original[6]))
    nominal_liquid = (D('.3')-nominal[1])/(nominal[0]-nominal[1])
    if abs(nominal_liquid-(1-original[7])) > D('1e-38'):
        raise AssertionError('frozen PR11 nominal reference mismatch')
    cases = [('outside', float(dew[2]-D('1e-10')*(bubble[2]-dew[2])), -1),
             ('endpoint', float(dew[2]), 0),
             ('target_minus_ulp', math.nextafter(p0, -math.inf), 1),
             ('target', p0, 1),
             ('target_plus_ulp', math.nextafter(p0, math.inf), 1)]
    cases += [(name, float(dew[2]+D(q)*(bubble[2]-dew[2])), 2)
              for name, q in [('above_2e-10','2e-10'), ('above_4e-10','4e-10'),
                              ('inside_1e-8','1e-8')]]
    rows = []
    for name, pf, region in cases:
        p = D.from_float(pf)
        x, y, zl, zv = equilibrium(True, p, (dew[4], dew[5]))
        liquid = (z-y)/(x-y)
        _, ml = phase_at(True, p, x, False)
        _, mv = phase_at(True, p, y, True)
        _, feed_mu = phase_at(True, p, z, True)
        distance = x*(ml[0]-feed_mu[0])+(1-x)*(ml[1]-feed_mu[1])
        right = (z*x/y+(1-z)*(1-x)/(1-y)).ln()
        h = D('1e-25')
        xp, xm = chemical_residual(True,p,x+h,y)[0], chemical_residual(True,p,x-h,y)[0]
        yp, ym = chemical_residual(True,p,x,y+h)[0], chemical_residual(True,p,x,y-h)[0]
        a,c = [(xp[i]-xm[i])/(2*h) for i in range(2)]
        b,d = [(yp[i]-ym[i])/(2*h) for i in range(2)]
        det = a*d-b*c
        dl_dx, dl_dy = -liquid/(x-y), (liquid-1)/(x-y)
        sensitivity = abs((dl_dx*d-dl_dy*c)/det)+abs((-dl_dx*b+dl_dy*a)/det)
        if max(abs(ml[i]-mv[i]) for i in range(2)) >= D('1e-60'):
            raise ArithmeticError('reference equilibrium residual')
        if abs(liquid*x+(1-liquid)*y-z) >= D('1e-65'):
            raise ArithmeticError('reference lever balance')
        if region == 1:
            if not (D('6e-11') < liquid < D('8e-11') and distance < D('-1.1e-10')
                    and right > D('1e-11')):
                raise AssertionError('target no longer separates TPD, RR and disappearance gates')
            if sensitivity*D('1e-11') >= (D('1e-10')-liquid)/4:
                raise AssertionError('local residual sensitivity too large for audit classification')
        if region == 2 and not liquid > D('1.2e-10'):
            raise AssertionError('upper neighbor no longer above the declared gate')
        if region == -1 and not liquid < 0:
            raise AssertionError('outside reference must stay a signed extension')
        rows.append((name, pf, region, x, y, liquid, zl, zv, distance, right, sensitivity))
    target = next(row for row in rows if row[0] == 'target')
    if abs(target[5]-nominal_liquid) > D('1e-15'):
        raise AssertionError('input quantization exceeds separate audit budget')
    if abs(rows[4][5]-rows[2][5]) > D('1e-15'):
        raise AssertionError('adjacent binary64 input uncertainty exceeds audit budget')
    return rows, nominal_liquid


def render(rows):
    text = ('// Generated only by reference_dew_limit_decimal.py; no production output.\n'
            '// Pressure is exact binary64; other fields use the independent decimal PR76 model.\n'
            '#ifndef MPMC_TEST_DEW_LIMIT_REFERENCES_HPP\n#define MPMC_TEST_DEW_LIMIT_REFERENCES_HPP\n'
            'namespace dew_limit_reference {\n'
            'struct State { const char* name; double p; int region;\n'
            '    long double x, y, liquid, zl, zv, tpd, rr_right, sensitivity; };\n'
            'inline constexpr State states[] = {\n')
    for name,p,region,*values in rows:
        literals = [format(v, '.30g')+'L' for v in values]
        text += '    {"'+name+'", '+p.hex()+', '+str(region)+', '+', '.join(literals)+'},\n'
    return text+'};\n} // namespace dew_limit_reference\n#endif\n'


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--precision', type=int, choices=(80,96), default=80)
    parser.add_argument('--write', action='store_true', help='Write only the NEW audit fixture')
    args = parser.parse_args()
    with localcontext() as context:
        context.prec = args.precision
        rows, nominal = compute()
        text = render(rows)
        path = Path(__file__).with_name('dew_limit_references.hpp')
        if args.write:
            path.write_text(text, encoding='utf-8')
        elif text != path.read_text(encoding='utf-8'):
            raise AssertionError('new audit fixture mismatch; do not fit C++ results')
        print(f'Decimal({args.precision}): {len(rows)} near-dew input states checked')
        print('Nominal PR11 liquid amount:', format(nominal,'.24g'))
        for name,p,region,x,y,liquid,zl,zv,tpd,right,sensitivity in rows:
            print(name, 'p=', repr(p), 'liquid=', format(liquid,'.20g'),
                  'TPD=',format(tpd,'.12g'), 'RR-right=',format(right,'.12g'),
                  'local residual sensitivity=',format(sensitivity,'.12g'))
        print('Local sensitivities are diagnostics, not nonlinear/global certificates.')
