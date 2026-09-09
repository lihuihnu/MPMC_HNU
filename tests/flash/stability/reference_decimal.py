"""Independent fixed-state PR76 TPD anchors; stdlib only, no production imports.

Raw fixtures are declared here independently. Original Z cubic, double analytic
seeds refined in Decimal(80), classical full double sum and direct fugacity
expression. All roots have checked residuals. This verifies algebra/model
regressions, not experiments or a global stability proof. See flash/README.md.
"""
from decimal import Decimal as D, localcontext
import math
from pathlib import Path
import re


def roots(a, b):
    c2, c1, c0 = b - 1, a - 3*b*b - 2*b, b*b*b + b*b - a*b
    q, r = (3*c1-c2*c2)/9, (9*c2*c1-27*c0-2*c2**3)/54
    disc = q**3+r*r
    if disc >= 0:
        sd = math.sqrt(float(disc))
        def cbrt(x):
            return math.copysign(abs(x)**(1/3), x)
        seeds = [cbrt(float(r)+sd)+cbrt(float(r)-sd)-float(c2)/3]
    else:
        theta = math.acos(float(r/(-q**3).sqrt()))
        seeds = [2*math.sqrt(float(-q))*math.cos((theta+2*k*math.pi)/3)-float(c2)/3
                 for k in range(3)]
    values = []
    for seed in seeds:
        z = D(str(seed))
        for _ in range(40):
            f = ((z+c2)*z+c1)*z+c0
            df = (3*z+2*c2)*z+c1
            z -= f/df
        f = ((z+c2)*z+c1)*z+c0
        if abs(f) > D('1e-70') * (1+abs(z)**3+abs(c2*z*z)+abs(c1*z)+abs(c0)):
            raise ArithmeticError('reference root did not converge')
        if z > b and (3*z+2*c2)*z+c1 > 0:
            values.append(z)
    return sorted(values)


def phase(p, t, composition, specs, interactions):
    p, t = D(str(p)), D(str(t))
    w = [D(str(x)) for x in composition]
    r = D('8.31446261815324')
    ap, bp = [], []
    for tc, pc, omega in specs:
        tc, pc, omega = (D(str(x)) for x in (tc, pc, omega))
        kappa = D('.37464')+D('1.54226')*omega-D('.26992')*omega*omega
        ap.append(D('.45724')*r*r*tc*tc/pc*(1+kappa*(1-(t/tc).sqrt()))**2)
        bp.append(D('.07780')*r*tc/pc)
    pairs = [[(1-D(str(interactions[i][j])))*(ap[i]*ap[j]).sqrt()
              for j in range(len(w))] for i in range(len(w))]
    row = [sum(w[j]*pairs[i][j] for j in range(len(w))) for i in range(len(w))]
    a = sum(w[i]*row[i] for i in range(len(w)))
    b = sum(w[i]*bp[i] for i in range(len(w)))
    aa, bb = a*p/(r*t)**2, b*p/(r*t)
    rt2 = D(2).sqrt()
    choices = []
    for z in roots(aa, bb):
        log_ratio = ((z+(1+rt2)*bb)/(z+(1-rt2)*bb)).ln()
        phi = [bp[i]/b*(z-1)-(z-bb).ln()-aa/(2*rt2*bb)*
               (2*row[i]/a-bp[i]/b)*log_ratio for i in range(len(w))]
        choices.append((sum(w[i]*phi[i] for i in range(len(w))), z, phi))
    return min(choices)[1:]


def tpd(p, t, feed, trial, specs, interactions):
    z, w = ([D(str(x)) for x in fractions] for fractions in (feed, trial))
    _, reference = phase(p, t, z, specs, interactions)
    _, candidate = phase(p, t, w, specs, interactions)
    return sum(w[i]*((w[i]/z[i]).ln()+candidate[i]-reference[i]) for i in range(len(w)))


def anchors():
    n2 = [(126.2, 3390000, .04), (305.4, 4880000, .098)]
    co2 = [(304.2, 7380000, .225), (190.6, 4600000, .008)]
    # One scalar each from published negative trial compositions (not complete tables).
    first = tpd(7600000, 270, ['.18','.82'], ['.4943','.5057'], n2, [[0,.08],[.08,0]])
    second = tpd(6080000, 220, ['.20','.80'], ['.4972','.5028'], co2, [[0,.095],[.095,0]])
    artificial = [(400,4000000,-.125),(500,3000000,.25),(600,5000000,.5)]
    pairs = [[0,.125,-.0625],[.125,0,.0625],[-.0625,.0625,0]]
    third = tpd(1000000,450,['.25','.5','.25'],['.4','.2','.4'], artificial,pairs)
    return [first, second, third]


if __name__ == '__main__':
    with localcontext() as context:
        context.prec = 80
        computed = anchors()
        path = Path(__file__).with_name('stability_test.cpp')
        text = path.read_text(encoding='utf-8')
        block = text.split('constexpr long double tpd_golden[] = {',1)[1].split('};',1)[0]
        literals = re.findall(r'([-+]?[0-9]+\.[0-9]*(?:[eE][-+]?[0-9]+)?)L', block)
        if len(literals) != len(computed):
            raise AssertionError('anchor count changed')
        for literal, value in zip(literals, computed):
            if abs(D(literal)-value) > D('1e-35')*abs(value):
                raise AssertionError(f'anchor mismatch: {literal} vs {value}')
            print(format(value, '.40g'))
        print('Independent Decimal(80) TPD anchors: 3/3 checked')
