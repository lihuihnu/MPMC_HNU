"""Independent binary equilibrium references: Decimal(80), stdlib only.

Original Z polynomial + refined analytic roots, direct PR formula, and Newton
in the two phase COMPOSITIONS (not production logK-SSI/RR). Numerical derivatives
are only for this independent reference; analytic fixed-K tests are separate.
Parameters: Hua/Brennecke/Stadtherr October 1997 author manuscript, sections 4.3/4.4.
These are PR numerical regressions, not experimental fluid validation.
"""
from decimal import Decimal as D, localcontext
import math
from pathlib import Path
import re


def roots(a, b):
    c2, c1, c0 = b-1, a-3*b*b-2*b, b*b*b+b*b-a*b
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
    found = []
    for seed in seeds:
        z = D(str(seed))
        for _ in range(40):
            z -= (((z+c2)*z+c1)*z+c0)/((3*z+2*c2)*z+c1)
        residual = ((z+c2)*z+c1)*z+c0
        if abs(residual) > D('1e-70')*(1+abs(z)**3+abs(c2*z*z)+abs(c1*z)+abs(c0)):
            raise ArithmeticError('reference root did not converge')
        if z > b and (3*z+2*c2)*z+c1 > 0:
            found.append(z)
    if not found:
        raise ArithmeticError('no admissible reference roots')
    return sorted(found)


def phase(x, specs, kij, p, t, vapor):
    w = [x, 1-x]
    r = D('8.31446261815324')
    ai, bi = [], []
    for tc, pc, omega in specs:
        k = D('.37464')+D('1.54226')*omega-D('.26992')*omega**2
        ai.append(D('.45724')*r*r*tc*tc/pc*(1+k*(1-(t/tc).sqrt()))**2)
        bi.append(D('.07780')*r*tc/pc)
    aij = [[(ai[i]*ai[j]).sqrt()*(1-(kij if i != j else 0)) for j in range(2)] for i in range(2)]
    sums = [sum(w[j]*aij[i][j] for j in range(2)) for i in range(2)]
    a = sum(w[i]*sums[i] for i in range(2))
    b = sum(w[i]*bi[i] for i in range(2))
    aa, bb = a*p/(r*t)**2, b*p/(r*t)
    z = roots(aa, bb)[-1 if vapor else 0]
    s = D(2).sqrt()
    lr = ((z+(1+s)*bb)/(z+(1-s)*bb)).ln()
    phi = [bi[i]/b*(z-1)-(z-bb).ln()-aa/(2*s*bb)*(2*sums[i]/a-bi[i]/b)*lr for i in range(2)]
    return z, [w[i].ln()+phi[i] for i in range(2)]


def binary_equilibrium(carbon):
    specs = [('304.2','7380000','.225'),('190.6','4600000','.008')] if carbon else [
             ('126.2','3390000','.04'),('305.4','4880000','.098')]
    specs = [[D(v) for v in row] for row in specs]
    kij, p, t = (D('.095'), D('6080000'), D('220')) if carbon else (D('.08'), D('7600000'), D('270'))
    x, y = (D('.48'), D('.18')) if carbon else (D('.16'), D('.50'))
    def residual(l, v):
        _, ml = phase(l, specs, kij, p, t, False)
        _, mv = phase(v, specs, kij, p, t, True)
        return [ml[i]-mv[i] for i in range(2)]
    h = D('1e-25')
    for _ in range(32):
        f = residual(x,y)
        norm = max(abs(v) for v in f)
        if norm < D('1e-60'):
            break
        fxp, fxm = residual(x+h,y), residual(x-h,y)
        fyp, fym = residual(x,y+h), residual(x,y-h)
        a,c = [(fxp[i]-fxm[i])/(2*h) for i in range(2)]
        b,d = [(fyp[i]-fym[i])/(2*h) for i in range(2)]
        determinant = a*d-b*c
        dx,dy = (-d*f[0]+b*f[1])/determinant, (c*f[0]-a*f[1])/determinant
        alpha = D(1)
        for _ in range(40):
            xn,yn = x+alpha*dx,y+alpha*dy
            if 0 < xn < 1 and 0 < yn < 1 and abs(xn-yn)>D('.05'):
                if max(abs(v) for v in residual(xn,yn)) < norm:
                    x,y = xn,yn
                    break
            alpha /= 2
        else:
            raise ArithmeticError('reference Newton backtracking failed')
    else:
        raise ArithmeticError('reference equilibrium did not converge')
    zl, ml = phase(x, specs, kij, p, t, False)
    zv, mv = phase(y, specs, kij, p, t, True)
    assert zl < zv and max(abs(ml[i]-mv[i]) for i in range(2)) < D('1e-60')
    return x,y,zl,zv


def anchors():
    values = []
    for carbon, feeds in [(False,['.18','.30','.44']), (True,['.20','.30','.43'])]:
        x,y,zl,zv = binary_equilibrium(carbon)
        for feed in feeds:
            z = D(feed)
            beta = (z-x)/(y-x)
            assert 0 < beta < 1 and abs((1-beta)*x+beta*y-z)<D('1e-70')
            values.append((x,y,beta,zl,zv))
    return values


if __name__ == '__main__':
    with localcontext() as context:
        context.prec = 80
        computed = anchors()
        for row in computed:
            print(', '.join(format(x, '.40g')+'L' for x in row))
        text = Path(__file__).with_name('split_test.cpp').read_text(encoding='utf-8')
        block = text.split('constexpr long double golden[6][5] = {',1)[1].split('};',1)[0]
        literals = re.findall(r'([-+]?[0-9]+\.[0-9]*(?:[eE][-+]?[0-9]+)?)L',block)
        flat = [x for row in computed for x in row]
        assert len(literals) == len(flat)
        for literal, value in zip(literals,flat):
            assert abs(D(literal)-value) < D('1e-35')*abs(value), (literal,value)
        print('Independent Decimal(80) binary split anchors: 30/30 checked')
