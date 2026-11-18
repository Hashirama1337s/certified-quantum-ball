"""build-1 QB1 positive controls and the sign-flip null, all through ajet / charts (the EC code path).

C0  Zahn Hayward xi = 0 (metric e^{-2u}(-2 dtau^2 + 2 du^2 + dOmega^2), u = 0): V_tautau = 1/80, V_uu = 17/240,
    V_thth = -7/480 (Zahn eq. 33).
C3  Reissner-Nordstrom M = 1, Q = 1/2 at r = 3 (design-1 rncheck closed forms): V(xi=0) = V(1/6) = closed form,
    Bach = 60 V, Route A (RicS := geometric Ric) = V.
C2  Minkowski patch (KHA N = e^{-x}, A = 1, om = V = 0  =>  F2 = 1, At = Wt = Vt = 0): B = 0, through polar AND cart.
C6  flat Friedmann k = 1/2 (engine/frw.py closed form, KHA chart; AMEND-LB-E-1), through polar AND cart (x-mode):
    r^4 V_uu(xi=0) = 2 z^4/243 at z = 1/2, 1, 3/2 (z the frw.py parameter), V_uu(1/6) = 0, Bach = 0, Einstein holds.
NULL-flip: sign of the R^{rs}R_{m r n s} term in T(2) flipped: Hayward must FAIL (V_thth = -1/32).
"""
from flint import arb, ctx
import ajet as J
from ajet import Jet, n4
import charts as CH
from aser import Ser, revert

ctx.cap = 64


def hayward(flip=False):
    J.setup(3, 4)
    d1, d2 = Jet.var(1), Jet.var(2)
    sn = (d2 + arb.pi()/2).sin()
    e = (d1*(-2)).exp()
    g = [[Jet() for _ in n4] for _ in n4]
    g[0][0] = e*(-2); g[1][1] = e*2; g[2][2] = e; g[3][3] = e*sn*sn
    d = J.curv(g); T1, T2 = J.zahn(d, flip=flip); V = J.Vxi(T1, T2, 0)
    got = [V[i][i].c0() for i in range(3)]
    exp_ = [arb(1)/80, arb(17)/240, arb(-7)/480]
    tr = sum((d['gi'][a][b].c0()*V[a][b].c0() for a in n4 for b in n4), arb(0))
    return dict(got=got, expected=exp_, contains=[g_.contains(e_) for g_, e_ in zip(got, exp_)], trace=tr)


def rn(R0=3):
    J.setup(3, 4)
    d1, d2 = Jet.var(1), Jet.var(2)
    sn = (d2 + arb.pi()/2).sin()
    r = d1 + R0
    f = 1 - 2/r + (r*r).inv()*arb('0.25')
    g = [[Jet() for _ in n4] for _ in n4]
    g[0][0] = -f; g[1][1] = f.inv(); g[2][2] = r*r; g[3][3] = r*r*sn*sn
    d = J.curv(g); T1, T2 = J.zahn(d); Bc = J.bach(T1, T2)
    R = arb(R0)
    ex = [(4*R**2 - 8*R + 1)**2/(640*R**10), -1/(120*R**6), (4*R**2 - 8*R + 1)/(240*R**6)]
    res = {}
    for nm, xi in (('xi0', 0), ('xi1/6', arb(1)/6)):
        V = J.Vxi(T1, T2, xi)
        res['V_' + nm] = [V[i][i].c0() for i in range(3)]
        res['V_' + nm + '_contains'] = [V[i][i].c0().contains(ex[i]) for i in range(3)]
    res['Bach'] = [Bc[i][i].c0() for i in range(3)]
    res['Bach_eq_60V'] = [Bc[i][i].c0().contains(60*ex[i]) for i in range(3)]
    VA = J.routeA(d, d['Ric'])
    res['routeA'] = [VA[i][i].c0() for i in range(3)]
    res['routeA_contains'] = [VA[i][i].c0().contains(ex[i]) for i in range(3)]
    res['expected_V'] = ex
    res['R'] = d['R'].c0()
    return res


def minkowski(x0='0.3'):
    x0 = arb(x0)
    one = [arb(1)] + [arb(0)]*5; zero = [arb(0)]*6
    op = CH.polar(x0, one, zero, zero, zero, arb(1)/3)
    oc = CH.cart(x0.exp(), 'x', x0, one, zero, zero, zero, arb(1)/3)
    return op, oc


def frw_series(z0, K=6):
    """frw.py k = 1/2 closed form -> x-series (about x(z0)) of F2, At, Wt, Vt, plus x0"""
    q = arb(4)/9; c = q*(1 - q)
    zs = Ser.var(arb(z0), K)
    Vz = zs*(-q)
    Nz = (1 + zs*zs*c)/zs
    Az = (1 - zs*zs*q*q).inv()
    oz = zs*zs*q*q*arb('1.5')/(1 - zs*zs*q*q)
    xz = zs.log() - (1 + zs*zs*c).log()/(2*(1 - q))
    ex = xz.exp()                       # e^x = z_KHA
    F2 = Az*Nz*Nz*ex*ex
    At = (Az - 1)/(ex*ex)
    Wt = oz/(ex*ex)
    Vt = Vz/ex
    x0 = xz.c[0]
    gi = revert(xz.shift0())            # h(dx)
    conv = lambda f: f.shift0().compose_into(gi) + f.c[0]
    return x0, [conv(f).c for f in (F2, At, Wt, Vt)]


def frw(z0):
    x0, (cF2, cAt, cWt, cVt) = frw_series(z0)
    op = CH.polar(x0, cF2, cAt, cWt, cVt, arb(1)/2, want_routeA=False)
    oc = CH.cart(x0.exp(), 'x', x0, cF2, cAt, cWt, cVt, arb(1)/2, want_routeA=False)
    zf = arb(z0)
    pred = 2*zf**4/243
    return op, oc, pred, -8*zf**2/27
