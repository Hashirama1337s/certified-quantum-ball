"""build-2 QB1b: 'direct s-chart' adapter for Taylor-model cells (no base-point factor e^{x0} anywhere).

Coordinates (s, r, theta, phi), s = KHA time, base point s0 = x0, r0 = 1, so delta x = ln(1 + dr) + ds exactly and
the cell centre x0 never enters the jets.  Metric (identical to KHA/build-1 charts):
    g_ss = -G r^2,  G = A N^2;   g_rr = A;   g_thth = r^2;  g_phph = r^2 sin^2 theta.
Inputs: x-Taylor lists (arb or TM) of G, A, om, V.  Fluid: rho = om/(4 pi r^2 A), u = Gamma (n + V e_r),
n = d_s/(r sqrt G), e_r = d_r/sqrt A.  Tensor algebra, Route A and projections: build-1 ajet / charts._fluid_and_out
(imported, unedited).  b_ab = r0^4 B_(ab) = B_(ab).
"""
import sys
sys.path.insert(0, '<LAB>/quantum/seats/build-1')
from flint import arb, ctx
import ajet as J
from ajet import Jet, n4
import charts as CH
ctx.cap = 64


def schart(cG, cA, co, cV, k, deg=4, want_routeA=True, want_dx=False):
    J.setup(3, deg)
    ds, dr, dth = Jet.var(0), Jet.var(1), Jet.var(2)
    one = Jet.const(1)
    dx = (one + dr).log() + ds
    G = J.horner(cG, dx); A = J.horner(cA, dx); om = J.horner(co, dx); V = J.horner(cV, dx)
    r = one + dr
    sn = (dth + arb.pi()/2).sin()
    g = [[Jet() for _ in n4] for _ in n4]
    g[0][0] = -(G*r*r); g[1][1] = A; g[2][2] = r*r; g[3][3] = r*r*sn*sn
    d = J.curv(g)
    Gm = (one - V*V).pow_half_neg()
    sG = G.sqrt(); sA = A.sqrt()
    uup = [Gm/(sG*r), Gm*V/sA, Jet(), Jet()]
    rho = om/(r*r*A*4*arb.pi())
    G0, A0 = G.c0(), A.c0()
    frames = dict(n=[1/G0.sqrt(), 0, 0, 0], e=[0, 1/A0.sqrt(), 0, 0], t=[0, 0, 1, 0], u=[uup[0].c0(), uup[1].c0(), 0, 0])
    out, Bc = CH._fluid_and_out(d, g, uup, rho, k, frames, arb(1), arb(1), want_routeA)
    out['chart'] = 'schart'
    if want_dx:
        Jn = r*r*r*r*Bc[0][0]/(G*r*r)          # b_nn along r at fixed s; d/dx = r d/dr there (r0 = 1)
        out['db_nn_dx'] = Jn.coef((0, 1, 0))
        out['J_dr'] = [Jn.coef((0, j, 0)) for j in range(deg - 3)]
    return out


def profiles_direct(cN, co, cV, K, Ser, k):
    """(N, om, V) x-lists -> G = A N^2, A = 1 + 2 om (1 + (k+1) N V + k V^2)/(1 - V^2), om, V (x-lists)"""
    N = Ser(cN, K); om = Ser(co, K); V = Ser(cV, K)
    one = Ser.const(1, K)
    A = one + om*2*(one + N*V*(k + 1) + V*V*k)/(one - V*V)
    G = A*N*N
    return [G.c, A.c, om.c, V.c]
