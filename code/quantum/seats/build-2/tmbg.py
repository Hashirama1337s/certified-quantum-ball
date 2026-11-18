"""build-2 QB1b: Taylor-model (cell variable eps) evaluation of b_nn on one closed cell, on top of the pinned build-1
modules (imported, never edited; tm.install() rebinds constructor names only).

Cell kinds
  record  x-cell [a, b] (exact rationals) inside the Cauchy disc of one x-step of the tube (E5) or outward (BG2) record:
          state TM y(x_m + eps) = shifted step polynomial (TM Horner) + record tail bound (order 0),
          x-jets from the ODE y' = f(y) in TM arithmetic (taylor.jets with TM entries), profiles (bg), chart.
  manifold  theta-cell of the centre manifold cm or the sonic manifold sm: theta-jets as TMs (shifted polynomial +
          per-order tail bound), profiles, reversion x(theta) (all in TM arithmetic), chart.  dx/dtheta > 0 is
          checked on the cell, so x(theta) is monotone there.
"""
import os, sys
from math import comb
from flint import arb, ctx
ctx.cap = 64
sys.path.insert(0, '<LAB>/quantum/seats/build-1')
import tm
from tm import TM
tm.install()
import bg
import charts as CH
from aser import Ser, revert

KX = 7


def kap():
    return arb(1)/3


def tail_bound(piece, i, j, pa):
    """bound of |j-th Taylor coefficient of the tail E_i| at every |p| <= pa (aser.Piece formula)"""
    if piece.kind is None:
        return arb(0)
    q = arb(pa)/piece.R
    K = piece.K
    s = arb(0); n = K + 1
    term = comb(n, j)*q**(n - j)
    while True:
        rq = q*(n + 1)/(n + 1 - j)
        if rq.upper() <= arb('0.9') or n > K + 20000:
            break
        s += term; n += 1; term = term*rq
    if not rq < 1:
        raise ArithmeticError('tail ratio not < 1')
    s = s + term/(1 - rq)
    return piece.M[i]/piece.R**j*s


def piece_jets_tm(piece, pm, J, comps=None):
    """TM (in eps, p = pm + eps, |eps| <= tm.H) of the Taylor coefficients 0..J of each component"""
    pa = arb(pm).abs_upper() + tm.H
    if not pa < piece.R:
        raise ArithmeticError('cell leaves the disc')
    ps = TM.var(pm)
    comps = range(len(piece.C)) if comps is None else comps
    out = []
    for i in comps:
        c = piece.C[i]; K = piece.K
        row = []
        for j in range(J + 1):
            acc = TM([arb(0)])
            for n in range(K, j - 1, -1):
                acc = acc*ps + c[n]*comb(n, j)
            t = tail_bound(piece, i, j, pa)
            row.append(TM(acc.c, acc.r + t.abs_upper()))
        out.append(row)
    return out


def xseries_tm(pj, xj, K):
    X = Ser(xj, K)
    if not X.c[1] > 0:
        raise ArithmeticError('dx/dp not certified positive on the cell')
    g = revert(X.shift0())
    return xj[0], [Ser(c, K).compose_into(g).c for c in pj]


def chart_eval(x0, prof, chart, deg=4, want_dx=False, routeA=False):
    if chart == 'cart':
        return CH.cart(x0.exp(), 'x', x0, *prof, kap(), deg=deg, want_routeA=routeA)
    return CH.polar(x0, *prof, kap(), deg=deg, want_dx=want_dx, want_routeA=routeA)


def eval_record_cell(Bg, which, a, b, chart, deg=4, want_dx=False):
    """x-cell [a, b] (exact) -> charts output (TM values)"""
    a, b = arb(a), arb(b)
    xm = (a + b)/2; h = (b - a)/2
    tm.H = h
    steps = Bg.tube_xsteps() if which == 'tube' else Bg.outw_xsteps()
    rec = Bg.find_step(steps, xm)
    piece = Bg.tube_piece(rec) if which == 'tube' else Bg.outw_piece(rec)
    st = piece_jets_tm(piece, xm - rec['x0'], 0)
    y = [st[i][0] for i in range(3)]
    cN, co, cV = Bg.ode_jets(y, KX)
    x0 = TM.var(xm)
    prof, chk = bg.profiles_from_NomV(x0, cN, co, cV, KX)
    o = chart_eval(x0, prof, chart, deg, want_dx)
    o['_chk'] = chk
    return o


def centre_profiles_tm(Bg, thm, K):
    """TM version of bg.Background.centre_profiles_theta (same formulas) on the theta-cell thm +- tm.H"""
    Pc, xc = Bg.cm_piece()
    jj = piece_jets_tm(Pc, thm, K)
    Pj, Wj, Uj, Xj = (Ser(c, K) for c in jj)
    th = Ser.var(TM.var(thm), K)
    Q = Pj/th
    E = (Xj + xc).exp()
    k = kap(); one = Ser.const(1, K)
    PU = Pj*Uj
    At = Wj*(Q/E)*(Q/E)*2*(one + Uj*(k + 1) + PU*PU*k)/(one - PU*PU)
    z = th*E
    A = one + At*z*z
    F2 = A*(E/Q)*(E/Q)
    Wt = Wj*(Q/E)*(Q/E)
    Vt = Uj*Q/E
    xs = th.log() + Xj + xc
    return xseries_tm([F2.c, At.c, Wt.c, Vt.c], xs.c, K)


def eval_centre_cell(Bg, tha, thb, chart='cart'):
    tha, thb = arb(tha), arb(thb)
    tm.H = (thb - tha)/2
    x0, prof = centre_profiles_tm(Bg, (tha + thb)/2, KX)
    return chart_eval(x0, prof, chart)


def eval_sonic_cell(Bg, tha, thb, which='tube', deg=4, want_dx=False):
    tha, thb = arb(tha), arb(thb)
    tm.H = (thb - tha)/2
    Ps = Bg.sm_piece(which)
    jj = piece_jets_tm(Ps, (tha + thb)/2, KX)
    x0, (cN, co, cV) = xseries_tm(jj[:3], jj[3], KX)
    prof, chk = bg.profiles_from_NomV(x0, cN, co, cV, KX)
    return chart_eval(x0, prof, 'polar', deg, want_dx)


def bnn(o):
    v = o['hat_b_nn'] if o.get('chart') == 'cart' else o['b_nn']
    return tm.rng(v)
