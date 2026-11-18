"""build-2 QB1b: b_nn on an exact rational x-cell [a, b] in Taylor-model arithmetic, eps = x - x_m, |eps| <= h.

x is EXACT in every cell (x = x_m + eps), so no base-point remainder enters (the build-1 polar/cart charts carry
e^{x0}; with a TM x0 its remainder is amplified by the cancellations, hence the direct s-chart, schart.py).
  record cells (tube E5 / outward BG2): state TM from the step polynomial + order-0 record tail, ODE x-jets in TM.
  manifold cells (centre cm, sonic sm): the parameter theta(x_m + eps) is solved as a TM from x(theta) = x by a
      contraction (Krawczyk-type bound: residual / (1 - L), L = sup|1 - x'(theta)/d0| on a widened interval), then
      theta-jets (shifted polynomial + per-order tail bound) and reversion of x(theta) give the x-jets.
"""
import sys
from math import comb
from flint import arb, ctx
ctx.cap = 64
sys.path.insert(0, '<LAB>/quantum/seats/build-1')
import tm
from tm import TM
tm.install()
import bg
from aser import Ser, revert
from schart import schart, profiles_direct
from tmbg import tail_bound

KX = 5


def kap():
    return arb(1)/3


def jets_at(piece, ps, J, comps=None):
    """Taylor coefficients 0..J of the components at the parameter ps (TM or arb), as TMs / arbs"""
    isT = isinstance(ps, TM)
    pa = ps.rng().abs_upper() if isT else arb(ps).abs_upper()
    if not pa < piece.R:
        raise ArithmeticError('parameter leaves the disc')
    comps = range(len(piece.C)) if comps is None else comps
    out = []
    for i in comps:
        c = piece.C[i]; K = piece.K; row = []
        for j in range(J + 1):
            acc = TM([arb(0)]) if isT else arb(0)
            for n in range(K, j - 1, -1):
                acc = acc*ps + c[n]*comb(n, j)
            t = tail_bound(piece, i, j, pa).abs_upper()
            row.append(TM(acc.c, acc.r + t) if isT else acc + arb(0, t))
        out.append(row)
    return out


def solve_theta(piece, xfun_jets, xm, th_guess):
    """theta TM with x(theta(eps)) = xm + eps.  xfun_jets(ps, J) -> list of x-jets at ps (TM or arb)."""
    xt = TM.var(xm)
    p0 = arb(th_guess)
    d0 = xfun_jets(p0, 1)[1]
    d0 = arb(d0.mid())
    th = TM([p0])
    for _ in range(8):
        th = th - (TM(xfun_jets(th, 0)[0].c) - xt)/d0 if False else th - (xfun_jets(th, 0)[0] - xt)/d0
    res = (xfun_jets(th, 0)[0] - xt)/d0
    Rr = res.rng().abs_upper()
    R0 = th.rng()
    delta = 4*Rr + arb('1e-60')
    I = R0 + arb(0, delta.upper())
    dI = xfun_jets(I, 1)[1]
    L = (1 - dI/d0).abs_upper()
    if not L < arb('0.5'):
        raise ArithmeticError('theta solve: no contraction (L = %s)' % L.str(5))
    e = Rr/(1 - L)
    if not e < delta:
        raise ArithmeticError('theta solve: enclosure not inside the contraction interval')
    if not dI > 0:
        raise ArithmeticError('theta solve: x(theta) not monotone on the interval')
    return TM(th.c, th.r + e.abs_upper()), dict(L=L, res=Rr)


def to_xseries(pj, xj, K):
    X = Ser([arb(0)] + list(xj[1:]), K)
    g = revert(X)
    return [Ser(c, K).compose_into(g).c for c in pj]


def centre_cell(Bg, xm, h, th_guess):
    Pc, xc = Bg.cm_piece()

    def xfun(ps, J):
        jj = jets_at(Pc, ps, J, comps=[3])[0]
        S_ = Ser(jj, J) if J > 0 else None
        if J == 0:
            return [ps.log() + jj[0] + xc]
        return [ps.log() + jj[0] + xc, 1/ps + jj[1]]
    th, info = solve_theta(Pc, xfun, xm, th_guess)
    K = KX
    jj = jets_at(Pc, th, K)
    Pj, Wj, Uj, Xj = (Ser(c, K) for c in jj)
    tser = Ser.var(th, K)
    xs = tser.log() + Xj + xc
    k = kap(); one = Ser.const(1, K)
    PU = Pj*Uj
    om = Wj*Pj*Pj
    A = one + om*2*(one + Uj*(k + 1) + PU*PU*k)/(one - PU*PU)
    G = A/(Pj*Pj)
    V = PU
    prof = to_xseries([G.c, A.c, om.c, V.c], xs.c, K)
    return prof, info


def sonic_cell(Bg, xm, h, th_guess, which='tube'):
    Ps = Bg.sm_piece(which)

    def xfun(ps, J):
        return jets_at(Ps, ps, J, comps=[3])[0]
    th, info = solve_theta(Ps, xfun, xm, th_guess)
    K = KX
    jj = jets_at(Ps, th, K)
    cN, co, cV = to_xseries(jj[:3], jj[3], K)
    return profiles_direct(cN, co, cV, K, Ser, kap()), info


def record_cell(Bg, which, xm, h):
    steps = Bg.tube_xsteps() if which == 'tube' else Bg.outw_xsteps()
    rec = Bg.find_step(steps, xm)
    piece = Bg.tube_piece(rec) if which == 'tube' else Bg.outw_piece(rec)
    st = jets_at(piece, TM.var(xm - rec['x0']), 0)
    y = [st[i][0] for i in range(3)]
    cN, co, cV = Bg.ode_jets(y, KX)
    return profiles_direct(cN, co, cV, KX, Ser, kap()), dict(step_x0=rec['x0'])


def cell(Bg, region, a, b, deg=4, want_dx=False, routeA=False, th_guess=None):
    a, b = arb(a), arb(b)
    xm = (a + b)/2; h = (b - a)/2
    tm.H = h
    if region == 'centre':
        prof, info = centre_cell(Bg, xm, h, th_guess)
    elif region == 'sonic':
        prof, info = sonic_cell(Bg, xm, h, th_guess)
    else:
        prof, info = record_cell(Bg, region, xm, h)
    o = schart(*prof, kap(), deg=deg, want_routeA=routeA, want_dx=want_dx)
    o['_info'] = info
    return o


def guess_theta(Bg, region, x):
    """float-level Newton (arb point arithmetic; only a starting value, not part of any enclosure)"""
    import math
    x = arb(x)
    if region == 'centre':
        Pc, xc = Bg.cm_piece()
        f = lambda t, J: [t.log() + jets_at(Pc, t, J, [3])[0][0] + xc] + ([1/t + jets_at(Pc, t, 1, [3])[0][1]] if J else [])
        t = arb(math.exp(float((x - xc).mid())))
    else:
        Ps = Bg.sm_piece('tube')
        f = lambda t, J: jets_at(Ps, t, J, [3])[0]
        t = arb(x.mid())
    for _ in range(30):
        v = f(t, 1)
        t = arb((t - (v[0] - x)/v[1]).mid())
    return t


def centre_cell_cart(Bg, a, b, th_guess, deg=4, routeA=False):
    """centre cell through the build-1 Cartesian chart (t = -1 scale, regular at the centre); x0 exact TM"""
    import charts as CH
    a, b = arb(a), arb(b)
    xm = (a + b)/2; tm.H = (b - a)/2
    Pc, xc = Bg.cm_piece()

    def xfun(ps, J):
        jj = jets_at(Pc, ps, J, comps=[3])[0]
        if J == 0:
            return [ps.log() + jj[0] + xc]
        return [ps.log() + jj[0] + xc, 1/ps + jj[1]]
    th, info = solve_theta(Pc, xfun, xm, th_guess)
    K = KX
    jj = jets_at(Pc, th, K)
    Pj, Wj, Uj, Xj = (Ser(c, K) for c in jj)
    tser = Ser.var(th, K)
    xs = tser.log() + Xj + xc
    k = kap(); one = Ser.const(1, K)
    E = (Xj + xc).exp()                 # z/theta
    Q = Pj/tser                         # P/theta
    PU = Pj*Uj
    At = Wj*(Q/E)*(Q/E)*2*(one + Uj*(k + 1) + PU*PU*k)/(one - PU*PU)
    z = tser*E
    A = one + At*z*z
    F2 = A*(E/Q)*(E/Q)
    Wt = Wj*(Q/E)*(Q/E)
    Vt = Uj*Q/E
    prof = to_xseries([F2.c, At.c, Wt.c, Vt.c], xs.c, K)
    x0 = TM.var(xm)
    o = CH.cart(x0.exp(), 'x', x0, *prof, kap(), deg=deg, want_routeA=routeA)
    o['_info'] = info
    return o


# ------------------------------------------------------------------ deep centre: w = z^2 series with a Cauchy tail
_WP = {}


def wpiece(Bg, Kw=26, rho_th='0.02'):
    """Piece over w: exact w-Taylor coefficients (n <= Kw) of F2, At, Wt, Vt at w = 0 (bg.centre_profiles_w, from the
    stored cm coefficients), tail |g_n| <= Mg/Rw^n for n > Kw (Cauchy on |w| <= Rw).
    Mg = sup |g(theta)| over the complex box |Re|,|Im| <= rho_th (contains the disc |theta| <= rho_th, corners inside
    the cm tail disc |theta| <= 0.03); Rw = rho_th^2 e^{2 xc}/4.  Rouche on |theta| = rho_th
    (|e^{2X(theta)} - 1| < 3/4 checked) gives exactly two roots +-theta(w) of theta^2 E(theta)^2 = w for |w| <= Rw,
    so G(w) = g(theta(w)) (g even) is analytic there and |G| <= Mg."""
    if 'p' in _WP:
        return _WP['p']
    from flint import acb
    from aser import Piece
    outs, oddmax, xc = Bg.centre_profiles_w(Kw)
    Pc, _ = Bg.cm_piece()
    rho = arb(rho_th)
    th = acb(arb(0, rho), arb(0, rho))
    nu = Pc.R
    assert (rho*arb(2).sqrt()) < nu
    eps = Pc.M[0]; epsX = Pc.M[3]
    def hz(cs, t):
        s = acb(0)
        for c in reversed(cs): s = s*t + c
        return s
    tail = acb(arb(0, eps.upper()), arb(0, eps.upper()))
    P_ = hz(Pc.C[0], th) + tail; W_ = hz(Pc.C[1], th) + tail; U_ = hz(Pc.C[2], th) + tail
    X_ = hz(Pc.C[3], th) + acb(arb(0, epsX.upper()), arb(0, epsX.upper()))
    Qe = arb(0, (eps/nu).upper())
    Q_ = hz(Pc.C[0][1:], th) + acb(Qe, Qe)          # P/theta; |E_P/theta| <= eps/nu (Schwarz)
    E_ = (X_ + xc).exp()
    k = kap()
    PU = P_*U_
    At = W_*(Q_/E_)**2*2*(1 + U_*(k + 1) + PU*PU*k)/(1 - PU*PU)
    z = th*E_
    A_ = 1 + At*z*z
    F2 = A_*(E_/Q_)**2
    Wt = W_*(Q_/E_)**2
    Vt = U_*Q_/E_
    rou = abs((2*X_).exp() - 1)
    assert rou < arb('0.75'), 'Rouche condition fails'
    Mg = [abs(v).upper() for v in (F2, At, Wt, Vt)]
    Rw = rho**2*(2*xc).exp()/4
    Rw = arb(Rw.lower())
    p = Piece([c[:Kw + 1] for c in outs], Mg, Rw, 'sch', K=Kw)
    _WP['p'] = p; _WP['info'] = dict(Mg=[m.str(5) for m in Mg], Rw=Rw.str(10), Kw=Kw, rouche=rou.str(5), oddmax=oddmax.str(3))
    return p


def wcell(Bg, a, b, deg=4, routeA=False):
    """x-cell [a, b] with b <= -11/2: w-cell [lower e^{2a}, upper e^{2b}], Cartesian w-mode (regular, O(1) scale)"""
    import charts as CH
    a, b = arb(a), arb(b)
    p = wpiece(Bg)
    wa = arb((2*a).exp().lower()); wb = arb((2*b).exp().upper())
    wm = (wa + wb)/2; tm.H = (wb - wa)/2
    w0 = TM.var(wm)
    jj = jets_at(p, w0, deg)
    zeta = w0.sqrt()
    o = CH.cart(zeta, 'w', w0, *jj, kap(), deg=deg, want_routeA=routeA)
    o['_info'] = dict(w=[wa.str(10), wb.str(10)])
    return o
