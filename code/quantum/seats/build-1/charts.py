"""build-1 QB1: the two metric-in / tensor-out adapters on top of ajet (one interval code path for EC and controls).

Both take the self-similar background as Taylor coefficient lists (arb balls) of four regular profile functions
in an argument variable (KHA x, or w = z^2 at the centre):
    F2 = A N^2 z^2     (-g_tt in polar-areal coordinates; lapse^2)
    At = (A - 1)/z^2
    Wt = om/z^2        (om = 4 pi r^2 A rho)
    Vt = V/z           (V = fluid 3-velocity relative to the Eulerian observer, u = Gamma (n + V e_r))
with z = r/(-t) = e^x.  The KHA line element r^2{A[(1-N^2)ds^2 - 2ds dx + dx^2] + dOmega^2} with s = -ln(-t),
x = ln r + s is exactly  -A N^2 z^2 dt^2 + A dr^2 + r^2 dOmega^2  (derivation in RUN.md), so the Eulerian frame
n (normal to t = const = s = const) and e_r agree with DESIGN-QB1's KHA frame.

polar(...)  coordinates (t, r, theta, phi), base point t0 = -e^{-x0}, r0 = 1: b_ab = r0^4 B_(ab) directly.
cart(...)   coordinates (t, X, Y, Z), base point t0 = -1, X0 = zeta = e^{x0}, Y = Z = 0, metric
            -F2 dt^2 + (delta_ij + At X_i X_j/t^2) dX^i dX^j; regular at the centre (zeta = 0 allowed in w-mode);
            returns Bhat_(ab) = B_(ab) at t = -1 = b_ab e^{-4x}, and b_ab = zeta^4 Bhat_(ab).
"""
from flint import arb, ctx
import ajet as J
from ajet import Jet, n4

ctx.cap = 64
PI = None


def _pi():
    return arb.pi()


def _fluid_and_out(d, g, uup, rho, k, frames, scale4, scale2, want_routeA=True, flip=False):
    """common tail: Zahn, Bach, fluid T, Route A, projections.  frames: dict name -> 4-vector at base point."""
    gi = d['gi']
    T1, T2 = J.zahn(d, flip=flip)
    Bc = J.bach(T1, T2)
    p = rho*k
    # u_a = g_ab u^b
    udn = [J._sum(g[a][b]*uup[b] for b in n4 if g[a][b].c and uup[b].c) for a in n4]
    Tdn = [[(rho + p)*udn[a]*udn[b] + p*g[a][b] for b in n4] for a in n4]
    trT = p*3 - rho
    eightpi = 8*_pi()
    RicS = [[(Tdn[a][b] - g[a][b]*trT*arb('0.5'))*eightpi for b in n4] for a in n4]
    P = J.proj
    nv, ev, tv, uv = frames['n'], frames['e'], frames['t'], frames['u']
    out = {}
    out['norm_n'] = P(g, nv, nv); out['norm_u'] = P(g, uv, uv); out['n_e'] = P(g, nv, ev)
    hat = dict(nn=P(Bc, nv, nv), nr=P(Bc, nv, ev), rr=P(Bc, ev, ev), T=P(Bc, tv, tv), uu=P(Bc, uv, uv))
    for kk, v in hat.items():
        out['b_' + kk] = v*scale4
    out['trB'] = sum((gi[a][b].c0()*Bc[a][b].c0() for a in n4 for b in n4), arb(0))*scale4
    out['trace_dev'] = out['b_T'] - (out['b_nn'] - out['b_rr'])/2
    out['R'] = d['R'].c0()*scale2
    out['T1_uu'] = P(T1, uv, uv)*scale4; out['T1_nn'] = P(T1, nv, nv)*scale4
    out['V0_uu'] = P(J.Vxi(T1, T2, 0), uv, uv)*scale4
    out['V6_uu'] = P(J.Vxi(T1, T2, arb(1)/6), uv, uv)*scale4
    out['xi_diff_uu'] = out['V0_uu'] - out['V6_uu']
    E = [[d['Ric'][a][b] - RicS[a][b] for b in n4] for a in n4]
    out.update(E_nn=P(E, nv, nv)*scale2, E_nr=P(E, nv, ev)*scale2, E_rr=P(E, ev, ev)*scale2, E_T=P(E, tv, tv)*scale2)
    out['rho8pi_r2'] = (eightpi*rho).c0()*scale2
    if want_routeA:
        VA = J.routeA(d, RicS)
        A60 = [[VA[a][b]*60 for b in n4] for a in n4]
        out.update(a_nn=P(A60, nv, nv)*scale4, a_nr=P(A60, nv, ev)*scale4, a_rr=P(A60, ev, ev)*scale4,
                   a_T=P(A60, tv, tv)*scale4, a_uu=P(A60, uv, uv)*scale4)
    return out, Bc


def polar(x0, cF2, cAt, cWt, cVt, k, deg=4, want_routeA=True, want_dx=False, flip=False, tchart='t'):
    """x-series inputs about x0 (x0 arb, point or interval).  deg 5 + want_dx: also d b_nn/dx.
    tchart 't': polar-areal time t (g_tt = -F2), base t0 = -e^{-x0}, r0 = 1;
    tchart 's': KHA time s = -ln(-t) (g_ss = -F2 e^{-2x} r^2 = -A N^2 r^2), base s0 = x0, r0 = 1.  Same slices
    (s is a function of t), so n, e_r, u and every b_ab are identical; 's' avoids the e^{x0} time scale at large x."""
    J.setup(3, deg)
    k = arb(k)
    x0 = arb(x0)
    dt, dr, dth = Jet.var(0), Jet.var(1), Jet.var(2)
    ex0 = x0.exp()
    one = Jet.const(1)
    if tchart == 't':
        dx = (one + dr).log() - (one - dt*ex0).log()
    else:
        dx = (one + dr).log() + dt
    F2 = J.horner(cF2, dx); At = J.horner(cAt, dx); Wt = J.horner(cWt, dx); Vt = J.horner(cVt, dx)
    ex = dx.exp()*ex0
    e2x = ex*ex
    A = one + At*e2x
    om = Wt*e2x
    V = Vt*ex
    r = one + dr
    sn = (dth + _pi()/2).sin()
    if tchart == 's':
        F2 = F2/e2x*r*r                  # -g_ss
    g = [[Jet() for _ in n4] for _ in n4]
    g[0][0] = -F2; g[1][1] = A; g[2][2] = r*r; g[3][3] = r*r*sn*sn
    d = J.curv(g)
    Gm = (one - V*V).pow_half_neg()
    sF = F2.sqrt(); sA = A.sqrt()
    uup = [Gm/sF, Gm*V/sA, Jet(), Jet()]
    rho = om/(r*r*A*4*_pi())
    F20, A0 = F2.c0(), A.c0()
    frames = dict(n=[1/F20.sqrt(), 0, 0, 0], e=[0, 1/A0.sqrt(), 0, 0], t=[0, 0, 1, 0], u=[uup[0].c0(), uup[1].c0(), 0, 0])
    out, Bc = _fluid_and_out(d, g, uup, rho, k, frames, arb(1), arb(1), want_routeA, flip)
    out['chart'] = 'polar'
    out['inputs0'] = dict(F2=F20, A=A0, om=om.c0(), V=V.c0())
    if want_dx:
        assert deg >= 5
        Jn = r*r*r*r*Bc[0][0]/F2
        out['db_nn_dx'] = Jn.coef((0, 1, 0))      # d/dx = r d/dr at fixed t, r0 = 1
        # J(dr) = b_nn(x0 + ln(1 + dr)) exactly to degree deg - 4 (Bach jets exact to degree deg - 4)
        out['J_dr'] = [Jn.coef((0, j, 0)) for j in range(deg - 3)]
    return out


def cart(zeta, mode, arg0, cF2, cAt, cWt, cVt, k, deg=4, want_routeA=True):
    """mode 'x': coefficient lists in (x - arg0), zeta = e^{arg0} > 0 (zeta, arg0 arb);
       mode 'w': coefficient lists in (w - arg0), w = z^2, arg0 = zeta^2 (zeta = 0 allowed)."""
    J.setup(4, deg)
    k = arb(k)
    zeta = arb(zeta)
    dt, dX, dY, dZ = (Jet.var(i) for i in range(4))
    one = Jet.const(1)
    X = dX + zeta; Y = dY; Z = dZ
    L = X*X + Y*Y + Z*Z
    mt = one - dt                                   # -t
    mt2 = mt*mt
    w = L/mt2
    if mode == 'x':
        # L/zeta^2 = 1 + (2 zeta dX + dX^2 + dY^2 + dZ^2)/zeta^2 with an exact constant 1 (zeta may be an interval)
        Lrel = one + (dX*(2*zeta) + dX*dX + dY*dY + dZ*dZ)*(1/(zeta*zeta))
        dx = Lrel.log()*arb('0.5') - mt.log()
        arg = dx
    else:
        arg = w.nil()
    F2 = J.horner(cF2, arg); At = J.horner(cAt, arg); Wt = J.horner(cWt, arg); Vt = J.horner(cVt, arg)
    A = one + At*w
    Xs = [X, Y, Z]
    g = [[Jet() for _ in n4] for _ in n4]
    g[0][0] = -F2
    imt2 = mt2.inv()
    for i in range(3):
        for j in range(3):
            v = At*Xs[i]*Xs[j]*imt2
            if i == j:
                v = v + 1
            g[i + 1][j + 1] = v
    d = J.curv(g)
    Gm = (one - Vt*Vt*w).pow_half_neg()
    sF = F2.sqrt(); sA = A.sqrt()
    fac = Gm*Vt/(mt*sA)
    uup = [Gm/sF, fac*X, fac*Y, fac*Z]
    rho = Wt/(mt2*A*4*_pi())
    F20, A0 = F2.c0(), A.c0()
    frames = dict(n=[1/F20.sqrt(), 0, 0, 0], e=[0, 1/A0.sqrt(), 0, 0], t=[0, 0, 1, 0],
                  u=[uup[0].c0(), uup[1].c0(), 0, 0])
    z4 = zeta**4; z2 = zeta**2
    out, Bc = _fluid_and_out(d, g, uup, rho, k, frames, arb(1), arb(1), want_routeA)
    # out values are at t = -1 scale (hat); b = zeta^4 hat, r^2-scaled Einstein = zeta^2 hat
    hatk = [kk for kk in out if kk.startswith(('b_', 'a_', 'T1_', 'V0_', 'V6_', 'xi_')) or kk in ('trB', 'trace_dev')]
    res = {'chart': 'cart', 'zeta': zeta}
    for kk in hatk:
        res['hat_' + kk] = out[kk]
        res[kk] = out[kk]*z4
    for kk in ('R', 'E_nn', 'E_nr', 'E_rr', 'E_T', 'rho8pi_r2'):
        res['hat_' + kk] = out[kk]
        res[kk] = out[kk]*z2
    for kk in ('norm_n', 'norm_u', 'n_e'):
        res[kk] = out[kk]
    res['inputs0'] = dict(F2=F20, A=A0)
    return res
