"""eqs_l.py - build seat build-1 (B7 "No Sideways Fall", SB0 + SB1), 2026-11-17.

Own symbolic derivation (sympy) of the linear non-spherical perturbation equations of Einstein-Euler,
p = k rho, about the CSS spherical background, in KHA self-similar variables (s, x), growth exp(kappa s).
Engine: lin_core.py (first-order variation of the Einstein tensor and of div T, background as x-jets).
Gundlach gr-qc/9906124 (OPENED, seats/audit-1/gr-qc_9906124.txt) is used ONLY as a check, never as a source.

Sectors and formulations (stated choices):
  AX-beta  axial fluid, every l >= 1: covariant du_phi = R E b(x) sin(th) Y'(th).
  AX-Pi    axial metric, l >= 2, RW gauge h_{ab}^{axial} = 0: h_{A phi} = R^2 E h_A(x) sin Y';
           gauge-invariant curl Q = kappa h_x - h_s' (Q = A N Pi up to the factor E R^-2), master ODE for Q.
  P1       polar l = 1, gauge h_{A theta} = 0 and K = 0 ('k = 0'), plus the alpha-gauge condition at the centre.
  P(l)     polar l >= 2, Regge-Wheeler gauge (h_{A theta} = 0, G = 0); unknowns h_AB = R^2 E H_AB(x) Y, K.
  TOY      massless scalar wave equation on the flat background (N = e^-x, A = 1), G App. D check.
Fluid perturbation: drho = rho E dl(x) Y, du_x = R E v(x) Y, du_theta = R E w(x) Y', du_s from u.u = -1.

Outputs (python main): eqs_l.pkl (sympy srepr strings only) and eqs_checks.json (every symbolic check).
usage: memcap 8: py -3.11 eqs_l.py
"""
import sys, os, json, time, pickle, hashlib, random
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, '<LAB>/engine')
import sympy as sp
import lin_core as LC
from lin_core import R, E, Y0, Y1, cs, sn, D

kap, k, L = LC.kap, LC.k, LC.L
N0, A0, om0, V0 = LC.N0, LC.A0, LC.om0, LC.V0
N1, A1, om1, V1 = LC.N1, LC.A1, LC.om1, LC.V1
P_, W_, U_ = sp.symbols('P W U')


# ======================================================================  axial beta
def beta_equation():
    B = LC.Lin()
    hs = LC.pert('hs'); hx = LC.pert('hx'); b = LC.pert('b')
    h = sp.zeros(4, 4)
    h[0, 3] = h[3, 0] = R**2*E*hs[0]*sn*Y1
    h[1, 3] = h[3, 1] = R**2*E*hx[0]*sn*Y1
    du = [0, 0, 0, R*E*b[0]*sn*Y1]
    fl = B.fluid(h, 0, du, comps=(3,))
    F3 = sp.cancel(LC.finalize(fl[3]/(sn*Y1)))
    F3b = sp.cancel(LC.finalize(fl[3]/(sn*Y1), triple=(sp.Rational(5, 13), sp.Rational(12, 13))))
    chk = {'angular_triple_consistent': sp.cancel(F3 - F3b) == 0,
           'metric_free': all(sp.diff(F3, s_) == 0 for s_ in hs[:3] + hx[:3]),
           'l_free': sp.diff(F3, L) == 0}
    Kf = om0*(1 + k)/(A0*N0*LC.gV*LC.sA)
    a = sp.cancel(sp.diff(F3, b[1])/Kf)
    kc = sp.cancel(sp.diff(sp.diff(F3, b[0]), kap)/Kf)
    c = sp.cancel(sp.diff(F3, b[0]).subs(kap, 0)/Kf)
    chk['form K (kappa b + a b_x + c b)'] = (kc == 1) and sp.cancel(F3 - Kf*(kap*b[0] + a*b[1] + c*b[0])) == 0
    chk['a == 1 + N V'] = sp.cancel(a - (1 + N0*V0)) == 0
    c_form = (1 + N0*V0)*(om1/om0 - A1/(2*A0)) + 2*N0*V0 + N1*V0 + V1*(N0 + V0)/(1 - V0**2)
    chk['c == (1+NV)(om_x/om - A_x/2A) + 2NV + N_x V + V_x (N+V)/(1-V^2)'] = sp.cancel(c - c_form) == 0
    return dict(a=a, c=sp.cancel(c_form), check=chk)


def centre_chart_subs():
    """background x-jets (orders 0, 1) in the blown-up centre chart (P, W, U) at general k:
    N = 1/P, om = W P^2, V = U P, A = Acons(N, om, V), d/dx = G (engine/system.py rhs_centre_blown)."""
    import system as S
    G = S.rhs_centre_blown()                       # u dZ/du = G(Z), u = e^x  ->  d/dx Z = G
    kk = S.k
    Zs = [S.P, S.W, S.U]
    cen = {S.N: 1/S.P, S.om: S.W*S.P**2, S.V: S.U*S.P}
    Acen = sp.cancel(S.Acons.subs(cen))
    def dx(f):
        return sum(sp.diff(f, z)*g for z, g in zip(Zs, G))
    vals = {N0: 1/S.P, om0: S.W*S.P**2, V0: S.U*S.P, A0: Acen}
    ders = {N1: dx(1/S.P), om1: dx(S.W*S.P**2), V1: dx(S.U*S.P), A1: dx(Acen)}
    return vals, ders, (S.P, S.W, S.U), G, kk, Acen


def beta_centre(beta):
    """symbolic centre limits a0, c0 at general k and the closed-form spectrum identity (PS3)."""
    vals, ders, (P, W, U), G, kk, Acen = centre_chart_subs()
    rep = {**vals, **ders}
    a = beta['a'].xreplace(rep).subs(k, kk); c = beta['c'].xreplace(rep).subs(k, kk)
    a = sp.cancel(a); c = sp.cancel(c)
    Ustar = -sp.Rational(2)/(3*(1 + kk))
    # centre equilibrium of the blown-up field at P = 0
    G0 = [sp.cancel(g.subs(P, 0)) for g in G]
    eqU = sp.factor(G0[2])
    a0 = sp.cancel(a.subs(P, 0).subs(U, Ustar)); c0 = sp.cancel(c.subs(P, 0).subs(U, Ustar))
    l_, n_ = sp.symbols('l n', integer=True, nonnegative=True)
    kap_spec = -c0 - a0*(l_ + 2*n_)
    G20 = (2*(1 - 3*kk) - (1 + 3*kk)*(l_ + 2*n_))/(3*(1 + kk))
    # parity: G(-P, W, U) = (-GP, GW, GU)  =>  centre manifold symmetric, a and c even in P
    par = [sp.cancel(G[0].subs(P, -P) + G[0]), sp.cancel(G[1].subs(P, -P) - G[1]), sp.cancel(G[2].subs(P, -P) - G[2])]
    a_even = sp.cancel(a.subs(P, -P) - a) == 0
    c_even = sp.cancel(c.subs(P, -P) - c) == 0
    out = dict(a_centre=a, c_centre=c, a0=a0, c0=c0, kappa_spectrum=sp.factor(kap_spec),
               check={'U* = -2/(3(1+k)) is an equilibrium of G at P=0 (G_U(0,W,U*) = 0)': sp.cancel(G0[2].subs(U, Ustar)) == 0,
                      'G_P(0,W,U) == 0 and G_W(0,W,U*) == 0 hold': sp.cancel(G0[0]) == 0 and sp.cancel(G0[1].subs(U, Ustar)) == 0,
                      'a0 == (1+3k)/(3(1+k)) == 1 + U*': sp.cancel(a0 - (1 + 3*kk)/(3*(1 + kk))) == 0 and sp.cancel(a0 - 1 - Ustar) == 0,
                      'c0 == -2(1-3k)/(3(1+k))': sp.cancel(c0 + 2*(1 - 3*kk)/(3*(1 + kk))) == 0,
                      'PS3 identity: -c0 - a0 (l+2n) - G(20) simplifies to 0': sp.simplify(kap_spec - G20) == 0,
                      'parity of the blown-up centre field (P odd; W, U even)': all(p == 0 for p in par),
                      'a even in P': a_even, 'c even in P': c_even,
                      'value kappa(1/3, 1, 0) == -1/2': kap_spec.subs({kk: sp.Rational(1, 3), l_: 1, n_: 0}) == -sp.Rational(1, 2),
                      'value kappa(1/3, l, n) == -(l+2n)/2': sp.simplify(kap_spec.subs(kk, sp.Rational(1, 3)) + (l_ + 2*n_)/2) == 0})
    out['eqU_factor'] = eqU
    return out


# ======================================================================  axial Pi
def axial_pi():
    B = LC.Lin()
    hs = LC.pert('hs'); hx = LC.pert('hx'); b = LC.pert('b')
    h = sp.zeros(4, 4)
    h[0, 3] = h[3, 0] = R**2*E*hs[0]*sn*Y1
    h[1, 3] = h[3, 1] = R**2*E*hx[0]*sn*Y1
    ein = B.einstein(h, 0, [0, 0, 0, 0], [(0, 3), (1, 3), (2, 3)])
    Ys = -cs/sn*Y1 - L*Y0
    rows = {}
    for key, fac in (((0, 3), sn*Y1), ((1, 3), sn*Y1), ((2, 3), sn*(Ys - cs/sn*Y1))):
        rows[key] = sp.cancel(LC.finalize(ein[key]/fac))
    Q, Q1, Q2 = sp.symbols('Q Q1 Q2')
    # substitute hs' = kappa hx - Q, hs'' = kappa hx' - Q'
    sub = {hs[1]: kap*hx[0] - Q, hs[2]: kap*hx[1] - Q1}
    Es = sp.cancel(rows[(0, 3)].subs(sub)); Ex = sp.cancel(rows[(1, 3)].subs(sub))
    chk = {}
    chk['E_s phi after substitution depends on hx, hx_x?'] = [str(s_) for s_ in (hx[0], hx[1], hx[2]) if sp.cancel(sp.diff(Es, s_)) != 0]
    chk['E_x phi after substitution depends on hs, hs_x?'] = [str(s_) for s_ in (hs[0],) if sp.cancel(sp.diff(Ex, s_)) != 0] + \
        [str(s_) for s_ in (hx[1], hx[2]) if sp.cancel(sp.diff(Ex, s_)) != 0]
    al_s = sp.cancel(sp.diff(Es, Q1)); be_s = sp.cancel(sp.diff(Es, Q)); ga_s = sp.cancel(sp.diff(Es, hs[0]))
    al_x = sp.cancel(sp.diff(Ex, Q1)); be_x = sp.cancel(sp.diff(Ex, Q)); ga_x = sp.cancel(sp.diff(Ex, hx[0]))
    chk['E_s phi == al_s Q_x + be_s Q + ga_s h_s'] = sp.cancel(Es - al_s*Q1 - be_s*Q - ga_s*hs[0]) == 0
    chk['E_x phi == al_x Q_x + be_x Q + ga_x h_x'] = sp.cancel(Ex - al_x*Q1 - be_x*Q - ga_x*hx[0]) == 0
    # master equation: h_s = -(al_s Q' + be_s Q)/ga_s, h_x = -(al_x Q' + be_x Q)/ga_x, Q = kappa h_x - h_s'
    a2 = al_s/ga_s
    b2 = LC.Dx(al_s/ga_s) + be_s/ga_s - kap*al_x/ga_x
    c2 = LC.Dx(be_s/ga_s) - kap*be_x/ga_x - 1
    # exact light-cone exponent (N = 1), using ga_s/ga_x = 1 (on-shell, checked numerically below)
    lam_lc = -N1
    expo = sp.cancel(-(LC.Dx(al_s) + be_s - kap*al_x)/LC.Dx(al_s))
    expo_lc = sp.cancel(expo.subs(N0, 1))
    chk['al_s == (N^2-1)/(2 A N^2)'] = sp.cancel(al_s - (N0**2 - 1)/(2*A0*N0**2)) == 0
    chk['first-order (Q_x level) light-cone exponent at N = 1 equals (1 - lambda_lc - kappa)/lambda_lc, lambda_lc = -N_x'] = \
        sp.cancel(expo_lc - (1 - lam_lc - kap)/lam_lc) == 0
    return dict(rows=rows, al_s=al_s, be_s=be_s, ga_s=ga_s, al_x=al_x, be_x=be_x, ga_x=ga_x, a2=a2, b2=b2, c2=c2,
                expo_Qx=expo, expo_Qx_lc=expo_lc, check=chk, beta_src=None)


def lc_frob2(a, b, at):
    """second-order light-cone Frobenius (shared by AX-Pi and TOY): a f'' + b f' + c f = 0 with a(x*) = 0:
    exponents 0 and sigma = 1 - b/a'; returns the first-order (f') exponent sigma - 1 = -b/a' at the point."""
    ap = LC.Dx(a)
    return sp.cancel((-b/ap).xreplace(at))


# ======================================================================  polar fluid rows
def polar_fluid():
    B = LC.Lin()
    names = ['Hss', 'Hsx', 'Hxx', 'K', 'dl', 'v', 'w']
    P = {n: LC.pert(n) for n in names}
    h = sp.zeros(4, 4)
    h[0, 0] = R**2*E*P['Hss'][0]*Y0; h[0, 1] = h[1, 0] = R**2*E*P['Hsx'][0]*Y0; h[1, 1] = R**2*E*P['Hxx'][0]*Y0
    h[2, 2] = R**2*E*P['K'][0]*Y0; h[3, 3] = R**2*E*P['K'][0]*Y0*sn**2
    du = [0, R*E*P['v'][0]*Y0, R*E*P['w'][0]*Y1, 0]
    du[0] = B.norm_du_s(h, du)
    fl = B.fluid(h, B.rho*E*P['dl'][0]*Y0, du, comps=(0, 1, 2))
    rows = [sp.cancel(LC.finalize(fl[0]/Y0)), sp.cancel(LC.finalize(fl[1]/Y0)), sp.cancel(LC.finalize(fl[2]/Y1))]
    f0 = [P[n][0] for n in ('dl', 'v', 'w')]; f1 = [P[n][1] for n in ('dl', 'v', 'w')]
    C = sp.Matrix([[sp.diff(r_, v) for v in f1] for r_ in rows])
    Bk = sp.Matrix([[sp.diff(sp.diff(r_, v), kap) for v in f0] for r_ in rows])
    Dm = sp.Matrix([[sp.diff(r_, v).subs(kap, 0) for v in f0] for r_ in rows])
    lam = sp.Symbol('lam')
    chi = sp.factor(LC.reduce_roots((C - lam*Bk).det()))
    sq = sp.sqrt(k)
    lp = 1 + N0*(V0 + sq)/(1 + sq*V0); lm = 1 + N0*(V0 - sq)/(1 - sq*V0)
    quad = sp.cancel(chi/(1 + N0*V0 - lam))
    chk = {'w_x enters only the theta row; rows s,x free of w_x': all(C[i, 2] == 0 for i in (0, 1)) and C[2, 0] == 0 and C[2, 1] == 0,
           'L enters only the w column': all(sp.diff(Dm[i, j], L) == 0 and sp.diff(Bk[i, j], L) == 0 for i in range(3) for j in (0, 1)) and sp.diff(C, L) == sp.zeros(3, 3),
           'det(C - lam B) has the transport factor (1 + N V - lam)': sp.cancel(sp.diff(quad, lam, 3)) == 0,
           'quadratic factor vanishes at lam = 1 + N(V +- sqrt k)/(1 +- sqrt k V)': sp.simplify(quad.subs(lam, lp)) == 0 and sp.simplify(quad.subs(lam, lm)) == 0,
           'det C_S vanishes on the sonic set N = (1 - sqrt(k) V)/(sqrt(k) - V)':
               sp.simplify(LC.reduce_roots(C[:2, :2].det()).subs(N0, (1 - sq*V0)/(sq - V0))) == 0}
    # x-flux structure (lemma hypothesis H-flux): C_ff columns == d(dT^x_nu)/d(f) up to the row normalisation
    return dict(C=C, B=Bk, D=Dm, rows=rows, charpoly=chi, check=chk)


def fluxcheck():
    """exact identity: the fluid principal part equals the x-flux d(delta T^x_nu)/d f (nu = s, x, theta), with
    metric perturbation zero; this is the hypothesis used by the sonic (metric-frozen) lemma."""
    B = LC.Lin()
    dl, v, w = sp.symbols('dl v w')
    du = [0, R*E*v*Y0, R*E*w*Y1, 0]
    du[0] = B.norm_du_s(sp.zeros(4, 4), du)
    drho = B.rho*E*dl*Y0
    g, gi, ud = B.g, B.gi, B.udn
    dT = sp.zeros(4, 4)
    for a in range(4):
        for b in range(4):
            dT[a, b] = (1 + k)*(drho*ud[a]*ud[b] + B.rho*(du[a]*ud[b] + ud[a]*du[b])) + k*drho*g[a, b]
    dTm = gi*dT
    # sqrt(-g) factor: dT^x_nu enters as d_x(sqrt(-g) dT^x_nu)/sqrt(-g): principal coefficient dT^x_nu itself
    flux = [LC.finalize(dTm[1, 0]/Y0), LC.finalize(dTm[1, 1]/Y0), LC.finalize(dTm[1, 2]/Y1)]
    return [[sp.cancel(sp.diff(fx, s_)) for s_ in (dl, v, w)] for fx in flux]


# ======================================================================  l = 1 residual gauge
def gauge_l1():
    B = LC.Lin()
    Tf, Xf = sp.symbols('Tf Xf')
    TH = LC.pert('TH')
    xi = [E*Tf*Y0, E*Xf*Y0, E*TH[0]*Y1, 0]
    g, u = B.g, B.udn
    h = sp.zeros(4, 4)
    for a in range(4):
        for b in range(4):
            h[a, b] = sum(xi[c]*D(g[a, b], c) + g[c, b]*D(xi[c], a) + g[a, c]*D(xi[c], b) for c in range(4))
    cond = [sp.cancel(LC.finalize(h[0, 2]/Y1)), sp.cancel(LC.finalize(h[1, 2]/Y1))]
    sol = sp.solve(cond, [Tf, Xf], dict=True)[0]
    tt = LC.finalize(h[2, 2]); pp = LC.finalize(h[3, 3]/sn**2)
    Kc = sp.cancel(((tt + pp)/2).subs(L, 2))
    Kc = sp.factor(sp.cancel(Kc.subs({Y1: 0}).subs(Y0, 1)).subs(sol))
    q = sp.cancel(sp.solve(Kc, TH[1])[0]/TH[0])
    dw = sp.factor(sp.cancel(LC.finalize((u[0]*xi[0] + u[1]*xi[1])/(R*Y0)).subs(sol)).subs(TH[1], q*TH[0]))
    # centre leading coefficient of w_gauge: Theta = exp(-int A) ~ c e^{-x}, 1/N = P ~ e^{x}
    vals, ders, (P, W, U), G, kk, Acen = centre_chart_subs()
    th_cen = sp.Symbol('thc')            # Theta * P  -> finite nonzero constant at the centre
    w_c = sp.cancel((dw/TH[0]).xreplace({N0: 1/P, V0: U*P, A0: Acen, LC.sA: sp.sqrt(Acen), LC.gV: sp.sqrt(1 - U**2*P**2)}))
    w_c = sp.cancel(w_c/P)               # w = Theta * P * (...) and Theta * P -> const
    Ustar = -sp.Rational(2)/(3*(1 + kk))
    lead = sp.factor(sp.limit(w_c.subs(U, Ustar).subs(k, kk), P, 0))
    A_c0 = sp.limit(Acen.subs(U, Ustar), P, 0)
    rate = sp.solve(lead, kap)
    chk = {"residual ODE Theta' = -A Theta (kappa-independent)": sp.cancel(q + A0) == 0,
           'T = (kappa - A) Theta/(A N^2)': sp.cancel(sol[Tf].subs(TH[1], -A0*TH[0]) - (kap - A0)*TH[0]/(A0*N0**2)) == 0,
           'X = (kappa + A(N^2 - 1)) Theta/(A N^2)': sp.cancel(sol[Xf].subs(TH[1], -A0*TH[0]) - (kap + A0*(N0**2 - 1))*TH[0]/(A0*N0**2)) == 0,
           'w_gauge = -(kappa - A(1 + N V)) Theta/(N gV sA)': sp.cancel(dw + (kap - A0*(1 + N0*V0))*TH[0]/(N0*LC.gV*LC.sA)) == 0,
           'A -> 1 at the centre (Theta ~ e^-x, xi^theta ~ 1/r smooth for every kappa)': A_c0 == 1,
           'alpha-gauge rate (unique zero of the centre coefficient)': [str(r_) for r_ in rate],
           'rate == (1+3k)/(3(1+k))': len(rate) == 1 and sp.cancel(rate[0] - (1 + 3*kk)/(3*(1 + kk))) == 0,
           'zero is simple (centre coefficient linear in kappa)': sp.degree(sp.numer(sp.together(lead)), kap) == 1,
           'translation (kappa = 1) excluded: coefficient nonzero': sp.cancel(lead.subs(kap, 1)) != 0,
           'boost (kappa = 0) excluded: coefficient nonzero': sp.cancel(lead.subs(kap, 0)) != 0}
    return dict(T=sol[Tf], X=sol[Xf], q=q, w_gauge=dw, centre_coeff=lead, rate=rate, check=chk)


# ======================================================================  toy
def toy():
    """massless scalar wave on the flat background in the same coordinates: Phi = E f(x) Y (R^0)."""
    B = LC.Lin()
    f = LC.pert('f')
    Phi = E*f[0]*Y0
    g, gi, G = B.g, B.gi, B.G
    box = sum(gi[a, b]*(D(D(Phi, b), a) - sum(G[c][a][b]*D(Phi, c) for c in range(4))) for a in range(4) for b in range(4))
    e = sp.cancel(LC.finalize(box/Y0))
    n0 = sp.Symbol('n0')
    flat = {N0: n0, A0: 1, om0: 0, V0: 0, LC.sA: 1, LC.gV: 1}
    flat.update({LC.BJ['N'][i]: (-1)**i*n0 for i in range(1, 5)})
    flat.update({LC.BJ[F][i]: 0 for F in ('A', 'om', 'V') for i in range(1, 5)})
    e = sp.cancel(e.xreplace(flat))
    a = sp.cancel(sp.diff(e, f[2])); b = sp.cancel(sp.diff(e, f[1])); c = sp.cancel(sp.diff(e, f[0]))
    # first-order exponent at the light cone n0 = 1 (x = 0): -b/a', a' = da/dx = -n0 da/dn0
    rho = sp.cancel((-b/(-n0*sp.diff(a, n0))).subs(n0, 1))
    lam, nn, l_, m_ = sp.symbols('lam n l m')
    # Gundlach normalisation (D4): kappa = lambda + n - l
    rho_G = sp.expand(rho.subs(kap, lam + nn - l_))
    weak = sp.solve(sp.Eq(rho_G.subs(nn, 0), m_ + 1), lam)
    chk = {'toy first-order light-cone exponent == -kappa': sp.cancel(rho + kap) == 0,
           'in G (D4) normalisation rho == l - n - lambda': sp.expand(rho_G - (l_ - nn - lam)) == 0,
           'weak rates (rho = m+1, n = 0) == l - 1 - m  [G (F2)]': len(weak) == 1 and sp.expand(weak[0] - (l_ - 1 - m_)) == 0}
    return dict(a=a, b=b, c=c, rho=rho, check=chk)


# ======================================================================  principal symbols (PL1)
def principal_symbols(ax, pr_rows):
    lam = sp.Symbol('lam')
    out = {}
    # axial master equation for Q: symbol c_ss lam^2 - c_sx lam + c_xx (kappa <-> d/ds)
    num = sp.together(ax['a2']*sp.Symbol('Q2') + ax['b2']*sp.Symbol('Q1') + ax['c2']*sp.Symbol('Q'))
    c_xx = ax['a2']
    c_sx = sp.cancel(sp.diff(ax['b2'], kap))
    c_ss = sp.cancel(sp.diff(ax['c2'], kap, 2)/2)
    sym = sp.factor(sp.cancel(c_ss*lam**2 - c_sx*lam + c_xx))
    out['axial Q symbol'] = str(sym)
    out['axial Q symbol roots 1 +- N'] = sp.cancel(sym.subs(lam, 1 + N0)) == 0 and sp.cancel(sym.subs(lam, 1 - N0)) == 0
    # polar l>=2, fluid eliminated: K obeys a SOUND-speed wave operator in Z1 (GS structure, G l.1165-1166)
    Z1 = pr_rows['Z'][0]; Kj = pr_rows['jets']['K']
    cxx = sp.diff(Z1, Kj[2]); csx = sp.diff(sp.diff(Z1, Kj[1]), kap); css = sp.diff(sp.diff(Z1, Kj[0]), kap, 2)/2
    sK = sp.cancel(css*lam**2 - csx*lam + cxx)
    sq = sp.sqrt(k)
    lp = 1 + N0*(V0 + sq)/(1 + sq*V0); lm = 1 + N0*(V0 - sq)/(1 - sq*V0)
    out['K symbol in Z1 vanishes at the sonic pair'] = sp.simplify(LC.reduce_roots(sK.subs(lam, lp))) == 0 and \
        sp.simplify(LC.reduce_roots(sK.subs(lam, lm))) == 0
    out['K symbol in Z1 is not identically zero'] = sp.cancel(sK) != 0
    return out


def sha(path):
    return hashlib.sha256(open(path, 'rb').read()).hexdigest()


if __name__ == '__main__':
    T0 = time.time()
    log = lambda *a: print('[%6.1fs]' % (time.time() - T0), *a, flush=True)
    import polar_sys as PS
    checks = {}
    be = beta_equation(); checks['beta equation'] = be['check']; log('beta', be['check'])
    bc = beta_centre(be); checks['beta centre / PS3'] = bc['check']; log('beta centre', bc['check'])
    ax = axial_pi(); checks['axial Pi'] = ax['check']; log('axial', ax['check'])
    pf = polar_fluid(); checks['polar fluid / PL1'] = pf['check']; log('polar fluid', pf['check'])
    fx = fluxcheck()
    Cf = pf['C']
    # rows are normalised by the same factor as the flux? check proportionality row by row
    flux_ok = [sp.cancel(LC.reduce_roots(Cf[i, j] - fx[i][j])) == 0 for i in range(3) for j in range(3)]
    checks['polar fluid / PL1']['principal part == x-flux d(dT^x_nu)/df exactly (lemma H-flux)'] = all(flux_ok)
    log('flux', flux_ok)
    g1 = gauge_l1(); checks['l=1 gauge / PL4'] = g1['check']; log('gauge', g1['check'])
    ty = toy(); checks['toy / PL5'] = ty['check']; log('toy', ty['check'])
    s1 = PS.build('l1'); log('polar l=1 system n=%d' % s1['n'])
    s2 = PS.build('lge2'); log('polar l>=2 system n=%d' % s2['n'])
    ps = principal_symbols(ax, dict(Z=s2['Z'], jets=s2['P'])); checks['principal symbols / PL1'] = ps; log('symbols', ps)
    # numeric on-shell checks (background point of the engine ODE; 90 digits): ga_s == ga_x == (L-2)/2
    import bgnum
    num = {}
    for pt in (('1.7', '0.21', '0.13'), ('0.6', '0.05', '0.4'), ('1', '0.0527', '0.644')):
        sub = bgnum.jets_at(*pt)
        for Lv in (6, 12, 20):
            s_ = dict(sub); s_[L] = Lv; s_[kap] = sp.Rational(3, 7)
            gs = sp.N(ax['ga_s'].xreplace(s_), 60); gx = sp.N(ax['ga_x'].xreplace(s_), 60)
            num['%s L=%d' % ('/'.join(pt), Lv)] = [str(sp.Float(gs - (Lv - 2)/sp.Integer(2), 3)), str(sp.Float(gx - (Lv - 2)/sp.Integer(2), 3))]
    checks['axial Pi']['on-shell ga_s - (L-2)/2, ga_x - (L-2)/2 (numeric, 60 digits)'] = num
    log('on-shell', num)
    out = dict(beta_a=sp.srepr(be['a']), beta_c=sp.srepr(be['c']),
               ax={kk_: sp.srepr(ax[kk_]) for kk_ in ('al_s', 'be_s', 'ga_s', 'al_x', 'be_x', 'ga_x', 'a2', 'b2', 'c2', 'expo_Qx')},
               pf={kk_: [sp.srepr(e) for e in pf[kk_]] for kk_ in ('C', 'B', 'D')},
               sys_l1={kk_: [sp.srepr(e) for e in s1[kk_]] for kk_ in ('C', 'Dfull')}, n_l1=s1['n'],
               sys_l2={kk_: [sp.srepr(e) for e in s2[kk_]] for kk_ in ('C', 'Dfull')}, n_l2=s2['n'],
               toy={kk_: sp.srepr(ty[kk_]) for kk_ in ('a', 'b', 'c', 'rho')},
               gauge={kk_: sp.srepr(g1[kk_]) for kk_ in ('T', 'X', 'q', 'w_gauge', 'centre_coeff')},
               beta_centre={kk_: sp.srepr(bc[kk_]) for kk_ in ('a0', 'c0', 'kappa_spectrum')},
               src_sha256={f: sha(os.path.join(HERE, f)) for f in ('lin_core.py', 'polar_sys.py', 'eqs_l.py')})
    pickle.dump(out, open(os.path.join(HERE, 'eqs_l.pkl'), 'wb'))
    def js(o):
        if isinstance(o, dict): return {str(k_): js(v) for k_, v in o.items()}
        if isinstance(o, (list, tuple)): return [js(v) for v in o]
        if isinstance(o, (bool, int, float, str)) or o is None: return o
        return str(o)
    json.dump(js(checks), open(os.path.join(HERE, 'eqs_checks.json'), 'w'), indent=1)
    allb = []
    def walk(o):
        if isinstance(o, dict):
            for v in o.values(): walk(v)
        elif isinstance(o, bool): allb.append(o)
    walk(checks)
    log('boolean checks: %d, all true: %s' % (len(allb), all(allb)))
    log('wrote eqs_l.pkl, eqs_checks.json')
