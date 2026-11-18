"""Build seat build-1 (B7, SB0/SB1), 2026-11-17.  Own symbolic linearisation of Einstein-Euler (p = k rho)
about a CSS spherical background in KHA self-similar variables (s, x):

    ds^2 = R^2 { A [ (1 - N^2) ds^2 - 2 ds dx + dx^2 ] + dOmega^2 },   R = r = e^{x - s},
    4 pi rho = om / (A R^2),   u^mu = ( G/(R N sA), G (1/N + V)/(R sA), 0, 0 ),  G = 1/sqrt(1 - V^2), sA = sqrt(A)

(the same background as engine/derive.py l.26-68, re-entered here, not imported).  Perturbations carry the
factor E = e^{kappa s} and the CSS weight that makes every equation homogeneous in R.

Representation: everything is a rational function of JET SYMBOLS.  Background fields F in (N, A, om, V) are
F0, F1, F2, F3, F4 (x-derivatives); sA = sqrt(A0), gV = sqrt(1 - V0^2) are extra symbols with their own
derivative rules; perturbation functions p carry jets p0..p3; angular dependence uses Y0 = Y, Y1 = dY/dtheta,
cs = cos theta, sn = sin theta, with Y'' = -cot(theta) Y' - L Y (L = l(l+1)).
d/ds acts on R (-> -R) and E (-> kappa E) only; d/dx on jets and R (-> R); d/dtheta on (Y0, Y1, cs, sn).

The linearised Einstein tensor is computed from the textbook first-order formulas
    dGam^a_bc = 1/2 g^ad (nabla_b h_dc + nabla_c h_db - nabla_d h_bc),
    dR_bc = nabla_a dGam^a_bc - nabla_c dGam^a_ab,
    dG_bc = dR_bc - 1/2 h_bc R - 1/2 g_bc (g^ad dR_ad - h^ad R_ad),
and the fluid equations as the first-order variation of nabla_mu T^mu_nu.  No background equation is used in
the derivation (coefficients keep F1, F2, ...), so every background derivative term is kept (failure mode 1).
"""
import sympy as sp

R, E = sp.symbols('R E', positive=True)
kap, k, L = sp.symbols('kappa k L')
cs, sn, Y0, Y1 = sp.symbols('cs sn Y0 Y1')
sA, gV = sp.symbols('sA gV', positive=True)
FN = ['N', 'A', 'om', 'V']
BJ = {f: sp.symbols('%s0:5' % f) for f in FN}          # background jets F0..F4
N0, A0, om0, V0 = BJ['N'][0], BJ['A'][0], BJ['om'][0], BJ['V'][0]
N1, A1, om1, V1 = BJ['N'][1], BJ['A'][1], BJ['om'][1], BJ['V'][1]
NEXT = {}
for f in FN:
    for i in range(4):
        NEXT[BJ[f][i]] = BJ[f][i + 1]
PJ = {}                                                  # perturbation jets


def pert(name, order=4):
    s = sp.symbols('%s0:%d' % (name, order + 1))
    PJ[name] = s
    for i in range(order):
        NEXT[s[i]] = s[i + 1]
    return s


def Dx(f):
    f = sp.sympify(f)
    out = 0
    for s_ in f.free_symbols:
        if s_ in NEXT:
            out += sp.diff(f, s_)*NEXT[s_]
    if f.has(R):
        out += sp.diff(f, R)*R
    if f.has(sA):
        out += sp.diff(f, sA)*A1/(2*sA)
    if f.has(gV):
        out += sp.diff(f, gV)*(-V0*V1/gV)
    return out


def Ds(f):
    f = sp.sympify(f)
    return sp.diff(f, R)*(-R) + sp.diff(f, E)*kap*E


def Dth(f):
    f = sp.sympify(f)
    return (sp.diff(f, Y0)*Y1 + sp.diff(f, Y1)*(-cs/sn*Y1 - L*Y0) + sp.diff(f, cs)*(-sn) + sp.diff(f, sn)*cs)


def D(f, mu):
    if mu == 0: return Ds(f)
    if mu == 1: return Dx(f)
    if mu == 2: return Dth(f)
    return sp.S(0)                                      # m = 0: no phi dependence


# ---------------- background
def background(Nb=None, Ab=None, omb=None, Vb=None, sAb=None, gVb=None):
    """metric, inverse, u_mu, u^mu, rhohat = 4 pi rho for the CSS background (jets as default)."""
    Nb = N0 if Nb is None else Nb; Ab = A0 if Ab is None else Ab
    omb = om0 if omb is None else omb; Vb = V0 if Vb is None else Vb
    sAb = sA if sAb is None else sAb; gVb = gV if gVb is None else gVb
    g = sp.zeros(4, 4)
    g[0, 0] = R**2*Ab*(1 - Nb**2); g[0, 1] = g[1, 0] = -R**2*Ab; g[1, 1] = R**2*Ab
    g[2, 2] = R**2; g[3, 3] = R**2*sn**2
    gi = sp.zeros(4, 4)
    # inverse of the (s,x) block  R^2 A [[1-N^2, -1], [-1, 1]]: det = -R^4 A^2 N^2
    gi[0, 0] = -1/(R**2*Ab*Nb**2); gi[0, 1] = gi[1, 0] = -1/(R**2*Ab*Nb**2)
    gi[1, 1] = -(1 - Nb**2)/(R**2*Ab*Nb**2)
    gi[2, 2] = 1/R**2; gi[3, 3] = 1/(R**2*sn**2)
    uup = [1/(gVb*R*Nb*sAb), (1/Nb + Vb)/(gVb*R*sAb), 0, 0]
    udn = [sum(g[a, b]*uup[b] for b in range(4)) for a in range(4)]
    rho = omb/(Ab*R**2)
    return g, gi, uup, udn, rho


def christoffel(g, gi):
    G = [[[0]*4 for _ in range(4)] for _ in range(4)]
    for a in range(4):
        for b in range(4):
            for c in range(b, 4):
                v = sum(gi[a, d]*(D(g[d, b], c) + D(g[d, c], b) - D(g[b, c], d)) for d in range(4))/2
                v = sp.cancel(sp.together(v))
                G[a][b][c] = G[a][c][b] = v
    return G


def ricci(G):
    Rc = sp.zeros(4, 4)
    for b in range(4):
        for c in range(b, 4):
            v = (sum(D(G[a][b][c], a) for a in range(4)) - sum(D(G[a][b][a], c) for a in range(4))
                 + sum(G[a][a][d]*G[d][b][c] for a in range(4) for d in range(4))
                 - sum(G[a][c][d]*G[d][b][a] for a in range(4) for d in range(4)))
            Rc[b, c] = Rc[c, b] = sp.cancel(sp.together(v))
    return Rc


def reduce_roots(e):
    """reduce sA^2 -> A0, gV^2 -> 1 - V0^2, cs^2 -> 1 - sn^2 in numerator and denominator; cancel."""
    e = sp.together(sp.sympify(e))
    n, d = sp.fraction(e)

    def red(p):
        p = sp.expand(p)
        if p == 0: return p
        P = sp.Poly(p, sA, gV, cs)
        out = 0
        for (i, j, c_), co in P.terms():
            out += co*A0**(i//2)*sA**(i % 2)*(1 - V0**2)**(j//2)*gV**(j % 2)*(1 - sn**2)**(c_//2)*cs**(c_ % 2)
        return sp.expand(out)
    return sp.cancel(red(n)/red(d))


class Lin:
    """first-order variation of Einstein and fluid equations for a given perturbation."""
    def __init__(self):
        self.g, self.gi, self.uup, self.udn, self.rho = background()
        self.G = christoffel(self.g, self.gi)
        self.Ric = ricci(self.G)
        self.Rs = sp.cancel(sp.together(sum(self.gi[a, b]*self.Ric[a, b] for a in range(4) for b in range(4))))
        g, gi, uu, ud, rho = self.g, self.gi, self.uup, self.udn, self.rho
        self.T = sp.zeros(4, 4)          # T-hat_{ab} = 4 pi T_ab
        for a in range(4):
            for b in range(4):
                self.T[a, b] = (1 + k)*rho*ud[a]*ud[b] + k*rho*g[a, b]
        self.Tm = (gi*self.T).applyfunc(lambda e: sp.cancel(sp.together(e)))   # T^a_b

    def norm_du_s(self, h, du):
        """solve u^mu du_mu = 1/2 h_{mu nu} u^mu u^nu for du_s (du given for mu = 1, 2, 3)."""
        uu = self.uup
        huu = sum(h[a, b]*uu[a]*uu[b] for a in range(4) for b in range(4))/2
        return (huu - sum(uu[m]*du[m] for m in (1, 2, 3)))/uu[0]

    def einstein(self, h, drho, du, comps):
        g, gi, G = self.g, self.gi, self.G
        def nab_h(c, a, b):   # (nabla_c h)_{ab}
            return D(h[a, b], c) - sum(G[e][c][a]*h[e, b] + G[e][c][b]*h[a, e] for e in range(4))
        NH = [[[nab_h(c, a, b) for b in range(4)] for a in range(4)] for c in range(4)]
        dG = [[[0]*4 for _ in range(4)] for _ in range(4)]
        for a in range(4):
            for b in range(4):
                for c in range(b, 4):
                    v = sum(gi[a, d]*(NH[b][d][c] + NH[c][d][b] - NH[d][b][c]) for d in range(4))/2
                    dG[a][b][c] = dG[a][c][b] = v
        def dRic(b, c):
            v = 0
            for a in range(4):
                v += D(dG[a][b][c], a) - D(dG[a][a][b], c)
                for f in range(4):
                    v += (G[a][a][f]*dG[f][b][c] - G[f][a][b]*dG[a][f][c] - G[f][a][c]*dG[a][b][f])
                    v -= (G[a][c][f]*dG[f][a][b] - G[f][c][a]*dG[a][f][b] - G[f][c][b]*dG[a][a][f])
            return v
        need = set(comps)
        allpairs = [(a, b) for a in range(4) for b in range(a, 4)]
        dR = {}
        for (a, b) in allpairs:
            dR[(a, b)] = dRic(a, b)
        dRm = lambda a, b: dR[(min(a, b), max(a, b))]
        hup = gi*h*gi
        dRs = sum(gi[a, b]*dRm(a, b) for a in range(4) for b in range(4)) - \
            sum(hup[a, b]*self.Ric[a, b] for a in range(4) for b in range(4))
        ud, rho = self.udn, self.rho
        out = {}
        for (a, b) in need:
            dGab = dRm(a, b) - h[a, b]*self.Rs/2 - g[a, b]*dRs/2
            dT = (1 + k)*(drho*ud[a]*ud[b] + rho*(du[a]*ud[b] + ud[a]*du[b])) + k*drho*g[a, b] + k*rho*h[a, b]
            out[(a, b)] = dGab - 2*dT
        self._dG = dG
        return out

    def fluid(self, h, drho, du, comps=(0, 1, 2, 3)):
        """first-order variation of nabla_mu T^mu_nu (T-hat), nu in comps."""
        g, gi, G, ud, rho = self.g, self.gi, self.G, self.udn, self.rho
        hup = gi*h*gi
        dT = sp.zeros(4, 4)
        for a in range(4):
            for b in range(4):
                dT[a, b] = (1 + k)*(drho*ud[a]*ud[b] + rho*(du[a]*ud[b] + ud[a]*du[b])) + k*drho*g[a, b] + k*rho*h[a, b]
        dTm = gi*dT - hup*self.T                    # delta(T^m_n)
        def nab_h(c, a, b):
            return D(h[a, b], c) - sum(G[e][c][a]*h[e, b] + G[e][c][b]*h[a, e] for e in range(4))
        dGam = lambda a, b, c: sum(gi[a, d]*(nab_h(b, d, c) + nab_h(c, d, b) - nab_h(d, b, c)) for d in range(4))/2
        out = {}
        for n in comps:
            v = sum(D(dTm[m, n], m) for m in range(4))
            v += sum(G[m][m][a]*dTm[a, n] - G[a][m][n]*dTm[m, a] for m in range(4) for a in range(4))
            v += sum(dGam(m, m, a)*self.Tm[a, n] - dGam(a, m, n)*self.Tm[m, a] for m in range(4) for a in range(4))
            out[n] = v
        return out


def finalize(e, triple=(sp.Rational(3, 5), sp.Rational(4, 5))):
    """drop overall E and R (homogeneous), substitute a Pythagorean (cs, sn), reduce roots."""
    e = sp.sympify(e).subs({E: 1, R: 1})
    e = e.subs({cs: triple[0], sn: triple[1]})
    return reduce_roots(e)


def lin_coeffs(e, Y):
    """e linear in Y0, Y1 -> (coef Y0, coef Y1)"""
    e = sp.expand(sp.numer(sp.together(e)))/sp.denom(sp.together(e))
    return sp.cancel(sp.diff(e, Y0)), sp.cancel(sp.diff(e, Y1))
