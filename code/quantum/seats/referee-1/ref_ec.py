"""referee-1: independent float/mpmath recomputation of the Evans-Coleman Bach profiles (sonic b_uu, centre c0).
Own geometry (ref_geom.Geom) on the generic spherical metric -F(t,r) dt^2 + G(t,r) dr^2 + r^2 dOmega^2.
Background input = certified series coefficients (mid strings) of the E5/BG2 manifolds (data, read-only):
  sonic: outward_C.pkl['sm']  (N, om, V)(theta), x = X(theta)      [BG2 ed61d09268]
  centre: lin_tube2_k1_3.pkl['cm'] (P, W, U)(theta), x = log theta + X(theta) + xc; N = 1/P, om = W P^2, V = U P
KHA metric in (t, r) (build-1 RUN.md identity, re-derived here): z = r/(-t) = e^x, F = A N^2 z^2, G = A.
A from the algebraic constraint A = 1 + 2 om (1 + (1+k) N V + k V^2)/(1 - V^2); the Einstein residual G_ab - 8 pi T_ab
(own code) is printed as the independent check that A, N, om, V are read correctly (incl. the sign of V).
"""
import sys, time, pickle, itertools
import mpmath as mp
import sympy as sp
from ref_geom import Geom

mp.mp.dps = 60
k = mp.mpf(1)/3
t, r, th, ph = sp.symbols('t r theta phi')
Ff = sp.Function('F')(t, r); Gf = sp.Function('G')(t, r)
OUTW = '<LAB>/sideways/build-5/phase2/outward_C.pkl'
TUBE = '<LAB>/engine/cert/lin_tube2_k1_3.pkl'
KX = 6   # x-Taylor order kept (need <= 4 for Bach)


# ---------------- truncated power series in mpmath
def smul(a, b):
    c = [mp.mpf(0)]*(KX + 1)
    for i in range(KX + 1):
        for j in range(KX + 1 - i):
            c[i + j] += a[i]*b[j]
    return c


def scomp(a, h):
    """a(h(s)) where h[0] == 0"""
    out = [mp.mpf(0)]*(KX + 1); p = [mp.mpf(1)] + [mp.mpf(0)]*KX
    for i in range(len(a)):
        if i > KX:
            break
        out = [o + a[i]*q for o, q in zip(out, p)]
        p = smul(p, h)
    return out


def srevert(h):
    """inverse series g with h(g(s)) = s, h[0] = 0, h[1] != 0 (fixed-point iteration)"""
    g = [mp.mpf(0), 1/h[1]] + [mp.mpf(0)]*(KX - 1)
    for _ in range(KX + 2):
        hg = scomp(h, g)
        # g <- g - (h(g) - s)/h1  (linearly convergent, exact after KX+1 sweeps on truncated series)
        g = [gi - (hgi - (1 if i == 1 else 0))/h[1] for i, (gi, hgi) in enumerate(zip(g, hg))]
    return g


def taylor_at(coeffs, th0):
    """re-expand polynomial sum c_n theta^n about th0: list of coefficients in (theta - th0)"""
    out = []
    n = len(coeffs)
    for j in range(KX + 1):
        s = mp.mpf(0)
        for i in range(j, n):
            s += coeffs[i]*mp.binomial(i, j)*th0**(i - j)
        out.append(s)
    return out


def sexp(a):
    out = [mp.exp(a[0])] + [mp.mpf(0)]*KX
    for n in range(1, KX + 1):
        out[n] = sum(j*a[j]*out[n - j] for j in range(1, n + 1))/n
    return out


def slog_1p(h0, h):
    """log(h0 + h(s)) - log(h0), h[0] == 0"""
    u = [hi/h0 for hi in h]
    out = [mp.mpf(0)]*(KX + 1); p = [mp.mpf(1)] + [mp.mpf(0)]*KX
    for m in range(1, KX + 1):
        p = smul(p, u)
        out = [o + (-1)**(m + 1)*q/m for o, q in zip(out, p)]
    return out


def sinv(a):
    out = [1/a[0]] + [mp.mpf(0)]*KX
    for n in range(1, KX + 1):
        out[n] = -sum(a[j]*out[n - j] for j in range(1, n + 1))/a[0]
    return out


# ---------------- generic geometry, built once
_GEO = {}


def geometry():
    if not _GEO:
        t0 = time.time()
        G = Geom(sp.diag(-Ff, Gf, r**2, r**2*sp.sin(th)**2), [t, r, th, ph], simp=lambda e: e)
        G.tensors()
        Gein = (G.Ric - G.g*G.R/2)
        _GEO.update(B=G.B, Ein=Gein, R=G.R, T1=G.T1, T2=G.T2)
        print('generic geometry built in %.1f s' % (time.time() - t0)); sys.stdout.flush()
    return _GEO


def evaluate(xjets, x0):
    """xjets: dict name -> x-Taylor list about x0 for N, om, V. Evaluate at t = -1, r = e^{x0}, theta = pi/2."""
    N, om, V = xjets['N'], xjets['om'], xjets['V']
    one = [mp.mpf(1)] + [mp.mpf(0)]*KX
    num = smul(om, [2*(1 + (1 + k)*a + k*b) if i == 0 else 2*((1 + k)*a + k*b) for i, (a, b) in
                    enumerate(zip(smul(N, V), smul(V, V)))])
    den = [o - vv for o, vv in zip(one, smul(V, V))]
    A = [o + q for o, q in zip(one, smul(num, sinv(den)))]
    ex = sexp([x0] + [mp.mpf(1)] + [mp.mpf(0)]*(KX - 1))      # e^{x0 + s}
    Fx = smul(smul(A, smul(N, N)), smul(ex, ex))
    Gx = A
    # partial derivatives of F(t,r) = f(log r - log(-t)) at t = -1, r = e^{x0}: via sympy on the truncated polynomial
    s = sp.Symbol('s')
    def poly(c):
        return sum(sp.Float(mp.nstr(ci, 55), 55)*s**i for i, ci in enumerate(c))
    xx = sp.log(r) - sp.log(-t) - sp.Float(mp.nstr(x0, 55), 55)
    fpol = poly(Fx).subs(s, xx); gpol = poly(Gx).subs(s, xx)
    r0 = sp.Float(mp.nstr(mp.exp(x0), 55), 55)
    pt = {t: -1, r: r0, th: sp.pi/2}
    geo = geometry()
    def dsub(expr):
        reps = {}
        for d in expr.atoms(sp.Derivative):
            fun = fpol if d.expr == Ff else gpol
            reps[d] = sp.diff(fun, *d.variable_count).subs(pt)
        reps[Ff] = fpol.subs(pt); reps[Gf] = gpol.subs(pt)
        return sp.N(expr.subs(reps).subs(pt), 50)
    Bt = {key: dsub(geo['B'][i, j]) for key, (i, j) in {'tt': (0, 0), 'tr': (0, 1), 'rr': (1, 1), 'hh': (2, 2)}.items()}
    Et = {key: dsub(geo['Ein'][i, j]) for key, (i, j) in {'tt': (0, 0), 'tr': (0, 1), 'rr': (1, 1), 'hh': (2, 2)}.items()}
    Rv = dsub(geo['R'])
    T1tt = dsub(geo['T1'][0, 0])
    F0, G0 = mp.mpf(str(fpol.subs(pt))), mp.mpf(str(gpol.subs(pt)))
    r0m = mp.exp(x0)
    m = lambda e: mp.mpf(str(e))
    b = {'nn': r0m**4*m(Bt['tt'])/F0, 'nr': r0m**4*m(Bt['tr'])/mp.sqrt(F0*G0), 'rr': r0m**4*m(Bt['rr'])/G0,
         'T': r0m**4*m(Bt['hh'])/r0m**2}
    V0 = V[0]
    b['uu'] = (b['nn'] + 2*V0*b['nr'] + V0**2*b['rr'])/(1 - V0**2)
    b['uu_flipV'] = (b['nn'] - 2*V0*b['nr'] + V0**2*b['rr'])/(1 - V0**2)
    # Einstein residual in the frame, scaled by r^2:  r^2 (G_(ab) - 8 pi T_(ab)),  rho = om/(4 pi r^2 A)
    rho = om[0]/(4*mp.pi*r0m**2*A[0]); gam2 = 1/(1 - V0**2)
    T = {'nn': rho*(mp.mpf(4)/3*gam2 - mp.mpf(1)/3), 'nr': -mp.mpf(4)/3*rho*gam2*V0,  # u.n = -gamma
         'rr': rho*(mp.mpf(4)/3*gam2*V0**2 + mp.mpf(1)/3), 'T': rho/3}
    E = {'nn': m(Et['tt'])/F0, 'nr': m(Et['tr'])/mp.sqrt(F0*G0), 'rr': m(Et['rr'])/G0, 'T': m(Et['hh'])/r0m**2}
    res = {a: r0m**2*(E[a] - 8*mp.pi*T[a]) for a in E}
    resflip = r0m**2*(E['nr'] + 8*mp.pi*T['nr'])  # residual if V had the opposite sign
    return b, res, resflip, r0m**2*m(Rv), A[0]


def sonic():
    d = pickle.load(open(OUTW, 'rb'))['sm']
    Y = [[mp.mpf(p[0]) for p in comp] for comp in d['Y']]
    X = [mp.mpf(p[0]) for p in d['X']]
    print('sonic X[0], X[1] =', mp.nstr(X[0], 5), mp.nstr(X[1], 20))
    h = taylor_at(X, mp.mpf(0)); g = srevert(h)
    jets = {nm: scomp(taylor_at(Y[i], mp.mpf(0)), g) for i, nm in enumerate(['N', 'om', 'V'])}
    print('sonic N0, om0, V0 =', [mp.nstr(jets[n][0], 22) for n in ['N', 'om', 'V']])
    N0, V0 = jets['N'][0], jets['V'][0]
    print('sonic condition (1+NV)^2 - k(N+V)^2 =', mp.nstr((1 + N0*V0)**2 - k*(N0 + V0)**2, 5))
    b, res, rf, R, A0 = evaluate(jets, mp.mpf(0))
    print('sonic A0 =', mp.nstr(A0, 22))
    for a in ['nn', 'nr', 'rr', 'T', 'uu', 'uu_flipV']:
        print('sonic b_%s = %s' % (a, mp.nstr(b[a], 25)))
    print('sonic traceless check b_T - (b_nn - b_rr)/2 =', mp.nstr(b['T'] - (b['nn'] - b['rr'])/2, 5))
    print('sonic Einstein residual r^2(G-8piT):', {a: mp.nstr(v, 5) for a, v in res.items()}, ' with V flipped (nr):', mp.nstr(rf, 5))
    print('sonic r^2 R =', mp.nstr(R, 5))


def centre(xs):
    d = pickle.load(open(TUBE, 'rb'))
    cm = d['cm']; xc = mp.mpf(d['xc'][0])
    Y = [[mp.mpf(p[0]) for p in comp] for comp in cm['Y']]
    X = [mp.mpf(p[0]) for p in cm['X']]
    out = []
    for x0 in xs:
        x0 = mp.mpf(x0)
        th0 = mp.findroot(lambda q: mp.log(q) + mp.polyval(X[::-1], q) + xc - x0, mp.exp(x0 - xc))
        # x - x0 as a series in h = theta - th0
        Xh = taylor_at(X, th0); Xh[0] = mp.mpf(0)
        xh = [a + b for a, b in zip(slog_1p(th0, [mp.mpf(0), mp.mpf(1)] + [mp.mpf(0)]*(KX - 1)), Xh)]
        g = srevert(xh)
        P, W, U = (scomp(taylor_at(Y[i], th0)[:1] + [mp.mpf(0)]*KX if False else taylor_at(Y[i], th0), g) for i in range(3))
        Nj = sinv(P); omj = smul(W, smul(P, P)); Vj = smul(U, P)
        b, res, rf, R, A0 = evaluate({'N': Nj, 'om': omj, 'V': Vj}, x0)
        hat = b['nn']*mp.exp(-4*x0)
        out.append((x0, hat, b, res))
        print('centre x0 = %s  b_nn e^{-4x} = %s  b_rr/b_nn = %s  b_T/b_nn = %s  max|res| = %s' % (
            mp.nstr(x0, 6), mp.nstr(hat, 22), mp.nstr(b['rr']/b['nn'], 12), mp.nstr(b['T']/b['nn'], 12),
            mp.nstr(max(abs(v) for v in res.values()), 3)))
        sys.stdout.flush()
    # Richardson in e^{2x}: hat(x) = c0 + c1 e^{2x} + c2 e^{4x}
    if len(out) >= 3:
        (x1, h1, _, _), (x2, h2, _, _), (x3, h3, _, _) = out[-3:]
        M = mp.matrix([[1, mp.exp(2*x1), mp.exp(4*x1)], [1, mp.exp(2*x2), mp.exp(4*x2)], [1, mp.exp(2*x3), mp.exp(4*x3)]])
        c = mp.lu_solve(M, mp.matrix([h1, h2, h3]))
        print('centre Richardson c0 =', mp.nstr(c[0], 22), ' c1 =', mp.nstr(c[1], 10))


if __name__ == '__main__':
    which = sys.argv[1]
    if which == 'sonic':
        sonic()
    else:
        centre([-7, -8, -9, -10, -11, -12])
