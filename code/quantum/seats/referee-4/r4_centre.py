"""referee-4: independent centre derivation of the Bach tensor of the k=1/3 self-similar perfect-fluid solution.

Nothing from engine/ or build-3/ is used. Metric (Schwarzschild-like, central proper time t<0, alpha(t,0)=1):
    ds^2 = -alpha^2 dt^2 + A dr^2 + r^2 dOmega^2,  alpha = alpha(z), A = A(z), z = r/(-t).
Fluid p = rho/3, Eulerian 3-velocity vf(z):  u = Gamma (n + vf e_r),  n = alpha^-1 d_t, e_r = A^-1/2 d_r.
8 pi rho t^2 = mu(z).  Einstein: G_ab = 8 pi T_ab (mostly plus; Wald conventions).
Bach: B_ab = nabla^c nabla^d C_acbd + (1/2) R^cd C_acbd  (DESIGN-QB1 l.~58 definition).
Outputs t^4 B_(nn)(t,0), etc., as functions of W0 := 4 pi rho_c t^2 = mu0/2, and the frame profiles in z.
Mode: python r4_centre.py [NUM]  (NUM: optional rational mu0 to run numerically instead of symbolically)
"""
import sys, time
import sympy as sp

T0 = time.time()
t, r, th, ph = sp.symbols('t r theta phi')
X = [t, r, th, ph]
m0 = sp.Symbol('m0')
ORD = 6  # metric to z^6
A2, A4, A6, c2, c4, c6, m2, m4, v1, v3, v5 = sp.symbols('A2 A4 A6 c2 c4 c6 m2 m4 v1 v3 v5')
z = sp.Symbol('z')

# ---------- step 1: Einstein series with generic functions ----------
al = sp.Function('al')(t, r); AA = sp.Function('AA')(t, r)
g = sp.diag(-al**2, AA, r**2, r**2*sp.sin(th)**2)
gi = sp.diag(-1/al**2, 1/AA, 1/r**2, 1/(r**2*sp.sin(th)**2))


def christ(g, gi):
    G = [[[0]*4 for _ in range(4)] for _ in range(4)]
    for a in range(4):
        for b in range(4):
            for c in range(4):
                G[a][b][c] = sum(gi[a, d]*(sp.diff(g[d, b], X[c]) + sp.diff(g[d, c], X[b]) - sp.diff(g[b, c], X[d]))
                                 for d in range(4))/2
    return G


def riemann(Gm):
    # R^a_{bcd} = d_c Gam^a_{db} - d_d Gam^a_{cb} + Gam^a_{ce} Gam^e_{db} - Gam^a_{de} Gam^e_{cb}   (Wald)
    R = {}
    for a in range(4):
        for b in range(4):
            for c in range(4):
                for d in range(4):
                    if c >= d:
                        continue
                    e_ = (sp.diff(Gm[a][d][b], X[c]) - sp.diff(Gm[a][c][b], X[d])
                          + sum(Gm[a][c][e]*Gm[e][d][b] - Gm[a][d][e]*Gm[e][c][b] for e in range(4)))
                    R[(a, b, c, d)] = e_
                    R[(a, b, d, c)] = -e_
                R[(a, b, c, c)] = 0 if (a, b, c, c) not in R else R[(a, b, c, c)]
    for a in range(4):
        for b in range(4):
            for c in range(4):
                R[(a, b, c, c)] = 0
    return R


Gm = christ(g, gi)
Rm = riemann(Gm)
Ric = sp.Matrix(4, 4, lambda b, d: sum(Rm[(a, b, a, d)] for a in range(4)))
Rs = sum(gi[a, b]*Ric[a, b] for a in range(4) for b in range(4))
Ein = Ric - Rs*g/2

zz = r/(-t)
alS = 1 + c2*zz**2 + c4*zz**4 + c6*zz**6
AS = 1 + A2*zz**2 + A4*zz**4 + A6*zz**6
muS = m0 + m2*zz**2 + m4*zz**4
vS = v1*zz + v3*zz**3 + v5*zz**5
Gam = 1/sp.sqrt(1 - vS**2)
uup = [Gam/alS, Gam*vS/sp.sqrt(AS), 0, 0]
gS = g.subs({al: alS, AA: AS}).doit()
udn = [sum(gS[a, b]*uup[b] for b in range(4)) for a in range(4)]
rho8pi = muS/t**2
Tm = sp.Matrix(4, 4, lambda a, b: rho8pi*(sp.Rational(4, 3)*udn[a]*udn[b] + gS[a, b]/3))


def subs_fun(e):
    e = e.subs({al: alS, AA: AS}).doit()
    return e


eqs = []
for (a, b) in [(0, 0), (0, 1), (1, 1), (2, 2)]:
    e_ = subs_fun(Ein[a, b]) - Tm[a, b]
    e_ = e_.subs({t: -1, th: sp.pi/2})
    ser = sp.series(e_, r, 0, 6).removeO()
    for n in range(0, 6):
        cf = sp.expand(ser.coeff(r, n))
        if cf != 0:
            eqs.append(((a, b), n, cf))
print('einstein series built (%.0fs)' % (time.time() - T0), flush=True)
unk = [A2, c2, v1, A4, c4, m2, v3, A6, c6, m4, v5]
sol = {}
# solve progressively by order
for n in range(0, 6):
    E = [sp.expand(cf.subs(sol)) for (_, nn, cf) in eqs if nn == n]
    E = [e for e in E if e != 0]
    if not E:
        continue
    fs = sorted(set().union(*[e.free_symbols for e in E]) - {m0}, key=str)
    s_ = sp.solve(E, fs, dict=True)
    print('order r^%d: unknowns %s -> %s' % (n, fs, s_), flush=True)
    assert len(s_) == 1, s_
    sol.update({k: sp.factor(v) for k, v in s_[0].items()})
    sol = {k: sp.factor(v.subs(sol)) for k, v in sol.items()}
# residual check
bad = [((a, b), n) for ((a, b), n, cf) in eqs if sp.simplify(cf.subs(sol)) != 0]
print('all Einstein series equations to r^5 hold:', not bad, bad)
print('SOL', sol, flush=True)
missing = [u for u in unk if u not in sol]
print('unsolved (beyond reach of r^5 eqs):', missing, flush=True)
