"""blind-qb2: symbolic derivation (own formulation).

Coordinates (u, x, th, ph), u = -ln(-t), x = r/(-t) with t the polar time normalised
so that alpha(0) = 1 (central proper time).  Metric
    g = e^{-2u} [ -a du^2 + 2 b du dx + c dx^2 + x^2 dOmega^2 ],
background (polar-areal written in (u,x)):  a = al^2 - A x^2, b = -A x, c = A.
Fluid (p = rho/3): u^mu = e^{u} (P, Q, 0, 0), 8 pi rho = e^{2u} w.
Background: P = G/al, Q = x G/al + G v/sqrt(A), G = 1/sqrt(1-v^2)
(v = 3-velocity w.r.t. the polar slicing normal observers, radial).

Perturbation gauge (NOT polar-areal): delta g_ux = 0 and delta g_thth = 0 in (u,x) coords;
unknowns da, dc, dw, dQ  (each times eps e^{lam u}); dP from normalisation.

Quantum source: V_ab at xi = 1/6 with R = 0, from Zahn eq. (19) with (12):
    V_ab = (1/120) Box R_ab - (1/240) g_ab R_cd R^cd + (1/60) R_acbd R^cd
(cross-checked against Zahn eq. (32) coefficients).  Bach B_ab = nabla^c nabla^d C_acbd
+ 1/2 R^cd C_acbd computed separately to test V = B/60.
"""
import sympy as sp, pickle, sys, time

t0 = time.time()
u, x, th, ph = sp.symbols('u x th ph', real=True)
lam, eps = sp.symbols('lam eps')
X = [u, x, th, ph]

def christoffel(g, ginv):
    n = 4
    Gam = [[[0]*n for _ in range(n)] for _ in range(n)]
    dg = [[[sp.diff(g[i, j], X[k]) for k in range(n)] for j in range(n)] for i in range(n)]
    for a in range(n):
        for b in range(n):
            for c in range(b, n):
                s = 0
                for d in range(n):
                    if ginv[a, d] != 0:
                        s += ginv[a, d]*(dg[d][b][c] + dg[d][c][b] - dg[b][c][d])
                s = sp.simplify(s/2) if False else s/2
                Gam[a][b][c] = s
                Gam[a][c][b] = s
    return Gam

def riemann(Gam):
    # R^a_{bcd} = d_c Gam^a_{db} - d_d Gam^a_{cb} + Gam^a_{ce} Gam^e_{db} - Gam^a_{de} Gam^e_{cb}
    n = 4
    R = {}
    for a in range(n):
        for b in range(n):
            for c in range(n):
                for d in range(c+1, n):
                    s = sp.diff(Gam[a][d][b], X[c]) - sp.diff(Gam[a][c][b], X[d])
                    for e in range(n):
                        s += Gam[a][c][e]*Gam[e][d][b] - Gam[a][d][e]*Gam[e][c][b]
                    R[(a, b, c, d)] = s
                    R[(a, b, d, c)] = -s
                R[(a, b, c, c)] = 0
    return R

def ricci(R):
    n = 4
    Ric = sp.zeros(4, 4)
    for b in range(n):
        for d in range(n):
            Ric[b, d] = sum(R[(a, b, a, d)] for a in range(n))
    return Ric

def metric(a, b, c):
    e = sp.exp(-2*u)
    g = sp.Matrix([[-a*e, b*e, 0, 0], [b*e, c*e, 0, 0], [0, 0, x**2*e, 0], [0, 0, 0, x**2*sp.sin(th)**2*e]])
    det2 = -a*c - b**2
    ginv = sp.Matrix([[c/det2/e, -b/det2/e, 0, 0], [-b/det2/e, -a/det2/e, 0, 0],
                      [0, 0, 1/(x**2*e), 0], [0, 0, 0, 1/(x**2*sp.sin(th)**2*e)]])
    return g, ginv

def einstein_and_fluid(a, b, c, w, P, Q):
    g, ginv = metric(a, b, c)
    Gam = christoffel(g, ginv)
    R = riemann(Gam)
    Ric = ricci(R)
    Rs = sum(ginv[i, j]*Ric[i, j] for i in range(4) for j in range(4))
    G = Ric - Rs*g/2
    U = [sp.exp(u)*P, sp.exp(u)*Q, 0, 0]
    Ul = [sum(g[i, j]*U[j] for j in range(4)) for i in range(4)]
    rho8 = sp.exp(2*u)*w          # 8 pi rho
    T8 = sp.Matrix(4, 4, lambda i, j: rho8*(sp.Rational(4, 3)*Ul[i]*Ul[j] + g[i, j]/3))
    # mixed T^a_b
    T8m = ginv*T8
    div = []
    for bb in range(2):
        s = 0
        for aa in range(4):
            s += sp.diff(T8m[aa, bb], X[aa])
            for e in range(4):
                s += Gam[aa][aa][e]*T8m[e, bb] - Gam[e][aa][bb]*T8m[aa, e]
        div.append(s)
    E = G - T8
    norm = sum(g[i, j]*U[i]*U[j] for i in range(4) for j in range(4)) + 1
    return dict(g=g, ginv=ginv, Gam=Gam, R=R, Ric=Ric, Rs=Rs, E=E, div=div, norm=norm, T8=T8, Ul=Ul)

# ---------------- background
al, A, w, v = [sp.Function(n)(x) for n in ('al', 'A', 'w', 'v')]
Gm = 1/sp.sqrt(1 - v**2)
a0 = al**2 - A*x**2; b0 = -A*x; c0 = A
P0 = Gm/al; Q0 = x*Gm/al + Gm*v/sp.sqrt(A)

def tosym(expr, funcs, names, order=4):
    rep = {}
    for f, nm in zip(funcs, names):
        for k in range(order, 0, -1):
            rep[sp.Derivative(f, (x, k))] = sp.Symbol(f'{nm}{k}')
        rep[f] = sp.Symbol(nm)
    return expr.subs(rep)

bgf = [al, A, w, v]; bgn = ['al', 'A', 'w', 'v']
out = {}
print('background...', flush=True)
B = einstein_and_fluid(a0, b0, c0, w, P0, Q0)
print('norm check', sp.simplify(B['norm']), flush=True)
Ebg = {k: sp.simplify(sp.expand(B['E'][i, j])) for k, (i, j) in
       {'uu': (0, 0), 'ux': (0, 1), 'xx': (1, 1), 'tt': (2, 2)}.items()}
divbg = [sp.simplify(B['div'][i]*sp.exp(-2*u)) for i in range(2)]
for k in Ebg:
    out['Ebg_'+k] = tosym(Ebg[k], bgf, bgn)
for i in range(2):
    out['divbg_%d' % i] = tosym(divbg[i], bgf, bgn)
print('bg done', time.time()-t0, flush=True)
for k, e in out.items():
    print(k, 'free:', sorted(str(s) for s in e.free_symbols), 'ops', sp.count_ops(e))

# 2m/r: 1 - 2m/r = g^{ab} d_a r d_b r, r = x e^{-u}
rr = x*sp.exp(-u)
dr = [sp.diff(rr, u), sp.diff(rr, x)]
one_m = sum(B['ginv'][i, j]*dr[i]*dr[j] for i in range(2) for j in range(2))
out['one_minus_2m_r_bg'] = tosym(sp.simplify(one_m), bgf, bgn)
print('1-2m/r bg:', out['one_minus_2m_r_bg'])

# ---------------- quantum V tensor on the background (R_ab = 8 pi T_ab on shell, R = 0)
g, ginv, Gam, Rm = B['g'], B['ginv'], B['Gam'], B['R']
Ric = B['T8']   # on-shell R_ab = 8 pi T_ab (traceless fluid)
# covariant derivative of symmetric 2-tensor: (nabla_c S)_{ab}
def cov2(S):
    D = {}
    for cc in range(4):
        for aa in range(4):
            for bb in range(4):
                s = sp.diff(S[aa, bb], X[cc])
                for e in range(4):
                    s -= Gam[e][cc][aa]*S[e, bb] + Gam[e][cc][bb]*S[aa, e]
                D[(cc, aa, bb)] = s
    return D
print('V tensor...', flush=True)
D1 = cov2(Ric)
# second derivative (nabla_d nabla_c R)_{ab} for the needed (a,b) only, contracted with g^{dc}
def boxS(aa, bb):
    s = 0
    for d in range(4):
        for cc in range(4):
            if ginv[d, cc] == 0:
                continue
            t = sp.diff(D1[(cc, aa, bb)], X[d])
            for e in range(4):
                t -= Gam[e][d][cc]*D1[(e, aa, bb)] + Gam[e][d][aa]*D1[(cc, e, bb)] + Gam[e][d][bb]*D1[(cc, aa, e)]
            s += ginv[d, cc]*t
    return s
Ricup = ginv*Ric*ginv
RR = sum(Ric[i, j]*Ricup[i, j] for i in range(4) for j in range(4))
# R_acbd with lowered first index
def Rlow(a_, b_, c_, d_):
    return sum(g[a_, e]*Rm[(e, b_, c_, d_)] for e in range(4))
Vc = {}
for k, (i, j) in {'uu': (0, 0), 'ux': (0, 1), 'xx': (1, 1), 'tt': (2, 2)}.items():
    s = boxS(i, j)/120 - g[i, j]*RR/240
    s += sum(Rlow(i, cc, j, d)*Ricup[cc, d] for cc in range(4) for d in range(4))/60
    Vc[k] = (s*sp.exp(-2*u)).subs({u: 0, th: sp.pi/3})   # V_ab = e^{2u} Vhat_ab
    out['Vhat_'+k] = tosym(Vc[k], bgf, bgn)
    print('Vhat', k, 'ops', sp.count_ops(out['Vhat_'+k]), time.time()-t0, flush=True)

# ---------------- perturbation (gauge delta b = 0, delta g_thth = 0)
da, dc, dw, dQ = [sp.Function(n)(x) for n in ('da', 'dc', 'dw', 'dQ')]
pf = [da, dc, dw, dQ]; pn = ['da', 'dc', 'dw', 'dQ']
E_ = sp.exp(lam*u)*eps
a1 = a0 + E_*da; c1 = c0 + E_*dc; w1 = w + E_*dw; Q1 = Q0 + E_*dQ
# normalisation fixes dP: -a P^2 + 2 b P Q + c Q^2 = -1 ; linearise
dPs = sp.Symbol('dPs')
nrm = -a1*(P0 + E_*dPs)**2 + 2*b0*(P0 + E_*dPs)*Q1 + c1*Q1**2 + 1
lin = sp.diff(nrm, eps).subs(eps, 0)
dPsol = sp.solve(lin, dPs)[0]
P1 = P0 + E_*dPsol
print('perturbation...', flush=True)
Pt = einstein_and_fluid(a1, b0, c1, w1, P1, Q1)
for k, (i, j) in {'uu': (0, 0), 'ux': (0, 1), 'xx': (1, 1), 'tt': (2, 2)}.items():
    e1 = sp.diff(Pt['E'][i, j], eps).subs(eps, 0)*sp.exp(-lam*u)
    e1 = e1.subs({u: 0, th: sp.pi/3})
    out['dE_'+k] = tosym(tosym(e1, pf, pn), bgf, bgn)
    print('dE', k, 'free', sorted(str(s) for s in out['dE_'+k].free_symbols), time.time()-t0, flush=True)
for i in range(2):
    e1 = sp.diff(Pt['div'][i], eps).subs(eps, 0)*sp.exp(-lam*u)*sp.exp(-2*u)
    e1 = e1.subs({u: 0, th: sp.pi/3})
    out['ddiv_%d' % i] = tosym(tosym(e1, pf, pn), bgf, bgn)
    print('ddiv', i, 'free', sorted(str(s) for s in out['ddiv_%d' % i].free_symbols), time.time()-t0, flush=True)
# perturbation of 1 - 2m/r
om1 = sum(Pt['ginv'][i, j]*dr[i]*dr[j] for i in range(2) for j in range(2))
d_om = (sp.diff(om1, eps).subs(eps, 0)*sp.exp(-lam*u)).subs(u, 0)
out['d_one_minus_2m_r'] = tosym(tosym(d_om, pf, pn), bgf, bgn)
print('d(1-2m/r):', out['d_one_minus_2m_r'])
out['dP'] = tosym(tosym(sp.simplify(dPsol), pf, pn), bgf, bgn)

with open('<LAB>/quantum/seats/blind-qb2/qb2_exprs.pkl', 'wb') as f:
    pickle.dump({k: sp.srepr(e) for k, e in out.items()}, f)
print('DONE', time.time()-t0)
