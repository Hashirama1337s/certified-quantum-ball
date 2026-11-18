"""blind-qb2: Bach tensor B_ab = nabla^c nabla^d C_acbd + 1/2 R^cd C_acbd on the background metric
(metric-only: al(x), A(x) and derivatives to order 4), and V_ab built from the metric Ricci
(V = (1/120) Box R_ab - (1/240) g R.R + (1/60) R_acbd R^cd), both as Vhat/Bhat = e^{-2u} (.)_ab.
Output: qb2_bach.pkl (srepr strings)."""
import sympy as sp, pickle, time
t0 = time.time()
u, x, th, ph = sp.symbols('u x th ph', real=True)
X = [u, x, th, ph]
al, A = sp.Function('al')(x), sp.Function('A')(x)
e = sp.exp(-2*u)
a0 = al**2 - A*x**2; b0 = -A*x; c0 = A
g = sp.Matrix([[-a0*e, b0*e, 0, 0], [b0*e, c0*e, 0, 0], [0, 0, x**2*e, 0], [0, 0, 0, x**2*sp.sin(th)**2*e]])
ginv = sp.simplify(g.inv())
n = 4
dg = [[[sp.diff(g[i, j], X[k]) for k in range(n)] for j in range(n)] for i in range(n)]
Gam = [[[sp.simplify(sum(ginv[a, d]*(dg[d][b][c] + dg[d][c][b] - dg[b][c][d]) for d in range(n))/2)
         for c in range(n)] for b in range(n)] for a in range(n)]
print('Gam', time.time()-t0, flush=True)
Rup = {}
for a in range(n):
    for b in range(n):
        for c in range(n):
            for d in range(n):
                s = sp.diff(Gam[a][d][b], X[c]) - sp.diff(Gam[a][c][b], X[d])
                for f in range(n):
                    s += Gam[a][c][f]*Gam[f][d][b] - Gam[a][d][f]*Gam[f][c][b]
                Rup[(a, b, c, d)] = s
Rl = {(a, b, c, d): sp.simplify(sum(g[a, f]*Rup[(f, b, c, d)] for f in range(n)))
      for a in range(n) for b in range(n) for c in range(n) for d in range(n)}
print('Riem', time.time()-t0, flush=True)
Ric = sp.Matrix(4, 4, lambda b, d: sp.simplify(sum(Rup[(a, b, a, d)] for a in range(n))))
Rs = sp.simplify(sum(ginv[i, j]*Ric[i, j] for i in range(n) for j in range(n)))
# Weyl (n=4)
def Cw(a, b, c, d):
    return (Rl[(a, b, c, d)] - (g[a, c]*Ric[d, b] - g[a, d]*Ric[c, b] - g[b, c]*Ric[d, a] + g[b, d]*Ric[c, a])/2
            + Rs*(g[a, c]*g[d, b] - g[a, d]*g[c, b])/6)
C = {(a, b, c, d): Cw(a, b, c, d) for a in range(n) for b in range(n) for c in range(n) for d in range(n)}
print('Weyl', time.time()-t0, flush=True)
# first covariant derivative DC[(f,a,b,c,d)] = nabla_f C_abcd
DC = {}
for f in range(n):
    for (a, b, c, d), val in C.items():
        s = sp.diff(val, X[f])
        for h in range(n):
            s -= (Gam[h][f][a]*C[(h, b, c, d)] + Gam[h][f][b]*C[(a, h, c, d)]
                  + Gam[h][f][c]*C[(a, b, h, d)] + Gam[h][f][d]*C[(a, b, c, h)])
        DC[(f, a, b, c, d)] = s
print('DC', time.time()-t0, flush=True)
def nnC(e_, f, a, b, c, d):
    # nabla_e nabla_f C_abcd
    s = sp.diff(DC[(f, a, b, c, d)], X[e_])
    for h in range(n):
        s -= (Gam[h][e_][f]*DC[(h, a, b, c, d)] + Gam[h][e_][a]*DC[(f, h, b, c, d)]
              + Gam[h][e_][b]*DC[(f, a, h, c, d)] + Gam[h][e_][c]*DC[(f, a, b, h, d)]
              + Gam[h][e_][d]*DC[(f, a, b, c, h)])
    return s
Ricup = ginv*Ric*ginv
out = {}
comps = {'uu': (0, 0), 'ux': (0, 1), 'xx': (1, 1), 'tt': (2, 2)}
for k, (a, b) in comps.items():
    s = 0
    for c in range(n):
        for d in range(n):
            for e_ in range(n):
                if ginv[c, e_] == 0:
                    continue
                for f in range(n):
                    if ginv[d, f] == 0:
                        continue
                    s += ginv[c, e_]*ginv[d, f]*nnC(e_, f, a, c, b, d)
            s += Ricup[c, d]*C[(a, c, b, d)]/2
    out['Bhat_'+k] = (s*sp.exp(-2*u)).subs({u: 0, th: sp.pi/3})
    print('B', k, time.time()-t0, flush=True)
# V from metric Ricci
def cov2(Sm):
    D = {}
    for cc in range(n):
        for aa in range(n):
            for bb in range(n):
                s = sp.diff(Sm[aa, bb], X[cc])
                for h in range(n):
                    s -= Gam[h][cc][aa]*Sm[h, bb] + Gam[h][cc][bb]*Sm[aa, h]
                D[(cc, aa, bb)] = s
    return D
D1 = cov2(Ric)
RR = sum(Ric[i, j]*Ricup[i, j] for i in range(n) for j in range(n))
for k, (i, j) in comps.items():
    s = 0
    for d in range(n):
        for cc in range(n):
            if ginv[d, cc] == 0:
                continue
            t = sp.diff(D1[(cc, i, j)], X[d])
            for h in range(n):
                t -= Gam[h][d][cc]*D1[(h, i, j)] + Gam[h][d][i]*D1[(cc, h, j)] + Gam[h][d][j]*D1[(cc, i, h)]
            s += ginv[d, cc]*t
    s = s/120 - g[i, j]*RR/240 + sum(Rl[(i, cc, j, d)]*Ricup[cc, d] for cc in range(n) for d in range(n))/60
    out['Vm_'+k] = (s*sp.exp(-2*u)).subs({u: 0, th: sp.pi/3})
out['Rs'] = (Rs*sp.exp(-2*u)).subs({u: 0})
rep = {}
for f, nm in [(al, 'al'), (A, 'A')]:
    for kk in range(6, 0, -1):
        rep[sp.Derivative(f, (x, kk))] = sp.Symbol(f'{nm}{kk}')
    rep[f] = sp.Symbol(nm)
out = {k: v.subs(rep) for k, v in out.items()}
for k, v in out.items():
    print(k, sorted(str(s) for s in v.free_symbols))
with open('<LAB>/quantum/seats/blind-qb2/qb2_bach.pkl', 'wb') as fh:
    pickle.dump({k: sp.srepr(v) for k, v in out.items()}, fh)
print('DONE', time.time()-t0)
