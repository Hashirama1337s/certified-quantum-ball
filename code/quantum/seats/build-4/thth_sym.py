"""build-4 C8 theta-theta residual, exact sympy part.
(1) R_thth := coefficient of eps e^{kap s} in E_thth = G_thth - 8 pi T_thth - Sigma_thth, for the HKA bar perturbation
    (same ansatz and metric as n1_rows.py / referee r2_forced.py), as a polynomial in
    Y = (a, n, w, v), Y' (a_x ...), Y'' (a_xx ...), sigma_TT, with background symbols Ab, Nb, ob, Vb and the
    x-derivatives Abx, Abxx, Nbx, Nbxx (kept as symbols: evaluated from the background ODE numerically).
(2) D_x := e^{2(x-s)} e^{-kap s} (div Sigma)_x  as a linear form in sigma_ss, sigma_sx, sigma_xx, sigma_TT and their
    x-derivatives; sigma_TT(D_x = 0) solved (coefficient checked == -2, dt_struct.log).
Writes thth_exprs.txt (srepr).  kappa = 2.
"""
import sys, time
sys.path.insert(0, '<LAB>/engine')
import sympy as sp
t0 = time.time()
kv = sp.Rational(1, 3)
s, x, th, ph = sp.symbols('s x theta phi', real=True)
X = [s, x, th, ph]
eps = sp.Symbol('epsilon'); kap = sp.Integer(2)
Ab, Nb, ob, Vb = [sp.Function(nm)(x) for nm in ('Ab', 'Nb', 'ob', 'Vb')]
a, n, w, v = [sp.Function(nm)(x) for nm in ('a', 'n', 'w', 'v')]
E = sp.exp(kap*s)
A = Ab*(1 + eps*a*E); N = Nb*(1 + eps*n*E); om = ob*(1 + eps*w*E); V = Vb + eps*v*E
r2 = sp.exp(2*(x - s))
g = sp.zeros(4, 4)
g[0, 0] = r2*A*(1 - N**2); g[0, 1] = g[1, 0] = -r2*A; g[1, 1] = r2*A; g[2, 2] = r2; g[3, 3] = r2*sp.sin(th)**2
lin = lambda e: sp.diff(e, eps).subs(eps, 0)
zer = lambda e: e.subs(eps, 0)
g0 = g.applyfunc(zer); g1 = g.applyfunc(lin)
gi0 = g0.inv().applyfunc(sp.simplify); gi1 = (-gi0*g1*gi0).applyfunc(sp.simplify)
gi = gi0 + eps*gi1
tr1 = lambda e: zer(sp.expand(e)) + eps*lin(sp.expand(e))
Gm = [[[tr1(sum(gi[a_, d]*(sp.diff(g[d, b], X[c]) + sp.diff(g[d, c], X[b]) - sp.diff(g[b, c], X[d])) for d in range(4))/2)
         for c in range(4)] for b in range(4)] for a_ in range(4)]
Ric = sp.zeros(4, 4)
for (b, c) in ((0, 0), (0, 1), (1, 1), (2, 2)):
    e = sum(sp.diff(Gm[a_][b][c], X[a_]) - sp.diff(Gm[a_][b][a_], X[c]) for a_ in range(4))
    e += sum(Gm[a_][a_][d]*Gm[d][b][c] - Gm[a_][c][d]*Gm[d][b][a_] for a_ in range(4) for d in range(4))
    Ric[b, c] = tr1(e)
Ric[1, 0] = Ric[0, 1]; Ric[3, 3] = Ric[2, 2]*sp.sin(th)**2
Rs = tr1(sum(gi[a_, b]*Ric[a_, b] for a_ in range(4) for b in range(4)))
rho = om/(4*sp.pi*A*r2); p = kv*rho
sTT = sp.Symbol('sg_TT')
Eth = lin(Ric[2, 2] - Rs*g[2, 2]/2 - 8*sp.pi*p*g[2, 2]) - E*sTT
Eth = sp.expand(Eth/E)
# plain symbols
Ys = sp.symbols('a n w v'); Y1 = sp.symbols('a_x n_x w_x v_x'); Y2 = sp.symbols('a_xx n_xx w_xx v_xx')
bs = sp.symbols('Ab Nb ob Vb'); b1 = sp.symbols('Abx Nbx obx Vbx'); b2 = sp.symbols('Abxx Nbxx obxx Vbxx')
rep = {}
for fn, s0, s1_, s2_ in zip((Ab, Nb, ob, Vb), bs, b1, b2):
    rep[sp.Derivative(fn, (x, 2))] = s2_; rep[sp.Derivative(fn, x)] = s1_
for fn, s0, s1_, s2_ in zip((a, n, w, v), Ys, Y1, Y2):
    rep[sp.Derivative(fn, (x, 2))] = s2_; rep[sp.Derivative(fn, x)] = s1_
val = dict(zip((Ab, Nb, ob, Vb), bs)); val.update(dict(zip((a, n, w, v), Ys)))
Rth = sp.simplify(Eth.subs(rep).subs(val))
assert not (Rth.free_symbols & {s, x, th, ph}), Rth.free_symbols
print('R_thth free symbols:', sorted(map(str, Rth.free_symbols)), ' (%.0fs)' % (time.time() - t0), flush=True)
# D_x with sigma functions
sf = [sp.Function(nm)(x) for nm in ('fss', 'fsx', 'fxx', 'fTT')]
Sig = sp.zeros(4, 4)
Sig[0, 0] = E*sf[0]; Sig[0, 1] = Sig[1, 0] = E*sf[1]; Sig[1, 1] = E*sf[2]; Sig[2, 2] = E*sf[3]; Sig[3, 3] = E*sf[3]*sp.sin(th)**2
gb = g0; gbi = gi0
Gb = [[[sp.simplify(sum(gbi[a_, d]*(sp.diff(gb[d, b], X[c]) + sp.diff(gb[d, c], X[b]) - sp.diff(gb[b, c], X[d])) for d in range(4))/2)
        for c in range(4)] for b in range(4)] for a_ in range(4)]
Smix = gbi*Sig
nu = 1
e = sum(sp.diff(Smix[m, nu], X[m]) for m in range(4))
e += sum(Gb[m][m][l]*Smix[l, nu] for m in range(4) for l in range(4))
e -= sum(Gb[l][m][nu]*Smix[m, l] for m in range(4) for l in range(4))
Dx = sp.simplify(e*sp.exp(2*(x - s))/E)
ss = sp.symbols('sg_ss sg_sx sg_xx sg_TT'); ssx = sp.symbols('sgx_ss sgx_sx sgx_xx sgx_TT')
rep2 = {sp.Derivative(sf[i], x): ssx[i] for i in range(4)}
rep2.update({sp.Derivative(Ab, x): b1[0], sp.Derivative(Nb, x): b1[1]})
Dx = sp.expand(Dx.subs(rep2).subs({sf[i]: ss[i] for i in range(4)}).subs({Ab: bs[0], Nb: bs[1]}))
cTT = sp.simplify(Dx.coeff(ss[3])); cTTx = sp.simplify(Dx.coeff(ssx[3]))
print('D_x: coefficient of sigma_TT =', cTT, ' of sigma_TT\' =', cTTx)
sTT_sol = sp.simplify(sp.solve(Dx, ss[3])[0])
print('sigma_TT(D_x = 0) =', sTT_sol)
open('<LAB>/quantum/seats/build-4/thth_exprs.txt', 'w').write(sp.srepr(Rth) + '\n' + sp.srepr(sTT_sol) + '\n' + sp.srepr(Dx) + '\n')
print('total %.0fs' % (time.time() - t0))
