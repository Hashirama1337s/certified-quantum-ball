"""build-3 QB2 node N1: the four HKA source rows and sigma_I, re-derived by the build (exact sympy).

Method (different from referee-2's r2_forced.py, which solves the KHA rows for (A_x, N_x, A_s), differences the
evolution rows against the HKA pencil and reads off F):  UNDETERMINED COEFFICIENTS.
  * linearised Einstein rows E_ss, E_sx, E_xx of  G - 8 pi T - Sigma  and the two fluid rows div T = 0, in HKA bar
    variables  A = Ab(1 + eps a e^{kap s}), N = Nb(1 + eps n e^{kap s}), om = ob(1 + eps w e^{kap s}),
    V = Vb + eps v e^{kap s},  Sigma_{mu nu} = eps e^{kap s} sigma_{mu nu}(x)  (ss, sx, xx);
  * the HKA pencil from the engine (engine/cert/lin_mc.py hka(): Q, P0, P1, identity row), background on shell
    (engine/system.py rhs_regular, Acons);
  * unknowns: Shat = Q^{-1} S (4x3) and sigma_I (3).  Substitute Y' = Q^{-1}(P0 + kap P1) Y + Shat sigma into each
    of the five rows; require  residual_r = lam_r (I.Y - sigma_I.sigma)  identically in (Y, sigma, kap).
    lam_r is read off from the Y part (and checked to make the whole Y part proportional to I);
    the sigma part gives 15 linear equations for the 15 unknowns.
Then S = Q Shat is compared EXACTLY with srcprop.log rows 0-1 (Ham, slice; HKA sign convention: srcprop l.13 / (-A)
see below) and with r2_forced_k13.log K5 rows 0-3 and sigma_I (K5 l.25), and the frame forms
  S_0 = A b_nn, S_1 = -(1/2) A (b_nn - b_rr), sigma_I = A(-b_nn + N b_nr)   (per unit e^{-2x}, i.e. beta) are checked.
"""
import sys, time
sys.path.insert(0, '<LAB>/engine'); sys.path.insert(0, '<LAB>/engine/cert')
import sympy as sp
import system as SY
t0 = time.time()
kv = sp.Rational(1, 3)
s, x, th, ph = sp.symbols('s x theta phi', real=True)
X = [s, x, th, ph]
eps, kap = sp.symbols('epsilon kappa')
Ab, Nb, ob, Vb = [sp.Function(nm)(x) for nm in ('Ab', 'Nb', 'ob', 'Vb')]
a, n, w, v = [sp.Function(nm)(x) for nm in ('a', 'n', 'w', 'v')]
sg = sp.symbols('sg_ss sg_sx sg_xx')
E = sp.exp(kap*s)
A = Ab*(1 + eps*a*E); N = Nb*(1 + eps*n*E); om = ob*(1 + eps*w*E); V = Vb + eps*v*E
r2 = sp.exp(2*(x - s))
g = sp.zeros(4, 4)
g[0, 0] = r2*A*(1 - N**2); g[0, 1] = g[1, 0] = -r2*A; g[1, 1] = r2*A; g[2, 2] = r2; g[3, 3] = r2*sp.sin(th)**2


def lin(e):
    """first-order coefficient in eps"""
    return sp.diff(e, eps).subs(eps, 0)


def zeroth(e):
    return e.subs(eps, 0)


# metric inverse to first order: g = g0 + eps g1  ->  gi = gi0 - eps gi0 g1 gi0
g0 = g.applyfunc(zeroth); g1 = g.applyfunc(lin)
gi0 = g0.inv().applyfunc(sp.simplify)
gi1 = (-gi0*g1*gi0).applyfunc(sp.simplify)
gi = gi0 + eps*gi1


def christ(g, gi):
    return [[[sum(gi[a_, d]*(sp.diff(g[d, b], X[c]) + sp.diff(g[d, c], X[b]) - sp.diff(g[b, c], X[d])) for d in range(4))/2
              for c in range(4)] for b in range(4)] for a_ in range(4)]


def trunc1(e):
    """keep orders eps^0, eps^1"""
    e = sp.expand(e)
    return zeroth(e) + eps*lin(e)


Gm = christ(g, gi)
Gm = [[[trunc1(Gm[a_][b][c]) for c in range(4)] for b in range(4)] for a_ in range(4)]
print('christoffel %.0fs' % (time.time() - t0), flush=True)
Ric = sp.zeros(4, 4)
for b in range(2):
    for c in range(b, 2):
        e = sum(sp.diff(Gm[a_][b][c], X[a_]) - sp.diff(Gm[a_][b][a_], X[c]) for a_ in range(4))
        e += sum(Gm[a_][a_][d]*Gm[d][b][c] - Gm[a_][c][d]*Gm[d][b][a_] for a_ in range(4) for d in range(4))
        Ric[b, c] = Ric[c, b] = trunc1(e)
# R_thth needed for the Ricci scalar
e = sum(sp.diff(Gm[a_][2][2], X[a_]) - sp.diff(Gm[a_][2][a_], X[2]) for a_ in range(4))
e += sum(Gm[a_][a_][d]*Gm[d][2][2] - Gm[a_][2][d]*Gm[d][2][a_] for a_ in range(4) for d in range(4))
Ric[2, 2] = trunc1(e); Ric[3, 3] = Ric[2, 2]*sp.sin(th)**2
Rs = trunc1(sum(gi[a_, b]*Ric[a_, b] for a_ in range(4) for b in range(4)))
print('ricci %.0fs' % (time.time() - t0), flush=True)
gam = 1/sp.sqrt(1 - V**2); es = sp.exp(s - x)
uup = [es*gam/(N*sp.sqrt(A)), es*gam/sp.sqrt(A)*(1/N + V), 0, 0]
rho = om/(4*sp.pi*A*r2); p = kv*rho
udn = [sum(g[m, l]*uup[l] for l in range(4)) for m in range(4)]
rows = []
for (m, l, si) in ((0, 0, 0), (0, 1, 1), (1, 1, 2)):
    Gml = Ric[m, l] - Rs*g[m, l]/2
    Tml = (rho + p)*udn[m]*udn[l] + p*g[m, l]
    rows.append(lin(Gml - 8*sp.pi*Tml) - E*sg[si])
Tup = [[(rho + p)*uup[m]*uup[l] + p*gi[m, l] for l in range(4)] for m in range(4)]
for nu in range(2):
    e = sum(sp.diff(Tup[m][nu], X[m]) for m in range(4))
    e += sum(Gm[m][m][l]*Tup[l][nu] for m in range(4) for l in range(4))
    e += sum(Gm[nu][m][l]*Tup[m][l] for m in range(4) for l in range(4))
    rows.append(lin(e)*r2**2)
print('rows built %.0fs' % (time.time() - t0), flush=True)
# plain symbols, on-shell background
Y = sp.symbols('a n w v'); Yx = sp.symbols('a_x n_x w_x v_x')
Nn, on_, Vn = SY.N, SY.om, SY.V
Acons = sp.cancel(SY.Acons.subs(SY.k, kv))
f = [sp.cancel(sp.together(e_.subs(SY.k, kv))) for e_ in SY.rhs_regular()]


def ddx(e_):
    return sp.cancel(sp.together(sum(sp.diff(e_, y)*d for y, d in zip((Nn, on_, Vn), f))))


dA = ddx(Acons)
d2 = {}
rep = {}
for fn, val, dval in ((Ab, Acons, dA), (Nb, Nn, f[0]), (ob, on_, f[1]), (Vb, Vn, f[2])):
    rep[sp.Derivative(fn, x)] = dval
for fn, sy, sx in zip((a, n, w, v), Y, Yx):
    rep[sp.Derivative(fn, x)] = sx
val0 = {Ab: Acons, Nb: Nn, ob: on_, Vb: Vn, a: Y[0], n: Y[1], w: Y[2], v: Y[3]}
lrows = []
for i, e_ in enumerate(rows):
    e_ = sp.expand(e_/E)
    e_ = e_.subs(rep).subs(val0)
    assert not (e_.free_symbols & {s, th, ph}), (i, e_.free_symbols & {s, th, ph})
    lrows.append(e_)
print('on shell %.0fs' % (time.time() - t0), flush=True)
# second derivatives must not appear (first-order rows)
for i, e_ in enumerate(lrows):
    assert not e_.has(sp.Derivative), i
# HKA pencil (engine lin_mc.hka, the M instrument's matrices)
import hka_copy as LM     # verbatim copy of engine/cert/lin_mc.py hka() (lin_mc writes lin_fields2.pkl on import)
Ms, Mx, P0, P1, Irow = LM.hka(LM.omx, LM.Vx)
bgx = {LM.omx: f[1], LM.Vx: f[2], SY.k: kv, LM.kap: kap}
Qm = sp.eye(4); Qm[2:, 2:] = Mx
Qm = Qm.subs(bgx); P = (P0 + kap*P1).subs(bgx); Ir = sp.Matrix([Irow]).subs(bgx)
MH = (Qm.inv()*P).applyfunc(lambda e_: sp.cancel(sp.together(e_)))
print('HKA pencil %.0fs' % (time.time() - t0), flush=True)
Sh = sp.Matrix(4, 3, lambda i, j: sp.Symbol('Sh%d%d' % (i, j)))
sI = sp.symbols('sI0:3')
Yv = sp.Matrix(Y); sgv = sp.Matrix(sg)
Ypr = MH*Yv + Sh*sgv
sub = {Yx[i]: Ypr[i] for i in range(4)}
eqs = []
lams = []
allv = list(Y) + list(Yx) + list(sg)
for i, e_ in enumerate(lrows):
    # the row is linear in (Y, Y', sigma): coefficients by differentiation (cheap, no expansion of the row)
    cY = [sp.diff(e_, y) for y in Y]; cYx = [sp.diff(e_, y) for y in Yx]; cs = [sp.diff(e_, y) for y in sg]
    z0 = e_.subs({v_: 0 for v_ in allv})
    assert sp.cancel(sp.together(z0)) == 0, ('constant part', i)
    for c_ in cY + cYx + cs:
        assert not (c_.free_symbols & set(allv)), ('nonlinear', i)
    cY = [sp.cancel(sp.together(c_)) for c_ in cY]; cYx = [sp.cancel(sp.together(c_)) for c_ in cYx]
    cs = [sp.cancel(sp.together(c_)) for c_ in cs]
    RY = [sp.cancel(sp.together(cY[j] + sum(cYx[k]*MH[k, j] for k in range(4)))) for j in range(4)]
    Rs_ = [cs[j] + sum(cYx[k]*Sh[k, j] for k in range(4)) for j in range(3)]
    lam = sp.cancel(sp.together(RY[0]/Ir[0]))
    okY = all(sp.cancel(sp.together(RY[j] - lam*Ir[j])) == 0 for j in range(4))
    print('row %d: Y-part = lam*I : %s ; lam = %s' % (i, okY, sp.factor(lam)), flush=True)
    assert okY
    lams.append(lam)
    for j in range(3):
        eqs.append(sp.together(Rs_[j] + lam*sI[j]))
unk = list(Sh) + list(sI)
sol = sp.solve(eqs, unk, dict=True)
print('solutions:', len(sol), '%.0fs' % (time.time() - t0), flush=True)
assert len(sol) == 1
sol = sol[0]
assert all(u in sol for u in unk), 'underdetermined'
Shs = Sh.subs(sol).applyfunc(lambda e_: sp.cancel(sp.together(e_)))
sIs = [sp.factor(sp.cancel(sol[u])) for u in sI]
S = (Qm*Shs).applyfunc(lambda e_: sp.factor(sp.cancel(sp.together(e_))))
for i in range(4):
    print('N1 S row %d:' % i, [S[i, j] for j in range(3)])
print('N1 sigma_I:', sIs)
# exact comparisons
NN, VV = Nn, Vn
K5 = [[1/NN**2, 2/NN**2, 1/NN**2],
      [-1/(2*NN**2), -1/NN**2, (NN - 1)*(NN + 1)/(2*NN**2)],
      [VV/NN, (NN + 6*VV)/(3*NN), -(2*NN**2*VV - NN - 3*VV)/(3*NN)],
      [1/(3*NN), -(NN*VV - 2)/(3*NN), -(2*NN**2 + NN*VV - 1)/(3*NN)]]
K5sI = [-1/NN**2, (NN**2 - 2)/NN**2, (NN - 1)*(NN + 1)/NN**2]
for i in range(4):
    print('N1 row %d == r2_forced_k13.log K5 row %d :' % (i, i), all(sp.cancel(S[i, j] - K5[i][j]) == 0 for j in range(3)))
print('N1 sigma_I == K5 l.25 :', all(sp.cancel(sIs[j] - K5sI[j]) == 0 for j in range(3)))
# srcprop.log l.13-14 (derive normalisation: Ham row = Acons-weighted; see its own print): row Ham = (8NVom+...)/(3N^2(V^2-1)) (1,2,1)
mA = sp.cancel(-(8*NN*VV*on_ + 2*VV**2*on_ - 3*VV**2 + 6*on_ + 3)/(3*(VV - 1)*(VV + 1)))   # = Acons (xv.log l.2)
print('check: srcprop prefactor == Acons:', sp.cancel(mA - Acons) == 0)
src13 = [(8*NN*VV*on_ + 2*VV**2*on_ - 3*VV**2 + 6*on_ + 3)/(3*NN**2*(VV - 1)*(VV + 1))*c for c in (1, 2, 1)]
src14 = [1/(2*NN), 1/NN, -(NN - 1)*(NN + 1)/(2*NN)]
r13 = [sp.cancel(S[0, j]/src13[j]) for j in range(3)]
r14 = [sp.cancel(S[1, j]/src14[j]) for j in range(3)]
print('N1 row 0 / srcprop l.13 =', r13, '(constant ratio -1/A expected: srcprop row Ham is -A x HKA row 0)')
print('N1 row 1 / srcprop l.14 =', r14)
# frame forms with beta: sigma_ss = A(N^2 bnn - 2N bnr + brr), sigma_sx = A(N bnr - brr), sigma_xx = A brr
bnn, bnr, brr = sp.symbols('beta_nn beta_nr beta_rr')
sgm = [Acons*(NN**2*bnn - 2*NN*bnr + brr), Acons*(NN*bnr - brr), Acons*brr]
Sf = [sp.factor(sp.cancel(sum(S[i, j]*sgm[j] for j in range(3)))) for i in range(4)]
sIf = sp.factor(sp.cancel(sum(sIs[j]*sgm[j] for j in range(3))))
print('N1 S_0 frame == A beta_nn :', sp.cancel(Sf[0] - Acons*bnn) == 0)
print('N1 S_1 frame == -(1/2) A (beta_nn - beta_rr) :', sp.cancel(Sf[1] + Acons*(bnn - brr)/2) == 0)
print('N1 sigma_I frame == A(-beta_nn + N beta_nr) :', sp.cancel(sIf - Acons*(-bnn + NN*bnr)) == 0)
print('N1 S_2 frame:', Sf[2]); print('N1 S_3 frame:', Sf[3])
open('<LAB>/quantum/seats/build-3/n1_rows.txt', 'w').write(
    '\n'.join(sp.srepr(e_) for e_ in [S[i, j] for i in range(4) for j in range(3)] + sIs) + '\n')
print('total %.0fs' % (time.time() - t0))
