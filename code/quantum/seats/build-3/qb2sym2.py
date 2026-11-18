"""build-3 QB2: exact sympy preparation of every expression the certified instrument evaluates (non-counting).
Inputs: n1_rows.txt (node N1: the build's own HKA source rows S (4x3, columns sigma_ss, sigma_sx, sigma_xx) and
sigma_I), bach_exprs.txt (b_nn, b_nr, b_rr, b_T of gt).  Output: qb2_exprs.pkl (srepr strings) + printed checks.

Per unit beta (beta_ab = e^{-2x} b_ab):  sigma_ss = A(N^2 b_nn - 2N b_nr + b_rr), sigma_sx = A(N b_nr - b_rr),
sigma_xx = A b_rr  (PROOF-QB2-LEMMAS-v2 sec. 1).  Q^{-1} S = (S_0, S_1, Mx^{-1}(S_2, S_3)).
  generic (tube / sonic):  rows 0,1 = S_0, S_1 ; rows 2,3 = adj(Mx).(S_2,S_3) and detMx separately (sonic shift).
  centre chart: u = 1/N = theta Q, q = N V (= U), V = q u ; the rows per unit beta are written as rational
      functions of (u, q, A) and checked analytic at u = 0 (denominator nonzero at u = 0).
"""
import sys, pickle
sys.path.insert(0, '<LAB>/engine')
import sympy as sp
import system as SY
kv = sp.Rational(1, 3); g = 1 + kv
N, om, V = SY.N, SY.om, SY.V
A = sp.Symbol('A')
L = [eval(l, sp.__dict__.copy()) for l in open('<LAB>/quantum/seats/build-3/n1_rows.txt')]
S = sp.Matrix(4, 3, L[:12]); sI = L[12:15]
assert not any(om in e.free_symbols for e in L), 'rows depend on om'
bnn, bnr, brr = sp.symbols('beta_nn beta_nr beta_rr')
sgm = [A*(N**2*bnn - 2*N*bnr + brr), A*(N*bnr - brr), A*brr]
Sf = [sp.expand(sum(S[i, j]*sgm[j] for j in range(3))) for i in range(4)]
sIf = sp.expand(sum(sI[j]*sgm[j] for j in range(3)))
q = 1 - V**2
Mx = sp.Matrix([[1 + N*V, g*(N + V)/q], [(g - 1)*(N + V), g*(1 + N*V)/q]])
adj = Mx.adjugate(); det = sp.factor(Mx.det())
num23 = [sp.expand(adj[i, 0]*Sf[2] + adj[i, 1]*Sf[3]) for i in range(2)]
out = {}
betas = (bnn, bnr, brr)
# generic per-unit-beta coefficients: R[row][c] for rows 0,1 and numerators of rows 2,3; detMx; sigma_I
gen = []
for e in (Sf[0], Sf[1], num23[0], num23[1]):
    gen.append([sp.cancel(sp.together(sp.diff(e, b))) for b in betas])
gen_sI = [sp.cancel(sp.diff(sIf, b)) for b in betas]
out['gen_rows'] = [[sp.srepr(e) for e in r] for r in gen]
out['gen_det'] = sp.srepr(det)
out['gen_sI'] = [sp.srepr(e) for e in gen_sI]
print('rows 0,1 per beta:', [[sp.factor(e) for e in r] for r in gen[:2]])
print('sigma_I per beta:', [sp.factor(e) for e in gen_sI])
print('detMx =', det)
# centre chart: N = 1/u, V = q u
u, qq = sp.symbols('u q')
cen = {N: 1/u, V: qq*u}
QiS = []
for i in range(2):
    QiS.append([sp.cancel(sp.together(e.subs(cen))) for e in gen[i]])
detc = sp.cancel(sp.together(det.subs(cen)))
for i in range(2):
    QiS.append([sp.cancel(sp.together(e.subs(cen)/detc)) for e in gen[2 + i]])
sIc = [sp.cancel(sp.together(e.subs(cen))) for e in gen_sI]
ok = True
for i in range(4):
    for j in range(3):
        e = QiS[i][j]
        d0 = sp.denom(sp.together(e)).subs(u, 0)
        nm = sp.numer(sp.together(e))
        lo = min(sp.Poly(nm, u).monoms())[0] if nm != 0 else None
        hd = sp.Poly(sp.denom(sp.together(e)), u)
        dlo = min(hd.monoms())[0]
        m = (lo if lo is not None else 0) - dlo
        hat = sp.cancel(sp.together(e*u**(-m))) if nm != 0 else sp.Integer(0)
        assert hat == 0 or sp.denom(sp.together(hat)).subs(u, 0) != 0, (i, j)
        print('centre Q^-1 S row %d per %s: u-order %d ; hat(u=0) = %s' % (i, betas[j], m, sp.factor(hat.subs(u, 0))))
        out.setdefault('cen_rows', {})[(i, j)] = sp.srepr(hat)
        out.setdefault('cen_rows_uorder', {})[(i, j)] = m
print('centre sigma_I per beta:', [sp.factor(e) for e in sIc])
out['cen_sI'] = [sp.srepr(e) for e in sIc]
# centre Acons in (P, W, U)
P, W, U = sp.symbols('P W U')
Ac = sp.cancel(sp.together(SY.Acons.subs(SY.k, kv).subs({N: 1/P, om: W*P**2, V: U*P})))
print('A centre =', Ac, '; A(P=0) =', sp.simplify(Ac.subs(P, 0)))
out['Acen'] = sp.srepr(Ac)
out['Acons'] = sp.srepr(sp.cancel(SY.Acons.subs(SY.k, kv)))
# Bach expressions: exp(l0) -> 1/u (u = 1/N)
Bx = [eval(l, sp.__dict__.copy()) for l in open('<LAB>/quantum/seats/build-3/bach_exprs.txt')]
l0 = sp.Symbol('l0'); uu = sp.Symbol('u')
Bu = []
for e in Bx:
    e2 = sp.expand(sp.simplify(e.subs(sp.exp(l0), 1/uu)))
    e2 = sp.cancel(sp.together(e2.subs(sp.exp(2*l0), 1/uu**2)))
    assert l0 not in e2.free_symbols, e2.free_symbols
    Bu.append(e2)
    print('bach in u: free', sorted(map(str, e2.free_symbols)), ' denom', sp.factor(sp.denom(e2)))
out['bach_u'] = [sp.srepr(e) for e in Bu]
pickle.dump(out, open('<LAB>/quantum/seats/build-3/qb2_exprs.pkl', 'wb'))
print('saved qb2_exprs.pkl')
