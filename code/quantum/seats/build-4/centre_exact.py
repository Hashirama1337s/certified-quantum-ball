"""build-4 QB2b item 2: EXACT centre orders of b_ab (sympy), and the certified order-4 coefficient.

Regular centre (engine/system.py rhs_centre_blown, k = 1/3):  Z = (P, W, U),  dZ/dx = G(Z),  equilibrium
P = 0, U = -1/2, W = W0 free.  Regular solution as a power series in v = e^{x - xc}:
    P = v + sum_{n>=2} p_n v^n,  W = W0 + sum w_n v^n,  U = -1/2 + sum u_n v^n,     v dZ/dv = G(Z)
(order n: (n - J0) Z_n = [lower]_n with J0 = DG(equilibrium); the order-1 data and the invertibility of (n - J0),
n = 2..NMAX, are checked exactly).  On it, b_ab = Bach_(ab)[gt] (build-3 bach_u.txt) with u = P, A = Acons,
A_j = D^j A, l_j = D^{j-1} l', l' = (ln N)' = -D ln P, D = d/dx = v d/dv.
Claims (exact, every coefficient a rational function of W0, checked == 0 by cancel):
    b_nn, b_rr, b_T = O(v^4),  b_nr = O(v^5);  b_rr_4 = b_T_4 = b_nn_4/3 (isotropy).
Any other parametrisation theta = v(1 + O(v)) (e.g. the cm record's) has the same orders.
Then c0 = lim b_nn e^{-4x} = b_nn_4(W0) e^{-4 xc}, certified on the E5 boxes of W0, xc (arb) -> overlap with
-173.10044439229957 (build-1 result_qb1.json centre c0).
"""
import sys, time, json
sys.path.insert(0, '<LAB>/engine')
import sympy as sp
import system as SY
T0 = time.time()
kv = sp.Rational(1, 3)
P, W, U = sp.symbols('P W U'); W0, v = sp.symbols('W0 v')
G = [sp.cancel(sp.together(e.subs(SY.k, kv))) for e in SY.rhs_centre_blown()]
NMAX = 6
eq = {P: 0, W: W0, U: sp.Rational(-1, 2)}
G0 = [sp.simplify(g.subs(eq)) for g in G]
print('G at equilibrium:', G0)
J0 = sp.Matrix(3, 3, lambda i, j: sp.simplify(sp.diff(G[i], [P, W, U][j]).subs(eq)))
print('J0 =', J0, ' eigenvalues', J0.eigenvals())
# series ansatz
pc = [0, 1] + list(sp.symbols('p2:%d' % (NMAX + 1)))
wc = [W0] + list(sp.symbols('w1:%d' % (NMAX + 1)))
uc = [sp.Rational(-1, 2)] + list(sp.symbols('u1:%d' % (NMAX + 1)))
sol = {}


def ser(c, n):
    return sum(c[i]*v**i for i in range(n + 1))


for n in range(1, NMAX + 1):
    Zs = [ser(pc, n).subs(sol), ser(wc, n).subs(sol), ser(uc, n).subs(sol)]
    eqs = []
    for i in range(3):
        lhs = sp.expand(v*sp.diff(Zs[i], v))
        rhs = sp.series(G[i].subs({P: Zs[0], W: Zs[1], U: Zs[2]}), v, 0, n + 1).removeO()
        eqs.append(sp.cancel(sp.expand(lhs - rhs).coeff(v, n)))
    unk = ([pc[n]] if n >= 2 else []) + [wc[n], uc[n]]
    s_ = sp.solve(eqs, unk, dict=True)
    if n == 1:
        # P_1 = 1 is the normalisation (defines xc); the order-1 equations must hold with it
        print('order 1 solutions:', s_)
    assert len(s_) == 1, (n, s_)
    for kk, vv in s_[0].items():
        sol[kk] = sp.cancel(vv)
    M = sp.Matrix([[sp.diff(e, u_) for u_ in unk] for e in eqs])
    if n == 1:
        res = [sp.simplify(e.subs(s_[0])) for e in eqs]
        print('order 1 (P_1 = 1 normalisation): all three equations hold:', all(r == 0 for r in res), flush=True)
        assert all(r == 0 for r in res)
    else:
        print('order %d: unique, det of the linear system = %s  (%.0fs)' % (n, sp.factor(M.det()), time.time() - T0), flush=True)
Pz = sp.expand(ser(pc, NMAX).subs(sol)); Wz = sp.expand(ser(wc, NMAX).subs(sol)); Uz = sp.expand(ser(uc, NMAX).subs(sol))
print('P =', sp.series(Pz, v, 0, 4)); print('U =', sp.series(Uz, v, 0, 3)); print('W =', sp.series(Wz, v, 0, 3))
N_ = NMAX


def tr(e):
    return sp.series(e, v, 0, N_ + 1).removeO()


def Dx(e):
    return sp.expand(v*sp.diff(e, v))


Acen = sp.cancel(sp.together(SY.Acons.subs(SY.k, kv).subs({SY.N: 1/P, SY.om: W*P**2, SY.V: U*P})))
A = tr(Acen.subs({P: Pz, W: Wz, U: Uz}))
Aj = [A]
for j in range(4):
    Aj.append(tr(Dx(Aj[-1])))
l1 = tr(-Dx(Pz)/Pz)
lj = [l1]
for j in range(3):
    lj.append(tr(Dx(lj[-1])))
print('inputs ready (%.0fs)' % (time.time() - T0), flush=True)
Bx = [eval(l, sp.__dict__.copy()) for l in open('<LAB>/quantum/seats/build-3/bach_u.txt')]
syms = [sp.Symbol('u')] + list(sp.symbols('A0:5')) + list(sp.symbols('l1:5'))
vals = [Pz] + Aj + lj
out = {}
coeffs = {}
for nm, e in zip(('nn', 'nr', 'rr', 'T'), Bx):
    num, den = sp.fraction(sp.together(e))
    numv = sp.expand(tr(num.subs(dict(zip(syms, vals)))))
    denv = sp.expand(tr(den.subs(dict(zip(syms, vals)))))
    q = tr(numv/denv) if False else sp.series(numv/denv, v, 0, N_).removeO()
    cs = [sp.cancel(sp.expand(q).coeff(v, i)) for i in range(N_)]
    coeffs[nm] = cs
    print('b_%s coefficients v^0..v^%d:' % (nm, N_ - 1), [('0' if c == 0 else 'nonzero') for c in cs], '(%.0fs)' % (time.time() - T0), flush=True)
claims = dict(
    nn_O4=all(coeffs['nn'][i] == 0 for i in range(4)), rr_O4=all(coeffs['rr'][i] == 0 for i in range(4)),
    T_O4=all(coeffs['T'][i] == 0 for i in range(4)), nr_O5=all(coeffs['nr'][i] == 0 for i in range(5)),
    isotropy_rr4=sp.cancel(coeffs['rr'][4] - coeffs['nn'][4]/3) == 0, isotropy_T4=sp.cancel(coeffs['T'][4] - coeffs['nn'][4]/3) == 0)
print('CLAIMS', claims)
bnn4 = sp.factor(coeffs['nn'][4]); bnr5 = sp.factor(coeffs['nr'][5])
print('b_nn_4(W0) =', bnn4); print('b_nr_5(W0) =', bnr5)
# certified evaluation on the E5 boxes (lin_tube_v2.ENC, the boxes of the cm record)
from flint import arb, ctx
ctx.prec = 200
W0b = arb('8.90132327775564352') + arb(0, arb('7.92e-18'))
xcb = arb('0.212365646758945385') + arb(0, arb('7.95e-19'))
f = sp.lambdify(W0, bnn4, modules=[{'sqrt': lambda z: z.sqrt()}])
num_, den_ = sp.fraction(sp.together(bnn4))
pn = sp.Poly(sp.expand(num_), W0); pd = sp.Poly(sp.expand(den_), W0)
ev = lambda p: sum((arb(int(sp.Rational(c).p))/arb(int(sp.Rational(c).q)))*W0b**m[0] for m, c in zip(p.monoms(), p.coeffs()))
c0 = ev(pn)/ev(pd)*(-4*xcb).exp()
ref = arb('-173.1004443922995663248427') + arb(0, arb('2.99e-13'))
print('c0 =', c0.str(25), ' overlaps QB1 c0:', c0.overlaps(ref))
json.dump(dict(claims={k_: bool(v_) for k_, v_ in claims.items()}, b_nn_4=str(bnn4), b_nr_5=str(bnr5),
               c0=[c0.mid().str(25, radius=False), c0.rad().str(3, radius=False)], c0_overlaps_QB1=bool(c0.overlaps(ref)),
               NMAX=NMAX, J0=str(J0), time_s=round(time.time() - T0, 1)),
          open('<LAB>/quantum/seats/build-4/centre_exact.json', 'w'), indent=1)
