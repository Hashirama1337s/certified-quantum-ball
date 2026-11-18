"""referee-4: Bach tensor at the regular centre, from the definition, with the centre series of r4_centre.log.

B_ab = nabla^c nabla^d C_acbd + (1/2) R^cd C_acbd  (DESIGN-QB1 sec. 1.2 definition), Wald conventions.
Metric: ds^2 = -alpha^2 dt^2 + A dr^2 + r^2 dOmega^2, alpha, A even polynomials in z = r/(-t) to z^6
(coefficients from r4_centre.log; A6, c6 keep the free symbol v5 as a control: it must drop out of
B_(nn), B_(rr), B_(TT) at r^0 and of B_(nr) at r^1).
All intermediates are Laurent polynomials in r, truncated at r^KMAX (error analysis in the report).
Output: frame components at t = -1 (central proper time, so v = z), W0 = m0/2.
Usage: r4_bach.py            -> symbolic m0
"""
import sys, time, json
import sympy as sp

T0 = time.time()
t, r, th, ph = sp.symbols('t r theta phi')
X = [t, r, th, ph]
m0, v5 = sp.symbols('m0 v5')
KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 12
MODE = sys.argv[2] if len(sys.argv) > 2 else "on"
sol = {'A2': m0/3, 'c2': m0/6, 'A4': -m0*(m0 - 12)/45, 'c4': -m0*(m0 - 42)/360,
       'A6': -m0*(-8*m0**2 - 18*m0 + 2160*v5 - 297)/2160, 'c6': -m0*(-10*m0**2 + 18*m0 + 1200*v5 - 321)/4320}
if MODE == "off":  # POST HOC off-shell null: break Einstein at orders z^2 and z^4
    sol["A2"] = sol["A2"] + sp.Rational(1, 10); sol["c4"] = sol["c4"] + sp.Rational(1, 7)
print("KMAX", KMAX, "MODE", MODE, flush=True)
zz = r/(-t)
alS = 1 + sol['c2']*zz**2 + sol['c4']*zz**4 + sol['c6']*zz**6
AS = 1 + sol['A2']*zz**2 + sol['A4']*zz**4 + sol['A6']*zz**6


def trunc(e, K=KMAX):
    e = sp.expand(e)
    out = 0
    for term in sp.Add.make_args(e):
        # degree() fails for negative powers; use as_coeff_exponent
        c, p = term.as_coeff_exponent(r)
        if p <= K:
            out += term
    return out


def inv_series(f, K=KMAX):
    # 1/f for f = 1 + O(r^2), polynomial in r
    d = sp.expand(f - 1)
    s, pw = 1, 1
    for _ in range(K//2 + 1):
        pw = trunc(-pw*d, K)
        s += pw
    return trunc(s, K)


def sqrt_series(f, K=KMAX):
    d = sp.expand(f - 1)
    s = 0
    for n in range(K//2 + 1):
        s += sp.binomial(sp.Rational(1, 2), n)*trunc(d**n, K)
    return trunc(s, K)


g = sp.diag(trunc(-alS**2), AS, r**2, r**2*sp.sin(th)**2)
gi = sp.diag(trunc(-inv_series(alS**2)), inv_series(AS), 1/r**2, 1/(r**2*sp.sin(th)**2))
Gm = [[[trunc(sum(gi[a, a]*(sp.diff(g[a, b], X[c]) + sp.diff(g[a, c], X[b]) - sp.diff(g[b, c], X[a]))/2 for _ in [0]))
        for c in range(4)] for b in range(4)] for a in range(4)]
print('christoffel (%.0fs)' % (time.time() - T0), flush=True)
# Riemann R^a_{bcd} (Wald): d_c Gam^a_{db} - d_d Gam^a_{cb} + Gam^a_{ce}Gam^e_{db} - Gam^a_{de}Gam^e_{cb}
Rud = {}
for a in range(4):
    for b in range(4):
        for c in range(4):
            for d in range(4):
                if c == d:
                    Rud[(a, b, c, d)] = 0
                elif c < d:
                    e_ = (sp.diff(Gm[a][d][b], X[c]) - sp.diff(Gm[a][c][b], X[d])
                          + sum(Gm[a][c][e]*Gm[e][d][b] - Gm[a][d][e]*Gm[e][c][b] for e in range(4)))
                    e_ = trunc(e_)
                    Rud[(a, b, c, d)] = e_
                    Rud[(a, b, d, c)] = -e_
Rdd = {k: trunc(g[k[0], k[0]]*v) for k, v in Rud.items()}
Ric = sp.Matrix(4, 4, lambda b, d: trunc(sum(Rud[(a, b, a, d)] for a in range(4))))
Rs = trunc(sum(gi[a, a]*Ric[a, a] for a in range(4)))
print('riemann (%.0fs); R scalar at t=-1, theta=pi/2, low orders:' % (time.time() - T0),
      sp.factor(trunc(Rs.subs(t, -1), 2)), flush=True)
C = {}
for a in range(4):
    for b in range(4):
        for c in range(4):
            for d in range(4):
                C[(a, b, c, d)] = trunc(Rdd[(a, b, c, d)]
                                        - (g[a, c]*Ric[b, d] - g[a, d]*Ric[b, c] - g[b, c]*Ric[a, d] + g[b, d]*Ric[a, c])/2
                                        + Rs*(g[a, c]*g[b, d] - g[a, d]*g[b, c])/6)
print('weyl (%.0fs)' % (time.time() - T0), flush=True)
# trace check of Weyl
trw = [trunc(sum(gi[a, a]*C[(a, b, a, d)] for a in range(4))) for b in range(4) for d in range(4)]
def minpow(e):
    e = sp.expand(e)
    if e == 0:
        return None
    return min(term.as_coeff_exponent(r)[1] for term in sp.Add.make_args(e))
print('Weyl trace residual: lowest r-power per component (None = exactly 0):', [minpow(x) for x in trw], flush=True)


def covC(e, a, c, b, d):
    # nabla_e C_acbd
    v = sp.diff(C[(a, c, b, d)], X[e])
    for f in range(4):
        v -= Gm[f][e][a]*C[(f, c, b, d)] + Gm[f][e][c]*C[(a, f, b, d)] + Gm[f][e][b]*C[(a, c, f, d)] + Gm[f][e][d]*C[(a, c, b, f)]
    return v


Kt = {}
for a in range(4):
    for c in range(4):
        for b in range(4):
            Kt[(a, c, b)] = trunc(sum(gi[d, d]*covC(d, a, c, b, d) for d in range(4)))
print('div weyl (%.0fs)' % (time.time() - T0), flush=True)
Rup = sp.Matrix(4, 4, lambda c, d: trunc(gi[c, c]*gi[d, d]*Ric[c, d]))
B = sp.zeros(4, 4)
for a in range(4):
    for b in range(a, 4):
        v = 0
        for c in range(4):
            w = sp.diff(Kt[(a, c, b)], X[c])
            for f in range(4):
                w -= Gm[f][c][a]*Kt[(f, c, b)] + Gm[f][c][c]*Kt[(a, f, b)] + Gm[f][c][b]*Kt[(a, c, f)]
            v += gi[c, c]*w
        v += sum(Rup[c, d]*C[(a, c, b, d)] for c in range(4) for d in range(4))/2
        B[a, b] = B[b, a] = trunc(v)
print('bach (%.0fs)' % (time.time() - T0), flush=True)
trB = trunc(sum(gi[a, a]*B[a, a] for a in range(4)))
sym = trunc(B[0, 1] - B[1, 0])
# frame components
Bnn = trunc(B[0, 0]*inv_series(alS**2))
Bnr = trunc(B[0, 1]*inv_series(trunc(alS*sqrt_series(AS))))
Brr = trunc(B[1, 1]*inv_series(AS))
BT = trunc(B[2, 2]/r**2)
res = {}
for nm, e in (('nn', Bnn), ('nr', Bnr), ('rr', Brr), ('T', BT)):
    e1 = sp.expand(e.subs(t, -1).subs(th, sp.pi/3))
    cs = {}
    for term in sp.Add.make_args(e1):
        c_, p_ = term.as_coeff_exponent(r)
        cs[p_] = cs.get(p_, 0) + c_
    lo = {int(p): sp.factor(cs[p]) for p in sorted(cs) if p <= 2}
    res[nm] = lo
    print('B_(%s) at t=-1: coefficients of r^p, p<=2:' % nm, lo, flush=True)
print('trace of B (orders <= 2) at t=-1:', sp.factor(trunc(trB.subs(t, -1).subs(th, sp.pi/3), 2)))
W0 = sp.Symbol('W0')
bnn4 = sp.factor(res['nn'].get(0, 0).subs(m0, 2*W0))
bnr5 = sp.factor(res['nr'].get(1, 0).subs(m0, 2*W0))
print('b_nn_4 = t^4 B_(nn)(t,0) =', bnn4, '  claim -2*W0*(8*W0-3)/3 ->', sp.simplify(bnn4 + sp.Rational(2, 3)*W0*(8*W0 - 3)) == 0)
print('b_nr_5 (coef of z^5 in r^4 B_(nr), t=-1) =', bnr5, '  claim -8*W0*(8*W0-3)/9 ->', sp.simplify(bnr5 + sp.Rational(8, 9)*W0*(8*W0 - 3)) == 0)
print('b_nr_5/b_nn_4 =', sp.simplify(res['nr'].get(1, 0)/res['nn'].get(0, 0)))
print('b_rr_4/b_nn_4 =', sp.simplify(res['rr'].get(0, 0)/res['nn'].get(0, 0)), ' b_T_4/b_nn_4 =', sp.simplify(res['T'].get(0, 0)/res['nn'].get(0, 0)))
print('v5 present in the leading coefficients:', any(sp.sympify(x).has(v5) for d_ in res.values() for x in d_.values()))
json.dump({k: {str(p): str(c) for p, c in d_.items()} for k, d_ in res.items()},
          open('<LAB>/quantum/seats/referee-4/r4_bach_%d_%s.json' % (KMAX, MODE), 'w'), indent=1)
print('done (%.0fs)' % (time.time() - T0))
