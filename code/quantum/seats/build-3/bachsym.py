"""build-3 QB2, exact sympy (non-counting derivation, output feeds the certified Tapes).

Conformal reduction (PAPER, checked numerically against the QB1 instrument in n0check.py):
  g = r^2 gt,  r = e^{x-s},  gt = -N^2 A ds^2 + A (dx - ds)^2 + dth^2 + sin^2 th dph^2   (s, x coordinates)
  (identical to referee-2 r2_forced.py metric(); A, N functions of x only on the background).
In 4D the Bach tensor has conformal weight -2:  B_ab[g] = r^{-2} B_ab[gt], so the frame components obey
  b_ab := r^4 B_(ab)[g] = Bt_(ab)[gt],
with the gt-frame n = (d_s + d_x)/(N sqrt A), e = d_x/sqrt A, e_th = d_th.  gt depends on x only.
Bach convention here: B_ab = nabla^c nabla^d C_acbd + (1/2) R^cd C_acbd (the QB1 normalisation is fixed by
comparison, n0check.py).

Outputs (bach_exprs.txt, srepr):  b_nn, b_nr, b_rr, b_T as polynomials/rational functions of
  A0..A4 (A and its x-derivatives) and l1..l4 (derivatives of l = ln N; l itself must not appear).
Checks printed: no l0; traceless b_nn = b_rr + 2 b_T; conservation D_t == 0 and D_x == 0 identically
(dt_struct.log formula for D_t; D_x from the covariant divergence) with sigma built from beta = e^{-2x} b.
"""
import time, sympy as sp
t0 = time.time()
s, x, th, ph = sp.symbols('s x theta phi', real=True)
X = [s, x, th, ph]
Af = sp.Function('A')(x); lf = sp.Function('l')(x)
Nf = sp.exp(lf)
g = sp.zeros(4, 4)
g[0, 0] = -Nf**2*Af + Af; g[0, 1] = g[1, 0] = -Af; g[1, 1] = Af
g[2, 2] = 1; g[3, 3] = sp.sin(th)**2
gi = g.inv().applyfunc(sp.simplify)
n = 4
Gam = [[[sp.simplify(sum(gi[a, d]*(sp.diff(g[d, b], X[c]) + sp.diff(g[d, c], X[b]) - sp.diff(g[b, c], X[d]))
                         for d in range(n))/2) for c in range(n)] for b in range(n)] for a in range(n)]
# Riemann R^a_{bcd} = d_c Gam^a_{db} - d_d Gam^a_{cb} + Gam^a_{ce} Gam^e_{db} - Gam^a_{de} Gam^e_{cb}
Rud = [[[[sp.simplify(sp.diff(Gam[a][d][b], X[c]) - sp.diff(Gam[a][c][b], X[d])
                      + sum(Gam[a][c][e]*Gam[e][d][b] - Gam[a][d][e]*Gam[e][c][b] for e in range(n)))
          for d in range(n)] for c in range(n)] for b in range(n)] for a in range(n)]
Rdd = [[[[sp.simplify(sum(g[a, e]*Rud[e][b][c][d] for e in range(n))) for d in range(n)] for c in range(n)] for b in range(n)] for a in range(n)]
Ric = sp.Matrix(n, n, lambda b, d: sp.simplify(sum(Rud[a][b][a][d] for a in range(n))))
Rs = sp.simplify(sum(gi[a, b]*Ric[a, b] for a in range(n) for b in range(n)))
print('curvature %.0fs' % (time.time() - t0), flush=True)
# Weyl (all lower)
def Cw(a, b, c, d):
    return (Rdd[a][b][c][d] - (g[a, c]*Ric[b, d] - g[a, d]*Ric[b, c] - g[b, c]*Ric[a, d] + g[b, d]*Ric[a, c])/2
            + Rs*(g[a, c]*g[b, d] - g[a, d]*g[b, c])/6)
C = [[[[sp.simplify(Cw(a, b, c, d)) for d in range(n)] for c in range(n)] for b in range(n)] for a in range(n)]
print('weyl %.0fs' % (time.time() - t0), flush=True)
# nabla_d C_{abce} etc.  Bach_{ab} = nabla^c nabla^d C_{acbd} + 1/2 R^{cd} C_{acbd}
def cov_C(T):
    """nabla_e T_{abcd} -> D[e][a][b][c][d]"""
    D = [[[[[None]*n for _ in range(n)] for _ in range(n)] for _ in range(n)] for _ in range(n)]
    for e in range(n):
        for a in range(n):
            for b in range(n):
                for c in range(n):
                    for d in range(n):
                        v = sp.diff(T[a][b][c][d], X[e])
                        v -= sum(Gam[m][e][a]*T[m][b][c][d] + Gam[m][e][b]*T[a][m][c][d]
                                 + Gam[m][e][c]*T[a][b][m][d] + Gam[m][e][d]*T[a][b][c][m] for m in range(n))
                        D[e][a][b][c][d] = v
    return D
# contract first: W_{abc} := nabla^d C_{acbd} = g^{de} nabla_e C_{acbd}
DC = cov_C(C)
W = [[[sp.simplify(sum(gi[d, e]*DC[e][a][c][b][d] for d in range(n) for e in range(n))) for c in range(n)] for b in range(n)] for a in range(n)]
print('div weyl %.0fs' % (time.time() - t0), flush=True)
# nabla^c W_{abc}
B = sp.zeros(n, n)
Ricu = gi*Ric*gi
for a in range(n):
    for b in range(a, n):
        v = 0
        for c in range(n):
            for f in range(n):
                if gi[c, f] == 0:
                    continue
                dv = sp.diff(W[a][b][c], X[f]) - sum(Gam[m][f][a]*W[m][b][c] + Gam[m][f][b]*W[a][m][c] + Gam[m][f][c]*W[a][b][m] for m in range(n))
                v += gi[c, f]*dv
        v += sum(Ricu[c, d]*C[a][c][b][d] for c in range(n) for d in range(n))/2
        B[a, b] = B[b, a] = sp.simplify(v)
print('bach %.0fs' % (time.time() - t0), flush=True)
tr = sp.simplify(sum(gi[a, b]*B[a, b] for a in range(n) for b in range(n)))
print('trace Bach == 0:', tr == 0)
# frame components
nN = 1/(Nf*sp.sqrt(Af)); eN = 1/sp.sqrt(Af)
bnn = sp.simplify(nN**2*(B[0, 0] + 2*B[0, 1] + B[1, 1]))
bnr = sp.simplify(nN*eN*(B[0, 1] + B[1, 1]))
brr = sp.simplify(eN**2*B[1, 1])
bT = sp.simplify(B[2, 2])
# plain symbols
As = sp.symbols('A0:6'); ls = sp.symbols('l0:6')
rep = {}
for j in range(5, 0, -1):
    rep[sp.Derivative(Af, (x, j))] = As[j]; rep[sp.Derivative(lf, (x, j))] = ls[j]
def plain(e):
    e = e.subs(rep).subs({Af: As[0], lf: ls[0]})
    return sp.simplify(e)
outs = [plain(e) for e in (bnn, bnr, brr, bT)]
for nm, e in zip(('b_nn', 'b_nr', 'b_rr', 'b_T'), outs):
    print(nm, 'free:', sorted(map(str, e.free_symbols)), ' ops', sp.count_ops(e))
print('traceless b_nn - b_rr - 2 b_T == 0:', sp.simplify(outs[0] - outs[2] - 2*outs[3]) == 0)
with open('<LAB>/quantum/seats/build-3/bach_exprs.txt', 'w') as fh:
    for e in outs:
        fh.write(sp.srepr(e) + '\n')
# ---- conservation identity (dt_struct.log D_t, kappa = 2) and D_x, sigma from beta = e^{-2x} b
Ab = Af; Nb = Nf
bf = [bnn, bnr, brr, bT]
beta = [sp.exp(-2*x)*e for e in bf]
fss = Ab*(Nb**2*beta[0] - 2*Nb*beta[1] + beta[2]); fsx = Ab*(Nb*beta[1] - beta[2]); fxx = Ab*beta[2]; fTT = beta[3]
kap = 2
Dt = ((1/Ab - 2/(Ab*Nb**2))*sp.diff(fsx, x) + (1/Ab - 1/(Ab*Nb**2))*sp.diff(fxx, x)
      + (-kap/(Ab*Nb**2) + 2*sp.diff(Nb, x)/(Ab*Nb**3) + sp.diff(Ab, x)/(2*Ab**2*Nb**2))*fss
      + (-2*kap/(Ab*Nb**2) + 2/Ab + sp.diff(Nb, x)/(Ab*Nb) + 4*sp.diff(Nb, x)/(Ab*Nb**3) + sp.diff(Ab, x)/(Ab**2*Nb**2))*fsx
      + (-kap/(Ab*Nb**2) + 2/Ab + sp.diff(Nb, x)/(Ab*Nb) + 2*sp.diff(Nb, x)/(Ab*Nb**3) - sp.diff(Ab, x)/(2*Ab**2) + sp.diff(Ab, x)/(2*Ab**2*Nb**2))*fxx
      - sp.diff(fss, x)/(Ab*Nb**2))
print('D_t == 0 identically:', sp.simplify(plain(Dt)) == 0, flush=True)
# D_x and D_s directly: covariant divergence of Sigma = e^{2s} sigma_{mu nu}(x) w.r.t. g = r^2 gt
gfull = sp.exp(2*(x - s))*g
gfi = gfull.inv().applyfunc(sp.simplify)
GamF = [[[sp.simplify(sum(gfi[a, d]*(sp.diff(gfull[d, b], X[c]) + sp.diff(gfull[d, c], X[b]) - sp.diff(gfull[b, c], X[d]))
                          for d in range(n))/2) for c in range(n)] for b in range(n)] for a in range(n)]
Sig = sp.zeros(4, 4)
ek = sp.exp(2*s)
Sig[0, 0] = ek*fss; Sig[0, 1] = Sig[1, 0] = ek*fsx; Sig[1, 1] = ek*fxx; Sig[2, 2] = ek*fTT; Sig[3, 3] = ek*fTT*sp.sin(th)**2
Smix = gfi*Sig
divs = []
for nu in range(2):
    e = sum(sp.diff(Smix[m, nu], X[m]) for m in range(4))
    e += sum(GamF[m][m][l]*Smix[l, nu] for m in range(4) for l in range(4))
    e -= sum(GamF[l][m][nu]*Smix[m, l] for m in range(4) for l in range(4))
    divs.append(sp.simplify(plain(e)))
print('(div Sigma)_s == 0:', divs[0] == 0, ' (div Sigma)_x == 0:', divs[1] == 0)
print('total %.0fs' % (time.time() - t0))
