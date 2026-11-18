"""Lohner variants for attempt 3 (new file).  Same set representation as lohner.py.
mode 'qrw' : QR of A = mid(J B) with columns sorted by ||A_j|| * |r_j|  (Lohner's weighted pivoting)
mode 'par' : parallelepiped: B' = mid(J B), r' = (B'^{-1} J B) r + B'^{-1}(p - c')
mode 'qr'  : as lohner.py (columns sorted by ||A_j|| only)
Extra: the step-start set and the step Jacobian enclosure are returned, so derivative sets can be
carried as their own Lohner sets (C^1-Lohner, see lohner_step_c1)."""
import math
from flint import arb, arb_mat
from core import Dual, sym
import taylor as TY
from param import val, _c
import lohner as LO

def qr_w(A, w):
    n = len(A)
    cols = [[float(A[i][j].mid()) for i in range(n)] for j in range(n)]
    order = sorted(range(n), key=lambda j: -math.sqrt(sum(x*x for x in cols[j]))*w[j])
    Q = []
    for j in order:
        v = cols[j][:]
        for qv in Q:
            d = sum(a*b for a, b in zip(v, qv)); v = [a - d*b for a, b in zip(v, qv)]
        nv = math.sqrt(sum(x*x for x in v))
        if nv == 0: raise ArithmeticError('degenerate QR')
        Q.append([x/nv for x in v])
    return [[arb(Q[j][i]) for j in range(n)] for i in range(n)]

def lohner_step2(F, S, h, K, q, mode='qrw'):
    n = len(S.c)
    p, _, _ = TY.step_cauchy2(F, S.c, h, K, q=q, conv=lambda c: _c(c, arb(0)))
    X = S.hull()
    z = arb(0); one = arb(1)
    Xd = [Dual(X[i], [one if j == i else z for j in range(n)]) for i in range(n)]
    yd, Bx, _ = TY.step_cauchy2(F, Xd, h, K, q=q, conv=lambda c: _c(c, Xd[0]))
    J = [[yd[i].d[j] for j in range(n)] for i in range(n)]
    M = LO.matmul(J, S.B)
    if mode == 'par':
        Bn = [[arb(M[i][j].mid()) for j in range(n)] for i in range(n)]
    else:
        w = [float(t.abs_upper()) + 1e-300 for t in S.r] if mode == 'qrw' else [1.0]*n
        Bn = qr_w(M, w)
    Bi = arb_mat(Bn).inv(); Bi = [[Bi[i, j] for j in range(n)] for i in range(n)]
    BiM = LO.matmul(Bi, M)
    c1 = [arb(v.mid()) for v in p]
    r1 = [a + b for a, b in zip(LO.matvec(BiM, S.r), LO.matvec(Bi, [p[i] - c1[i] for i in range(n)]))]
    vh1 = []; t1 = []
    for vh, t in zip(S.vh, S.t):
        Jv = LO.matvec(J, vh); nv = [arb(x.mid()) for x in Jv]
        vh1.append(nv)
        t1.append([a + b for a, b in zip(LO.matvec(BiM, t), LO.matvec(Bi, [Jv[i] - nv[i] for i in range(n)]))])
    S.c, S.B, S.r, S.vh, S.t = c1, Bn, r1, vh1, t1
    return Bx, J
