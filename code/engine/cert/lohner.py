"""C^0/C^1 Lohner integrator with QR wrapping control (certified, arb balls).

Set representation:   y in  c + B r,   c thin centre, B thin (near-orthogonal), r interval vector.
Tangent vectors:      v_j in  vh_j + B t_j   (shared basis), propagated by the same step Jacobians.
One step of length h (x -> x+h, sign arbitrary), using taylor.step_cauchy2 (Cauchy remainder,
complex centred Picard a priori enclosure):
   p  = phi_h(c)                     (thin step, interval result)
   J  = D phi_h over X = hull(c + B r)   (Dual step over the box; includes derivative remainder)
   mean value:  phi_h(c + B r) in p + (J B) r
   A = mid(J B) = Q R (Gram-Schmidt in floats), Q exact, Qi = Q^{-1} rigorous
   c' = mid(p),  B' = Q,  r' = (Qi J B) r + Qi (p - c')
   v_j' in J v_j:  vh_j' = mid(J vh_j),  t_j' = (Qi J B) t_j + Qi (J vh_j - vh_j')
Rigorous because every true trajectory from the set stays in X during the step's start
(X contains the set) and J encloses D phi_h on X (Dual a priori enclosure over X).
"""
import time, math
from flint import arb, arb_mat
from core import Dual, sym
import taylor as TY
from param import val, _c

def _mid(x): return arb(x.mid())

def matvec(M, v):
    return [sum((M[i][j]*v[j] for j in range(len(v))), arb(0)) for i in range(len(M))]

def matmul(A, B):
    n, m, p = len(A), len(B), len(B[0])
    return [[sum((A[i][l]*B[l][j] for l in range(m)), arb(0)) for j in range(p)] for i in range(n)]

def qr_float(A):
    """Gram-Schmidt on float midpoints, columns ordered by norm (Lohner's pivoting)."""
    n = len(A)
    cols = [[float(A[i][j].mid()) for i in range(n)] for j in range(n)]
    order = sorted(range(n), key=lambda j: -math.sqrt(sum(x*x for x in cols[j])))
    Q = []
    for j in order:
        v = cols[j][:]
        for qv in Q:
            d = sum(a*b for a, b in zip(v, qv))
            v = [a - d*b for a, b in zip(v, qv)]
        nv = math.sqrt(sum(x*x for x in v))
        if nv == 0: raise ArithmeticError('degenerate QR')
        Q.append([x/nv for x in v])
    return [[arb(Q[j][i]) for j in range(n)] for i in range(n)]   # columns

class LohnerSet:
    def __init__(self, y, vecs=()):
        self.c = [_mid(v) for v in y]
        n = len(y)
        self.B = [[arb(1) if i == j else arb(0) for j in range(n)] for i in range(n)]
        self.r = [y[i] - self.c[i] for i in range(n)]
        self.vh = [[_mid(x) for x in v] for v in vecs]
        self.t = [[v[i] - self.vh[k][i] for i in range(n)] for k, v in enumerate(vecs)]
    def hull(self):
        return [self.c[i] + sum((self.B[i][j]*self.r[j] for j in range(len(self.r))), arb(0)) for i in range(len(self.c))]
    def vecs(self):
        return [[vh[i] + sum((self.B[i][j]*t[j] for j in range(len(t))), arb(0)) for i in range(len(vh))] for vh, t in zip(self.vh, self.t)]

def lohner_step(F, S, h, K, q):
    n = len(S.c)
    cv = lambda c: _c(c, arb(0))
    p, _, _ = TY.step_cauchy2(F, S.c, h, K, q=q, conv=cv)
    X = S.hull()
    z = arb(0); one = arb(1)
    Xd = [Dual(X[i], [one if j == i else z for j in range(n)]) for i in range(n)]
    yd, Bx, _ = TY.step_cauchy2(F, Xd, h, K, q=q, conv=lambda c: _c(c, Xd[0]))
    J = [[yd[i].d[j] for j in range(n)] for i in range(n)]
    M = matmul(J, S.B)
    Q = qr_float(M)
    Qi = arb_mat(Q).inv()
    Qi = [[Qi[i, j] for j in range(n)] for i in range(n)]
    QiM = matmul(Qi, M)
    c1 = [_mid(v) for v in p]
    r1 = [a + b for a, b in zip(matvec(QiM, S.r), matvec(Qi, [p[i] - c1[i] for i in range(n)]))]
    vh1 = []; t1 = []
    for vh, t in zip(S.vh, S.t):
        Jv = matvec(J, vh)
        nv = [_mid(x) for x in Jv]
        vh1.append(nv)
        t1.append([a + b for a, b in zip(matvec(QiM, t), matvec(Qi, [Jv[i] - nv[i] for i in range(n)]))])
    S.c, S.B, S.r, S.vh, S.t = c1, Q, r1, vh1, t1
    return Bx, J

def integrate(F, S, x0, x1, K=30, q='0.25', h0=0.002, grow=1.25, log=None, tmax=None):
    """advance LohnerSet S from x0 to x1 (x1 < x0 or > x0). Returns info (steps, boxes, time)."""
    T0 = time.time()
    q = arb(q)
    x = arb(x0); x1 = arb(x1)
    sgn = 1 if float(x1.mid()) > float(x.mid()) else -1
    h = h0; n = 0; boxes = []
    while (float(x1.mid()) - float(x.mid()))*sgn > 0:
        rem = abs(float((x1 - x).mid()))
        hh = min(h, rem); last = hh >= rem - 1e-15
        ha = (x1 - x) if last else arb(str(sgn*hh))
        try:
            Bx, J = lohner_step(F, S, ha, K, q)
        except ArithmeticError:
            h = hh/2
            if h < 1e-7: raise ArithmeticError('Lohner step collapse at x=%s' % x)
            continue
        boxes.append((float(x.mid()), float(ha.mid()), Bx))
        x = x1 if last else x + ha
        n += 1; h = hh*grow
        if log and n % 50 == 0:
            w = max(v.rad() for v in S.hull())
            log('  lohner step %d x=%.4f h=%.4g width=%s t=%.1fs' % (n, float(x.mid()), hh, w.str(3), time.time() - T0))
        if tmax and time.time() - T0 > tmax:
            return dict(steps=n, boxes=boxes, time=time.time() - T0, x=x, done=False)
    return dict(steps=n, boxes=boxes, time=time.time() - T0, x=x, done=True)
