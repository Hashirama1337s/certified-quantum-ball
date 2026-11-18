"""B1: regular-centre series (blown-up chart) + 3x3 matching for the EC background.
Centre chart Z = (P, W, U), P = 1/N, om = W P^2, V = U P, u = e^x, u dZ/du = G(Z).
Regular centre: Z0 = (0, W0, -2/(3(1+k))), P1 = p (scaling), J0 eigenvalues {1, 0, -3}.
"""
import sys
sys.path.insert(0, '<LAB>/engine')
import sympy as sp
from mpmath import mp, mpf, matrix, lu_solve, findroot, exp, log, nstr
import system as S
from jet import Tape
from bg import Background

class Centre:
    def __init__(self, B, K=40):
        self.B = B; self.K = K
        G = S.rhs_centre_blown()
        self.t = Tape(G, [S.P, S.W, S.U], {S.k: B.k}, K)

    def series(self, W0, K=None):
        """coefficients of Z(v) with P = v + ..., v = p u"""
        K = K or self.K
        t = self.t; t.reset(K + 1)
        Z = [[mpf(0)]*(K + 1) for _ in range(3)]
        Z[1][0] = W0; Z[2][0] = self.B.Ufix
        def put(n):
            for i in range(3): t.set_input(i, n, Z[i][n])
            t.order(n)
        put(0)
        Z[0][1] = mpf(1)
        # order 1: W1, U1 from (1 - J0) restricted; get affine map numerically
        for n in range(1, K + 1):
            keep0 = Z[0][n] if n == 1 else None
            base = [Z[i][n] for i in range(3)]
            # known part with Z_n = (P1 fixed if n==1, else 0)
            for i in range(3): Z[i][n] = mpf(0)
            if n == 1: Z[0][n] = keep0
            put(n); g0 = [t.out(i, n) for i in range(3)]
            J = [[None]*3 for _ in range(3)]
            for j in range(3):
                if n == 1 and j == 0: continue
                Z[j][n] += 1; put(n)
                for i in range(3): J[i][j] = t.out(i, n) - g0[i]
                Z[j][n] -= 1
            # solve n Z_n = g0 + J Z_n for free components
            idx = [1, 2] if n == 1 else [0, 1, 2]
            Mx = matrix(len(idx), len(idx)); rh = matrix(len(idx), 1)
            for a, i in enumerate(idx):
                rh[a] = g0[i]
                for b_, j in enumerate(idx):
                    Mx[a, b_] = (n if i == j else 0) - J[i][j]
            if n == 1:
                # P equation at order 1 must be consistent: 1*P1 = g0[0]
                self.check1 = g0[0] - Z[0][1]
            sol = lu_solve(Mx, rh)
            for a, i in enumerate(idx): Z[i][n] = sol[a]
            put(n)
        return Z

    def eval(self, W0, p, x, Z=None):
        Z = Z or self.series(W0)
        v = p*exp(x)
        P, W, U = Background.evalser(Z, v)
        return [1/P, W*P**2, U*P], Z

def match(B, C, V0, W0, lp, xm, frac):
    N0, om0 = B.sonic_point(V0, om_guess=B._om)
    B._om = om0
    co, info = B.sonic_series(V0, om0, 0)
    h = frac*B.radius(co)
    y = B.evalser(co, -h)
    x, y = B.integrate(-h, y, xm, frac=frac)
    yc, Z = C.eval(W0, exp(lp), xm)
    return [(y[0] - yc[0])/yc[0], y[1]/yc[1] - 1, y[2] - yc[2]]
