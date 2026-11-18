"""B1 background numerics (non-rigorous, mpmath): Evans-Coleman CSS solution.

Chart: y = (N, om, V), A = Acons (Maison (8)); sonic point fixed at x = 0
(KHA sec. 4 gauge).  Steps:
  sonic_point(V0)          -> (N0, om0) on the sonic line
  sonic_series(V0, branch) -> Taylor coefficients of the analytic branch
  taylor_integrate         -> regular-point Taylor integration
  residual(V0)             -> centre-regularity residual (blown-up U = N V)
"""
import sys, time
import sympy as sp
from mpmath import mp, mpf, sqrt, findroot, matrix, polyroots, log, fabs, exp
sys.path.insert(0, '<LAB>/engine')
import system as S
from jet import Tape

class Background:
    def __init__(self, k=sp.Rational(1, 3), dps=30, K=30):
        mp.dps = dps
        self.k = k
        self.kf = mpf(k.p)/k.q
        self.c = sqrt(self.kf)
        self.K = K
        ins = [S.N, S.om, S.V]
        pk = {S.k: k}
        self.tS = Tape(list(S.Msym) + list(S.bsym) + [S.gN], ins, pk, K)   # M00 M01 M10 M11 b0 b1 gN
        self.tR = Tape(S.rhs_regular(), ins, pk, K)
        self.Ufix = -mpf(2)/(3*(1 + self.kf))
        lam = {s: sp.lambdify([S.N, S.om, S.V], e.subs(S.k, k), 'mpmath')
               for s, e in [('b0', S.bsym[0]), ('b1', S.bsym[1]), ('M00', S.Msym[0, 0]),
                            ('M01', S.Msym[0, 1]), ('M10', S.Msym[1, 0]), ('M11', S.Msym[1, 1]),
                            ('A', S.Acons.subs(S.k, k))]}
        self.lam = lam

    # ---------- sonic line
    def sonic_point(self, V0, om_guess=None):
        c = self.c
        N0 = (1 - V0*c)/(c - V0)
        L = self.lam
        def ell_b(o):
            M00, M10 = L['M00'](N0, o, V0), L['M10'](N0, o, V0)
            return M10*L['b0'](N0, o, V0) - M00*L['b1'](N0, o, V0)
        if om_guess is None:
            # scan for positive roots
            roots = []
            grid = [mpf(i)/200 for i in range(1, 400)]
            vals = [ell_b(o) for o in grid]
            for i in range(len(grid) - 1):
                if vals[i]*vals[i + 1] < 0:
                    roots.append(findroot(ell_b, (grid[i], grid[i + 1]), solver='anderson'))
            return N0, roots
        return N0, findroot(ell_b, om_guess)

    # ---------- analytic branch through the sonic point
    def sonic_series(self, V0, om0, branch=0, K=None):
        K = K or self.K
        t = self.tS; t.reset(K + 1)
        N0 = (1 - V0*self.c)/(self.c - V0)
        y = [[mpf(0)]*(K + 2) for _ in range(3)]
        y[0][0], y[1][0], y[2][0] = N0, om0, V0
        def put(n):
            for i in range(3):
                t.set_input(i, n, y[i][n])
            t.order(n)
        put(0)
        M0 = [[t.out(0, 0), t.out(1, 0)], [t.out(2, 0), t.out(3, 0)]]
        ell = (M0[1][0], -M0[0][0])
        r = (M0[0][1], -M0[0][0])
        nr = sqrt(r[0]**2 + r[1]**2); r = (r[0]/nr, r[1]/nr)
        row = 0 if abs(M0[0][0]) + abs(M0[0][1]) > abs(M0[1][0]) + abs(M0[1][1]) else 1
        a1, a2 = M0[row]
        def minnorm(rhs):  # M0 w = rhs (rank one)
            bi = rhs[row]
            d = a1**2 + a2**2
            return (bi*a1/d, bi*a2/d)
        def Mn(m):
            return [[t.out(0, m), t.out(1, m)], [t.out(2, m), t.out(3, m)]]
        def Rn(n):  # residual coefficient n with w_{n+1} = 0
            R = [-t.out(4, n), -t.out(5, n)]
            for j in range(n):  # j+1 <= n
                Mm = Mn(n - j)
                wj = (y[1][j + 1], y[2][j + 1])
                R[0] += Mm[0][0]*(j + 1)*wj[0] + Mm[0][1]*(j + 1)*wj[1]
                R[1] += Mm[1][0]*(j + 1)*wj[0] + Mm[1][1]*(j + 1)*wj[1]
            return R
        b0 = (t.out(4, 0), t.out(5, 0))
        wp = minnorm(b0)
        y[0][1] = t.out(6, 0)
        info = {'slopes': [], 'roots1': None, 'ell': ell, 'r': r}
        for n in range(1, K + 1):
            def lR(beta):
                y[1][n] = wp[0] + beta*r[0]; y[2][n] = wp[1] + beta*r[1]
                put(n)
                R = Rn(n)
                return ell[0]*R[0] + ell[1]*R[1]
            if n == 1:
                qm, q0, qp = lR(-1), lR(0), lR(1)
                q2 = (qp + qm)/2 - q0; q1 = (qp - qm)/2
                disc = q1**2 - 4*q2*q0
                rts = sorted([(-q1 + sqrt(disc))/(2*q2), (-q1 - sqrt(disc))/(2*q2)], key=lambda z: mp.re(z))
                info['roots1'] = rts
                beta = rts[branch]
                slope = 2*q2*beta + q1
            else:
                l0 = lR(0); l1 = lR(1)
                slope = l1 - l0
                beta = -l0/slope
            info['slopes'].append(slope)
            lR(beta)
            R = Rn(n)
            info.setdefault('lres', []).append(ell[0]*R[0] + ell[1]*R[1])
            if n < K:
                wp = minnorm((-R[0]/(n + 1), -R[1]/(n + 1)))
                y[0][n + 1] = t.out(6, n)/(n + 1)
        # resonance index from slopes (affine in n for n >= 2)
        s2, s3 = info['slopes'][1], info['slopes'][2]
        b_ = s3 - s2; a_ = s2 - 2*b_
        info['n_star'] = -a_/b_
        # truncate last (beta_K free part handled): drop index K
        coeffs = [y[i][:K] for i in range(3)]
        return coeffs, info

    # ---------- regular Taylor integrator
    def taylor_step(self, y0, K=None):
        K = K or self.K
        t = self.tR; t.reset(K + 1)
        Y = [[mpf(0)]*(K + 2) for _ in range(3)]
        for i in range(3):
            Y[i][0] = y0[i]
        for n in range(K):
            for i in range(3):
                t.set_input(i, n, Y[i][n])
            t.order(n)
            for i in range(3):
                Y[i][n + 1] = t.out(i, n)/(n + 1)
        return [Y[i][:K + 1] for i in range(3)]

    @staticmethod
    def radius(coeffs):
        K = len(coeffs[0]) - 1
        est = []
        for n in range(K//2, K + 1):
            m = max(abs(c[n]) for c in coeffs)
            if m > 0:
                est.append(m**(-mpf(1)/n))
        return min(est) if est else mpf(1)

    @staticmethod
    def evalser(coeffs, h):
        out = []
        for c in coeffs:
            v = mpf(0)
            for a in reversed(c):
                v = v*h + a
            out.append(v)
        return out

    def integrate(self, x0, y0, x1, frac=mpf('0.35'), K=None, maxsteps=400, record=False):
        x, y = mpf(x0), list(y0)
        sgn = 1 if x1 > x0 else -1
        path = [(x, y)]
        for it in range(maxsteps):
            co = self.taylor_step(y, K)
            h = frac*self.radius(co)
            if abs(x1 - x) <= h:
                h = abs(x1 - x)
                y = self.evalser(co, sgn*h); x = mpf(x1)
                path.append((x, y))
                return (x, y, path) if record else (x, y)
            y = self.evalser(co, sgn*h); x += sgn*h
            if record: path.append((x, y))
        raise RuntimeError('maxsteps at x=%s' % x)

    def leave_sonic(self, V0, om0, branch, h0=mpf('0.25')):
        co, info = self.sonic_series(V0, om0, branch)
        rho = self.radius(co)
        h = min(h0, rho/3)
        return co, info, h

    def residual(self, V0, om0=None, branch=0, L=10, return_all=False):
        if om0 is None:
            _, rts = self.sonic_point(V0); om0 = rts[0]
        co, info, h = self.leave_sonic(V0, om0, branch)
        ym = self.evalser(co, -h)
        x, y = self.integrate(-h, ym, -L)
        U = y[0]*y[2]
        res = U - self.Ufix
        if return_all:
            return res, dict(info=info, y=y, h=h, coeffs=co)
        return res
