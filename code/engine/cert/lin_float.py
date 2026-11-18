"""B4 proof seat, new file.  NON-CERTIFIED reference implementation (mpmath + scipy double) of the
characteristic function of the linearised Evans-Coleman problem.  It fixes the exact algorithm that
lin_cert.py must reproduce in ball arithmetic; nothing here counts.

Unknowns Y = (Abar_p, Nbar_p, wbar_p, V_p)  (HKA bar variables; system checked in lin_sym.py).
  sonic side : Frobenius series at x = 0 of the analytic solution s_kappa with
               (Abar_p, Nbar_p)(0) = (a, b) [default (0,1): no gauge zero], identity imposed at x = 0;
               Taylor sum at x = -h, then DOP853 (joint background + Y) to the matching point x_m.
  centre side: regular solutions rho_1, rho_2 = power series in v = p e^x of  dY/dx = Mc(v) Y,
               y_0 in ker Mc(0) = span(e_Nbar, e_wbar), (n - Mc_0) y_n = sum_{j<n} Mc_{n-j} y_j.
  D(kappa) = det[ s_kappa(x_m), rho_1(x_m), rho_2(x_m), e_V ].
Zero of D <=> s_kappa centre-regular <=> physical mode (PROOF-B4.md, sec. 2).
"""
import sys, pickle, math, cmath
sys.path.insert(0, '<LAB>/engine'); sys.path.insert(0, '<LAB>/engine/cert')
import numpy as np
import sympy as sp
from mpmath import mp, mpf, mpc, matrix, lu_solve, sqrt as msqrt, exp as mexp
from scipy.integrate import solve_ivp
import system as S
from bg import Background
from bgc import Centre
from core import Tape

mp.dps = 40
k = sp.Rational(1, 3); g = k + 1
N, om, V = S.N, S.om, S.V
omx_s, Vx_s = sp.symbols('omx Vx')
V0 = mpf('0.11243940138809833622')          # B3 certified enclosure midpoint (BAR-LB-E3 HOLD)
W0 = mpf('8.901323277755643519117')
PSCALE = mpf('0.808668956318442829477')     # centre scale p (DESIGN.md table)

# ---------------- symbolic HKA pencil with omx, Vx as free symbols (sonic series) -------------
def hka(A, omx, Vx):
    q = 1 - V**2; wbx = omx/om
    Ms = sp.Matrix([[1, g*V/q], [(g - 1)*V, g/q]])
    Mx = sp.Matrix([[1 + N*V, g*(N + V)/q], [(g - 1)*(N + V), g*(1 + N*V)/q]])
    E1 = -(g + 2)/2*A*N*V
    E2 = (6 - 3*g)/2*N*V - (2 + g)/2*A*N*V + (2 - g)*om*N*V - N*V*wbx - g*N*Vx/q
    E3 = (2 - g)*om*N*V
    E4 = (6 - 3*g)/2*N - (2 + g)/2*A*N + (2 - g)*om*N - N*wbx - g*(1 + 2*N*V + V**2)*Vx/q**2
    F1 = (2 - 3*g)/2*A*N
    F2 = (2 - g)*(g - 1)*om*N + (7*g - 6)/2*N + (2 - 3*g)/2*A*N - (g - 1)*N*wbx - g*N*V*Vx/q
    F3 = (2 - g)*(g - 1)*om*N
    F4 = -(g - 1)*wbx - g*(N + 2*V + N*V**2)*Vx/q**2
    P0 = sp.Matrix([[-A, 0, 2*(1 + (g - 1)*V**2)*om/q, 4*g*om*V/q**2], [A, 0, (g - 2)*om, 0],
                    [E1, E2, E3, E4], [F1, F2, F3, F4]])
    I0 = [-A, 2*g*N*V*om/q, 2*om*(1 + (g - 1)*V**2 + g*N*V)/q, 2*g*om*(N*(1 + V**2) + 2*V)/q**2]
    return Ms, Mx, P0, I0

A_s = S.Acons.subs(S.k, k)
Ms_s, Mx_s, P0_s, I0_s = hka(A_s, omx_s, Vx_s)
SON_OUT = list(Mx_s) + list(P0_s) + list(Ms_s) + I0_s           # 4 + 16 + 4 + 4 = 28

def mpconv(c):
    c = sp.Rational(c); return mpf(int(c.p))/int(c.q)

class Sonic:
    """kappa-independent background jets at the sonic point; kappa-dependent Frobenius recursion."""
    def __init__(self, K=60):
        self.K = K
        B = Background(dps=40, K=K + 2); self.B = B
        N0, om0 = B.sonic_point(V0, om_guess=mpf('0.3256888894230684'))
        co, info = B.sonic_series(V0, om0, 0, K=K + 2)
        self.co = co; self.info = info
        der = [[(n + 1)*c[n + 1] for n in range(len(c) - 1)] + [mpf(0)] for c in co]
        t = Tape(SON_OUT, [N, om, V, omx_s, Vx_s])
        t.jet_init(K + 1, mpf(0), conv=mpconv)
        ins = [co[0], co[1], co[2], der[1], der[2]]
        for n in range(K + 1):
            for i in range(5):
                t.C[i][n] = ins[i][n] if n < len(ins[i]) else mpf(0)
            if n == 0:
                for i, (op, a, b) in enumerate(t.nodes):
                    if op == 'const': t.C[i][0] = mpconv(a)
            t.jet_order(n)
        o = lambda j, n: t.C[t.outs[j]][n]
        self.Mx = [matrix([[o(0, n), o(1, n)], [o(2, n), o(3, n)]]) for n in range(K + 1)]
        self.P0 = [matrix([[o(4 + 4*i + j, n) for j in range(4)] for i in range(4)]) for n in range(K + 1)]
        self.Ms = [matrix([[o(20, n), o(21, n)], [o(22, n), o(23, n)]]) for n in range(K + 1)]
        self.I0 = [o(24 + j, 0) for j in range(4)]
        M = self.Mx[0]
        self.ell = matrix([[M[1, 0], -M[0, 0]]])
        self.r = matrix([[M[0, 1]], [-M[0, 0]]])
        self.Abar_x0 = co and None
        # background derivative data at x = 0 (A, N): Abar' = A_x/A, Nbar' = N_x/N
        Nv, ov, Vv = co[0][0], co[1][0], co[2][0]
        Av = A_s.subs({N: 0, om: 0, V: 0}) if False else None
        self.bg0 = (Nv, ov, Vv, co[0][1], co[1][1], co[2][1])
        # exponent rho(kappa) = (l E_u r - kappa l Ms r)/(l Mx' r)
        lr = lambda Mt: (self.ell*Mt*self.r)[0]
        Eu = matrix([[self.P0[0][2, 2], self.P0[0][2, 3]], [self.P0[0][3, 2], self.P0[0][3, 3]]])
        self.lMxr = lr(self.Mx[1]); self.lEr = lr(Eu); self.lMsr = lr(self.Ms[0])
        self.lam1 = self.lMxr/self.lMsr        # lambda_-'(0)
        self.d = self.lEr/self.lMsr

    def rho(self, kap):
        return (self.d - kap)/self.lam1

    def Pk(self, n, kap):
        P = self.P0[n].copy()
        for i in range(2):
            for j in range(2): P[2 + i, 2 + j] -= kap*self.Ms[n][i, j]
        return P

    def Qn(self, n):
        Q = matrix(4, 4)
        if n == 0: Q[0, 0] = 1; Q[1, 1] = 1
        for i in range(2):
            for j in range(2): Q[2 + i, 2 + j] = self.Mx[n][i, j]
        return Q

    def delta2(self, kap):
        """2x2 sonic determinant (wbar, V columns of [l.P_0 ; identity]): s_kappa has a pole where it vanishes"""
        P = self.Pk(0, kap)
        row1 = [(self.ell[0, 0]*P[2, j] + self.ell[0, 1]*P[3, j]) for j in range(4)]
        row2 = [self.I0[0] + kap] + self.I0[1:]
        return row1[2]*row2[3] - row1[3]*row2[2]

    def series(self, kap, a=0, b=1, K=None):
        K = K or self.K
        kap = mpc(kap)
        P = [self.Pk(n, kap) for n in range(K + 1)]
        Q = [self.Qn(n) for n in range(K + 1)]
        l4 = matrix([[0, 0, self.ell[0, 0], self.ell[0, 1]]])
        r4 = matrix([[0], [0], [self.r[0, 0]], [self.r[1, 0]]])
        Mc = matrix(4, 4); rc = matrix(4, 1)
        lp = l4*P[0]
        for j in range(4):
            Mc[0, j] = lp[j]; Mc[1, j] = (self.I0[j] + (kap if j == 0 else 0))
            Mc[2, j] = 1 if j == 0 else 0; Mc[3, j] = 1 if j == 1 else 0
        rc[2] = a; rc[3] = b
        y = [lu_solve(Mc, rc)]
        M0 = self.Mx[0]
        row = 0 if abs(M0[0, 0]) + abs(M0[0, 1]) > abs(M0[1, 0]) + abs(M0[1, 1]) else 1
        a1, a2 = M0[row, 0], M0[row, 1]; dd = a1*a1 + a2*a2
        def part(R, n1):
            w = matrix(4, 1)
            w[0] = R[0]/n1; w[1] = R[1]/n1
            bi = R[2 + row]/n1
            w[2] = bi*a1/dd; w[3] = bi*a2/dd
            return w
        def Rn(n, yn):
            ys = y[:n] + [yn]
            R = matrix(4, 1)
            for j in range(n + 1): R += P[n - j]*ys[j]
            for j in range(n): R -= Q[n - j]*((j + 1)*ys[j + 1])
            return R
        wp = part(P[0]*y[0], 1)
        for n in range(1, K):
            l0 = (l4*Rn(n, wp))[0]; l1 = (l4*Rn(n, wp + r4))[0]
            beta = -l0/(l1 - l0)
            yn = wp + beta*r4
            y.append(yn)
            wp = part(Rn(n, yn), n + 1)
        return y

    def eval(self, y, x):
        out = [mpc(0)]*4
        for c in reversed(y):
            out = [out[i]*x + c[i] for i in range(4)]
        return out

# ---------------- double-precision joint ODE --------------------------------------------------
f_bg = [e.subs(S.k, k) for e in S.rhs_regular()]
_fbg = sp.lambdify([N, om, V], f_bg, 'numpy')
Ms_r, Mx_r, P0_r, I0_r = hka(A_s, f_bg[1], f_bg[2])
_Mx = sp.lambdify([N, om, V], Mx_r, 'numpy'); _P0 = sp.lambdify([N, om, V], P0_r, 'numpy')
_Ms = sp.lambdify([N, om, V], Ms_r, 'numpy')
cS = sp.sqrt(k)
lam_p = 1 + N*(V + cS)/(1 + cS*V); lam_m = 1 + N*(V - cS)/(1 - cS*V)
lamd = [sp.diff(l, N)*f_bg[0] + sp.diff(l, V)*f_bg[2] for l in (lam_p, lam_m)]
_lam = sp.lambdify([N, om, V], [lam_p, lam_m, lamd[0], lamd[1]], 'numpy')
_A = sp.lambdify([N, om, V], A_s, 'numpy')

def joint_rhs(x, z, kap):
    n_, o_, v_ = z[0], z[1], z[2]
    fb = _fbg(n_, o_, v_)
    Y = z[3:7] + 1j*z[7:11]
    Mx = np.array(_Mx(n_, o_, v_), dtype=float); P0 = np.array(_P0(n_, o_, v_), dtype=complex)
    Ms = np.array(_Ms(n_, o_, v_), dtype=float)
    P0[2:, 2:] -= kap*Ms
    top = P0[:2] @ Y
    fl = np.linalg.solve(Mx, P0[2:] @ Y)
    dY = np.concatenate([top, fl])
    return np.concatenate([fb, dY.real, dY.imag])

def bg_rhs(x, z):
    return _fbg(z[0], z[1], z[2])

class Problem:
    def __init__(self, K=60, xm=-2.0, Kc=40):
        self.son = Sonic(K)
        self.xm = xm
        self.cen = Centre(self.son.B, K=Kc)
        self.Z = self.cen.series(W0, K=Kc)
        self.Kc = Kc
        self._centre_jets()

    def bg_at(self, x):
        """background (N, om, V) from the sonic series (|x| small) or the inward DOP853 path"""
        co = self.son.co
        return [float(sum(c[n]*mpf(x)**n for n in range(len(c)))) for c in co]

    def _centre_jets(self):
        d = pickle.load(open('<LAB>/engine/cert/lin_fields.pkl', 'rb'))
        Mc = [sp.sympify(eval('sp.' + s if False else s, {**sp.__dict__})) for s in d['Mc']]
        kap = S.kappa
        Mc0 = [e.subs(kap, 0) for e in Mc]; Mc1 = [sp.diff(e, kap) for e in Mc]
        t = Tape(Mc0 + Mc1, [S.P, S.W, S.U], {S.k: k})
        Kc = self.Kc
        t.jet_init(Kc + 1, mpf(0), conv=mpconv)
        for n in range(Kc + 1):
            for i in range(3):
                t.C[i][n] = self.Z[i][n] if n < len(self.Z[i]) else mpf(0)
            t.jet_order(n)
        self.C0 = [matrix([[t.C[t.outs[4*i + j]][n] for j in range(4)] for i in range(4)]) for n in range(Kc + 1)]
        self.C1 = [matrix([[t.C[t.outs[16 + 4*i + j]][n] for j in range(4)] for i in range(4)]) for n in range(Kc + 1)]

    def centre_regular(self, kap, x):
        kap = mpc(kap); Kc = self.Kc
        M = [self.C0[n] + kap*self.C1[n] for n in range(Kc + 1)]
        v = PSCALE*mexp(mpf(x))
        sols = []
        for e in (1, 2):
            y = [matrix(4, 1)]; y[0][e] = 1
            chk = M[0]*y[0]
            for n in range(1, Kc + 1):
                rhs = matrix(4, 1)
                for j in range(n): rhs += M[n - j]*y[j]
                A_ = matrix(4, 4)
                for i in range(4): A_[i, i] = n
                y.append(lu_solve(A_ - M[0], rhs))
            val = [sum(y[n][i]*v**n for n in range(Kc + 1)) for i in range(4)]
            sols.append(val)
        return sols, float(abs(chk[0]) + abs(chk[1]) + abs(chk[2]) + abs(chk[3]))

    def sonic_to(self, kap, a=0, b=1, h=None, x1=None, rtol=1e-12):
        x1 = self.xm if x1 is None else x1
        kap = complex(kap)
        if h is None: h = min(0.15, 4.0/max(abs(kap), 1.0))
        y = self.son.series(kap, a, b)
        Yh = self.son.eval(y, -mpf(h))
        z0 = self.bg_at(-h) + [float(c.real) for c in Yh] + [float(c.imag) for c in Yh]
        sol = solve_ivp(joint_rhs, [-h, x1], z0, args=(kap,), method='DOP853', rtol=rtol, atol=1e-14)
        z = sol.y[:, -1]
        return z[:3], z[3:7] + 1j*z[7:11], y

    def D(self, kap, a=0, b=1):
        bgm, Ys, _ = self.sonic_to(kap, a, b)
        rho, chk = self.centre_regular(kap, self.xm)
        Mt = np.zeros((4, 4), dtype=complex)
        Mt[:, 0] = Ys
        Mt[:, 1] = [complex(c) for c in rho[0]]; Mt[:, 2] = [complex(c) for c in rho[1]]
        Mt[3, 3] = 1
        return np.linalg.det(Mt), dict(bg=bgm, Ys=Ys, chk=chk)

if __name__ == '__main__':
    import time
    t0 = time.time()
    pr = Problem()
    son = pr.son
    print('[%.1fs] sonic jets ready; n* (bg) =' % (time.time() - t0), son.info['n_star'])
    print('lambda_-\'(0) =', son.lam1, '  d =', son.d, '  d + lambda\' - 1 =', son.d + son.lam1 - 1)
    print('rho(kappa) = (d - kappa)/lambda\' ; rho(0) =', son.rho(0), ' rho(1) =', son.rho(1))
    # background at x_m from both sides
    xm = pr.xm
    sol = solve_ivp(bg_rhs, [-0.15, xm], pr.bg_at(-0.15), method='DOP853', rtol=1e-13, atol=1e-15)
    yc, _ = pr.cen.eval(W0, PSCALE, mpf(xm), Z=pr.Z)
    print('background at x_m: inward', sol.y[:, -1], ' centre series', [float(c) for c in yc])
    for kap in (2.81055255, 0.3556992037, 1.0, 0.0, 2.0 + 3j, 1e-3j + 10j):
        d1, info = pr.D(kap)
        print('kappa=%s  D=%s  |Mc0 y0|=%.1e  Delta2=%s' % (kap, d1, info['chk'], complex(son.delta2(kap))))
    print('[%.1fs] done' % (time.time() - t0))
