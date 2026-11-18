"""blind-qb2: Chebyshev collocation in zeta = x^2 on [0, X^2].

Background: Gauss-Newton on 4 polynomials (al, A, w, nu), v = x nu, residuals
  al' - F_al, A' - F_A, M (w', v') - b   (non-singular form; analytic branch at the sonic point)
  plus gauge al(0) = 1.
Perturbation (rate lam): unknowns da, dc, dw (even), dQ = x q (odd); 6 equations
  dE_uu, dE_ux, dE_xx, dE_thth (Einstein), ddiv_0, ddiv_1 (fluid), least squares, plus a
  gauge row (default da(0) = 0).
"""
import sys, numpy as np, sympy as sp, warnings
from numpy.polynomial import chebyshev as C
sys.path.insert(0, '<LAB>/quantum/seats/blind-qb2')
from qb2_common import *
warnings.filterwarnings('ignore')

sal1, sA1, Mbg, bbg = bg_system()
_args = (x, al, A, w, v)
F_al = sp.lambdify(_args, sal1, 'numpy')
F_A = sp.lambdify(_args, sA1, 'numpy')
F_M = [[sp.lambdify(_args, Mbg[i, j], 'numpy') for j in range(2)] for i in range(2)]
F_b = [sp.lambdify(_args, bbg[i], 'numpy') for i in range(2)]

class Grid:
    def __init__(self, N, X):
        self.N, self.X, self.Z = N, X, X*X
        k = np.arange(N)
        s = np.cos(np.pi*(2*k + 1)/(2*N))[::-1]        # Gauss-Chebyshev nodes (interior)
        self.s = s
        self.zeta = (s + 1)*self.Z/2
        self.x = np.sqrt(self.zeta)
        I = np.eye(N)
        self.T0 = np.array([C.chebval(s, I[j]) for j in range(N)]).T
        self.T1 = np.array([C.chebval(s, C.chebder(I[j])) for j in range(N)]).T*(2/self.Z)
        self.T2 = np.array([C.chebval(s, C.chebder(I[j], 2)) for j in range(N)]).T*(2/self.Z)**2
        self.T3 = np.array([C.chebval(s, C.chebder(I[j], 3)) for j in range(N)]).T*(2/self.Z)**3
        self.T4 = np.array([C.chebval(s, C.chebder(I[j], 4)) for j in range(N)]).T*(2/self.Z)**4
        self.e0 = np.array([C.chebval(-1.0, I[j]) for j in range(N)])   # value at zeta = 0

    def basis_at(self, xx):
        s = 2*np.asarray(xx)**2/self.Z - 1
        I = np.eye(self.N)
        return np.array([C.chebval(s, I[j]) for j in range(self.N)]).T

def xderivs(g, c, even=True, order=2, at=None):
    """values and x-derivatives of f(x) = P(zeta) (even) or x P(zeta) (odd) at the nodes."""
    x_ = g.x
    P = [g.T0 @ c, g.T1 @ c, g.T2 @ c, g.T3 @ c, g.T4 @ c]
    # d/dx P(x^2): P' 2x ; d2: 2P' + 4x^2 P'' ; d3: 12 x P'' + 8 x^3 P''' ; d4: 12P''+48x^2P'''+16x^4P''''
    f = [P[0], 2*x_*P[1], 2*P[1] + 4*x_**2*P[2], 12*x_*P[2] + 8*x_**3*P[3],
         12*P[2] + 48*x_**2*P[3] + 16*x_**4*P[4]]
    if even:
        return f[:order + 1]
    # odd: x * f
    o = [x_*f[0], f[0] + x_*f[1], 2*f[1] + x_*f[2], 3*f[2] + x_*f[3], 4*f[3] + x_*f[4]]
    return o[:order + 1]

def bg_unpack(g, cvec, order=2):
    N = g.N
    ca, cA, cw, cn = cvec[:N], cvec[N:2*N], cvec[2*N:3*N], cvec[3*N:]
    return (xderivs(g, ca, True, order), xderivs(g, cA, True, order),
            xderivs(g, cw, True, order), xderivs(g, cn, False, order))

def bg_residual(g, cvec):
    (a_, A_, w_, v_) = bg_unpack(g, cvec, 1)
    xx = g.x
    args = (xx, a_[0], A_[0], w_[0], v_[0])
    r1 = a_[1] - F_al(*args)
    r2 = A_[1] - F_A(*args)
    M00, M01 = F_M[0][0](*args), F_M[0][1](*args)
    M10, M11 = F_M[1][0](*args), F_M[1][1](*args)
    r3 = M00*w_[1] + M01*v_[1] - F_b[0](*args)
    r4 = M10*w_[1] + M11*v_[1] - F_b[1](*args)
    N = g.N
    gauge = np.array([g.e0 @ cvec[:N] - 1.0])
    return np.concatenate([r1, r2, r3, r4, gauge])

def bg_newton(g, cvec, iters=30, tol=1e-13, verbose=True):
    for it in range(iters):
        r = bg_residual(g, cvec)
        J = np.empty((r.size, cvec.size))
        h = 1e-7
        for j in range(cvec.size):
            cp = cvec.copy(); cp[j] += h
            cm = cvec.copy(); cm[j] -= h
            J[:, j] = (bg_residual(g, cp) - bg_residual(g, cm))/(2*h)
        dc = np.linalg.lstsq(J, -r, rcond=None)[0]
        cvec = cvec + dc
        rn = np.abs(bg_residual(g, cvec)).max()
        if verbose:
            print('  newton', it, 'res', rn, 'step', np.abs(dc).max(), flush=True)
        if np.abs(dc).max() < tol:
            break
    return cvec, rn

def fit_to_grid(g, xs, ys):
    """least-squares Chebyshev fit of sampled (x, [al, A, w, v]) data into coefficient vector."""
    B = g.basis_at(xs)
    out = []
    for i in range(4):
        yy = ys[i] if i < 3 else ys[i]/xs
        out.append(np.linalg.lstsq(B, yy, rcond=None)[0])
    return np.concatenate(out)

def bg_eval(g, cvec, xx):
    B = g.basis_at(xx)
    N = g.N
    return [B @ cvec[:N], B @ cvec[N:2*N], B @ cvec[2*N:3*N], xx*(B @ cvec[3*N:])]
