import sys, json, pickle; sys.path.insert(0, '<LAB>/quantum/seats/blind-qb2')
from qb2_coll import *
from qb2_bg import rhs, detM
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
w0v = 17.802642822265625
x0 = 1e-4
y0 = [1 + w0v/6*x0**2, 1 + w0v/3*x0**2, w0v, -0.5*x0]
r = solve_ivp(rhs, [x0, 0.80], y0, method='DOP853', rtol=1e-12, atol=1e-14, dense_output=True)
X = float(sys.argv[2]) if len(sys.argv) > 2 else 0.95
res = {}
prev = None
for N in [int(a) for a in sys.argv[1].split(',')]:
    g = Grid(N, X)
    if prev is None:
        xs = np.linspace(1e-3, 0.80, 400)
        ys = r.sol(xs)
        # extend beyond 0.80 by low-order polynomial extrapolation in x of the last stretch
        xe = np.linspace(0.80, X, 60)[1:]
        ext = []
        for i in range(4):
            p = np.polyfit(xs[-80:], ys[i][-80:], 3)
            ext.append(np.polyval(p, xe))
        xs2 = np.concatenate([xs, xe]); ys2 = [np.concatenate([ys[i], ext[i]]) for i in range(4)]
        cvec = fit_to_grid(g, xs2, ys2)
    else:
        gp, cp = prev
        xs = np.linspace(1e-3, X, 600)
        cvec = fit_to_grid(g, xs, bg_eval(gp, cp, xs))
    cvec, rn = bg_newton(g, cvec)
    xx = np.linspace(0.01, X, 2000)
    a_, A_, w_, v_ = bg_eval(g, cvec, xx)
    det = np.array([detM(xi, [a_[i], A_[i], w_[i], v_[i]]) for i, xi in enumerate(xx)])
    i = np.where(np.sign(det[:-1]) != np.sign(det[1:]))[0]
    fdet = lambda t: detM(t, bg_eval(g, cvec, np.array([t]))[0:1] and [q[0] for q in bg_eval(g, cvec, np.array([t]))])
    xs_ = brentq(fdet, xx[i[0]], xx[i[0]+1], xtol=1e-15) if len(i) else np.nan
    fv = lambda t: bg_eval(g, cvec, np.array([t]))[3][0]
    j = np.where(np.sign(v_[:-1]) != np.sign(v_[1:]))[0]
    xp = brentq(fv, xx[j[0]], xx[j[0]+1], xtol=1e-15)
    Ap = bg_eval(g, cvec, np.array([xp]))[1][0]
    w0 = bg_eval(g, cvec, np.array([0.0]))[2][0]
    res[N] = dict(N=N, X=X, res=rn, w0=w0, x_sonic=xs_, v_sonic=fv(xs_), x_peak=xp, A_peak=Ap, two_m_r_peak=1 - 1/Ap)
    print(res[N], flush=True)
    prev = (g, cvec)
    pickle.dump({'N': N, 'X': X, 'cvec': cvec.tolist()}, open('bg_N%d_X%.3f.pkl' % (N, X), 'wb'))
json.dump({str(k): {kk: float(vv) for kk, vv in d.items()} for k, d in res.items()}, open('bg_X%.3f.json' % X, 'w'), indent=1)
