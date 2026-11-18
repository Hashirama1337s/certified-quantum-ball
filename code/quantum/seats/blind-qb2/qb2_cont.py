import sys; sys.path.insert(0, '<LAB>/quantum/seats/blind-qb2')
from qb2_bg import *
import numpy as np
def through(w0v, th=1e-4, X=6.0):
    x0 = 1e-4
    y0 = [1 + w0v/6*x0**2, 1 + w0v/3*x0**2, w0v, -0.5*x0]
    f = lambda xx, y: abs(detM(xx, y)) - th
    f.terminal = True; f.direction = -1
    r1 = solve_ivp(rhs, [x0, 2.0], y0, method='DOP853', rtol=1e-12, atol=1e-14, events=f, dense_output=True)
    xe = r1.t[-1]; ye = r1.y[:, -1]
    # linear extrapolation across using slope at xe, find point where |det| = th on the far side
    d = np.array(rhs(xe, ye))
    for jump in np.linspace(0.001, 0.05, 50):
        yj = ye + d*jump
        if abs(detM(xe + jump, yj)) > th:
            break
    g = lambda xx, y: abs(detM(xx, y)) - th/10
    g.terminal = True; g.direction = -1
    r2 = solve_ivp(rhs, [xe + jump, X], yj, method='DOP853', rtol=1e-12, atol=1e-14, events=g, dense_output=True)
    return r1, r2, xe, jump
if __name__ == '__main__':
    for w0v in [0.3, 0.5, 0.7]:
        r1, r2, xe, jump = through(w0v)
        print('w0', w0v, 'xs~', xe, 'jump', jump, 'r2 end', r2.t[-1], r2.status, flush=True)
        for xx in np.linspace(r2.t[0], r2.t[-1], 12):
            y = r2.sol(xx)
            print('   x %.4f al %.5f A %.6f w %.5f v %.5f  2m/r %.6f' % (xx, *y, 1 - 1/y[1]))
