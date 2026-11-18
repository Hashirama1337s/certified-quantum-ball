"""blind-qb2 background: Evans-Coleman CSS radiation fluid, own derivation (qb2_derive.py).

Step 1: centre series + shooting on w0 = 8 pi rho(0) t^2 (alpha(0)=1) to the sonic point.
Step 2 (qb2_coll.py): Newton Chebyshev collocation through the sonic point.
"""
import sys, json, numpy as np, sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
sys.path.insert(0, '<LAB>/quantum/seats/blind-qb2')
from qb2_common import *

sal1, sA1, M, b = bg_system()
f_al1 = sp.lambdify((x, al, A, w, v), sal1, 'numpy')
f_A1 = sp.lambdify((x, al, A, w, v), sA1, 'numpy')
f_M = sp.lambdify((x, al, A, w, v), M, 'numpy')
f_b = sp.lambdify((x, al, A, w, v), b, 'numpy')

def rhs(xx, y):
    Mm = np.array(f_M(xx, *y), dtype=float)
    bb = np.array(f_b(xx, *y), dtype=float).ravel()
    wv = np.linalg.solve(Mm, bb)
    return [f_al1(xx, *y), f_A1(xx, *y), wv[0], wv[1]]

def detM(xx, y):
    Mm = np.array(f_M(xx, *y), dtype=float)
    return np.linalg.det(Mm) / (np.abs(Mm).max()**2)

def numer(xx, y):
    Mm = np.array(f_M(xx, *y), dtype=float)
    bb = np.array(f_b(xx, *y), dtype=float).ravel()
    Mw = Mm.copy(); Mw[:, 0] = bb
    return np.linalg.det(Mw) / (np.abs(Mm).max()*max(np.abs(bb).max(), 1e-300))

p2, q2, r2, s1, s3 = sp.symbols('p2 q2 r2 s1 s3')
def centre_solution():
    w0s = sp.Symbol('w0s', positive=True)
    ser = {al: 1 + p2*x**2, A: 1 + q2*x**2, w: w0s + r2*x**2, v: s1*x + s3*x**3}
    dser = {al1: sp.diff(ser[al], x), A1: sp.diff(ser[A], x), w1: sp.diff(ser[w], x), v1: sp.diff(ser[v], x)}
    eqs = []
    for k in ['Ebg_xx', 'Ebg_ux', 'divbg_0', 'divbg_1']:
        e_ = EX[k].subs(dser).subs(ser)
        s = sp.series(e_, x, 0, 4).removeO()
        eqs += [sp.factor(sp.expand(s.coeff(x, n))) for n in range(4)]
    eqs = [e_ for e_ in eqs if e_ != 0]
    return w0s, eqs

if __name__ == '__main__':
    w0s, eqs = centre_solution()
    for e_ in eqs:
        print('centre eq:', e_)

def centre_coeffs(w0v):
    w0s, eqs = centre_solution()
    eqv = [e_.subs(w0s, w0v) for e_ in eqs]
    sol = sp.solve(eqv, [p2, q2, r2, s1, s3], dict=True)
    return sol

def shoot(w0v, x0=1e-4, xmax=20.0, dense=False):
    # leading-order regular centre data (eqs printed by __main__): s1=-1/2, q2=w0/3, p2=w0/6
    y0 = [1 + w0v/6*x0**2, 1 + w0v/3*x0**2, w0v, -0.5*x0]
    ev = lambda xx, y: detM(xx, y)
    ev.terminal = True
    ev2 = lambda xx, y: 1e4 - np.abs(rhs(xx, y)).max() if (abs(y[3]) < 1 and y[1] > 0) else -1.0
    ev2.terminal = True
    r = solve_ivp(rhs, [x0, xmax], y0, method='DOP853', rtol=1e-12, atol=1e-14, events=[ev, ev2], dense_output=dense)
    if r.t_events[1].size:
        xe = r.t_events[1][0]; ye = r.y_events[1][0]
        return xe, numer(xe, ye), r
    if r.t_events[0].size:
        xe = r.t_events[0][0]; ye = r.y_events[0][0]
        return xe, numer(xe, ye), r
    return None, None, r
