"""blind-qb2 common: load symbolic expressions, background ODE."""
import pickle, numpy as np, sympy as sp

HERE = '<LAB>/quantum/seats/blind-qb2/'
with open(HERE + 'qb2_exprs.pkl', 'rb') as f:
    _raw = pickle.load(f)
_ns = dict(sp.__dict__)
EX = {k: eval(s, _ns) for k, s in _raw.items()}

S = sp.Symbol
x = S('x', real=True)
lam = S('lam')
al, A, w, v = S('al'), S('A'), S('w'), S('v')
al1, A1, w1, v1 = S('al1'), S('A1'), S('w1'), S('v1')
al2, A2, w2, v2 = S('al2'), S('A2'), S('w2'), S('v2')

def bg_system():
    """(al1, A1) solved from E_xx, E_ux; fluid system M (w1, v1) = b."""
    sol_al1 = sp.solve(EX['Ebg_xx'], al1)[0]
    sol_A1 = sp.solve(EX['Ebg_ux'].subs(al1, sol_al1), A1)[0]
    d = [EX['divbg_0'].subs({al1: sol_al1, A1: sol_A1}), EX['divbg_1'].subs({al1: sol_al1, A1: sol_A1})]
    M = sp.Matrix(2, 2, lambda i, j: sp.diff(d[i], [w1, v1][j]))
    b = sp.Matrix([-(d[i].subs({w1: 0, v1: 0})) for i in range(2)])
    return sol_al1, sol_A1, M, b
