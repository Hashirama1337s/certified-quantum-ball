"""blind-qb2: linear perturbation operator at rate lam (collocation), quantum source, observable."""
import sys, pickle, numpy as np, sympy as sp, time
sys.path.insert(0, '<LAB>/quantum/seats/blind-qb2')
from qb2_coll import *

_ns = dict(sp.__dict__)
with open(HERE + 'qb2_bach.pkl', 'rb') as f:
    BX = {k: eval(s, _ns) for k, s in pickle.load(f).items()}

BGS = [x, al, al1, al2, A, A1, A2, w, w1, w2, v, v1, v2]
BGS_B = [x, al, S('al1'), S('al2'), S('al3'), S('al4'), A, S('A1'), S('A2'), S('A3'), S('A4')]
PS = ['da', 'da1', 'da2', 'dc', 'dc1', 'dc2', 'dw', 'dw1', 'dw2', 'dQ', 'dQ1', 'dQ2']
PSY = [S(p) for p in PS]
EQS = ['dE_uu', 'dE_ux', 'dE_xx', 'dE_tt', 'ddiv_0', 'ddiv_1']
COMP = ['uu', 'ux', 'xx', 'tt']

t0 = time.time()
COEF = {}
for e in EQS:
    ex = EX[e]
    # check linearity in perturbation symbols
    COEF[e] = {}
    for p in PSY:
        if p in ex.free_symbols:
            c = sp.diff(ex, p)
            COEF[e][p.name] = sp.lambdify(BGS + [lam], c, 'numpy', cse=True)
    rest = ex.subs({p: 0 for p in PSY})
    assert rest == 0 or sp.simplify(rest) == 0, (e, 'inhomogeneous part')
F_V = {k: sp.lambdify(BGS, EX['Vhat_' + k], 'numpy', cse=True) for k in COMP}
F_B = {k: sp.lambdify(BGS_B, BX['Bhat_' + k], 'numpy', cse=True) for k in COMP}
F_Vm = {k: sp.lambdify(BGS_B, BX['Vm_' + k], 'numpy', cse=True) for k in COMP}
F_Rs = sp.lambdify(BGS_B, BX['Rs'], 'numpy', cse=True)
dom = EX['d_one_minus_2m_r']
F_dom = {p.name: sp.lambdify([x, al, A], sp.diff(dom, p), 'numpy') for p in PSY if p in dom.free_symbols}
F_dP = None
print('lambdified in %.1fs' % (time.time() - t0), flush=True)

def bg_at_nodes(g, cvec, order=2):
    a_, A_, w_, v_ = bg_unpack(g, cvec, order)
    return a_, A_, w_, v_

def bgvals(g, cvec):
    a_, A_, w_, v_ = bg_unpack(g, cvec, 2)
    return [g.x, a_[0], a_[1], a_[2], A_[0], A_[1], A_[2], w_[0], w_[1], w_[2], v_[0], v_[1], v_[2]]

def bgvals_B(g, cvec):
    a_, A_, w_, v_ = bg_unpack(g, cvec, 4)
    return [g.x] + a_ + A_

def pert_basis(g):
    """dict symbol -> (N x N) matrix mapping coefficients of the corresponding unknown to values."""
    x_ = g.x[:, None]
    T0, T1, T2 = g.T0, g.T1, g.T2
    ev = [T0, 2*x_*T1, 2*T1 + 4*x_**2*T2]
    od = [x_*T0, T0 + 2*x_**2*T1, 6*x_*T1 + 4*x_**3*T2]
    return {'da': (0, ev[0]), 'da1': (0, ev[1]), 'da2': (0, ev[2]),
            'dc': (1, ev[0]), 'dc1': (1, ev[1]), 'dc2': (1, ev[2]),
            'dw': (2, ev[0]), 'dw1': (2, ev[1]), 'dw2': (2, ev[2]),
            'dQ': (3, od[0]), 'dQ1': (3, od[1]), 'dQ2': (3, od[2])}

def operator(g, cvec, lamv):
    N = g.N
    bv = bgvals(g, cvec)
    PB = pert_basis(g)
    L = np.zeros((6*N, 4*N))
    for i, e in enumerate(EQS):
        for pn, f in COEF[e].items():
            coef = np.broadcast_to(np.asarray(f(*bv, lamv), dtype=float), (N,))
            blk, Bm = PB[pn]
            L[i*N:(i+1)*N, blk*N:(blk+1)*N] += coef[:, None]*Bm
    return L

def source(g, cvec, which='B'):
    """RHS for the u e^{2u} level per unit K eps: dE_ab = K Bhat_ab/(30 pi); fluid rows 0."""
    N = g.N
    rhs_ = np.zeros(6*N)
    if which == 'B':
        bb = bgvals_B(g, cvec)
        for i, k in enumerate(COMP):
            rhs_[i*N:(i+1)*N] = np.asarray(F_B[k](*bb), dtype=float)/(30*np.pi)
    else:   # Zahn V (on-shell fluid form) normalised as (8 pi) u V/(4 pi^2) = 2 V/pi
        bv = bgvals(g, cvec)
        for i, k in enumerate(COMP):
            rhs_[i*N:(i+1)*N] = np.asarray(F_V[k](*bv), dtype=float)*2/np.pi
    return rhs_

def gauge_row(g, which='da0'):
    N = g.N
    r = np.zeros(4*N)
    if which == 'da0':
        r[:N] = g.e0
    return r

def solve_forced(g, cvec, lamv, rhs_, gauge_val=0.0, weights=None):
    L = operator(g, cvec, lamv)
    G = gauge_row(g)[None, :]
    Afull = np.vstack([L, G])
    bfull = np.concatenate([rhs_, [gauge_val]])
    if weights is not None:
        Afull = Afull*weights[:, None]; bfull = bfull*weights
    sol, *_ = np.linalg.lstsq(Afull, bfull, rcond=None)
    resid = Afull @ sol - bfull
    return sol, resid, L

def pert_eval(g, sol, xx):
    N = g.N
    B = g.basis_at(np.atleast_1d(xx))
    return (B @ sol[:N], B @ sol[N:2*N], B @ sol[2*N:3*N], np.atleast_1d(xx)*(B @ sol[3*N:]))

def Psi_at(g, cvec, sol, xx):
    """delta(2m/r)/(1-2m/r) = -A * delta(1 - 2m/r)  (1 - 2m/r = 1/A)."""
    a_, A_, w_, v_ = bg_eval(g, cvec, np.atleast_1d(xx))
    da, dc, dw, dQ = pert_eval(g, sol, xx)
    vals = {'da': da, 'dc': dc}
    d = sum(F_dom[p](np.atleast_1d(xx), a_, A_)*vals[p] for p in F_dom)
    return -A_*d
