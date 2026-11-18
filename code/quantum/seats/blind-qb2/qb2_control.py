"""Manufactured-source positive control: pick analytic regular (da, dc, dw, dQ) with da(0)=0,
build the source by evaluating the symbolic operator on their EXACT derivatives, solve with the
collocation solver, compare recovered fields and Psi with the known ones."""
import sys, json, pickle; sys.path.insert(0, '<LAB>/quantum/seats/blind-qb2')
from qb2_pert import *
from scipy.optimize import brentq
X = float(sys.argv[2]); Ns = [int(a) for a in sys.argv[1].split(',')]
xs_ = sp.Symbol('xs')
man = {'da': 0.7*xs_**2*sp.exp(-xs_**2), 'dc': 0.4*xs_**2/(1 + xs_**2), 'dw': 2 + sp.cos(1.3*xs_**2),
       'dQ': xs_*(0.5 - 0.2*xs_**2)}
F_man = {}
for k, e in man.items():
    for o in range(3):
        F_man[k + ('' if o == 0 else str(o))] = sp.lambdify(xs_, sp.diff(e, xs_, o), 'numpy')
out = {}
for N in Ns:
    d = pickle.load(open('bg_N%d_X%.3f.pkl' % (N, X), 'rb'))
    g = Grid(N, X); cvec = np.array(d['cvec'])
    bv = bgvals(g, cvec)
    rhs_ = np.zeros(6*N)
    for i, e in enumerate(EQS):
        for pn, f in COEF[e].items():
            rhs_[i*N:(i+1)*N] += np.asarray(f(*bv, 2.0), float)*F_man[pn](g.x)
    sol, resid, L = solve_forced(g, cvec, 2.0, rhs_, 0.0)
    xt = np.array([0.2, 0.5, 0.628039165607, 0.85])
    got = pert_eval(g, sol, xt)
    err = [float(np.abs(got[i] - F_man[k](xt)).max()) for i, k in enumerate(['da', 'dc', 'dw', 'dQ'])]
    xp = 0.6280391656072097
    a_, A_, w_, v_ = bg_eval(g, cvec, np.array([xp]))
    psi_exact = float(-A_[0]*sum(F_dom[p](np.array([xp]), a_, A_)[0]*F_man[p](xp) for p in F_dom))
    psi_got = float(Psi_at(g, cvec, sol, xp)[0])
    out[N] = dict(field_err=err, psi_exact=psi_exact, psi_solver=psi_got, resid=float(np.abs(resid).max()))
    print(N, out[N], flush=True)
json.dump(out, open('control_manufactured_X%.3f.json' % X, 'w'), indent=1)
