import sys, json, pickle; sys.path.insert(0, '<LAB>/quantum/seats/blind-qb2')
from qb2_pert import *
from scipy.optimize import minimize_scalar, brentq
X = float(sys.argv[2]); Ns = [int(a) for a in sys.argv[1].split(',')]
mode = sys.argv[3] if len(sys.argv) > 3 else 'all'
out = {}
for N in Ns:
    d = pickle.load(open('bg_N%d_X%.3f.pkl' % (N, X), 'rb'))
    g = Grid(N, X); cvec = np.array(d['cvec'])
    xp = brentq(lambda t: bg_eval(g, cvec, np.array([t]))[3][0], 0.5, 0.7, xtol=1e-15)
    Ap = bg_eval(g, cvec, np.array([xp]))[1][0]
    r = {'N': N, 'X': X, 'x_peak': xp, 'A_peak': Ap}
    # --- source checks
    bb = bgvals_B(g, cvec); bv = bgvals(g, cvec)
    ratio = []; vmdiff = []
    for k in COMP:
        Bv = np.asarray(F_B[k](*bb), float); Vv = np.asarray(F_V[k](*bv), float); Vm = np.asarray(F_Vm[k](*bb), float)
        sel = np.abs(Vv) > 1e-3*np.abs(Vv).max()
        ratio.append([float(np.min(Bv[sel]/Vv[sel])), float(np.max(Bv[sel]/Vv[sel]))])
        vmdiff.append(float(np.abs(Vm - Vv).max()/np.abs(Vv).max()))
    r['B_over_V_minmax'] = ratio; r['Vmetric_vs_Vfluid_rel'] = vmdiff
    r['Rs_max'] = float(np.abs(np.asarray(F_Rs(*bb), float)).max())
    print('source checks', ratio, vmdiff, r['Rs_max'], flush=True)
    G = gauge_row(g)[None, :]
    def smin(lv):
        L = operator(g, cvec, lv)
        sv = np.linalg.svd(np.vstack([L, G]), compute_uv=False)
        return sv[-1]/sv[0]
    if mode in ('all', 'kappa'):
        lams = np.linspace(2.5, 3.1, 13)
        sm = [smin(l) for l in lams]
        print('scan', list(zip(lams.round(3), np.array(sm))), flush=True)
        i = int(np.argmin(sm))
        mres = minimize_scalar(smin, bracket=(lams[max(i-1,0)], lams[i], lams[min(i+1,len(lams)-1)]), tol=1e-12)
        r['kappa'] = float(mres.x); r['kappa_smin'] = float(mres.fun); r['smin_at_2'] = float(smin(2.0))
        print('kappa', mres.x, mres.fun, 'smin(2)', r['smin_at_2'], flush=True)
    if mode in ('all', 'forced'):
        for which in ['B', 'V']:
            rhs_ = source(g, cvec, which)
            sol, resid, L = solve_forced(g, cvec, 2.0, rhs_, 0.0)
            psi = Psi_at(g, cvec, sol, xp)[0]
            sol1, resid1, _ = solve_forced(g, cvec, 2.0, rhs_, 1.0)     # add gauge
            psi1 = Psi_at(g, cvec, sol1, xp)[0]
            da, dc, dw, dQ = pert_eval(g, sol, np.array([0.0, xp]))
            xs_test = np.array([0.3, xp, 0.75])
            r['forced_' + which] = dict(Psi=float(psi), Psi_gauge1=float(psi1), resid_max=float(np.abs(resid).max()),
                                       rhs_max=float(np.abs(rhs_).max()), dc0=float(dc[0]),
                                       Psi_profile_g0=[float(t) for t in Psi_at(g, cvec, sol, xs_test)],
                                       Psi_profile_g1=[float(t) for t in Psi_at(g, cvec, sol1, xs_test)])
            print(which, r['forced_' + which], flush=True)
            if which == 'B':
                np.save('sol_N%d_X%.3f.npy' % (N, X), sol)
    out[N] = r
json.dump(out, open('main_%s_X%.3f_%s.json' % (mode, X, sys.argv[1].replace(',', '-')), 'w'), indent=1, default=float)
