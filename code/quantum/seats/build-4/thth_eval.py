"""build-4 C8 theta-theta residual (and N-b), certified, evaluated DIRECTLY on the computed solutions:
   R_thth(Y, Y', Y'', sigma_TT; background)  (thth_sym.py, exact linearisation of E_thth = G - 8 pi T - Sigma)
at x = -0.01 (forced sonic series) and x = -6 (forced centre series), for
   * C8 (sigma_TT from D_x = 0, thth_sym.py): must enclose 0;
   * the real source (sigma_TT = beta_T = e^{-2x} b_T): must enclose 0;
   * N-b: real source with b_T x (1 + 1e-6): must exclude 0.
Y = particular + homogeneous parts with the matching coefficients of final.json (b on the sonic side; c, d on the centre
side; Phi0 for C8).  Derivatives of every series with certified tails (coefficient majorant |y_n| <= B/nu'^n, the
bound behind the RegSing tail; nu' optimised).  usage (memcap 8): py -3.11 thth_eval.py -> thth.json"""
import sys, os, json, time, pickle, math
sys.path.insert(0, '<LAB>/quantum/seats/build-4')
import qb2b as B
from qb2b import ST, C, Q, sym, arb, acb, arb_mat, ctx
from qb2core import NA, IP, ZC, ser, des
import sympy as sp
from core import Tape
T0 = time.time()
de = lambda p: arb(p[0]) + sym(arb(p[1]))
SER = []
_orig = C.RegSingAug.evaluate


def _grab(self, y, t):
    SER.append((self, y))
    return _orig(self, y, t)


C.RegSingAug.evaluate = _grab


def derivs(y, t, S, nu, cn, K=None):
    """value, d/dt, d2/dt2 of sum y_n t^n (y_n column matrices, real) with certified tails; nu' optimised"""
    K = len(y) - 1 if K is None else K
    n_ = y[0].nrows()
    at = abs(t).upper() if not isinstance(t, arb) else t.abs_upper()
    best = None
    for f in [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
        nup = at + (nu - at)*arb(f)
        rat = nup/nu; mup = S*rat/(1 - rat)
        n0 = K; ok = True
        while True:
            c = cn(n0 + 1)
            if c is not None and c*mup <= 1: break
            n0 += 1
            if n0 > K: ok = False; break
        if not ok: continue
        Bm = arb(0)
        for n in range(K + 1):
            Bm = Bm.max(max(abs(y[n][i, 0]).upper() for i in range(n_))*nup**n)
        q = at/nup
        t0_ = Bm*q**(K + 1)/(1 - q)
        t1_ = Bm/nup*(K + 1)*q**K/(1 - q)**2
        t2_ = Bm/nup**2*(K + 1)**2*q**(K - 1)/(1 - q)**3
        if best is None or t2_ < best[2]:
            best = (t0_, t1_, t2_)
    assert best is not None
    out = []
    for i in range(n_):
        v0 = v1 = v2 = 0*t
        for n in range(K, -1, -1):
            v0 = v0*t + y[n][i, 0]
        for n in range(K, 0, -1):
            v1 = v1*t + n*y[n][i, 0]
        for n in range(K, 1, -1):
            v2 = v2*t + n*(n - 1)*y[n][i, 0]
        out.append((v0 + sym(best[0]), v1 + sym(best[1]), v2 + sym(best[2])))
    return out


def eng_series(rs, kap, y0):
    ys, M = rs.series(kap, y0)
    yr = [arb_mat([[ys[n][i, 0].real] for i in range(4)]) for n in range(len(ys))]
    M0n = max(sum((abs(M[0][i, j]).upper() for j in range(4)), arb(0)) for i in range(4))
    S = rs.supC0 + arb(abs(kap).upper())*rs.supC1
    if rs.rank1:
        trM = sum((M[0][i, i] for i in range(4)), acb(0))
        cn = lambda n: ((1 + M0n/abs(acb(n) - trM).lower())/n) if abs(acb(n) - trM).lower() > 0 else None
    else:
        cn = lambda n: (1/(arb(n) - M0n)) if arb(n) > M0n else None
    return yr, S, rs.nu, cn


RTH = [eval(l, sp.__dict__.copy()) for l in open(B.B4 + 'thth_exprs.txt')]
names = ['Ab', 'Abx', 'Abxx', 'Nb', 'Nbx', 'Nbxx', 'ob', 'a', 'a_x', 'a_xx', 'n', 'n_x', 'n_xx', 'w', 'sg_TT']
T_R = Tape([RTH[0]], [sp.Symbol(s_) for s_ in names], {})
TT_names = ['Ab', 'Abx', 'Nb', 'Nbx', 'sg_ss', 'sg_sx', 'sg_xx', 'sgx_ss', 'sgx_sx', 'sgx_xx']
T_TT = Tape([RTH[1]], [sp.Symbol(s_) for s_ in TT_names], {})


def background(state, x):
    """local x-series (length 7) of N, om, V at the state; A series; beta profiles; b (Bach) value"""
    loc = Q.ode_xseries(state, 8)
    A = Q.series_out(Q.T_ACONS, loc, 8, arb(0), Q.conv_for(arb(0)))[0]
    ex = Q.s_exp([x, arb(1)] + [arb(0)]*6, 8)
    bt = Q.bach_from_series(loc[0], loc[1], loc[2], [arb(1)] + [arb(0)]*7, 8)
    return loc, A, ex, bt


def residual(Yv, loc, A, ex, sTT):
    """Yv: list of (val, d/dx, d2/dx2) for a, n, w"""
    bg = [A[0], A[1], 2*A[2], loc[0][0], loc[0][1], 2*loc[0][2], loc[1][0]]
    ins = bg + [Yv[0][0], Yv[0][1], Yv[0][2], Yv[1][0], Yv[1][1], Yv[1][2], Yv[2][0], sTT]
    return T_R.eval(ins)[0]


def sigma_series(A, N, bnn, bnr, brr, m):
    sss = Q.s_mul(A, Q.s_add(Q.s_sub(Q.s_mul(Q.s_mul(N, N, m), bnn, m), Q.s_scal(Q.s_mul(N, bnr, m), 2)), brr), m)
    ssx = Q.s_mul(A, Q.s_sub(Q.s_mul(N, bnr, m), brr), m)
    sxx = Q.s_mul(A, brr, m)
    return sss, ssx, sxx


old = ctx.prec
out = {}
try:
    L = B.engine_b(log=lambda *a: None)
    tb = pickle.load(open(B.NEW_TUBE, 'rb'))
    DATA = os.environ.get('QB2B_DATA', B.B4)
    F = json.load(open(DATA + 'final.json'))
    stg = json.load(open(DATA + 'stages.json'))
    mt = F['matching']
    bR = arb(mt['R']['b'][0]) + sym(arb(mt['R']['b'][1])); b8 = arb(mt['C8']['b'][0]) + sym(arb(mt['C8']['b'][1]))
    Phi0 = arb(F['C8_Phi0'][0]) + sym(arb(F['C8_Phi0'][1]))
    # sonic: re-run the stage to grab the series objects (the same computation as stages_b.py)
    import stages_b_helpers as H
    th01 = H.theta_sonic(tb, '-0.01')
    SER.clear()
    sres, sinfo = ST.sonic_stage(L, tb)
    lam = B._STASH['lams'][-1]
    order = ['R', 'N0', 'C', 'C1', 'C8', 'E', 'Ph']
    sm = tb['sm']; X = [de(p) for p in sm['X']]
    Xv = H.horner_d(X, th01, de(sm['epsX']), arb(tb['nu_s']))
    state = [Q.horner([de(p) for p in sm['Y'][i]], th01) + sym(de(sm['eps'])) for i in range(3)]
    loc, A, ex, bt = background(state, arb(-1)/100)
    xp, xpp = Xv[1], Xv[2]

    def to_x(d):     # (v, d/dtheta, d2/dtheta2) -> (v, d/dx, d2/dx2)
        v, d1, d2 = d
        y1 = d1/xp
        return (v, y1, (d2 - y1*xpp)/xp**2)
    sol = {}
    for nm, (rs, y) in zip(order, SER[:7]):
        dv = derivs(y, th01, rs.S, rs.nu, rs.cn)
        sol[nm] = [to_x((dv[i][0]*lam[i], dv[i][1]*lam[i], dv[i][2]*lam[i])) for i in range(NA)]
    kap = acb(2)
    M0 = [[acb(L.sonic.M0[0][i][j]) + kap*acb(L.sonic.M1[0][i][j]) for j in range(4)] for i in range(4)]
    w = max((M0[2], M0[3]), key=lambda r: float(sum(abs(x_).mid() for x_ in r)))
    I0 = [acb(L.I0[0]) + kap] + [acb(x_) for x_ in L.I0[1:]]
    from flint import acb_mat
    y0 = acb_mat([w, I0, [1, 0, 0, 0], [0, 1, 0, 0]]).solve(acb_mat([[0], [0], [1], [0]]))
    yr, S_, nu_, cn_ = eng_series(L.sonic, kap, [y0[i, 0] for i in range(4)])
    dh = derivs(yr, th01, S_, nu_, cn_)
    hka = [to_x(dh[i]) for i in range(4)]

    def comb(parts):
        return [tuple(sum((c*p[i][k] for c, p in parts), arb(0)) for k in range(3)) for i in range(4)]
    yR = comb([(arb(1), sol['R']), (bR, hka)])
    y8 = comb([(arb(1), sol['C8']), (Phi0, sol['Ph']), (b8, hka)])
    Phi8 = tuple(sol['C8'][IP][k] + Phi0*sol['Ph'][IP][k] for k in range(3))
    xv = arb(-1)/100
    e2x = (-2*xv).exp()
    # real source: sigma_TT = beta_T
    sTT_R = e2x*bt['b'][3][0]
    rR = residual(yR, loc, A, ex, sTT_R)
    rNb = residual(yR, loc, A, ex, sTT_R*(1 + arb('1e-6')))
    # C8: beta_nn = A^{-1/2} Phi (series in dx from (Phi, Phi', Phi'')), beta_nr, beta_rr explicit
    m = 3
    Phis = [Phi8[0], Phi8[1], Phi8[2]/2]
    irA = Q.s_div([arb(1), arb(0), arb(0)], Q.s_sqrt(A[:m], m), m)
    bnn = Q.s_mul(irA, Phis, m)
    bnr, brr = Q.c8_profiles(ex[:m + 1], m)
    sss, ssx, sxx = sigma_series(A[:m], loc[0][:m], bnn, bnr, brr, m)
    sTT_8 = T_TT.eval([A[0], A[1], loc[0][0], loc[0][1], sss[0], ssx[0], sxx[0], sss[1], ssx[1], sxx[1]])[0]
    r8 = residual(y8, loc, A, ex, sTT_8)
    S = lambda a_: [a_.mid().str(12, radius=False), a_.rad().str(3, radius=False)]
    out['x=-0.01'] = dict(C8=S(r8), C8_contains_0=bool(r8.contains(0)), real=S(rR), real_contains_0=bool(rR.contains(0)),
                          N_b=S(rNb), N_b_excludes_0=bool(not rNb.contains(0)), sigma_TT_C8=S(sTT_8), sigma_TT_real=S(sTT_R))
    print(json.dumps(out, indent=1), flush=True)
finally:
    ctx.prec = old
json.dump(out, open(os.environ.get('QB2B_DATA', B.B4) + 'thth.json', 'w'), indent=1)
