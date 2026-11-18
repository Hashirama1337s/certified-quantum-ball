"""build-5 PQ2c-2: linearised E_thth residual (build-4 thth_exprs.txt, exact) on the REAL source and the N-b variant
(b_T x (1 + 1e-6)) at x = -6, -4.6, -2, -0.25, -0.01, on the QB2b run-1 record (build-4 final.json, pass_ckpt.json,
stages.json; read only).
Y = full real-source solution at the point:
   x = -6     centre side   R_c(theta6) - c rho1(theta6) - d rho2(theta6)          (stages.json x6, final matching)
   x = -4.6   sonic side    R(x_m) + b s^HKA(x_m)                                   (pass_ckpt vectors)
   x = -2, -0.25            R + b s^HKA at the pass check points                     (pass_ckpt checks)
   x = -0.01                R + b s^HKA on the forced sonic series (sonic stage re-run here; s^HKA from stages.json)
Y', Y'' from the ODE at the point:  Y' = M Y + F,  Y'' = M_x Y + M Y' + F_x,  with M = [Ptop; (F0 + 2 F1)/detMx] from the
engine HKA step tape (lin_fields2.pkl, LF.step) on the local x-series of the background (ODE jets), F the real-source
column (build-3 qb2core.generic_columns) on the same local series.  Background inputs of R_thth from the same series.
usage (memcap 8): py -3.11 thth_c.py -> thth_c.json"""
import os, sys, json, pickle
os.environ['QB2B_RUN1'] = '1'
sys.path.insert(0, '<LAB>/quantum/seats/build-4')
import qb2b as B
from qb2b import ST, C, Q, sym, arb, acb, arb_mat, ctx
from qb2core import NA, IP, des
import sympy as sp
from core import Tape
import stages_b_helpers as H
B4 = B.B4; B5 = '<LAB>/quantum/seats/build-5/'
de = lambda p: arb(p[0]) + sym(arb(p[1]))
RTH = [eval(l, sp.__dict__.copy()) for l in open(B4 + 'thth_exprs.txt')]
names = ['Ab', 'Abx', 'Abxx', 'Nb', 'Nbx', 'Nbxx', 'ob', 'a', 'a_x', 'a_xx', 'n', 'n_x', 'n_xx', 'w', 'sg_TT']
T_R = Tape([RTH[0]], [sp.Symbol(s_) for s_ in names], {})
S = lambda a_: [a_.mid().str(12, radius=False), a_.rad().str(3, radius=False)]


def at_point(L, state, x, Y):
    """residuals (real, N-b) at a point from the state ball, x, and the 4-vector Y"""
    n = 8
    loc = Q.ode_xseries(state, n)
    one = [arb(1)] + [arb(0)]*(n - 1)
    LF = L.LF
    dom = Q.s_deriv(loc[1]); dV = Q.s_deriv(loc[2])
    outs = LF.series_out(LF.step, [loc[0], loc[1], loc[2], dom, dV, one], 2)
    Ptop, F0, F1, dM = outs[0:8], outs[8:16], outs[16:24], outs[24]
    M = [[None]*4 for _ in range(4)]
    for i in range(2):
        for j in range(4):
            M[i][j] = Ptop[4*i + j][:2]
    for i in range(2):
        for j in range(4):
            s_ = [F0[4*i + j][m] + 2*F1[4*i + j][m] for m in range(3)]
            M[2 + i][j] = Q.s_div(s_, dM[:3], 2)
    ex = Q.s_exp([x, arb(1)] + [arb(0)]*(n - 2), n)
    cols = C.generic_columns(loc[0], loc[1], loc[2], one, ex, n)
    Ex = C.x_form_matrix_entries(cols, cols['n'], lambda num, d: Q.s_div(num, d, cols['n']))
    F = [Ex[(i, 6)] for i in range(4)]
    Y1 = [sum((M[i][j][0]*Y[j] for j in range(4)), arb(0)) + F[i][0] for i in range(4)]
    Y2 = [sum((M[i][j][1]*Y[j] + M[i][j][0]*Y1[j] for j in range(4)), arb(0)) + F[i][1] for i in range(4)]
    A = Q.series_out(Q.T_ACONS, loc, 3, arb(0), Q.conv_for(arb(0)))[0]
    bg = [A[0], A[1], 2*A[2], loc[0][0], loc[0][1], 2*loc[0][2], loc[1][0]]
    sTT = (-2*x).exp()*cols['b'][3][0]
    base = bg + [Y[0], Y1[0], Y2[0], Y[1], Y1[1], Y2[1], Y[2]]
    r = T_R.eval(base + [sTT])[0]
    rnb = T_R.eval(base + [sTT*(1 + arb('1e-6'))])[0]
    return dict(real=S(r), real_contains_0=bool(r.contains(0)), N_b=S(rnb), N_b_excludes_0=bool(not rnb.contains(0)),
                sigma_TT=S(sTT), predicted_N_b=S(-arb('1e-6')*sTT))


old = ctx.prec
out = {}
try:
    L = B.engine_b(log=lambda *a: None)
    tb = pickle.load(open(B.NEW_TUBE, 'rb'))
    F = json.load(open(B4 + 'final.json')); stg = json.load(open(B4 + 'stages.json'))
    import qb2pass as PS
    ck = json.load(open(B4 + 'pass_ckpt.json'))
    mR = [arb(v[0]) + sym(arb(v[1])) for v in (F['matching']['R'][k] for k in ('b', 'c', 'd', 'e'))]
    # x = -0.01 : sonic series (re-run of the sonic stage, R stashed at theta(-0.01))
    th01 = H.theta_sonic(tb, '-0.01')
    B._STASH['extra_t'] = th01
    ST.sonic_stage(L, tb)
    lam = B._STASH['lams'][-1]; Rv = [B._STASH['vals'][0][0][i]*lam[i] for i in range(NA)]
    B._STASH['extra_t'] = None
    hka01 = [des(t) for t in stg['C1_x001']['hka']]
    st01 = [des(t) for t in stg['C1_x001']['state']]
    Y = [Rv[i] + mR[0]*hka01[i] for i in range(4)]
    out['x=-0.01'] = at_point(L, st01, arb(-1)/100, Y)
    # x = -0.25, -2 : pass check points
    for nm, key in (('x=-0.25', 'm025'), ('x=-2', 'm2')):
        rec = ck['extra']['checks'][key]
        Y = [des(rec['R'][i]) + mR[0]*des(rec['hka'][i]) for i in range(4)]
        out[nm] = at_point(L, [des(t) for t in rec['state']], des(rec['x']), Y)
    # x = -4.6 : end of the transport
    v = {nm: PS.Lohner.load(ck['vec'][nm]).value() for nm in ('R', 'hka')}
    Y = [v['R'][i] + mR[0]*v['hka'][i] for i in range(4)]
    out['x=-4.6'] = at_point(L, [des(t) for t in ck['extra']['end']['state']], des(ck['extra']['end']['x']), Y)
    # x = -6 : centre side
    cen = stg['centre']
    r6 = [des(t) for t in cen['R']['x6']['val']]; p1 = [des(t) for t in cen['rho1']['x6']]; p2 = [des(t) for t in cen['rho2']['x6']]
    Y = [r6[i] - mR[1]*p1[i] - mR[2]*p2[i] for i in range(4)]
    cm = tb['cm']; th6 = ST.theta_at(tb, -6)
    P6, W6, U6 = [Q.horner([de(p) for p in cm['Y'][i]], th6) + sym(de(cm['eps'])) for i in range(3)]
    out['x=-6'] = at_point(L, [1/P6, W6*P6**2, U6*P6], arb(-6), Y)
    # x = -4.6 again, on the CENTRE-side representation (same solution; e encloses 0): small, sharp values
    rm = [des(t) for t in cen['R']['m']['val']]; q1 = [des(t) for t in cen['rho1']['m']]; q2 = [des(t) for t in cen['rho2']['m']]
    Y = [rm[i] - mR[1]*q1[i] - mR[2]*q2[i] for i in range(4)]
    thm = L.theta_m
    Pm, Wm, Um = [Q.horner([de(p) for p in cm['Y'][i]], thm) + sym(de(cm['eps'])) for i in range(3)]
    out['x=-4.6_centre_side'] = at_point(L, [1/Pm, Wm*Pm**2, Um*Pm], L.xbar, Y)
    out['x=-4.6']['note'] = 'sonic-side transported representation (|Y| ~ 600): width too large for N-b; see x=-4.6_centre_side'
    out['all_real_contain_0'] = all(v_['real_contains_0'] for k_, v_ in out.items() if k_.startswith('x='))
    out['all_N_b_exclude_0'] = all(v_['N_b_excludes_0'] for k_, v_ in out.items() if k_.startswith('x=') and k_ != 'x=-4.6')
    out['N_b_points_used'] = 'x=-6, x=-4.6_centre_side, x=-2, x=-0.25, x=-0.01'
finally:
    ctx.prec = old
json.dump(out, open(B5 + 'thth_c.json', 'w'), indent=1)
print(json.dumps(out, indent=1))
