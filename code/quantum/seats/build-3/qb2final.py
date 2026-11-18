"""build-3 QB2 final: matching (N6), controls C1 C2 C6 C8 C8', nulls N-a N-c, then (only if N1, C1, C8 pass)
N7 + N8 Psi2.  Inputs: stages.json, pass_ckpt.json (k = 1440), n1_rows.log, n0check.json.  Outputs: final.json
(everything), then make_outputs.py writes result_qb2.json / controls_qb2.json.
usage (memcap 8): py -3.11 qb2final.py
"""
import sys, json, time, math, pickle
from flint import ctx as _ctx
_ctx.cap = 32
sys.path.insert(0, '<LAB>/engine/cert'); sys.path.insert(0, '<LAB>/engine')
sys.path.insert(0, '<LAB>/quantum/seats/build-3')
from flint import arb, acb, arb_mat, ctx
from core import sym
import qb2stages as ST
import qb2lib as Q
import qb2core as C
from qb2core import NA, IJ, IP, ZC, ser, des
HERE = '<LAB>/quantum/seats/build-3/'
DELTA = arb('1e-6')
T0 = time.time()
S = lambda a: [a.mid().str(20, radius=False), a.rad().str(3, radius=False)]


def vals(lst):
    return [des(t) for t in lst]


def main():
    old = ctx.prec
    out = {}
    try:
        L = ST.engine(log=lambda *a: None)
        ctx.prec = 128
        tb = pickle.load(open('<LAB>/engine/cert/lin_tube2_k1_3.pkl', 'rb'))
        stg = json.load(open(HERE + 'stages.json'))
        ck = json.load(open(HERE + 'pass_ckpt.json'))
        assert ck['k'] == len(L.steps), 'pass not complete'
        import qb2pass as PS
        vec = {nm: PS.Lohner.load(ck['vec'][nm]).value() for nm in ck['vec']}
        ex = ck['extra']
        xe = des(ex['end']['x']); st_e = vals(ex['end']['state'])
        out['x_m'] = S(xe); out['x_m_overlaps_XBAR'] = bool(xe.overlaps(L.xbar))
        # identity row at x_m (engine Tape, kappa-free part + kappa)
        I0 = L.LF.irow.eval(st_e, conv=Q.sym_const)
        Irow = [I0[0] + 2] + list(I0[1:])
        cV = Irow[3]
        out['c_V_xm'] = S(cV); out['c_V_overlaps_verdict_M'] = bool(cV.overlaps(arb('0.19288155165634') + sym(arb('1.41e-15'))))
        cen = stg['centre']
        rho1 = vals(cen['rho1']['m']); rho2 = vals(cen['rho2']['m'])
        hka = vec['hka'][:4]
        M = arb_mat([[hka[i], rho1[i], rho2[i], arb(1) if i == 3 else arb(0)] for i in range(4)])
        detM = M.det()
        out['N6_matrix_det_D_HKA2'] = S(detM)
        out['N6_det_overlaps_pq7b'] = bool(detM.overlaps(arb('-19.942295267987575083') + sym(arb('1.12e-8'))))
        cs = {nm: vals(cen[nm]['m']['val']) for nm in ('R', 'N0', 'C', 'C1', 'C8', 'E')}
        ys = {nm: vec[nm] for nm in ('R', 'N0', 'C', 'C1', 'C8', 'E', 'Ph')}
        # C8: Phi continuity at x_m fixes the sonic homogeneous amount Phi0
        Phi0 = (cs['C8'][IP] - ys['C8'][IP])/ys['Ph'][IP]
        y8 = [ys['C8'][i] + Phi0*ys['Ph'][i] for i in range(NA)]
        out['C8_Phi0'] = S(Phi0)
        match = {}
        for nm, ysn in (('R', ys['R']), ('N0', ys['N0']), ('C', ys['C']), ('C1', ys['C1']), ('C8', y8), ('E', ys['E'])):
            rhs = arb_mat([[cs[nm][i] - ysn[i]] for i in range(4)])
            sol = M.solve(rhs)
            match[nm] = [sol[i, 0] for i in range(4)]
        mt = {nm: dict(b=S(v[0]), c=S(v[1]), d=S(v[2]), e=S(v[3]), e_width=(2*v[3].rad()).str(3)) for nm, v in match.items()}
        out['matching'] = mt
        eR = match['R'][3]
        out['N6'] = dict(e=S(eR), e_contains_0=bool(eR.contains(0)), e_width=(2*eR.rad()).str(3),
                         width_le_1em10=bool(2*eR.rad() <= arb('1e-10')))
        # ---- nulls
        eN0 = match['N0'][3]
        Js = ys['N0'][IJ]; Jc = cs['N0'][IJ]
        pred_full = (Jc - Js)/cV
        pred_design = Js.abs_upper()/abs(cV)          # DESIGN-QB2-v2 sec. 4 formula (assumes J_c = 0)
        eNa = eR + DELTA*eN0
        out['N-a'] = dict(e_Na=S(eNa), excludes_0=bool(not eNa.contains(0)),
                          e_rowzero_per_unit_delta=S(eN0), J_s_xm=S(Js), J_c_xm=S(Jc),
                          predicted_e_per_unit_delta_full=S(pred_full), predicted_abs_e_design_formula_per_unit_delta=S(abs(Js)/abs(cV)),
                          e_rowzero_encloses_full_prediction=bool(eN0.overlaps(pred_full)),
                          predicted_e_Na=S(DELTA*pred_full),
                          e_Na_encloses_predicted=bool(eNa.overlaps(eR.mid() + DELTA*pred_full)),
                          design_formula_matches=bool(abs(eN0).overlaps(abs(Js)/abs(cV))))
        eC = match['C'][3]
        out['N-c'] = dict(e=S(eC), refuses=bool(not eC.contains(0)), width=(2*eC.rad()).str(3))
        # ---- controls
        e8 = match['C8'][3]; eE = match['E'][3]
        e8p = (1 + DELTA)*e8 - DELTA*eE
        out['C8'] = dict(e=S(e8), e_contains_0=bool(e8.contains(0)), e_width=(2*e8.rad()).str(3),
                         width_le_1em10=bool(2*e8.rad() <= arb('1e-10')), e_explicit_only=S(eE))
        out["C8'"] = dict(e=S(e8p), refuses=bool(not e8p.contains(0)), predicted=S(-DELTA*eE))
        # C1: phi G_2 reproduced (gauge Nbar(0) = 0): b, e contain 0.  Sanity (CORRECTED 2026-11-18 after the first final run:
        # the first version expected c = -1/2, d = -1, a sign slip): phi G_2 -> -(1/2) rho1 - rho2 at the centre means
        # phi G_2 = y_c^p + alpha rho1 + beta rho2 with alpha = -1/2, beta = -1, so y_c^p - y_s^p = -alpha rho1 - beta rho2,
        # i.e. the matching coefficients are c = +1/2, d = +1.  The design criterion (point equality) is unchanged.
        m1 = match['C1']
        c1 = dict(b=S(m1[0]), e=S(m1[3]), b_contains_0=bool(m1[0].contains(0)), e_contains_0=bool(m1[3].contains(0)),
                  c_overlaps_plus_half=bool(m1[1].overlaps(arb(1)/2)), d_overlaps_plus_one=bool(m1[2].overlaps(arb(1))),
                  first_run_sign_slip='expected c=-1/2, d=-1 in the first final run (wrong sign); corrected, see comment')
        pts = {}

        def g2phi(state, x):
            loc = Q.ode_xseries(state, 3)
            Aser = Q.series_out(Q.T_ACONS, loc, 3, arb(0), Q.conv_for(arb(0)))[0]
            G2 = [Aser[1]/Aser[0], loc[0][1]/loc[0][0] + 2, loc[1][1]/loc[1][0], loc[2][1]]
            ph = Q.phi_val(x)
            return [ph*g for g in G2]
        for nm in ('m2', 'm025'):
            rec = ex['checks'][nm]
            x_ = des(rec['x']); stt = vals(rec['state'])
            yv = [des(rec['C1'][i]) + m1[0]*des(rec['hka'][i]) for i in range(4)]
            ex_ = g2phi(stt, x_)
            pts[nm] = dict(x=S(x_), solver=[S(v) for v in yv], exact=[S(v) for v in ex_],
                           overlap=bool(all(yv[i].overlaps(ex_[i]) for i in range(4))))
        yv = [ys['C1'][i] + m1[0]*hka[i] for i in range(4)]
        ex_ = g2phi(st_e, xe)
        pts['m46'] = dict(x=S(xe), solver=[S(v) for v in yv], exact=[S(v) for v in ex_],
                          overlap=bool(all(yv[i].overlaps(ex_[i]) for i in range(4))))
        # x = -6: centre side y_c = y_c^p - c rho1 - d rho2 at theta(-6)
        th6 = arb(stg['theta_x6'].split(' +/- ')[0].strip('[')) if False else None
        r16 = vals(cen['rho1']['x6']); r26 = vals(cen['rho2']['x6']); c16 = vals(cen['C1']['x6']['val'])
        yv = [c16[i] - m1[1]*r16[i] - m1[2]*r26[i] for i in range(4)]
        cm = tb['cm']; de = lambda p: arb(p[0]) + sym(arb(p[1]))
        Z = [[de(p) for p in cm['Y'][i]] for i in range(3)]
        th6 = ST.theta_at(tb, -6)
        P6, W6, U6 = [Q.horner(Z[i], th6) + sym(de(cm['eps'])) for i in range(3)]
        st6 = [1/P6, W6*P6**2, U6*P6]
        ex_ = g2phi(st6, arb(-6))
        pts['m6'] = dict(x='-6', solver=[S(v) for v in yv], exact=[S(v) for v in ex_],
                         overlap=bool(all(yv[i].overlaps(ex_[i]) for i in range(4))))
        c1['points'] = pts
        c1['outward_points_0.5_3'] = 'not computed (domain cut, lead lever 1): no outward transport'
        c1['pass'] = bool(c1['b_contains_0'] and c1['e_contains_0'] and c1['c_overlaps_plus_half'] and c1['d_overlaps_plus_one']
                          and all(p['overlap'] for p in pts.values()))
        out['C1'] = c1
        # C2: resonance null at kappa* (engine D at the kappa* ball; same transport and centre stage)
        ks = arb('2.810552549') + sym(arb('3.83e-10'))
        d_ks = L.D(acb(ks), 1, 0); dn_ks = L.D(acb(ks), 0, 1)
        out['C2'] = dict(D_HKA_kstar=S(d_ks.real), D_new_kstar=S(dn_ks.real), refused=bool(d_ks.real.contains(0) and dn_ks.real.contains(0)),
                         D_HKA_2_own_transport=S(detM), D_HKA_2_excludes_0=bool(not detM.contains(0)))
        gate = out['C1']['pass'] and out['C8']['e_contains_0']
        n1ok = 'N1 row 3 == r2_forced_k13.log K5 row 3 : True' in open(HERE + 'n1_rows.log').read()
        out['gate_N1_C1_C8'] = dict(N1=n1ok, C1=out['C1']['pass'], C8=out['C8']['e_contains_0'], open=bool(gate and n1ok))
        # ---- N7
        n7 = ex['N7']; cells = n7['cells']
        cov_lo = min(arb(c['lo'].split(' +/-')[0].strip('[')) for c in cells)
        out['N7'] = dict(cells=len(cells), all_Vprime_pos=all(c['Vprime_pos'] for c in cells),
                         zero_cells=[c for c in cells if 'zero' in c], undecided=[c for c in cells if c.get('undecided')],
                         x_V=S(des(n7['xV'])) if 'xV' in n7 else None)
        if not (gate and n1ok):
            out['N8'] = 'NOT COMPUTED: gate N1/C1/C8 failed'
            return out
        # ---- N8 (route A) and C6 (route B)
        mR = match['R']
        rec = ex['checks']['xV']
        AV = des(rec['A']); xV = des(rec['x'])
        yA = [des(rec['R'][i]) + mR[0]*des(rec['hka'][i]) for i in range(4)]
        psiA = yA[0]/AV
        out['N8'] = dict(x_V=S(xV), A_xV=S(AV), Abar2_xV=S(yA[0]), Psi2=S(psiA), sign=(1 if psiA > 0 else (-1 if psiA < 0 else 0)),
                         rel_width=(psiA.rad()/abs(psiA)).str(3) if not psiA.contains(0) else 'contains 0',
                         rel_width_le_1em3=bool((not psiA.contains(0)) and psiA.rad()/abs(psiA) <= arb('1e-3')),
                         units='per unit K eps s e^{2s}, K = N_c/(30 pi); Psi2 = Abar_2(x_V)/A(x_V)')
        RB = ck['RB']
        G = arb_mat([[des(t) for t in r] for r in RB['G']]); Sv = arb_mat([[des(t)] for t in RB['S']])
        ycm = arb_mat([[cs['R'][i] - mR[1]*rho1[i] - mR[2]*rho2[i]] for i in range(4)])
        yend = G*ycm - Sv
        vs = RB['Vstep']
        psiH = arb_mat([[des(t) for t in r] for r in vs['psiH']]); phiV = arb_mat([[des(t) for t in r] for r in vs['phiV']])
        wH = arb_mat([[des(t)] for t in vs['wH']]); wV = arb_mat([[des(t)] for t in vs['wV']])
        yB = phiV*(psiH*yend - (wH - wV))
        psiB = yB[0, 0]/AV
        out['C6'] = dict(Psi2_routeB=S(psiB), Psi2_routeA=S(psiA), overlap=bool(psiA.overlaps(psiB)),
                         yB_xV=[S(yB[i, 0]) for i in range(4)], yA_xV=[S(v) for v in yA])
        return out
    finally:
        ctx.prec = old
        out['time_s'] = round(time.time() - T0, 1)


if __name__ == '__main__':
    o = main()
    json.dump(o, open(HERE + 'final.json', 'w'), indent=1, default=str)
    print(json.dumps(o, indent=1, default=str)[:6000])
