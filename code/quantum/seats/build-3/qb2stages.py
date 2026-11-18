"""build-3 QB2 nodes N3 (forced sonic series), N4a + N4 (centre orders, forced centre series).
usage (memcap 8):  py -3.11 qb2stages.py   ->  stages.json  (vectors at theta0 and theta_m, as exact (man,exp) ints)
The engine (count_b4v3.install(), CachedLinear2, sha256-pinned cache) supplies the homogeneous sonic and centre
matrices M0[n] + kappa M1[n], their circle sups, theta0, theta_m, s^HKA and rho1, rho2; this file adds the
forcing columns (qb2core) and solves the augmented series.
"""
import sys, os, time, json, pickle, math
from flint import ctx as _ctx
_ctx.cap = 32
sys.path.insert(0, '<LAB>/engine/cert'); sys.path.insert(0, '<LAB>/engine')
sys.path.insert(0, '<LAB>/quantum/seats/build-3')
import count_b4v3 as CV
import count_b4 as CB
import lin_stages_v2 as LS2
import sympy as sp
from flint import arb, acb, arb_mat, acb_mat, ctx
from core import sym, csym, Tape, sym_const
import qb2lib as Q
import qb2core as C
from qb2core import NA, IJ, IP, ZC, ser, des
from qb2lib import s_mul, s_div, s_add, s_sub, s_scal, s_up, s_deriv, s_exp, s_sqrt, s_trunc, zlike

HERE = '<LAB>/quantum/seats/build-3/'
KAP = arb(2)
T0 = time.time()


def log(*a):
    print('[%7.1fs]' % (time.time() - T0), *a, flush=True)


def engine(log=log):
    CV.install()
    L = LS2.Linear2(sp.Rational(1, 3), log=log, nsub=CB.NSUB, tol=CB.TOL)
    return L


de = lambda p: arb(p[0]) + sym(arb(p[1]))
VECS = ['R', 'N0', 'C', 'C1', 'C8', 'E']           # particular solutions (z = 1)


def to_mat(rows):
    return arb_mat([[x for x in r] for r in rows])


def theta_form_sonic(cols, xp, K):
    """x-form -> theta dY/dtheta form on the sonic manifold, series of length K+1 per entry"""
    z = zlike(xp[0])
    det = cols['det']
    assert det[0].contains(0), 'detMx(theta=0) must contain 0'
    dsh = det[1:]
    tx = lambda e: s_up(s_mul(xp, e, K))[:K + 1]          # theta x' e
    out = {}
    E = C.x_form_matrix_entries(cols, K + 1, lambda num, d: ('num', num))
    for key, e in E.items():
        if e is None:
            continue
        if isinstance(e, tuple):
            out[key] = s_div(s_mul(xp, e[1], K + 1), dsh, K + 1)  # x' adj S / (det/theta)
        else:
            out[key] = tx(e)
    return out


def sonic_stage(L, tb):
    sm = tb['sm']
    Y = [[de(p) for p in sm['Y'][i]] for i in range(3)]; X = [de(p) for p in sm['X']]
    eps, epsX = de(sm['eps']), de(sm['epsX'])
    nu_s = arb(tb['nu_s']); nuc = arb('0.95')*nu_s
    n_in = len(X) - 1                         # 80
    xp = s_deriv(X)
    ex = s_exp(X, n_in)
    cols = C.generic_columns(Y[0][:n_in], Y[1][:n_in], Y[2][:n_in], xp[:n_in], ex, n_in)
    n = cols['n']; K = n - 2
    E = theta_form_sonic(cols, xp, K)
    Ms = L.sonic
    assert Ms.K >= K
    Mh = []
    for m in range(K + 1):
        M = arb_mat(NA, NA)
        for i in range(4):
            for j in range(4):
                M[i, j] = Ms.M0[m][i][j] + KAP*Ms.M1[m][i][j]
        for (i, j), s in E.items():
            M[i, j] = s[m]
        Mh.append(M)
    # circle sup of the extras (theta-form) minus their theta = 0 values
    narc = 128; pi = arb.pi(); dev = {}
    rad_av = nu_s - nuc*arb('1.0001')
    for a in range(narc):
        th = (arb(a) + arb('0.5') + sym(arb('0.5')))*2*pi/narc
        tb_ = acb(th.cos(), th.sin())*nuc
        loc = [C.local_series([acb(c) for c in Y[i]], tb_, 7, eps, rad_av) for i in range(3)]
        lX = C.local_series([acb(c) for c in X], tb_, 8, epsX, rad_av)
        lxp = s_deriv(lX)
        lex = s_exp(lX, 7)
        cb = C.generic_columns(loc[0], loc[1], loc[2], lxp[:7], lex, 7)
        Eb = C.x_form_matrix_entries(cb, 1, lambda num, d: ('num', num))
        rowsum = {}
        for (i, j), e in Eb.items():
            if e is None:
                continue
            if isinstance(e, tuple):
                v = lxp[0]*e[1][0]/(cb['det'][0]/tb_)
            else:
                v = tb_*lxp[0]*e[0]
            d = abs(v - acb(E[(i, j)][0])).upper()
            dev[(i, j)] = dev.get((i, j), arb(0)).max(d)
    s0, s1 = L.son_info['supC0'], L.son_info['supC1']
    Mh_u = Mh
    Mh, Snorm, lam = C.make_scaled(Mh, dev, s0 + KAP*s1)
    Xs = Snorm
    M0 = [[Mh[0][i, j] for j in range(4)] for i in range(4)]
    trM = sum((M0[i][i] for i in range(4)), arb(0))
    M0n = max(sum((abs(x).upper() for x in r), arb(0)) for r in M0)
    F0n = arb(0)
    for i in range(4):
        F0n = F0n.max(sum((Mh[0][i, j].abs_upper() for j in range(4, NA)), arb(0)))

    def cn(nn):
        gap = abs(arb(nn) - trM).lower()
        if not gap > 0:
            return None
        cM = (1 + M0n/gap)/nn
        return (cM*(1 + F0n/nn)).max(arb(1)/nn)
    RS = C.RegSingAug(Mh, Snorm, nuc, cn)
    # order-0 data
    w_i = 2 if float(sum(abs(Mh_u[0][2, j]).mid() for j in range(4))) >= float(sum(abs(Mh[0][3, j]).mid() for j in range(4))) else 3
    o_i = 5 - w_i
    I0 = [arb(L.I0[0]) + KAP] + [arb(x) for x in L.I0[1:]]
    sI0 = {nm: (cols[nm][2][0] if cols[nm][0] == 'beta' else arb(0)) for nm in ('R', 'C', 'C8', 'E', 'P')}
    sI0['N0'] = arb(0); sI0['C1'] = arb(0)
    res = {}; info = dict(K=K, Snorm=Snorm.str(5), lam=[l.str(3) for l in lam], w_row=w_i, det0=cols['det'][0].str(5), dev_max={str(k): v.str(3) for k, v in dev.items()})
    th0 = L.th0
    drivers = {nm: {ZC[nm]: arb(1)} for nm in VECS}
    drivers['Ph'] = {IP: arb(1)}
    for nm, d in drivers.items():
        rhs_w = -sum((Mh_u[0][w_i, j]*v for j, v in d.items()), arb(0))
        sig = sum(((sI0['P'] if j == IP else sI0[[k for k, z in ZC.items() if z == j][0]])*v for j, v in d.items()), arb(0))
        Mat = arb_mat([[Mh_u[0][w_i, j] for j in range(4)], I0, [1, 0, 0, 0], [0, 1, 0, 0]])
        y4 = Mat.solve(arb_mat([[rhs_w], [sig], [0], [0]]))
        y0 = arb_mat(NA, 1)
        for i in range(4): y0[i, 0] = y4[i, 0]
        for j, v in d.items(): y0[j, 0] = v
        other = sum((Mh_u[0][o_i, j]*y0[j, 0] for j in range(NA)), arb(0))
        y0s = arb_mat(NA, 1)
        for i in range(NA): y0s[i, 0] = y0[i, 0]/lam[i]
        ys = RS.series(y0s)
        val, tinf = RS.evaluate(ys, th0)
        val = [val[i]*lam[i] for i in range(NA)]
        res[nm] = dict(val=[ser(v) for v in val], other_row=other.str(5), other_row_contains_0=bool(other.contains(0)),
                       tail=tinf['tail'].str(3), n0=tinf['n0'])
        log('sonic %s: Y(theta0) = %s  other-row %s tail %s' % (nm, [v.str(8) for v in val[:4]], other.str(3), tinf['tail'].str(3)))
    Yh, hinfo = L.sonic_start(acb(2), 1, 0)
    res['hka'] = dict(val=[ser(v.real) for v in Yh] + [ser(arb(0))]*(NA - 4), tail=hinfo['tail'].str(3))
    # sonic-point values for the record (Psi2 gauge: Abar(0) = Nbar(0) = 0 on every particular)
    info['b_sonic'] = [cols['b'][c][0].str(15) for c in range(4)]
    info['trace_sonic'] = cols['trace'][0].str(5)
    info['Dt_sonic_coeffs_max_abs'] = max(abs(x).upper() for x in cols['Dt'][:40]).str(3)
    info['Dt_sonic_coeffs_contain_0'] = all(x.contains(0) for x in cols['Dt'][:40])
    return res, info


def centre_stage(L, tb, thetas):
    cm = tb['cm']
    Z = [[de(p) for p in cm['Y'][i]] for i in range(3)]; chi = [de(p) for p in cm['X']]
    eps, epsX = de(cm['eps']), de(cm['epsX'])
    nu = arb(cm['nu']); nuc = arb('0.95')*nu
    xc = de(tb['xc'])
    n = len(chi)                              # 61
    z = arb(0)
    P, W, U = [Z[i][:n] for i in range(3)]
    assert P[0].is_zero()
    Qs = P[1:] + [z]
    one = C.one_like(z, n)
    fac = s_add(one, s_up(s_deriv(chi))[:n])
    Dc = lambda s: s_div(s_up(s_deriv(s))[:n], fac, n)
    cv = Q.conv_for(z)
    A = Q.series_out(Q.T_ACEN, [P, W, U], n, z, cv)[0]
    Aw = Q.series_out(Q.T_ACEN, [P, s_scal(W, 1 + C.DELTA_NC), U], n, z, cv)[0]
    tq = s_up(s_deriv(Qs))[:n]
    l1 = s_scal(s_div(s_add(one, s_div(tq, Qs, n)), fac, n), -1)
    lj = [l1]
    for j in range(3): lj.append(Dc(lj[-1]))
    info = {}

    def bach(Ab):
        Aj = [Ab]
        for j in range(4): Aj.append(Dc(Aj[-1]))
        return Q.series_out(Q.T_BACH, [P] + Aj + lj, n, z, cv), Aj
    b, Aj = bach(A)
    bw, _ = bach(Aw)
    # N4a: the structurally vanishing low coefficients (b = O(theta^4), b_nr = O(theta^5); PAPER premise, L3s)
    low = {}
    bh = []
    for nm, bb in (('true', b), ('null_Nc', bw)):
        hats = []
        for c, p in ((0, 4), (1, 5), (2, 4), (3, 4)):
            ok = all(bb[c][k].contains(0) for k in range(p))
            low['%s_b%d' % (nm, c)] = dict(orders_below=p, contain_0=ok, max_abs=max(abs(bb[c][k]).upper() for k in range(p)).str(3))
            if not ok:
                raise ArithmeticError('centre: low coefficient does not contain 0 (%s, %d)' % (nm, c))
            hats.append(bb[c][p:] + [z]*p)
        bh.append(hats)
    info['N4a_b_low'] = low
    E1 = s_exp(s_add(chi[:n], [xc] + [z]*(n - 1)), n)
    E2 = s_div(one, s_mul(E1, E1, n), n)
    pw = {0: 2, 1: 3, 2: 2}
    betaR = [s_up(s_mul(E2, bh[0][c], n), pw[c])[:n] for c in range(3)]
    betaC = [s_up(s_mul(E2, bh[1][c], n), pw[c])[:n] for c in range(3)]
    # centre rows: R_ic = u^{m_ic} Rhat_ic(u, q, A),  u = P, q = U
    Rh = Q.series_out(Q.T_CROWS, [P, U, A], n, z, cv)
    mo = Q.EX['cen_rows_uorder']

    def rows_of(beta):
        out = []
        for i in range(4):
            acc = [z]*n
            for c in range(3):
                if beta[c] is None: continue
                t = s_mul(Rh[3*i + c], beta[c], n)
                for _ in range(mo[(i, c)]): t = s_mul(P, t, n)
                acc = s_add(acc, t)
            out.append(acc)
        return out
    E = {}
    FR = rows_of(betaR)
    # N4a on the rows: Q^{-1}S (x-form) orders 0,1 must be exactly 0 here (structural zeros)
    info['N4a_rows_orders01'] = [[FR[i][k].str(3) for k in range(3)] for i in range(4)]
    E1sq = s_mul(E1, E1, n)
    den = s_add(one, s_up(E1sq, 2)[:n]); den2 = s_mul(den, den, n)
    php = s_up(s_scal(s_div(E1sq, den2, n), 2), 2)[:n]           # phi' = theta^2 2E1^2/(1+theta^2 E1^2)^2
    G2 = [s_div(Aj[1], A, n), s_add(l1, [z + 2] + [z]*(n - 1)), s_sub(s_div(Dc(W), W, n), s_scal(l1, 2)),
          Dc(s_mul(U, P, n))]
    bnr = s_up(s_div(s_mul(E1sq, E1, n), den2, n), 3)[:n]          # theta^3 E1^3/(1+theta^2E1^2)^2
    brr = s_up(s_div(E1sq, den2, n), 2)[:n]
    FE = rows_of([None, bnr, brr])
    rA = s_sqrt(A, n); irA = s_div(one, rA, n)
    FP = rows_of([irA, None, None])
    # g = N bnr' + (2N + 2N l1 + N A'/A) bnr - (A'/(2A)) brr  with N = 1/(theta Q):  N * (theta^k h) exact
    AA = s_div(Aj[1], A, n)
    dbnr = Dc(bnr)
    Ninto = lambda s: s_div(s[1:] + [z], Qs, n)                  # N s for s with exact zero constant term
    assert all(x.is_zero() for x in dbnr[:3]) and all(x.is_zero() for x in bnr[:3])
    g = s_sub(s_add(Ninto(dbnr), s_mul(s_add(s_add([z + 2] + [z]*(n - 1), s_scal(l1, 2)), AA), Ninto(bnr), n)),
              s_scal(s_mul(AA, brr, n), arb(1)/2))
    for i in range(4):
        E[(i, IP)] = FP[i]; E[(i, 6)] = FR[i]; E[(i, 8)] = rows_of(betaC)[i]
        E[(i, 9)] = s_mul(php, G2[i], n); E[(i, 10)] = FE[i]; E[(i, 11)] = FE[i]
        E[(i, 7)] = FR[0] if i == 0 else [z]*n
    E[(IJ, IJ)] = s_scal(A, -1); E[(IJ, 7)] = s_mul(s_sub([z + 2] + [z]*(n - 1), A), FR[0], n)
    E[(IP, IP)] = [z - 2] + [z]*(n - 1); E[(IP, 10)] = s_mul(rA, g, n)
    # theta form: multiply by fac = 1 + theta chi'
    Et = {k: s_mul(fac, v, n) for k, v in E.items()}
    # check: engine centre fac Dc/lam == 1 + theta chi'
    import pipeline as PL
    Dtape = Tape([PL.load_expr(PL.FIELDS['pol_c'][3])], [sp.Symbol('P'), sp.Symbol('W'), sp.Symbol('U')], {sp.Symbol('k'): sp.Rational(1, 3)})
    dc = Q.series_out(Dtape, [P, W, U], n, z, cv)[0]
    lam = de(cm['lam'])
    fe = [x/lam for x in dc]
    info['fac_engine_overlaps_1_plus_theta_chi'] = all(fe[k].overlaps(fac[k]) for k in range(n))
    Mc = L.centre
    K = min(Mc.K, n - 1)
    Mh = []
    for m in range(K + 1):
        M = arb_mat(NA, NA)
        for i in range(4):
            for j in range(4):
                M[i, j] = Mc.M0[m][i][j] + KAP*Mc.M1[m][i][j]
        for (i, j), s in Et.items():
            M[i, j] = s[m]
        Mh.append(M)
    # circle sup of the extras
    narc = 128; pi = arb.pi(); dev = {}
    rad_av = nu - nuc*arb('1.0001')
    for a in range(narc):
        th = (arb(a) + arb('0.5') + sym(arb('0.5')))*2*pi/narc
        tb_ = acb(th.cos(), th.sin())*nuc
        L7 = 7
        lP, lW, lU = [C.local_series([acb(c) for c in Z[i]], tb_, L7, eps, rad_av) for i in range(3)]
        lchi = C.local_series([acb(c) for c in chi], tb_, L7 + 1, epsX, rad_av)
        tser = [tb_, acb(1)] + [acb(0)]*(L7 - 2)
        Ns = s_div(C.one_like(acb(0), L7), lP, L7)
        oms = s_mul(lW, s_mul(lP, lP, L7), L7); Vs = s_mul(lU, lP, L7)
        xp = s_add(s_div(C.one_like(acb(0), L7), tser, L7), s_deriv(lchi)[:L7])
        lex = s_mul(tser, s_exp(s_add(lchi[:L7], [acb(xc)] + [acb(0)]*(L7 - 1)), L7), L7)
        cb = C.generic_columns(Ns, oms, Vs, xp, lex, L7)
        Eb = C.x_form_matrix_entries(cb, 1, lambda num, d: num[0]/d[0])
        facb = tb_*xp[0]
        rowsum = {}
        for (i, j), e in Eb.items():
            if e is None: continue
            v = facb*(e if not isinstance(e, list) else e[0])
            d = abs(v - acb(Et[(i, j)][0])).upper()
            dev[(i, j)] = dev.get((i, j), arb(0)).max(d)
    Mh_u = Mh
    Mh, Snorm, lam = C.make_scaled(Mh, dev, L.cen_info['supC0'] + KAP*L.cen_info['supC1'])
    Xs = Snorm
    M0n = C.rowsum_mat(Mh[0])
    cn = lambda nn: (1/(arb(nn) - M0n)) if arb(nn) > M0n else None
    RS = C.RegSingAug(Mh, Snorm, nuc, cn)
    info.update(K=K, Snorm=Snorm.str(5), lam=[l.str(3) for l in lam], M0norm=M0n.str(5), dev_max={str(k): v.str(3) for k, v in dev.items()})
    res = {}
    for nm in VECS:
        y0 = arb_mat(NA, 1); y0[ZC[nm], 0] = 1
        chk = Mh_u[0]*y0
        assert all(chk[i, 0].contains(0) for i in range(NA)), ('centre order 0', nm)
        y0s = arb_mat(NA, 1)
        for i in range(NA): y0s[i, 0] = y0[i, 0]/lam[i]
        ys = RS.series(y0s)
        res[nm] = {}
        for tn, th in thetas.items():
            val, tinf = RS.evaluate(ys, th)
            val = [val[i]*lam[i] for i in range(NA)]
            res[nm][tn] = dict(val=[ser(v) for v in val], tail=tinf['tail'].str(3))
        log('centre %s at theta_m: %s' % (nm, [des(t).str(8) for t in res[nm]['m']['val'][:6]]))
    for tn, th in thetas.items():
        r1, r2 = [], []
        for e_, store in ((1, r1), (2, r2)):
            y0 = [acb(1) if i == e_ else acb(0) for i in range(4)]
            Yv, _ = Mc.evaluate(acb(2), y0, acb(th))
            store.extend(ser(v.real) for v in Yv)
        res.setdefault('rho1', {})[tn] = r1; res.setdefault('rho2', {})[tn] = r2
    return res, info


def theta_at(tb, x):
    """theta with ln theta + chi(theta) + xc = x (interval Newton, contraction checked)"""
    cm = tb['cm']; chi = [de(p) for p in cm['X']]; epsX = de(cm['epsX']); xc = de(tb['xc'])
    dchi = s_deriv(chi)
    x = arb(x)
    th = arb(((x - xc).exp()).mid())
    for _ in range(60):
        g = th.log() + Q.horner(chi, th) + xc - x
        th = arb((th - g/(1/th + Q.horner(dchi, th))).mid())
    r = arb('1e-40')
    for _ in range(80):
        Th = th + sym(r)
        g = th.log() + Q.horner(chi, th) + sym(epsX) + xc - x
        Nt = th - g/(1/Th + Q.horner(dchi, Th) + sym(epsX/(arb(cm['nu']) - Th.abs_upper())))
        if Th.contains_interior(Nt):
            return Nt
        r *= 2
    raise ArithmeticError('theta_at failed')


if __name__ == '__main__':
    old = ctx.prec
    try:
        L = engine()
        ctx.prec = 128
        tb = pickle.load(open('<LAB>/engine/cert/lin_tube2_k1_3.pkl', 'rb'))
        out = dict(prec=ctx.prec, theta0=L.th0.str(30), theta_m=L.theta_m.str(30))
        sres, sinfo = sonic_stage(L, tb)
        log('sonic done', json.dumps(sinfo)[:400])
        th6 = theta_at(tb, -6)
        cres, cinfo = centre_stage(L, tb, dict(m=L.theta_m, x6=th6))
        log('centre done', json.dumps(cinfo, default=str)[:600])
        out.update(sonic=sres, sonic_info=sinfo, centre=cres, centre_info=cinfo, theta_x6=th6.str(30),
                   time_s=round(time.time() - T0, 1))
        json.dump(out, open(HERE + 'stages.json', 'w'), indent=1, default=str)
        log('wrote stages.json')
    finally:
        ctx.prec = old
