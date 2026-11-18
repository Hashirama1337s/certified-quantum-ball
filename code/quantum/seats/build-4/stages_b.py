"""build-4: N3 / N4 on the BG2 tube (build-3 qb2stages functions, patched by qb2b), plus the C1 point x = -0.01 on the
sonic series.  usage (memcap 8): py -3.11 stages_b.py -> build-4/stages.json"""
import sys, json, time, pickle
sys.path.insert(0, '<LAB>/quantum/seats/build-4')
import qb2b as B
from qb2b import ST, C, Q, sym, arb, acb, arb_mat, ctx
from qb2core import NA, ser, des
T0 = time.time()
log = lambda *a: print('[%7.1fs]' % (time.time() - T0), *a, flush=True)
de = lambda p: arb(p[0]) + sym(arb(p[1]))


def theta_sonic(tb, x):
    """theta on the sonic manifold with X(theta) = x (interval Newton, contraction checked)"""
    sm = tb['sm']; X = [de(p) for p in sm['X']]; epsX = de(sm['epsX']); nu = arb(tb['nu_s'])
    dX = Q.s_deriv(X); x = arb(x)
    th = arb(x.mid())
    for _ in range(60):
        th = arb((th - (Q.horner(X, th) - x)/Q.horner(dX, th)).mid())
    r = arb('1e-40')
    for _ in range(100):
        Th = th + sym(r)
        g = Q.horner(X, th) + sym(epsX) - x
        Nt = th - g/(Q.horner(dX, Th) + sym(epsX/(nu - Th.abs_upper())))
        if Th.contains_interior(Nt):
            return Nt
        r *= 2
    raise ArithmeticError('theta_sonic failed')


old = ctx.prec
try:
    L = B.engine_b(log=log)
    log('engine ready, %d substeps, depths %s' % (len(L.steps), sorted(set(s.get('depth', 0) for s in L.steps))))
    tb = pickle.load(open(B.NEW_TUBE, 'rb'))
    out = dict(prec=ctx.prec, theta0=L.th0.str(30), theta_m=L.theta_m.str(30), tube=B.NEW_TUBE)
    th01 = theta_sonic(tb, '-0.01')
    B._STASH['extra_t'] = th01
    sres, sinfo = ST.sonic_stage(L, tb)
    lam = B._STASH['lams'][0]
    vals = B._STASH['vals']
    assert len(vals) == 7
    c1v = [vals[3][0][i]*lam[i] for i in range(NA)]
    B._STASH['extra_t'] = None
    # s^HKA at theta(-0.01): engine sonic stage (lin_stages.Linear.sonic_start with the evaluation point changed)
    kap = acb(2)
    M0 = [[acb(L.sonic.M0[0][i][j]) + kap*acb(L.sonic.M1[0][i][j]) for j in range(4)] for i in range(4)]
    w = max((M0[2], M0[3]), key=lambda r: float(sum(abs(x).mid() for x in r)))
    I0 = [acb(L.I0[0]) + kap] + [acb(x) for x in L.I0[1:]]
    from flint import acb_mat
    y0 = acb_mat([w, I0, [1, 0, 0, 0], [0, 1, 0, 0]]).solve(acb_mat([[0], [0], [1], [0]]))
    Yh, _ = L.sonic.evaluate(kap, [y0[i, 0] for i in range(4)], acb(th01))
    sm = tb['sm']
    st01 = [Q.horner([de(p) for p in sm['Y'][i]], th01) + sym(de(sm['eps'])) for i in range(3)]
    out['C1_x001'] = dict(theta=th01.str(30), C1=[ser(v) for v in c1v], hka=[ser(v.real) for v in Yh],
                          state=[ser(v) for v in st01], tail=vals[3][1]['tail'].str(3))
    log('sonic done', json.dumps({k: v for k, v in sinfo.items() if k != 'dev_max'})[:300])
    th6 = ST.theta_at(tb, -6)
    cres, cinfo = ST.centre_stage(L, tb, dict(m=L.theta_m, x6=th6))
    log('centre done; tails', {k: v['m']['tail'] for k, v in cres.items() if k not in ('rho1', 'rho2')})
    out.update(sonic=sres, sonic_info=sinfo, centre=cres, centre_info=cinfo, theta_x6=th6.str(30), time_s=round(time.time() - T0, 1))
    json.dump(out, open(B.B4 + 'stages.json', 'w'), indent=1, default=str)
    log('wrote stages.json')
finally:
    ctx.prec = old
