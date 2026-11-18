"""b2_outward.py - build-5 phase 2, B2: outward leg rebuilt from the B1 box, with x_lc and the state at x_lc.

Control flow = build-1/outward_rec.py (imported READ-ONLY, never edited), which is driver_v2.outward_leg_v2
(engine/cert/driver_v2.py:262-317) plus step records.  The only parameter changes: the V0 box (from b1_<tag>.json),
PREC and K_t (Lohner order AND record jet order; outward_rec.KT is rebound on the imported module object, the file is
untouched).  x_lc and its step / offset: build-1/l0_cert.lc_root (imported read-only).  Records are strings only.
Reports: relative half-widths on every step whose span meets d = x_lc - x in [0.012, 0.063] (hull and Yc[.][0]),
the state at x_lc, x_lc; overlap with the old record build-1/outward_rec.pkl (old state evaluated at every new step
start by the old step polynomial + Cauchy tail; x_lc; m/r; end state).

usage (memcap 8): py -3.11 b2_outward.py <tag> <b1json> <prec> <Kt> [--v0shift 1e-22]  -> outward_<tag>.pkl, b2_<tag>.json
"""
import sys, os, time, json, pickle, hashlib
HERE = os.path.dirname(os.path.abspath(__file__))
B1DIR = '<LAB>/sideways/build-1'
sys.path.insert(0, B1DIR); sys.path.insert(0, '<LAB>/engine/cert'); sys.path.insert(0, '<LAB>/engine')
import outward_rec as OR                 # pinned build-1 module (sets ctx.prec 200, ctx.cap 32 at import)
import l0_cert as L0
import sympy as sp
from flint import arb, ctx
import driver_v31 as V31
import driver_v32 as V32
import driver_v2 as V2
from core import sym
ctx.cap = 64                             # module top, after the imports that set 32 (F3); asserted in main
OLD = os.path.join(B1DIR, 'outward_rec.pkl')


def rel(a):
    a = arb(a)
    return float(a.rad()/abs(a)) if not a.contains(0) else float('inf')


def old_state_at(old, x):
    """old record: state at x from the old x-step containing x (polynomial + Cauchy tail)"""
    KT = old['KT']
    for st in old['steps']:
        if st['kind'] != 'x':
            continue
        x0 = L0.deser(st['x0']); h = L0.deser(st['h'])
        t = x - x0
        if t.mid() >= 0 and t.mid() <= h.mid():
            Rc = L0.deser(st['Rc']); r = abs(t).upper()/Rc
            out = []
            for i in range(3):
                Y = [L0.deser(p) for p in st['Yc'][i]][:KT + 1]
                tot = L0.deser(st['Bx'][i]).rad()
                out.append(L0.horner(Y, t) + sym(tot*r**(KT + 1)/(1 - r)))
            return out
    return None


def main():
    tag, b1json, PREC, KT = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
    shift = sys.argv[6] if len(sys.argv) > 6 and sys.argv[5] == '--v0shift' else None
    assert ctx.cap >= 64
    T0 = time.time()
    log = lambda *a: print('[%7.1fs]' % (time.time() - T0), *a, flush=True)
    old_prec = ctx.prec
    try:
        ctx.prec = PREC
        OR.KT = KT                                       # record jet order (module attribute; file untouched)
        b1 = json.load(open(b1json))['B1']
        V0box = arb(b1['V0'])
        if shift:
            V0box = V0box + arb(shift)
        k13 = sp.Rational(1, 3)
        Pb = V31.Problem2(k13, prec=PREC, K_sonic=OR.KS, K_centre=OR.KC, K_ty=KT)
        nu, _ = V2.choose_nu(Pb, arb(OR.V0_MID) + sym(OR.NU_AT), OR.OM_T, OR.NU_LADDER, OR.KS, log=log)
        M = V32.MatcherV32(Pb, nu, OR.NUC, Ks=OR.KS, x_switch='-0.05', Kt=KT)
        log('V0 box %s (width %s), nu_s %s, prec %d, Kt %d' % (V0box.str(35, radius=True), (2*V0box.rad()).str(3), nu, PREC, KT))
        r = OR.outward_rec(Pb, M, V0box, OR.OM_T, X1=12.0, log=log)
        log('outward: steps=%d D<0 %s V>0 %s lemma %s' % (r['steps'], r['ok_D'], r['ok_V'], r.get('lemma')))
        ser = OR.ser
        out = dict(k='1/3', V0=(V0box.mid().str(60, radius=False), V0box.rad().str(10, radius=False)), nu_s=nu, X1=12.0, KT=KT,
                   Q='0.25', prec=PREC, steps=r['records'], sm=r['sm'], ok_D=r['ok_D'], ok_V=r['ok_V'],
                   ok_D_manifold_piece=r['ok_D_manifold_piece'], end=[ser(t) for t in r['y1']], lemma=r.get('lemma'),
                   a=ser(r['a']), source='build-5 phase2 b2_outward.py from %s%s' % (os.path.basename(b1json), (' SHIFTED %s' % shift) if shift else ''))
        if r.get('lemma'):
            out.update(mu=ser(r['mu']), om_inf=ser(r['om_inf']), V_inf=ser(r['V_inf']), A_inf=ser(r['A_inf']),
                       m_over_r=ser(r['m_over_r']), Q=[ser(t) for t in r['Q']], DQ=ser(r['DQ']))
        out['instrument_sha256'] = hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()
        pk = os.path.join(HERE, 'outward_%s.pkl' % tag)
        pickle.dump(out, open(pk, 'wb'))                 # strings only (OR.ser), no flint objects
        log('saved %s (%d steps)' % (pk, len(out['steps'])))
        # ---------------- x_lc, state at x_lc, widths near the cone
        rep = dict(tag=tag, V0=V0box.str(40, radius=True), prec=PREC, Kt=KT, steps=len(out['steps']), lemma=r.get('lemma'),
                   ok_D=r['ok_D'], ok_V=r['ok_V'])
        lc = L0.lc_root(out, lambda *a: None)
        xlc = lc['x_lc']
        rep['x_lc'] = xlc.str(40, radius=True); rep['x_lc_width'] = (2*xlc.rad()).str(5); rep['root_step'] = lc['step']
        st = out['steps'][lc['step']]; T = lc['t']; Rc = L0.deser(st['Rc'])
        ylc = []
        for i in range(3):
            Y = [L0.deser(p) for p in st['Yc'][i]][:KT + 1]
            tot = L0.deser(st['Bx'][i]).rad(); q = T.abs_upper()/Rc
            ylc.append(L0.horner(Y, T) + sym(tot*q**(KT + 1)/(1 - q)))
        rep['state_at_xlc'] = [y.str(30, radius=True) for y in ylc]
        rep['state_at_xlc_rel'] = [rel(y) for y in ylc[1:]] + ['N := 1 exactly (root)']
        near = []
        for j, s_ in enumerate(out['steps']):
            if s_['kind'] != 'x':
                continue
            x0 = L0.deser(s_['x0']); h = L0.deser(s_['h'])
            d_hi = xlc - x0; d_lo = xlc - (x0 + h)
            if float(d_lo.mid()) <= 0.063 and float(d_hi.mid()) >= 0.012:
                hull = [L0.deser(p) for p in s_['hull']]
                y0 = [L0.deser(s_['Yc'][i][0]) for i in range(3)]
                near.append(dict(step=j, x0=x0.str(15), hull_rel=[rel(v) for v in hull], Yc0_rel=[rel(v) for v in y0]))
        rep['near_cone_steps'] = near
        rep['near_cone_max_rel'] = max(max(n_['hull_rel'] + n_['Yc0_rel']) for n_ in near)
        rep['state_at_xlc_max_rel'] = max(rep['state_at_xlc_rel'][:2])
        # ---------------- overlap with the old record
        old = pickle.load(open(OLD, 'rb'))
        ov_all = True; nchk = 0; worst = None
        for j, s_ in enumerate(out['steps']):
            if s_['kind'] != 'x':
                continue
            x0 = L0.deser(s_['x0'])
            os_ = old_state_at(old, x0)
            if os_ is None:
                continue
            hull = [L0.deser(p) for p in s_['hull']]
            ok = all(hull[i].overlaps(os_[i]) for i in range(3))
            nchk += 1
            if not ok:
                ov_all = False; worst = worst or j
        lc_old = L0.lc_root(old, lambda *a: None)
        rep['overlap_old'] = dict(x_steps_checked=nchk, all_overlap=ov_all, first_failing_step=worst,
                                  x_lc_old=lc_old['x_lc'].str(25, radius=True), x_lc_overlaps_old=bool(xlc.overlaps(lc_old['x_lc'])),
                                  x_lc_inside_old=bool(lc_old['x_lc'].contains(xlc)),
                                  end_overlaps_old=all(L0.deser(a).overlaps(L0.deser(b)) for a, b in zip(out['end'], old['end'])))
        if r.get('lemma'):
            mr = r['m_over_r']; mr_old = L0.deser(old['m_over_r'])
            rep['m_over_r_t0'] = mr.str(25, radius=True); rep['m_over_r_overlaps_old'] = bool(mr.overlaps(mr_old))
            rep['m_over_r_old'] = mr_old.str(10, radius=True)
        rep['time_s'] = round(time.time() - T0, 1)
        json.dump(rep, open(os.path.join(HERE, 'b2_%s.json' % tag), 'w'), indent=1)
        log('x_lc %s (width %s); near-cone max rel %.3e; state at x_lc max rel %.3e; overlap old %s; x_lc overlaps old %s'
            % (xlc.str(30, radius=True), rep['x_lc_width'], rep['near_cone_max_rel'], rep['state_at_xlc_max_rel'], ov_all,
               rep['overlap_old']['x_lc_overlaps_old']))
        return 0
    finally:
        ctx.prec = old_prec


if __name__ == '__main__':
    sys.exit(main())
