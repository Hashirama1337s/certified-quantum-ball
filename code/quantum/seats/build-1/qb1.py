"""build-1 QB1 instrument (rigorous, python-flint arb): b_ab = r^4 B_(ab) on the certified Evans-Coleman k = 1/3
background, Route A (on-shell) and Route B (metric jets, Zahn eq. 19 / Bach), controls, null, sign map, centre.

usage (memcap 8, from this folder):
  qb1.py controls            -> controls_qb1.json
  qb1.py main                -> result_qb1.json  (sonic, routes, sign_map, centre)
  qb1.py null                -> null_qb1.json    (needs null_bg.pkl from null_bg.py)
  qb1.py probe               -> timing probe only
Background objects and their HOLD seals: see bg.py.  Charts: see charts.py.  Nothing is sealed or pinned here.
"""
import sys, os, time, json, hashlib, math
from flint import arb, ctx
ctx.cap = 64
import bg
import charts as CH
import ctl
from aser import Ser

HERE = os.path.dirname(os.path.abspath(__file__))
K = 7                      # x-series order (>= 5 needed for deg-5 derivative jets)
def KAP():
    return arb(1)/3            # at the current precision (never at import time)
COMP = ('nn', 'nr', 'rr', 'T', 'uu')
INSTR = ['ajet.py', 'aser.py', 'charts.py', 'bg.py', 'ctl.py', 'qb1.py']


def S(a):
    """ball -> [mid, rad] strings (20 significant digits for mid, rad rounded up)"""
    a = arb(a)
    return [a.mid().str(25, radius=False), a.rad().str(3, radius=False)]


def lo_hi(a):
    a = arb(a)
    return [a.lower().str(20, radius=False), a.upper().str(20, radius=False)]


def sha(fn):
    return hashlib.sha256(open(os.path.join(HERE, fn), 'rb').read()).hexdigest()


def sgn(a):
    return 1 if a > 0 else (-1 if a < 0 else 0)


# ------------------------------------------------------------------ evaluators
def ev_centre(Bg, th, chart='cart', deg=4):
    x0, prof = Bg.centre_profiles_theta(arb(th), K)
    if chart == 'cart':
        return x0, CH.cart(x0.exp(), 'x', x0, *prof, KAP(), deg=deg)
    return x0, CH.polar(x0, *prof, KAP(), deg=deg)


def ev_record_point(Bg, which, x, chart='polar', deg=4, want_dx=False):
    x = arb(x)
    ybox, (cN, co, cV) = Bg.point_series(which, x, K)
    prof, chk = bg.profiles_from_NomV(x, cN, co, cV, K)
    if chart == 'cart':
        o = CH.cart(x.exp(), 'x', x, *prof, KAP(), deg=deg)
    else:
        o = CH.polar(x, *prof, KAP(), deg=deg, want_dx=want_dx)
    o['A_identity_check'] = chk
    o['state'] = ybox
    return x, o


def ev_record_interval(Bg, which, xa, xb, chart='polar', deg=4, want_dx=False):
    xa, xb = arb(xa), arb(xb)
    box = Bg.state_box(which, xa, xb)
    cN, co, cV = Bg.ode_jets(box, K)
    x0 = xa.union(xb)
    prof, chk = bg.profiles_from_NomV(x0, cN, co, cV, K)
    if chart == 'cart':
        o = CH.cart(x0.exp(), 'x', x0, *prof, KAP(), deg=deg, want_routeA=False)
    else:
        o = CH.polar(x0, *prof, KAP(), deg=deg, want_dx=want_dx, want_routeA=False)
    return x0, o


def ev_sonic(Bg, which, th, chart='polar', deg=4, want_dx=False, routeA=True):
    x0, cN, co, cV = Bg.sonic_series(which, arb(th), K)
    prof, chk = bg.profiles_from_NomV(x0, cN, co, cV, K)
    o = CH.polar(x0, *prof, KAP(), deg=deg, want_dx=want_dx, want_routeA=routeA)
    o['A_identity_check'] = chk
    return x0, o


def bnn_key(o):
    return o['hat_b_nn'] if o.get('chart') == 'cart' else o['b_nn']


# ------------------------------------------------------------------ sign map
def cover(evalf, a, b, n0, minw, log, tag):
    """tile [a, b] (exact arb endpoints) by intervals on which b_nn has a certified sign; bisect down to minw.
    evalf(lo, hi) -> (xrange, o).  Returns (cells, unresolved)."""
    cells = []; unres = []
    work = [(a + (b - a)*i/n0, a + (b - a)*(i + 1)/n0) for i in range(n0)]
    while work:
        lo, hi = work.pop(0)
        try:
            xr, o = evalf(lo, hi)
            v = bnn_key(o)
            s = sgn(v)
        except ArithmeticError as e:
            xr, v, s = None, None, 0
        if s != 0:
            cells.append(dict(region=tag, p=[S(lo), S(hi)], x=lo_hi(xr), sign=s, b_nn_enc=S(v),
                              scale=('b_nn e^{-4x} (cart, t=-1)' if o.get('chart') == 'cart' else 'b_nn (polar, r=1)')))
        elif (hi - lo).abs_upper() > minw:
            m = (lo + hi)/2
            work[0:0] = [(lo, m), (m, hi)]
        else:
            unres.append((lo, hi))
    log('  cover %s: %d certified cells, %d unresolved' % (tag, len(cells), len(unres)))
    return cells, unres


def sign_map(Bg, log):
    """Certified signs of b_nn at grid POINTS (point evaluations are sharp).  Interval covering was attempted (cover())
    and refused: naive and second- to fourth-order Taylor-remainder interval evaluations wrap by 1e4..1e20 per unit
    relative input width (see RUN.md), so no interval cell certifies at affordable widths.  Points only."""
    T0 = time.time()
    pts = []
    for i in range(0, 151):                         # [-12, -4.5] centre chart, step 0.05
        pts.append(('centre', arb(-12) + arb(i)/20))
    x = arb('-4.45')
    while x < arb('-1.5'):
        pts.append(('tube', x)); x = x + arb('0.05')
    x = arb('-1.5')
    while x < arb('-0.05'):
        pts.append(('tube', x)); x = x + arb('0.01')
    pts.append(('tube', arb('-0.05')))
    for i in range(-9, 10):
        pts.append(('sonic', arb(i)/200))           # sonic manifold, theta in [-0.045, 0.045]
    x = arb('0.05')
    while x < arb(1):
        pts.append(('outw', x)); x = x + arb('0.01')
    x = arb(1)
    while not x > arb(12):
        pts.append(('outw', x)); x = x + arb('0.05')
    rows = []; nfail = 0
    for reg, p in pts:
        try:
            if reg == 'centre':
                th = arb(math.exp(float(p.mid()) - 0.212365646758945385))
                x0, o = ev_centre(Bg, th, 'cart'); v = o['hat_b_nn']; scale = 'hat'
            elif reg == 'sonic':
                x0, o = ev_sonic(Bg, 'tube', p, routeA=False); v = o['b_nn']; scale = 'b'
            else:
                x0, o = ev_record_point(Bg, reg, p, 'polar'); v = o['b_nn']; scale = 'b'
            s = sgn(v)
        except ArithmeticError as e:
            x0, v, s, scale = p, arb(0, 1e300), 0, 'err'
        if s == 0:
            nfail += 1
        rows.append(dict(region=reg, x=S(x0), sign=s, value=S(v), scale=scale))
    signs = [r['sign'] for r in rows]
    changes = [i for i in range(len(signs) - 1) if signs[i]*signs[i + 1] == -1]
    out = dict(kind='certified signs at grid points (NOT interval cells)', n_points=len(rows), n_uncertified=nfail,
               sign_changes_between=[[rows[i]['x'], rows[i + 1]['x']] for i in changes],
               negative_points=sum(1 for s in signs if s < 0), positive_points=sum(1 for s in signs if s > 0),
               first_positive_x=next((r['x'] for r in rows if r['sign'] > 0), None),
               last_negative_x=[r['x'] for r in rows if r['sign'] < 0][-1] if any(s < 0 for s in signs) else None,
               rows=rows, time_s=round(time.time() - T0, 1),
               scale_note="'hat' = b_nn e^{-4x} (Cartesian chart at t = -1; same sign as b_nn), 'b' = b_nn (polar r = 1)")
    log('  sign map: %d points, %d uncertified, sign changes between %s' % (len(rows), nfail, out['sign_changes_between']))
    return out


def zero_analysis(Bg, log, a0='-0.86', b0='-0.82'):
    res = {}
    a, b = arb(a0), arb(b0)
    _, oa = ev_record_point(Bg, 'tube', a); _, ob = ev_record_point(Bg, 'tube', b)
    sa, sb = sgn(oa['b_nn']), sgn(ob['b_nn'])
    res['bracket_start'] = dict(a=S(a), b=S(b), b_nn_a=S(oa['b_nn']), b_nn_b=S(ob['b_nn']))
    res['sign_change_certified'] = (sa*sb == -1)
    if sa*sb != -1:
        res['error'] = 'no certified sign change on the start bracket'
        return res
    lo, hi = a, b; it = 0
    while (hi - lo).abs_upper() > arb('2e-9') and it < 80:
        m = (lo + hi)/2
        _, om_ = ev_record_point(Bg, 'tube', m)
        s = sgn(om_['b_nn'])
        if s == 0:
            res['stopped'] = 'midpoint sign not certified at %s' % m.str(15); break
        if s == sa: lo = m
        else: hi = m
        it += 1
    res['zero_enclosure'] = [lo.str(20, radius=False), hi.str(20, radius=False)]
    res['zero_enclosure_width'] = (hi - lo).str(5, radius=False)
    res['bisection_steps'] = it
    res['existence'] = 'b_nn continuous (real-analytic background), certified opposite signs at both ends: a zero exists in the enclosure'
    # uniqueness: certified d b_nn/dx sign on [c - d, c + d] by the mean-value form
    #   b'(I) in b'(c) + b''(I) [-d, d],  b'' = 2 J_2 + J_1  (J(dr) = b(x + ln(1 + dr)), deg-6 interval jets on I),
    # d doubled while it certifies.  b'(c) from a sharp point evaluation (deg 5).
    c = (lo + hi)/2
    d = (hi - lo)
    _, oc5 = ev_record_point(Bg, 'tube', c, 'polar', deg=5, want_dx=True)
    dc = oc5['db_nn_dx']
    best = None; tries = []
    for _ in range(40):
        try:
            x0, o = ev_record_interval(Bg, 'tube', c - d, c + d, 'polar', deg=6, want_dx=True)
            J1, J2 = o['J_dr'][1], o['J_dr'][2]
            b2 = 2*J2 + J1
            dv = dc + b2*arb(0, d.abs_upper()); ok = bool(dv > 0 or dv < 0)
        except ArithmeticError:
            dv, ok = None, False
        tries.append(dict(half_width=d.str(5, radius=False), db_nn_dx=S(dv) if dv is not None else None, certified=ok))
        if not ok:
            break
        best = (d, dv); d = d*2
    res['derivative_tries'] = tries
    if best is not None:
        res['uniqueness_interval'] = [(c - best[0]).str(20, radius=False), (c + best[0]).str(20, radius=False)]
        res['uniqueness_db_nn_dx'] = S(best[1])
        res['unique_zero_in_uniqueness_interval'] = True
    else:
        res['unique_zero_in_uniqueness_interval'] = False
    _, oc = ev_record_point(Bg, 'tube', c, 'polar', deg=5, want_dx=True)
    res['db_nn_dx_at_zero_point'] = S(oc['db_nn_dx'])
    log('  zero in [%s, %s] width %s; uniqueness interval %s' % (lo.str(15), hi.str(15), (hi - lo).str(3), res.get('uniqueness_interval')))
    return res


# ------------------------------------------------------------------ sonic, routes, centre
def interval_attempts(Bg):
    """documentation of the refused interval covering: naive interval b_nn on 0.01 cells"""
    out = []
    for a, which, ch in (('-3.75', 'tube', 'cart'), ('-2', 'tube', 'cart'), ('-1', 'tube', 'polar'), ('-0.3', 'tube', 'polar'),
                         ('1', 'outw', 'polar'), ('11', 'outw', 'polar')):
        for h in ('1e-2', '1e-4'):
            xa = arb(a); xb = xa + arb(h)
            try:
                _, o = ev_record_interval(Bg, which, xa, xb, ch); v = bnn_key(o)
                out.append(dict(cell=[a, h], chart=ch, b_nn_enclosure=S(v), sign_certified=bool(v > 0 or v < 0)))
            except ArithmeticError as e:
                out.append(dict(cell=[a, h], chart=ch, error=str(e)[:80]))
    return out


def pack(o, x0, tag, src):
    d = dict(tag=tag, source=src, x=S(x0), chart=o.get('chart'))
    for c in COMP:
        d['B_' + c] = S(o['b_' + c])
        d['A_' + c] = S(o['a_' + c])
        d['overlap_' + c] = bool(o['b_' + c].overlaps(o['a_' + c]))
    d['overlap_all'] = all(d['overlap_' + c] for c in COMP)
    E = [o['E_nn'], o['E_nr'], o['E_rr'], o['E_T']]
    d['einstein_residual_r2'] = dict(nn=S(E[0]), nr=S(E[1]), rr=S(E[2]), TT=S(E[3]))
    d['einstein_residual_max_abs_upper'] = max(e.abs_upper() for e in E).str(3, radius=False)
    d['einstein_contains_0_all'] = all(e.contains(0) for e in E)
    d['R_r2'] = S(o['R']); d['trace_dev'] = S(o['trace_dev']); d['xi_diff_uu'] = S(o['xi_diff_uu'])
    d['T1_uu'] = S(o['T1_uu'])
    if o.get('chart') == 'cart':
        d['hat_b_nn'] = S(o['hat_b_nn'])
    d['norms'] = dict(n=S(o['norm_n']), u=S(o['norm_u']), n_e=S(o['n_e']))
    if 'A_identity_check' in o:
        d['A_identity_check'] = o['A_identity_check']
    return d


def route_points(Bg, log):
    rows = []
    for xt in (-12, -9, -6, '-4.5'):
        th = arb(math.exp(float(arb(xt).mid()) - 0.212365646758945385))
        for ch in ('polar', 'cart'):
            x0, o = ev_centre(Bg, th, ch)
            rows.append(pack(o, x0, 'centre chart theta=%s' % th.str(10), 'tube.cm ' + ch))
    for xt in ('-4.5', '-3', '-2', '-1.5', '-1', '-0.84', '-0.5', '-0.2489'):
        x0, o = ev_record_point(Bg, 'tube', arb(xt), 'polar')
        rows.append(pack(o, x0, 'tube x=%s' % xt, 'tube x-step, method b, polar'))
        if xt in ('-3', '-1'):
            x0, o = ev_record_point(Bg, 'tube', arb(xt), 'cart')
            rows.append(pack(o, x0, 'tube x=%s' % xt, 'tube x-step, method b, cart'))
    x0, o = ev_sonic(Bg, 'bg2', arb(0))
    rows.append(pack(o, x0, 'sonic x=0', 'BG2 sonic manifold, polar'))
    xlc = arb('1.357941265765732')
    for xt in ('0.2489', '0.5', '1', xlc, '2', '3', '6', '9', '12'):
        x0, o = ev_record_point(Bg, 'outw', arb(xt), 'polar')
        rows.append(pack(o, x0, 'outward x=%s' % (xt if isinstance(xt, str) else 'x_lc'), 'BG2 x-step, method b, polar'))
    for r in rows:
        log('  %-28s x=%-24s b_nn %-30s A/B overlap %s Eres %s' % (r['tag'][:28], r['x'][0][:14], r['B_nn'][0][:16] + ' +/- ' + r['B_nn'][1], r['overlap_all'], r['einstein_residual_max_abs_upper']))
    return rows


def sonic_block(Bg, log):
    x0, o = ev_sonic(Bg, 'bg2', arb(0))
    out = dict(source='BG2 sonic manifold (outward_C.pkl sm, V0 box width 7.6e-29), exact Taylor coefficients at theta=0, polar chart')
    for c in COMP:
        out['b_' + c] = S(o['b_' + c]); out['width_' + c] = (2*o['b_' + c].rad()).str(3, radius=False)
        out['routeA_' + c] = S(o['a_' + c])
    out['max_width'] = max(2*o['b_' + c].rad() for c in COMP).str(3, radius=False)
    out['target_width_1e-15_met'] = bool(max(2*o['b_' + c].rad() for c in COMP) <= arb('1e-15'))
    out['b_T_equals_(b_nn-b_rr)/2'] = bool(o['trace_dev'].contains(0))
    x1, o1 = ev_sonic(Bg, 'tube', arb(0))
    out['overlap_with_tube_sm_(E5)'] = all(o['b_' + c].overlaps(o1['b_' + c]) for c in COMP)
    out['tube_sm_values'] = {c: S(o1['b_' + c]) for c in COMP}
    # shells: sonic manifold vs x-step records at the junctions x = -0.05 and +0.05
    sh = []
    for th in (arb('-0.048'), arb('0.0565')):
        xs, os_ = ev_sonic(Bg, 'tube', th)
        which = 'tube' if th < 0 else 'outw'
        try:
            xr, orr = ev_record_interval(Bg, which, xs, xs, 'polar') if False else ev_record_point(Bg, which, xs, 'polar')
            ok = all(os_['b_' + c].overlaps(orr['b_' + c]) for c in COMP)
            sh.append(dict(theta=S(th), x=S(xs), sm=S(os_['b_uu']), record=S(orr['b_uu']), overlap_all=ok))
        except ArithmeticError as e:
            sh.append(dict(theta=S(th), x=S(xs), error=str(e)))
    out['shell_overlaps'] = sh
    log('  sonic: ' + ', '.join('b_%s %s' % (c, o['b_' + c].str(15)) for c in COMP))
    return out


def centre_block(Bg, log):
    w, oddmax, xc = Bg.centre_profiles_w(4)
    o = CH.cart(arb(0), 'w', arb(0), *w, KAP())
    out = dict(source='tube.cm stored coefficients at theta=0 (exact low-order Taylor coefficients), Cartesian chart at r=0, t=-1',
               c0=S(o['hat_b_nn']), sign=sgn(o['hat_b_nn']),
               hat_b_rr=S(o['hat_b_rr']), hat_b_T=S(o['hat_b_T']), hat_b_nr=S(o['hat_b_nr']), hat_b_uu=S(o['hat_b_uu']),
               routeA_hat_nn=S(o['hat_a_nn']), routeA_overlap=bool(o['hat_a_nn'].overlaps(o['hat_b_nn'])),
               isotropy_rr_eq_c0_over_3=bool((o['hat_b_rr'] - o['hat_b_nn']/3).contains(0)),
               isotropy_T_eq_rr=bool((o['hat_b_T'] - o['hat_b_rr']).contains(0)),
               hat_einstein=dict(nn=S(o['hat_E_nn']), rr=S(o['hat_E_rr']), TT=S(o['hat_E_T'])),
               w_series_odd_part_max=oddmax.str(3, radius=False),
               note='c0 = lim_{x -> -inf} b_nn e^{-4x} = B_(nn) at the centre at t = -1 (KHA time units, z = r/(-t) = e^x)')
    # approach check: hat_b_nn at theta points (centre chart) must approach c0
    appr = []
    for xt in (-12, -9):
        th = arb(math.exp(xt - 0.212365646758945385))
        x0, oc = ev_centre(Bg, th, 'cart')
        appr.append(dict(x=S(x0), hat_b_nn=S(oc['hat_b_nn'])))
    out['approach'] = appr
    log('  centre: c0 = %s (sign %d)' % (o['hat_b_nn'].str(20), out['sign']))
    return out


# ------------------------------------------------------------------ controls
def controls(log):
    out = {}
    h = ctl.hayward()
    out['C0_hayward_xi0'] = dict(got=[S(v) for v in h['got']], expected=['1/80', '17/240', '-7/480'], contains=h['contains'],
                                 trace=S(h['trace']), pass_=all(h['contains']) and h['trace'].contains(0))
    hf = ctl.hayward(flip=True)
    out['NULL_signflip_hayward'] = dict(got=[S(v) for v in hf['got']], contains=hf['contains'], refuses=not all(hf['contains']))
    r = ctl.rn()
    out['C3_RN_M1_Q0.5_r3'] = dict(V_xi0=[S(v) for v in r['V_xi0']], V_xi16=[S(v) for v in r['V_xi1/6']],
                                   expected_V=[S(v) for v in r['expected_V']], Bach=[S(v) for v in r['Bach']],
                                   Bach_eq_60V=r['Bach_eq_60V'], routeA=[S(v) for v in r['routeA']], R=S(r['R']),
                                   pass_=all(r['V_xi0_contains'] + r['V_xi1/6_contains'] + r['Bach_eq_60V'] + r['routeA_contains']))
    op, oc = ctl.minkowski()
    zero_ok = lambda o, keys: all(o[k].contains(0) for k in keys)
    kp = ['b_' + c for c in COMP] + ['a_' + c for c in COMP] + ['R']
    out['C2_minkowski_KHA_N=e^-x'] = dict(polar={k: S(op[k]) for k in kp}, cart={k: S(oc[k]) for k in kp},
                                         pass_=zero_ok(op, kp) and zero_ok(oc, kp))
    fr = []
    allp = True
    for z in ('0.5', '1', '1.5'):
        op, oc, pred, rp = ctl.frw(z)
        row = dict(z=z, predicted_r4Vuu_xi0=S(pred), predicted_r2R=S(rp))
        for nm, o in (('polar', op), ('cart', oc)):
            ok = (o['V0_uu'].overlaps(pred) and o['V6_uu'].contains(0) and o['R'].overlaps(rp)
                  and all(o['b_' + c].contains(0) for c in COMP) and all(o[e].contains(0) for e in ('E_nn', 'E_nr', 'E_rr', 'E_T')))
            row[nm] = dict(V0_uu=S(o['V0_uu']), V6_uu=S(o['V6_uu']), r2R=S(o['R']), bach={c: S(o['b_' + c]) for c in COMP},
                           einstein={e: S(o[e]) for e in ('E_nn', 'E_nr', 'E_rr', 'E_T')}, xi_null_breaks=bool(not o['xi_diff_uu'].contains(0)),
                           pass_=ok)
            allp &= ok
        fr.append(row)
    out['C6_friedmann_k1/2'] = dict(rows=fr, pass_=allp)
    out['all_positive_controls_pass'] = (out['C0_hayward_xi0']['pass_'] and out['C3_RN_M1_Q0.5_r3']['pass_'] and
                                         out['C2_minkowski_KHA_N=e^-x']['pass_'] and out['C6_friedmann_k1/2']['pass_'])
    out['signflip_null_refuses'] = out['NULL_signflip_hayward']['refuses']
    log('  controls: positives pass %s, sign-flip null refuses %s' % (out['all_positive_controls_pass'], out['signflip_null_refuses']))
    return out


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else 'main'
    T0 = time.time()
    log = lambda *a: print('[%7.1fs]' % (time.time() - T0), *a, flush=True)
    old = ctx.prec
    try:
        ctx.prec = 256
        meta = dict(instrument_sha256={f: sha(f) for f in INSTR}, prec=256, cap=ctx.cap, K_xseries=K,
                    inputs=dict(tube=bg.TUBE, outward=bg.OUTW), python_flint=__import__('flint').__version__)
        if mode == 'controls':
            out = dict(meta=meta, controls=controls(log))
            out['time_s'] = round(time.time() - T0, 1)
            json.dump(out, open(os.path.join(HERE, 'controls_qb1.json'), 'w'), indent=1)
            log('wrote controls_qb1.json')
            return 0
        Bg = bg.Background()
        if mode == 'probe':
            for xt in ('-2', '0.5'):
                t = time.time(); ev_record_interval(Bg, 'tube' if xt[0] == '-' else 'outw', arb(xt), arb(xt) + arb('0.01'), 'polar')
                log('interval eval %s: %.2fs' % (xt, time.time() - t))
            return 0
        if mode == 'main':
            res = dict(meta=meta)
            log('sonic'); res['sonic'] = sonic_block(Bg, log)
            log('centre'); res['centre'] = centre_block(Bg, log)
            log('routes'); res['routes'] = route_points(Bg, log)
            res['routes_summary'] = dict(n_points=len(res['routes']), all_overlap=all(r['overlap_all'] for r in res['routes']),
                                         all_einstein_contain_0=all(r['einstein_contains_0_all'] for r in res['routes']))
            log('sign map'); sm = sign_map(Bg, log)
            log('zero'); sm['zero'] = zero_analysis(Bg, log)
            log('interval attempts'); sm['interval_cell_attempts'] = interval_attempts(Bg)
            res['sign_map'] = sm
            res['time_s'] = round(time.time() - T0, 1)
            json.dump(res, open(os.path.join(HERE, 'result_qb1.json'), 'w'), indent=1)
            log('wrote result_qb1.json')
            return 0
    finally:
        ctx.prec = old


if __name__ == '__main__':
    sys.exit(main())
