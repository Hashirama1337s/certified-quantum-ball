"""l2_cert.py - build seat build-4, instrument for BAR-SF-L2 "The Light Cone" (seal 9077a35520), 2026-11-17.

Certified argument-principle counts of zeros of the axial l = 2 characteristic function D_2(kappa) (l2_lib.py) on the
E5 Evans-Coleman solution, plus the flat toy and the class-teeth null.

  PL2-1  half-disc {Re kappa >= 0, |kappa| <= 10}: shared engine gravity_engine/argument_count.count (upper path of the
         circumscribed half-disc polygon, exact-Fraction arc check, plan check, 0.25 rule, parity, JSONL ledger).
  PL2-2  flat toy (G App. D, l = 2), regions T-a |k + 3/2| <= 3/4, T-b |k + 1| <= 0.3, T-c |k + 2| <= 0.3
         (answers 1, 0, 1 stated in RUN.md sec. 0 before any toy count ran).
  PL2-3  control disc |kappa - (-1.8 + 6.2i)| <= 0.5 (target N = 1) + an enclosure of its zero (count 1 on a small disc).
  PL2-4  null on |kappa - kappa*| <= 0.1, kappa* = 1 - lambda_lc: analytic class (X1 level, S = {1}) N = 0;
         'bounded' class (value at the cone, S = {}) winding -1 (the pole the collision branch creates).
Closed-disc counts (PL2-2, PL2-3, PL2-4): circumscribed regular polygon (16 vertices, vertex radius rho/cos(pi/16)
(1 + 1e-12)); every edge cut into pieces of length <= 1.6 r; per piece the engine's segment_arg + sharpen; total/2 pi
within 0.25 of an integer (BAR-LB-M rule); JSONL ledger, resumable, header-checked.  The counted region is the polygon,
which contains the disc.

usage (memcap 8, from build-4):
  py -3.11 l2_cert.py --controls   toy end to end, null, resume test, failure-mode tests -> controls_l2.json
  py -3.11 l2_cert.py --probe      setup + timed TMs on a few main-plan discs (no count) -> probe_l2.json
  py -3.11 l2_cert.py --main       LEAD ONLY: PL2-1..PL2-4 -> result_l2.json, ledgers/
"""
import os, sys, json, time, math, hashlib, copy
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, '<TOOLS>')
import l2_lib as LL                       # sets ctx.cap = 64 at its module top (once)
from flint import arb, acb, ctx
from gravity_engine import argument_count as AC

PREC = 200
SU = 5
J_TM = 24
R_MAIN = 10
N_ARC = 8
R_RULE = [(1.5, 0.2), (100.0, 0.4)]       # TM radius r(|kappa|): 0.2 near 0 keeps the discs off kappa* and sigma = 2 (probe)
CONTROL_C = acb('-1.8', '6.2'); CONTROL_R = '0.5'
TOY_REGIONS = {'T-a': ('-1.5', '0', '0.75', 1), 'T-b': ('-1', '0', '0.3', 0), 'T-c': ('-2', '0', '0.3', 1)}
LEDGER_DIR = os.path.join(HERE, 'ledgers')
FILES = ['l2_cert.py', 'l2_lib.py', 'ax_sys.py']


def sha(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()


def instrument_sha():
    d = {f: sha(os.path.join(HERE, f)) for f in FILES}
    d['gravity_engine/argument_count.py'] = sha('<TOOLS>/gravity_engine/argument_count.py')
    d['build-1/eqs_l.pkl'] = sha(os.path.join(LL.B1, 'eqs_l.pkl'))
    d['build-1/outward_rec.pkl'] = sha(LL.OUTW)
    d['engine/cert/lin_tube2_k1_3.pkl'] = sha(LL.TUBE)
    d['build-1/l0_cert.py'] = sha(os.path.join(LL.B1, 'l0_cert.py'))
    return d


# =========================================================================== closed-polygon count
def polygon_plan(c, rho, rfun_len, nv=16):
    """pieces (a, b, r) of the circumscribed nv-gon of the disc |kappa - c| <= rho (counter-clockwise, closed)"""
    c = complex(c)
    Rv = rho/math.cos(math.pi/nv)*(1 + 1e-12)
    V = [c + Rv*complex(math.cos(2*math.pi*j/nv), math.sin(2*math.pi*j/nv)) for j in range(nv)]
    plan = []
    for j in range(nv):
        a, b = V[j], V[(j + 1) % nv]
        L = abs(b - a)
        n = max(1, int(math.ceil(L/rfun_len)))
        pts = [a + (b - a)*k/n for k in range(n)] + [b]
        for k in range(n):
            plan.append((pts[k], pts[k + 1], abs(pts[k + 1] - pts[k])*0.5*(1 + 1e-6)))
    return plan, V


def check_closed_plan(plan):
    for i in range(len(plan)):
        a2 = plan[(i + 1) % len(plan)][0]
        assert plan[i][1] == a2, 'contour gap between pieces %d and %d' % (i, (i + 1) % len(plan))
        assert abs(plan[i][1] - plan[i][0]) <= 2*plan[i][2], 'piece longer than its disc diameter'
    return True


def count_closed(plan, tm_at, point, ledger=None, header=None, refine=3, log=print, stop_after=None):
    """N = (1/2 pi) Delta arg D around the closed polygon; engine segment_arg + sharpen; 0.25 rule"""
    check_closed_plan(plan)
    old = ctx.prec
    ctx.prec = PREC
    try:
        done = {}
        if ledger and os.path.exists(ledger):
            lines = open(ledger).read().splitlines()
            h0 = json.loads(lines[0]) if lines else {}
            if header is not None and h0 != json.loads(json.dumps(header)):
                raise SystemExit('ledger %s exists and its header differs (refusing to resume)' % ledger)
            for line in lines[1:]:
                d = json.loads(line); done[(d['a'], d['b'])] = d
        elif ledger:
            with open(ledger, 'w') as f:
                f.write(json.dumps(header) + '\n')
        total = arb(0); nref = 0; nnew = 0
        work = [(a, b, r, 0) for (a, b, r) in plan]
        cache = {}

        def pc(z):
            k_ = (z.real.mid().str(30), z.imag.mid().str(30))
            if k_ not in cache:
                cache[k_] = point(z)
            return cache[k_]
        while work:
            a, b, r, lev = work.pop(0)
            key = (repr(a), repr(b))
            if key in done:
                total += arb(done[key]['arg']); continue
            if stop_after is not None and nnew >= stop_after:
                return None
            pa, pb = acb(a.real, a.imag), acb(b.real, b.imag); kc = (pa + pb)/2
            try:
                tm = tm_at(kc, r)
                inc = AC.segment_arg(tm, pa, pb, kc)
            except ArithmeticError:
                inc = None
            if inc is None:
                if lev >= refine:
                    raise ArithmeticError('INCONCLUSIVE: segment %s -> %s not certified at refinement %d' % (a, b, lev))
                m_ = (a + b)/2
                work[0:0] = [(a, m_, r/2, lev + 1), (m_, b, r/2, lev + 1)]
                nref += 1
                continue
            inc = AC.sharpen(inc, lambda: pc(pa), lambda: pc(pb))
            total += inc; nnew += 1
            if ledger:
                with open(ledger, 'a') as f:
                    f.write(json.dumps(dict(a=key[0], b=key[1], r=r, lev=lev, R=arb(tm.R).str(5), arg=inc.str(40))) + '\n')
        w = total/(2*arb.pi())
        n = round(float(w.mid()))
        if not (w - n).abs_upper() < arb('0.25'):
            raise ArithmeticError('Delta/2pi not near an integer: %s' % w)
        return dict(N=n, winding=w.str(15), refinements=nref, pieces=len(plan))
    finally:
        ctx.prec = old


# =========================================================================== builders
def build_axial(log):
    SY = LL.Sys('axial', su=SU)
    P = LL.E5Path(SY, log=log)
    return SY, P


def build_toy(log):
    SYt = LL.Sys('toy', su=1)
    return SYt, LL.ToyPath(SYt)


def toy_counts(tp, log, ledger_tag='toy'):
    toy = LL.Toy(tp)
    tm_at, point = LL.make_callables(toy, J=J_TM)
    out = {}
    for name, (cr, ci, rho, expect) in TOY_REGIONS.items():
        plan, _ = polygon_plan(complex(float(cr), float(ci)), float(rho), 0.6*float(rho)*0.4)
        T0 = time.time()
        res = count_closed(plan, tm_at, point, ledger=os.path.join(LEDGER_DIR, '%s_%s.jsonl' % (ledger_tag, name)),
                           header=dict(obj='flat toy l=2', region=name, c=[cr, ci], rho=rho, J=J_TM), log=log)
        res['expected (RUN.md sec. 0)'] = expect
        res['pass'] = (res['N'] == expect)
        res['time_s'] = round(time.time() - T0, 1)
        out[name] = res
        log('toy %s: N = %d (expected %d) %s' % (name, res['N'], expect, res['winding']))
    return out


def null_counts(P, log, tag='null'):
    lam = P.light['lam']
    kst = 1 - lam
    c = float(kst.mid())
    out = dict(kappa_star=kst.str(15))
    for name, S_set in (('analytic (X1 level, S = {1})', [1]), ('bounded (value at the cone, S = {})', [])):
        ax = LL.Axial(P, S_set=S_set)
        tm_at, point = LL.make_callables(ax, J=J_TM)
        plan, _ = polygon_plan(c, 0.1, 0.03)
        T0 = time.time()
        res = count_closed(plan, tm_at, point, ledger=os.path.join(LEDGER_DIR, '%s_%s.jsonl' % (tag, 'an' if S_set else 'bdd')),
                           header=dict(obj='axial l=2 null', S=S_set, c=c, rho=0.1, J=J_TM), log=log)
        res['time_s'] = round(time.time() - T0, 1)
        out[name] = res
        log('null %s: N = %d %s' % (name, res['N'], res['winding']))
    # control: D_an / D_bdd = 1 - sigma(kappa) at one point of the contour
    z = acb(c + 0.1, 0.05)
    da, cx = LL.Axial(P, S_set=[1]).D_tm(z, arb(0), 0)
    db, _ = LL.Axial(P, S_set=[]).D_tm(z, arb(0), 0)
    sig = (1 - z)/lam
    ratio = cx.evalb(da, acb(0))/cx.evalb(db, acb(0))
    out['control D_an/D_bdd == 1 - sigma at kappa = %s' % z.str(5)] = dict(ratio=ratio.str(12), one_minus_sigma=(1 - sig).str(12),
                                                                          overlaps=bool(ratio.overlaps(1 - sig)))
    an = out['analytic (X1 level, S = {1})']['N']; bd = out['bounded (value at the cone, S = {})']['N']
    out['pass'] = (an == 0) and (bd == -1) and out['control D_an/D_bdd == 1 - sigma at kappa = %s' % z.str(5)]['overlaps']
    return out


# =========================================================================== controls
def failure_mode_tests(P, tp, log):
    out = {}
    # (1) cap truncation: the instrument's arb_series path (L0.lc_root jets) must refuse below the needed cap
    old = ctx.cap
    try:
        ctx.cap = 4
        try:
            LL.L0.lc_root(P.ow, lambda *a: None); out['cap truncation refused (ctx.cap = 4)'] = False
        except ArithmeticError as e:
            out['cap truncation refused (ctx.cap = 4)'] = True
    finally:
        ctx.cap = old
    # (2) contour gap: a plan with one piece removed must be refused (closed counter and the engine's check_plan)
    plan, _ = polygon_plan(-1.5, 0.75, 0.18)
    bad = plan[:3] + plan[4:]
    try:
        check_closed_plan(bad); out['contour gap refused (closed counter)'] = False
    except AssertionError:
        out['contour gap refused (closed counter)'] = True
    rfun = AC.radius_rule(R_RULE, r_edge_re=1e9)
    mp = AC.disc_plan(R_MAIN, rfun, N_ARC, 'halfdisc')
    gap = mp[:5] + mp[6:]
    try:
        AC.check_plan(gap, AC.upper_path(R_MAIN, N_ARC)[0], log=lambda *a: None); out['contour gap refused (engine check_plan)'] = False
    except AssertionError:
        out['contour gap refused (engine check_plan)'] = True
    # (3) pole of the normalisation: the bounded (S = {}) TM on a disc containing kappa* must refuse; analytic is finite
    kst = acb(1 - P.light['lam'])
    try:
        LL.Axial(P, S_set=[]).D_tm(kst, arb('0.05'), 8); out['pole of the normalisation refused (bounded TM on kappa*)'] = False
    except ArithmeticError:
        out['pole of the normalisation refused (bounded TM on kappa*)'] = True
    tm, cx = LL.Axial(P, S_set=[1]).D_tm(kst, arb('0.05'), 8)
    out['analytic TM on kappa* finite'] = bool(cx.sup(tm).is_finite())
    out['analytic D(kappa*) (value)'] = cx.evalb(tm, acb(0)).str(10)
    # (4) zero ON the contour: toy disc whose boundary passes through kappa = -2 must be INCONCLUSIVE
    toy = LL.Toy(tp)
    tm_at, point = LL.make_callables(toy, J=J_TM)
    zplan, V = polygon_plan(-1.5, 0.5*math.cos(math.pi/16)/(1 + 1e-12), 0.12)   # vertex radius 0.5: vertex at -2
    try:
        count_closed(zplan, tm_at, point, refine=2); out['zero on the contour refused (toy vertex at -2)'] = False
    except ArithmeticError:
        out['zero on the contour refused (toy vertex at -2)'] = True
    # (5) light-cone residue: sigma(kappa) = (1 - kappa)/lam and v0 != 0 (normalisation proof inputs)
    cx = LL.Ctx(acb(0), arb(0), 2)
    R0 = P.light['M']
    out['residue R0_12 (kappa-free, nonzero)'] = R0[0][0][0][1].str(12)
    out['residue kappa^1, kappa^2 parts of R0_12 exactly 0'] = all(R0[p][0][0][1] == 0 for p in (1, 2))
    s0 = R0[0][0][0][0] + R0[0][0][1][1]; s1 = R0[1][0][0][0] + R0[1][0][1][1]
    lam = P.light['lam']
    out['tr R0 = sigma: s0 overlaps 1/lam, s1 overlaps -1/lam'] = bool(s0.overlaps(1/lam)) and bool(s1.overlaps(-1/lam))
    out['tr R0 kappa^2 part exactly 0'] = (R0[2][0][0][0] + R0[2][0][1][1]) == 0
    # (6) at kappa* the normalised start z_1 = adj(1 - R0) R_1 v0 != 0  (X1: c2 != 0 there)
    cxs = LL.Ctx(kst, arb(0), 0)
    Y, info = LL.regsing_tm(cxs, P.light, None, acb('1e-30'), S_set=[1])
    v = [cxs.evalb(c, acb(0)) for c in Y]
    out['yhat(t -> 0) at kappa* (= z_0 (1 - sigma) = 0 there, so the solution starts at z_1 t): |yhat(1e-30)|'] = [x.str(5) for x in v]
    ynorm = LL.regsing_tm(cxs, P.light, None, acb('0.001'), S_set=[1])[0]
    vv = [cxs.evalb(c, acb(0)) for c in ynorm]
    out['normalised solution at kappa*, t = 0.001: certified nonzero'] = bool(any(not x.contains(0) for x in vv))
    return out


def resume_test(tp, log):
    toy = LL.Toy(tp)
    tm_at, point = LL.make_callables(toy, J=J_TM)
    plan, _ = polygon_plan(-2.0, 0.3, 0.6*0.3*0.4)
    led = os.path.join(LEDGER_DIR, 'resume_test.jsonl')
    if os.path.exists(led):
        os.remove(led)
    hdr = dict(obj='resume test (toy T-c)', J=J_TM)
    part = count_closed(plan, tm_at, point, ledger=led, header=hdr, stop_after=7)
    n_written = len(open(led).read().splitlines()) - 1
    res = count_closed(plan, tm_at, point, ledger=led, header=hdr)
    led2 = os.path.join(LEDGER_DIR, 'resume_ref.jsonl')
    if os.path.exists(led2):
        os.remove(led2)
    ref = count_closed(plan, tm_at, point, ledger=led2, header=hdr)
    a = [json.loads(l) for l in open(led).read().splitlines()[1:]]
    b = [json.loads(l) for l in open(led2).read().splitlines()[1:]]
    same = len(a) == len(b) and all(x['a'] == y['a'] and x['b'] == y['b'] and arb(x['arg']).overlaps(arb(y['arg'])) for x, y in zip(a, b))
    refused = False
    try:
        count_closed(plan, tm_at, point, ledger=led, header=dict(obj='different header'))
    except SystemExit:
        refused = True
    return dict(stopped_after=n_written, resumed=res, uninterrupted=ref, ledgers_identical=same, header_mismatch_refused=refused,
                resume_tested=bool(part is None and n_written == 7 and res['N'] == ref['N'] and same and refused))


# =========================================================================== main pieces
def main_plan():
    AC.check_arc(R_MAIN, N_ARC, log=lambda *a: None)
    rfun = AC.radius_rule(R_RULE, r_edge_re=1e9)
    plan = AC.disc_plan(R_MAIN, rfun, N_ARC, 'halfdisc')
    AC.check_plan(plan, AC.upper_path(R_MAIN, N_ARC)[0], log=lambda *a: None)
    return plan


def pl2_1(P, log):
    plan = main_plan()
    ax = LL.Axial(P, S_set=[1])
    tm_at, point = LL.make_callables(ax, J=J_TM)
    D0 = point(acb(0)); DR = point(acb(plan[0][1].real, plan[0][1].imag))
    hdr = dict(obj='BAR-SF-L2 PL2-1 axial l=2 half-disc', R=R_MAIN, n_arc=N_ARC, rule=R_RULE, J=J_TM, su=SU,
               lib=sha(os.path.join(HERE, 'l2_lib.py')))
    T0 = time.time()
    n, w = AC.count(plan, D0, DR, tm_at, point, ledger=os.path.join(LEDGER_DIR, 'pl2_1_halfdisc.jsonl'), header=hdr,
                    prec=PREC, log=log)
    return dict(N=n, Delta_over_pi=w.str(15), D0=D0.str(15), D_R=DR.str(15), discs=len(plan), time_s=round(time.time() - T0, 1),
                **{'pass (N == 0, boundary certified zero-free)': n == 0})


def pl2_3(P, log):
    ax = LL.Axial(P, S_set=[1])
    tm_at, point = LL.make_callables(ax, J=J_TM)
    plan, _ = polygon_plan(complex(-1.8, 6.2), 0.5, 0.6*0.5*0.4)
    T0 = time.time()
    res = count_closed(plan, tm_at, point, ledger=os.path.join(LEDGER_DIR, 'pl2_3_control.jsonl'),
                       header=dict(obj='BAR-SF-L2 PL2-3 control disc', c=[-1.8, 6.2], rho=0.5, J=J_TM, su=SU), log=log)
    res['time_s'] = round(time.time() - T0, 1)
    res['pass (N == 1)'] = res['N'] == 1
    if res['N'] == 1:
        # enclosure: Newton on thin point values (float, only to locate), then a certified count 1 on a small disc
        z = complex(-1.8, 6.2)
        for it in range(30):
            tm, cx = ax.D_tm(acb(z.real, z.imag), arb('1e-6'), 2)
            c0 = tm.p.coeffs()[0]; c1 = tm.p.coeffs()[1]
            dz = complex(float((c0/c1).real.mid()), float((c0/c1).imag.mid()))
            z = z - dz
            if abs(dz) < 1e-12:
                break
        rr = 1e-4
        pl2, _ = polygon_plan(z, rr, 0.6*rr)
        enc = count_closed(pl2, tm_at, point, ledger=os.path.join(LEDGER_DIR, 'pl2_3_enclosure.jsonl'),
                           header=dict(obj='PL2-3 enclosure disc', c=[z.real, z.imag], rho=rr, J=J_TM), log=log)
        res['enclosure'] = dict(centre=[repr(z.real), repr(z.imag)], radius='%g (circumscribed 16-gon)' % rr, count=enc,
                                certified=enc['N'] == 1)
    return res


def run(mode, log):
    os.makedirs(LEDGER_DIR, exist_ok=True)
    old = ctx.prec
    try:
        ctx.prec = PREC
        T0 = time.time()
        res = dict(bar='BAR-SF-L2 seal 9077a35520', mode=mode, instrument_sha256=instrument_sha(), prec=PREC, su=SU, J=J_TM)
        SYt, tp = build_toy(log)
        res['toy_symbolic_checks'] = {k: bool(v) for k, v in SYt.check.items()}
        if mode == '--controls':
            for f in os.listdir(LEDGER_DIR):          # controls start from empty ledgers (own files only)
                if f.startswith(('toy_', 'null_', 'resume_')):
                    os.remove(os.path.join(LEDGER_DIR, f))
            log('controls: toy counts (PL2-2 regions)')
            res['PL2-2 toy'] = toy_counts(tp, log)
            res['resume'] = resume_test(tp, log)
            res['resume_tested'] = res['resume']['resume_tested']
            SY, P = build_axial(log)
            res['symbolic_checks'] = {k: bool(v) for k, v in SY.check.items()}
            res['path_info'] = P.info
            res['failure_modes'] = failure_mode_tests(P, tp, log)
            res['PL2-4 null'] = null_counts(P, log)
            res['time_s'] = round(time.time() - T0, 1)
            out = os.path.join(HERE, 'controls_l2.json')
        elif mode == '--probe':
            SY, P = build_axial(log)
            res['setup_s'] = P.info['setup_s']
            plan = main_plan()
            ax = LL.Axial(P, S_set=[1])
            tm_at, point = LL.make_callables(ax, J=J_TM)
            idx = [0, len(plan)//3, len(plan)//2, len(plan) - 8, len(plan) - 1]
            tim = []
            for i in idx:
                e, a, b, r = plan[i]
                kc = acb((a + b).real/2, (a + b).imag/2)
                t1 = time.time(); tm = tm_at(kc, r); dt = time.time() - t1
                t1 = time.time(); p_ = point(kc); dp = time.time() - t1
                tim.append(dict(piece=i, kc=[round(kc.real.mid().__float__(), 4), round(kc.imag.mid().__float__(), 4)], r=r,
                                tm_s=round(dt, 2), point_s=round(dp, 2), R=tm.R.str(4), rel_R=(tm.R/abs(tm.eval(acb(0)))).str(3) if not tm.eval(acb(0)).contains(0) else 'n/a'))
                log('probe piece %d: TM %.2fs point %.2fs R %s' % (i, dt, dp, tm.R.str(4)))
            res['main_plan_discs'] = len(plan)
            res['timings'] = tim
            per = sum(t['tm_s'] + 2*t['point_s'] for t in tim)/len(tim)
            res['projected_main_count_s'] = round(len(plan)*per*1.3, 0)
            res['projected_total_main_s (setup + PL2-1 + PL2-3 (2 x 32 pieces) + toy 3 x 32 + null 2 x 32)'] = round(P.info['setup_s'] + (len(plan) + 64 + 64)*per*1.3 + 130, 0)
            res['time_s'] = round(time.time() - T0, 1)
            out = os.path.join(HERE, 'probe_l2.json')
        elif mode == '--main':
            SY, P = build_axial(log)
            res['symbolic_checks'] = {k: bool(v) for k, v in SY.check.items()}
            res['path_info'] = P.info
            res['PL2-1'] = pl2_1(P, log)
            if res['PL2-1']['N'] >= 1:
                res['KILL'] = 'N >= 1 in the half-disc: possible growing axial mode; STOP and report'
            res['PL2-2'] = toy_counts(tp, log, ledger_tag='main_toy')
            res['PL2-2']['pass'] = all(v['pass'] for v in res['PL2-2'].values() if isinstance(v, dict))
            res['PL2-3'] = pl2_3(P, log)
            res['PL2-4'] = null_counts(P, log, tag='main_null')
            res['verdict_inputs'] = {'PL2-1': res['PL2-1']['pass (N == 0, boundary certified zero-free)'],
                                     'PL2-2': res['PL2-2']['pass'], 'PL2-3 (N == 1)': res['PL2-3']['pass (N == 1)'],
                                     'PL2-4': res['PL2-4']['pass']}
            res['time_s'] = round(time.time() - T0, 1)
            out = os.path.join(HERE, 'result_l2.json')
        else:
            raise SystemExit(__doc__)

        def js(o):
            if isinstance(o, dict):
                return {str(k): js(v) for k, v in o.items()}
            if isinstance(o, (list, tuple)):
                return [js(v) for v in o]
            if isinstance(o, (arb, acb)):
                return o.str(20)
            return o
        json.dump(js(res), open(out, 'w'), indent=1)
        log('wrote %s (%.0fs)' % (out, time.time() - T0))
        return res
    finally:
        ctx.prec = old


if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else ''
    T00 = time.time()
    run(mode, lambda *a: print('[%7.1fs]' % (time.time() - T00), *a, flush=True))
