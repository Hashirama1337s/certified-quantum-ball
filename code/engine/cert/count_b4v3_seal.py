"""Sealing fixes for the v3 count (speedup seat, NEW file; count_b4.py stays unedited, its functions are
replaced in the running process only).  Implements REFEREE-B4FW2.md sec. 2 / "Blocks sealing" 1, 2, 3, 5:

  F1  arc polygon: vertex radius Rp = R1/cos(dphi/2)*(1 + 1e-12); asserted in EXACT rational arithmetic
      (fractions.Fraction of the float vertices actually used) that every chord's distance from 0 is >= R1
      and every vertex has |kappa| <= KMAX_WALL (55).  Plan continuity asserted exactly: plan[i].b == plan[i+1].a,
      first a == the first vertex (Rp, 0), last b == 0.  The disc plan is count_b4.disc_plan_path's rule with
      the last piece of every edge snapped to the edge end (so no sub-1e-12 piece is left uncovered).
  F2  parity: (-1)**N != sign D(0) sign D(Rp) RAISES.
  F3  certified preconditions, ball evaluations, printed:
        c_V(x_m) != 0      (Thm 1.6, PROOF-B4.md:103-106, :274): Irow[3] on the certified tube end box at xbar;
        detA root < 0      (detA affine in kappa: quadratic coefficient det[a, e_A-col..] is structurally 0;
                            root = -dA(0)/(dA(1)-dA(0)), also the second difference is printed);
        iw == 2            at every plan point (centres and endpoints) AND at run time in every sonic_tm,
                            detA and sonic_start call (guards; a 3 raises);
        Abar_ss'(0) != 0   (L.abar_x0 ball).
  F5  ledger: a header line (engine sha256s, rule, settings) is written to a FRESH ledger and asserted on
      resume; LB_R_FIXED and LB_MAXDISCS are refused in count mode; everything is printed.
The counting loop itself is count_b4.count's (copied, with only F2/F5 added); segment_arg, sharpen and
REFINE are count_b4's objects.
"""
import os, sys, json, math, time, hashlib
from fractions import Fraction
import count_b4 as CB

KMAX_WALL = 55
ENGINE_FILES = ['engine/cert/count_b4.py', 'engine/cert/count_b4v3.py', 'engine/cert/count_b4v3_seal.py',
                'engine/cert/lin_tm_v3.py', 'engine/cert/lin_tm_v2.py', 'engine/cert/lin_tm.py',
                'engine/cert/lin_stages_v2.py', 'engine/cert/lin_stages.py', 'engine/cert/core.py',
                'engine/cert/lin_fields2.pkl', 'engine/cert/lin_tube2_k1_3.pkl', 'b4speed/cacheL.py', 'b4speed/L13_nsub2.pkl']
ROOT = '<LAB>/'


def sha(fn):
    h = hashlib.sha256()
    with open(fn, 'rb') as f:
        for blk in iter(lambda: f.read(1 << 20), b''): h.update(blk)
    return h.hexdigest()


# ------------------------------------------------------------------ F1 contour
def upper_path_v3(R1, shape=None):
    shape = shape or CB.SHAPE
    if shape == 'box':
        return [complex(R1, 0), complex(R1, R1), complex(0, R1), complex(0, 0)]
    N = CB.N_ARC
    dphi = (math.pi/2)/N
    Rp = R1/math.cos(dphi/2)*(1 + 1e-12)
    pts = [Rp*complex(math.cos(j*dphi), math.sin(j*dphi)) for j in range(N + 1)]
    pts[0] = complex(Rp, 0); pts[-1] = complex(0, Rp)
    return pts + [complex(0, 0)]


def _seg_dist2(a, b):
    """exact squared distance from 0 to the segment [a, b] (a, b pairs of Fractions)"""
    dx, dy = b[0] - a[0], b[1] - a[1]
    L2 = dx*dx + dy*dy
    t = -(a[0]*dx + a[1]*dy)/L2
    if t <= 0: return a[0]**2 + a[1]**2
    if t >= 1: return b[0]**2 + b[1]**2
    px, py = a[0] + t*dx, a[1] + t*dy
    return px*px + py*py


def check_arc(R1, log):
    pts = upper_path_v3(R1)
    if CB.SHAPE == 'box':
        log('[seal F1] shape box: arc checks not applicable'); return pts
    F = [(Fraction(p.real), Fraction(p.imag)) for p in pts[:-1]]     # vertices (exact values of the floats used)
    R2 = Fraction(R1)**2; W2 = Fraction(KMAX_WALL)**2
    dmin = min(_seg_dist2(F[i], F[i + 1]) for i in range(len(F) - 1))
    vmax = max(x*x + y*y for x, y in F)
    assert dmin >= R2, 'F1: a chord dips below R1 (min |k|^2 - R1^2 = %s)' % float(dmin - R2)
    assert vmax <= W2, 'F1: a vertex lies beyond the Far Wall band'
    log('[seal F1] arc: %d chords, vertex radius %.15g; EXACT: min chord distance^2 - R1^2 = %.3e >= 0, max vertex |k| = %.10f <= %d' % (
        len(F) - 1, abs(pts[0]), float(dmin - R2), math.sqrt(float(vmax)), KMAX_WALL))
    return pts


def disc_plan_v3(R1, shape=None):
    pts = upper_path_v3(R1, shape)
    rfun = CB.r_of_k
    plan = []
    for e in range(len(pts) - 1):
        a, b = pts[e], pts[e + 1]
        L = abs(b - a); t = 0.0
        while True:
            z = a + (b - a)*(t/L)
            r = min(rfun(z), rfun(a + (b - a)*(min(L, t + 1.6*rfun(z))/L)))
            dt = 1.6*r
            if t + dt >= L - 1e-9:                      # last piece of the edge: snap to the edge end exactly
                plan.append((e, z, b, r)); break
            plan.append((e, z, a + (b - a)*((t + dt)/L), r))
            t += dt
    return plan


def check_plan(plan, R1, log):
    pts = upper_path_v3(R1)
    for i in range(len(plan) - 1):
        assert plan[i][2] == plan[i + 1][1], 'F1: plan gap between pieces %d and %d' % (i, i + 1)
    assert plan[0][1] == pts[0], 'F1: plan does not start at the first vertex'
    assert plan[-1][2] == 0j, 'F1: plan does not end at 0'
    for (e, a, b, r) in plan:
        assert abs(b - a) <= 2*r*(1 - 1e-9), 'F1: piece longer than the disc diameter'
    log('[seal F1] plan: %d pieces, exactly continuous, starts at the vertex %r, ends at 0; every piece <= 2r' % (len(plan), pts[0]))


# ------------------------------------------------------------------ F3 preconditions
def iw_of(L, kc):
    from flint import acb
    RS = L.sonic; kc = acb(kc)
    M0c = [[acb(RS.M0[0][i][j]) + kc*acb(RS.M1[0][i][j]) for j in range(4)] for i in range(4)]
    return max((2, 3), key=lambda i: float(sum(abs(x).mid() for x in M0c[i])))


def install_iw_guards(L, TMe, log):
    import lin_tm as LT
    stats = {'calls': 0}
    orig_tm = TMe.sonic_tm; orig_dA = TMe.detA; orig_ss = L.sonic_start
    def g(kc):
        stats['calls'] += 1
        w = iw_of(L, kc)
        if w != 2: raise AssertionError('F3: iw = %d at kappa = %s' % (w, kc))
    def sonic_tm(kc, r, a=0, b=1):
        g(kc); return orig_tm(kc, r, a, b)
    def detA(kc):
        g(kc); return orig_dA(kc)
    def sonic_start(kap, a=0, b=1):
        g(kap); return orig_ss(kap, a, b)
    TMe.sonic_tm = sonic_tm; TMe.detA = detA; L.sonic_start = sonic_start
    return stats


def certify_preconditions(L, TMe, plan, log):
    from flint import arb, acb
    from lin_stages import deser
    # c_V(x_m): Irow[3] on the certified tube end box (N, om, V) at x = xbar (lin_tube_v2: last step lands on XBAR)
    end = [deser(p) for p in L.tb['end']]
    xb = deser(L.tb['XBAR'])
    from core import sym_const
    I = L.LF.irow.eval(end, conv=sym_const)
    cV = I[3]
    log('[seal F3] c_V(x_m) on the tube end box at xbar = %s: %s  (nonzero: %s)' % (xb.str(10), cV.str(15), not cV.contains(0)))
    assert not cV.contains(0), 'F3: c_V(x_m) not certified nonzero'
    # Abar_ss'(0)
    log('[seal F3] Abar_ss\'(0) = %s  (nonzero: %s)' % (L.abar_x0.str(15), not L.abar_x0.contains(0)))
    assert not L.abar_x0.contains(0)
    # detA root
    d0 = TMe.detA(acb(0)); d1 = TMe.detA(acb(1)); d2 = TMe.detA(acb(2))
    sd = d2 - 2*d1 + d0
    root = -d0/(d1 - d0)
    log('[seal F3] detA(0) = %s, detA(1) = %s, second difference %s; affine root = %s' % (d0.str(10), d1.str(10), sd.str(3), root.str(12)))
    assert sd.real.contains(0) and sd.imag.contains(0), 'F3: detA not affine'
    assert root.real < 0 and root.imag.contains(0), 'F3: detA root not certified negative'
    # iw at every plan point
    pts = set()
    for (e, a, b, r) in plan:
        pts.add(a); pts.add(b); pts.add((a + b)/2)
    bad = [z for z in pts if iw_of(L, z) != 2]
    log('[seal F3] iw == 2 at all %d plan points (endpoints + centres): %s; run-time guards installed on sonic_tm, detA, sonic_start' % (len(pts), not bad))
    assert not bad, 'F3: iw != 2 at %s' % bad[:3]
    return dict(cV=cV.str(15), root=root.str(12), abar=L.abar_x0.str(15))


# ------------------------------------------------------------------ count (count_b4.count + F2, F5)
def count_v3(R1, ledger, log, header_extra):
    from flint import arb, acb
    import sympy as sp
    import lin_stages_v2 as LS2, lin_tm_v2 as LT2
    for v in ('LB_R_FIXED', 'LB_MAXDISCS'):
        if os.environ.get(v): raise SystemExit('F5: %s is set; refused in count mode' % v)
    check_arc(R1, log)
    plan = disc_plan_v3(R1)
    check_plan(plan, R1, log)
    header = dict(header=True, R1=R1, shape=CB.SHAPE, N_ARC=CB.N_ARC, R_RULE=CB.R_RULE, R_RIGHT=CB.R_RIGHT, R_EDGE_RE=CB.R_EDGE_RE,
                  M=CB.M_DEG, J=CB.J_DEG, NSUB=CB.NSUB, TOL=CB.TOL, REFINE=CB.REFINE, ndisc_plan=len(plan),
                  engine={f: sha(ROOT + f) for f in ENGINE_FILES})
    header.update(header_extra)
    log('[seal F5] header: %s' % json.dumps(header))
    done = {}
    if os.path.exists(ledger):
        lines = open(ledger).read().splitlines()
        h0 = json.loads(lines[0]) if lines else {}
        if h0 != json.loads(json.dumps(header)):
            raise SystemExit('F5: ledger %s exists and its header differs from this run (refusing to resume)' % ledger)
        for line in lines[1:]:
            d = json.loads(line); done[(d['a'], d['b'])] = d
        log('[seal F5] resuming %s: header matches, %d discs on record' % (ledger, len(done)))
    else:
        with open(ledger, 'w') as f: f.write(json.dumps(header) + '\n')
        log('[seal F5] fresh ledger %s with header' % ledger)
    L = LS2.Linear2(sp.Rational(1, 3), log=log, nsub=CB.NSUB, tol=CB.TOL)
    TMe = LT2.LinearTM2(L, m=CB.M_DEG, rem_method='c1', J=CB.J_DEG)
    pre = certify_preconditions(L, TMe, plan, log)
    guards = install_iw_guards(L, TMe, log)
    Rs = plan[0][1].real
    D0 = L.D(acb(0)); DR = L.D(acb(arb(repr(Rs))))
    if not ((D0.real > 0 or D0.real < 0) and (DR.real > 0 or DR.real < 0)):
        raise ArithmeticError('real endpoint sign not certified')
    log('D(0) = %s   D(Rp=%r) = %s' % (D0.str(6), Rs, DR.str(6)))
    total = arb(0)
    _pc = {}
    def ptcache(z):
        k_ = (z.real.mid().str(30), z.imag.mid().str(30))
        if k_ not in _pc: _pc[k_] = TMe.detA(z)*L.D(z)
        return _pc[k_]
    work = [(a, b, r, 0) for (_, a, b, r) in plan]
    nref = 0
    while work:
        a, b, r, lev = work.pop(0)
        key = (repr(a), repr(b))
        if key in done:
            total += arb(done[key]['arg']); continue
        pa, pb = acb(a.real, a.imag), acb(b.real, b.imag); kc = (pa + pb)/2
        t0 = time.time()
        try:
            tm = TMe.D_tm(kc, r)
            inc = CB.segment_arg(tm, pa, pb, kc)
        except ArithmeticError:
            inc = None
        if inc is None:
            if lev >= CB.REFINE:
                raise ArithmeticError('INCONCLUSIVE: segment %s -> %s not certified at refinement %d' % (a, b, lev))
            m_ = (a + b)/2
            work[0:0] = [(a, m_, r/2, lev + 1), (m_, b, r/2, lev + 1)]
            nref += 1
            log('refine %s -> %s to r=%g' % (a, b, r/2)); continue
        inc = CB.sharpen(inc, lambda: ptcache(pa), lambda: ptcache(pb))
        total += inc
        rec = dict(a=key[0], b=key[1], r=r, lev=lev, R=tm.R.str(5), arg=inc.str(40), sec=round(time.time() - t0, 1))
        with open(ledger, 'a') as f: f.write(json.dumps(rec) + '\n')
        log('disc %s r=%g  R=%s  arg %s  (%.0fs)' % (kc.str(5), r, tm.R.str(3), inc.str(8), time.time() - t0))
    w = total/arb.pi()
    n = round(float(w.mid()))
    if not (w - n).abs_upper() < arb('0.25'):
        raise ArithmeticError('Delta/pi not near an integer: %s' % w)
    parity = (1 if (D0.real > 0) == (DR.real > 0) else -1)
    log('N(R1=%g) = %d   (Delta/pi = %s; parity (-1)^N = sign D(0) sign D(Rp): %s; refinements %d; iw guard calls %d)' % (
        R1, n, w.str(10), (-1)**n == parity, nref, guards['calls']))
    if (-1)**n != parity:
        raise ArithmeticError('F2: parity check FAILED: N = %d but sign D(0) sign D(Rp) = %d' % (n, parity))
    log('[seal] preconditions: c_V(x_m) = %s, detA root = %s, Abar_ss\'(0) = %s' % (pre['cV'], pre['root'], pre['abar']))
    return n
