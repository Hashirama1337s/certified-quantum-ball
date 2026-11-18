"""Certified zero count of an analytic function by the argument principle (contour discs + half-plane test).

In plain words: give it a function D(kappa) with real coefficients (D(conj k) = conj D(k)) and, for any small disc,
an enclosure of D on that disc (a Taylor model). It walks the upper half of a right half-disc contour, proves the
argument change of D piece by piece, and returns the EXACT number of zeros inside, or refuses (INCONCLUSIVE). For a
stability problem, D is the mode determinant and the count is the number of growing modes.

Source (verbatim, entry rule doctrine #121): LIGHTBALL Raid 3 "Ball of Light", published as certified-evans-coleman
v1.0 (DOI 10.5281/zenodo.24303877).
- `sym`, `csym`: engine/cert/core.py l.156-165;
- `segment_arg`, `sharpen`: engine/cert/count_b4.py l.101-140 (pinned, unedited in the run);
- `upper_path`, `seg_dist2`, `check_arc`, `disc_plan`, `check_plan`, the counting loop of `count`:
  engine/cert/count_b4v3_seal.py l.42-108, l.167-238 (REFEREE-B4FW2 fixes F1, F2, F5);
- `radius_rule`: count_b4.py r_of / r_of_k l.45-49, l.72-73.
Changes: module globals (N_ARC, R_RULE, R_RIGHT, R_EDGE_RE, KMAX_WALL, REFINE) became parameters; the
problem-specific engine (Linear2, LinearTM2, F3 precondition checks, iw guards) became the callables `tm_at`, `point`.
Seal: BAR-LB-M b1b0ee1cde, settled HOLD (exactly one growing spherical mode, N = 1 on the half-disc of radius 40);
amendments AMEND-LB-M-1 b92e091b0c. Referees: REFEREE-B4FW, REFEREE-B4FW2, REFEREE-B4TAIL (LIGHTBALL seats/).
CERTIFIED node backed: M count N(R1 = 40) = 1.
Tests: tests/test_argument_count.py (rebuilds the published 99-piece plan key for key and replays the ledger
results/m/runs7/count_R40.jsonl: N = 1, Delta/pi ball containing the logged one; exact-polynomial positive controls,
a left-half-plane null, a zero-on-contour refusal, a plan-gap refusal).

Precision: arb sums run inside `count` at `prec` bits (default 128, the M run's); the caller's ctx.prec is restored
in `finally` (doctrine #115). No python-flint series are used, so ctx.cap does not apply.
"""
import json
import math
import os
from fractions import Fraction

from flint import arb, acb, ctx

MODULE_VERSION = 1
UNIT = arb(0, 1)          # the ball [-1, 1]


def sym(r):
    """the interval [-|r|_up, |r|_up] (rigorous, r any arb/number)"""
    return UNIT*arb(r).abs_upper()


def csym(r):
    """complex box containing the disc |z| <= r"""
    s = sym(r)
    return acb(s, s)


def segment_arg(tm, pa, pb, kc, max_depth=30):
    """certified arg(D(pb)/D(pa)) inside one TM disc (half-plane test, lin_tm.winding_tm logic)"""
    total = None
    stack = [(pa - kc, pb - kc, 0)]
    while stack:
        u, v, dep = stack.pop()
        Du = tm.eval(u); Dv = tm.eval(v)
        cmid = (u + v)/2; rad = abs(v - u).upper()/2
        E = tm.eval(cmid + csym(rad))
        if (E*Du.conjugate()).real > 0 and (E*Dv.conjugate()).real > 0:
            inc = (Dv/Du).arg()
            total = inc if total is None else total + inc
        else:
            if dep >= max_depth: return None
            mid = (u + v)/2
            stack.append((mid, v, dep + 1)); stack.append((u, mid, dep + 1))
    return total


def sharpen(inc, Da_f, Db_f):
    """TM arg increments are wide balls (N2 smoke: Delta/pi = 2 +- 0.67 after 14 discs).  The true
    increment is Arg(D(pb)/D(pa)) + 2 pi k for the unique integer k that lands in the TM ball
    (unique when rad(inc) < 0.9 pi); D(pa), D(pb) are thin point values in the TM's normalisation.
    Returns the sharpened ball, or inc unchanged if any condition fails (never wider than inc)."""
    pi = arb.pi()
    if not inc.rad() < pi*arb('0.9'):
        return inc
    p = (Db_f()/Da_f()).arg()
    if not p.rad() < arb('0.1'):
        return inc
    kk = round(float(((inc - p)/(2*pi)).mid()))
    cand = p + 2*pi*kk
    if cand.overlaps(inc) and (cand - inc).abs_upper() < pi and cand.rad() < inc.rad():
        return cand
    return inc


def radius_rule(r_rule, r_edge_re=5.0, r_right=0.3):
    """r(kappa): r_rule = [(|kappa| limit, r), ...] left of Re kappa = r_edge_re, r_right on the right edge"""
    def r_of(absk):
        for lim, r in r_rule:
            if absk <= lim: return r
        return r_rule[-1][1]
    return lambda z: r_of(abs(z)) if z.real < r_edge_re else r_right


def upper_path(R1, n_arc=8, shape='halfdisc'):
    """upper half of the contour, from the positive real axis to 0 (box: R1 -> R1+iR1 -> iR1 -> 0; halfdisc:
    circumscribed polygon of the quarter circle, vertex radius R1/cos(dphi/2)*(1 + 1e-12), then iR' -> 0)"""
    if shape == 'box':
        return [complex(R1, 0), complex(R1, R1), complex(0, R1), complex(0, 0)]
    N = n_arc
    dphi = (math.pi/2)/N
    Rp = R1/math.cos(dphi/2)*(1 + 1e-12)
    pts = [Rp*complex(math.cos(j*dphi), math.sin(j*dphi)) for j in range(N + 1)]
    pts[0] = complex(Rp, 0); pts[-1] = complex(0, Rp)
    return pts + [complex(0, 0)]


def seg_dist2(a, b):
    """exact squared distance from 0 to the segment [a, b] (a, b pairs of Fractions)"""
    dx, dy = b[0] - a[0], b[1] - a[1]
    L2 = dx*dx + dy*dy
    t = -(a[0]*dx + a[1]*dy)/L2
    if t <= 0: return a[0]**2 + a[1]**2
    if t >= 1: return b[0]**2 + b[1]**2
    px, py = a[0] + t*dx, a[1] + t*dy
    return px*px + py*py


def check_arc(R1, n_arc=8, kmax_wall=None, log=print):
    """F1: EXACT (Fraction) check that every chord stays at |kappa| >= R1 and (optionally) every vertex <= kmax_wall"""
    pts = upper_path(R1, n_arc)
    F = [(Fraction(p.real), Fraction(p.imag)) for p in pts[:-1]]     # vertices (exact values of the floats used)
    R2 = Fraction(R1)**2
    dmin = min(seg_dist2(F[i], F[i + 1]) for i in range(len(F) - 1))
    vmax = max(x*x + y*y for x, y in F)
    assert dmin >= R2, 'F1: a chord dips below R1 (min |k|^2 - R1^2 = %s)' % float(dmin - R2)
    assert kmax_wall is None or vmax <= Fraction(kmax_wall)**2, 'F1: a vertex lies beyond the outer wall'
    log('[F1] arc: %d chords, vertex radius %.15g; EXACT: min chord distance^2 - R1^2 = %.3e >= 0, max vertex |k| = %.10f' % (
        len(F) - 1, abs(pts[0]), float(dmin - R2), math.sqrt(float(vmax))))
    return pts


def disc_plan(R1, rfun, n_arc=8, shape='halfdisc'):
    """[(edge, a, b, r)]: pieces of the upper path, each covered by one disc of radius r; last piece of an edge snapped"""
    pts = upper_path(R1, n_arc, shape)
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


def check_plan(plan, start, log=print):
    """F1: the plan is exactly continuous from `start` to 0 and every piece fits inside its disc"""
    for i in range(len(plan) - 1):
        assert plan[i][2] == plan[i + 1][1], 'F1: plan gap between pieces %d and %d' % (i, i + 1)
    assert plan[0][1] == start, 'F1: plan does not start at the first vertex'
    assert plan[-1][2] == 0j, 'F1: plan does not end at 0'
    for (e, a, b, r) in plan:
        assert abs(b - a) <= 2*r*(1 - 1e-9), 'F1: piece longer than the disc diameter'
    log('[F1] plan: %d pieces, exactly continuous, starts at the vertex %r, ends at 0; every piece <= 2r' % (len(plan), start))


def count(plan, D0, DR, tm_at, point, ledger=None, header=None, refine=3, prec=128, log=print):
    """N = (1/pi) Delta_arg D along the upper path = number of zeros in the right half-disc (D real on the real axis).

    tm_at(kc: acb, r) -> object with .eval(u: acb) enclosing D(kc + u) for |u| <= r and .R (remainder, arb);
    point(z: acb) -> thin value of D in the TM's normalisation (used only to sharpen arg balls);
    D0, DR: certified values D(0), D(plan[0].a) (acb or arb; their real parts must have certified signs).
    ledger: JSONL path; line 1 is `header` (resume is refused if it differs), then one line per disc.
    Raises ArithmeticError for INCONCLUSIVE (a piece not certified after `refine` halvings), a Delta/pi not within
    0.25 of an integer, or a failed parity check (-1)^N = sign D(0) sign D(R)."""
    old = ctx.prec
    ctx.prec = prec
    try:
        done = {}
        if ledger and os.path.exists(ledger):
            lines = open(ledger).read().splitlines()
            h0 = json.loads(lines[0]) if lines else {}
            if header is not None and h0 != json.loads(json.dumps(header)):
                raise SystemExit('F5: ledger %s exists and its header differs from this run (refusing to resume)' % ledger)
            for line in lines[1:]:
                d = json.loads(line); done[(d['a'], d['b'])] = d
            log('[F5] resuming %s: %d discs on record' % (ledger, len(done)))
        elif ledger:
            with open(ledger, 'w') as f: f.write(json.dumps(header if header is not None else {'header': True}) + '\n')
        re0, reR = acb(D0).real, acb(DR).real
        if not ((re0 > 0 or re0 < 0) and (reR > 0 or reR < 0)):
            raise ArithmeticError('real endpoint sign not certified')
        total = arb(0)
        _pc = {}
        def ptcache(z):
            k_ = (z.real.mid().str(30), z.imag.mid().str(30))
            if k_ not in _pc: _pc[k_] = point(z)
            return _pc[k_]
        work = [(a, b, r, 0) for (_, a, b, r) in plan]
        nref = 0
        while work:
            a, b, r, lev = work.pop(0)
            key = (repr(a), repr(b))
            if key in done:
                total += arb(done[key]['arg']); continue
            pa, pb = acb(a.real, a.imag), acb(b.real, b.imag); kc = (pa + pb)/2
            try:
                tm = tm_at(kc, r)
                inc = segment_arg(tm, pa, pb, kc)
            except ArithmeticError:
                inc = None
            if inc is None:
                if lev >= refine:
                    raise ArithmeticError('INCONCLUSIVE: segment %s -> %s not certified at refinement %d' % (a, b, lev))
                m_ = (a + b)/2
                work[0:0] = [(a, m_, r/2, lev + 1), (m_, b, r/2, lev + 1)]
                nref += 1
                continue
            inc = sharpen(inc, lambda: ptcache(pa), lambda: ptcache(pb))
            total += inc
            if ledger:
                rec = dict(a=key[0], b=key[1], r=r, lev=lev, R=arb(tm.R).str(5), arg=inc.str(40))
                with open(ledger, 'a') as f: f.write(json.dumps(rec) + '\n')
        w = total/arb.pi()
        n = round(float(w.mid()))
        if not (w - n).abs_upper() < arb('0.25'):
            raise ArithmeticError('Delta/pi not near an integer: %s' % w)
        parity = (1 if (re0 > 0) == (reR > 0) else -1)
        log('N = %d   (Delta/pi = %s; parity (-1)^N = sign D(0) sign D(R): %s; refinements %d)' % (
            n, w.str(10), (-1)**n == parity, nref))
        if (-1)**n != parity:
            raise ArithmeticError('F2: parity check FAILED: N = %d but sign D(0) sign D(R) = %d' % (n, parity))
        return n, w
    finally:
        ctx.prec = old
