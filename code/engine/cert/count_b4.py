"""B4 counting driver (build seat 2).  WRITTEN, NOT RUN IN FULL.  Proposal input for BAR-LB-M.

Counts the zeros of D_new(kappa) = det[s_kappa(xbar), rho1, rho2, e_V] (data (Abar,Nbar)(0) = (0,1),
gauge mode exactly deflated, PROOF-B4 Lemma 1.5) in the box  B(R1) = [0, R1] x [-R1, R1].

Contour reduction (exact, no numerics):  the HKA pencil has real coefficients, so D(conj k) = conj D(k)
(lin_tm T0: identical values at +-40i, 40+-40i, b4t0b.txt).  The counter-clockwise boundary of B(R1) is
the upper path G+ = [R1 -> R1+iR1 -> iR1 -> 0] plus its mirror image traversed in the same sense, so
   N(R1) = (1/2pi) Delta_arg(dB) = (1/pi) Delta_arg(G+),
with D(0), D(R1) real and certified nonzero (their signs fix Delta_arg(G+) mod 2pi; consistency check).
A single cross-check:  (-1)^N = sign D(0) sign D(R1)  (real zeros in (0,R1) counted with parity).

Certification of Delta_arg(G+): G+ is covered by discs (centres on G+, spacing 1.6 r, r = r(|kappa|)
from a pre-registered radius rule).  On each disc a kappa-Taylor model (lin_tm_v2.LinearTM2, M1 + M2)
encloses D on the whole disc; each sub-segment passes the half-plane test (lin_tm.winding_tm) or is bisected
(depth <= 30).  A disc whose TM remainder makes the test fail is replaced by two discs of radius r/2 on the
same edge piece (at most REFINE levels, pre-registered); still failing => INCONCLUSIVE, no claim.
Every disc is appended to a JSONL ledger (centre, r, remainder, arg increment as a ball string), so the
run is resumable and auditable.  The count is the nearest integer to Delta/pi, required within 0.25.

kappa*: real bisection of D_new on [2.7, 2.9] (point evaluations; sign change certified), then gamma = 1/kappa*.
If N(R1) = 1, the zero in B(R1) is unique, simple, and real (non-real zeros come in conjugate pairs), so it
is the bisected kappa*.

usage (from <LAB>, every run through memcap):
  py -3.11 <MEMCAP> 8 py -3.11 -u engine/cert/count_b4.py count  <R1> <ledger.jsonl>
  py -3.11 <MEMCAP> 8 py -3.11 -u engine/cert/count_b4.py kstar
  py -3.11 <MEMCAP> 8 py -3.11 -u engine/cert/count_b4.py plan   <R1>     (disc plan + cost only)
"""
import sys, os, json, math, time
sys.path.insert(0, '<LAB>/engine/cert'); sys.path.insert(0, '<LAB>/engine')

M_DEG = 16; J_DEG = 4; NSUB = 2; TOL = 1e-20     # the settings of the RT smokes (b4build2/s3_rt_fast.log)
REFINE = 3
SEC_PER_TM = {5: 72.0, 10: 67.0, 20: 79.0, 40: 74.0}      # measured, k=1/3, one process (b4build2/s3_rt_fast.log)

# pre-registered radius rule r(|kappa|) (to be frozen in BAR-LB-M from the N1-N3 / RT smokes)
# near the imaginary axis |D| ~ 1.7/|kappa|^2 is small and the remainder ratio grows with |kappa|
# (RT smokes: 5i r=0.2 3.4e-4; 20i r=0.1 1.7e-2; 40i r=0.1 2.4); away from it D grows and r can be larger
# (40+40i r=0.1: 8e-8).  Values beyond 40 are EXTRAPOLATED.
R_RULE = [(5.0, 0.25), (10.0, 0.15), (20.0, 0.1), (40.0, 0.05), (90.0, 0.025)]
R_EDGE_RE = 5.0      # Re kappa >= R_EDGE_RE: right/top edge rule
R_RIGHT = 0.3

def r_of(absk):
    for lim, r in R_RULE:
        if absk <= lim: return r
    return R_RULE[-1][1]

SHAPE = os.environ.get('LB_SHAPE', 'halfdisc')    # 'halfdisc' (default, matches the Far Wall) or 'box'
N_ARC = 64                                        # chords of the circumscribed arc polygon

def upper_path(R1, shape=None):
    """upper half of the contour, from the positive real axis to 0.
    box:      R1 -> R1+iR1 -> iR1 -> 0.
    halfdisc: CIRCUMSCRIBED polygon of the quarter circle |kappa| = R1 (vertices at radius
              R1/cos(pi/(4 N_ARC)), so every chord stays at |kappa| >= R1, inside the Far Wall's
              certified zero-free annulus R0 = 55 <= |kappa| <= 2000, PROOF-B4-v4 / lead note), then iR' -> 0."""
    shape = shape or SHAPE
    if shape == 'box':
        return [complex(R1, 0), complex(R1, R1), complex(0, R1), complex(0, 0)]
    dphi = (math.pi/2)/N_ARC
    Rp = R1/math.cos(dphi/2)
    pts = [Rp*complex(math.cos(j*dphi), math.sin(j*dphi)) for j in range(N_ARC + 1)]
    pts[0] = complex(Rp, 0); pts[-1] = complex(0, Rp)
    return pts + [complex(0, 0)]

def disc_plan(R1, shape=None):
    """list of (edge, a, b, r) sub-segments of the upper path, each covered by one disc of radius r"""
    return disc_plan_path(upper_path(R1, shape), lambda z: r_of_k(z))

def r_of_k(z):
    return r_of(abs(z)) if z.real < R_EDGE_RE else R_RIGHT

def disc_plan_path(pts, rfun):
    plan = []
    for e in range(len(pts) - 1):
        a, b = pts[e], pts[e + 1]
        L = abs(b - a); t = 0.0
        while t < L - 1e-12:
            z = a + (b - a)*(t/L)
            r = min(rfun(z), rfun(a + (b - a)*(min(L, t + 1.6*rfun(z))/L)))
            dt = min(1.6*r, L - t)
            plan.append((e, a + (b - a)*(t/L), a + (b - a)*((t + dt)/L), r))
            t += dt
    return plan

def _old_disc_plan(R1):
    pts = upper_path(R1); plan = []
    for e in range(3):
        a, b = pts[e], pts[e + 1]
        L = abs(b - a); t = 0.0
        while t < L - 1e-12:
            z = a + (b - a)*(t/L)
            r = r_of(abs(z) + 1.6*r_of(abs(z)))
            dt = min(1.6*r, L - t)
            plan.append((e, a + (b - a)*(t/L), a + (b - a)*((t + dt)/L), r))
            t += dt
    return plan

def segment_arg(tm, pa, pb, kc, max_depth=30):
    """certified arg(D(pb)/D(pa)) inside one TM disc (half-plane test, lin_tm.winding_tm logic)"""
    from flint import acb
    from core import csym
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

R_FIXED = float(os.environ['LB_R_FIXED']) if os.environ.get('LB_R_FIXED') else None   # smokes only
MAXDISCS = int(os.environ.get('LB_MAXDISCS', '0'))                                     # smokes only (split a run)

def sharpen(inc, Da_f, Db_f):
    """TM arg increments are wide balls (N2 smoke: Delta/pi = 2 +- 0.67 after 14 discs).  The true
    increment is Arg(D(pb)/D(pa)) + 2 pi k for the unique integer k that lands in the TM ball
    (unique when rad(inc) < 0.9 pi); D(pa), D(pb) are thin point values in the TM's normalisation.
    Returns the sharpened ball, or inc unchanged if any condition fails (never wider than inc)."""
    from flint import arb
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

def count(R1, ledger, log):
    from flint import arb, acb
    import sympy as sp
    import lin_stages_v2 as LS2, lin_tm_v2 as LT2
    L = LS2.Linear2(sp.Rational(1, 3), log=log, nsub=NSUB, tol=TOL)
    TMe = LT2.LinearTM2(L, m=M_DEG, rem_method='c1', J=J_DEG)
    done = {}
    if os.path.exists(ledger):
        for line in open(ledger):
            d = json.loads(line); done[(d['a'], d['b'])] = d
    # real-axis endpoints
    Rs = upper_path(R1)[0].real
    D0 = L.D(acb(0)); DR = L.D(acb(arb(repr(Rs))))
    if not ((D0.real > 0 or D0.real < 0) and (DR.real > 0 or DR.real < 0)):
        raise ArithmeticError('real endpoint sign not certified')
    log('D(0) = %s   D(R1) = %s' % (D0.str(6), DR.str(6)))
    total = arb(0)
    _pc = {}
    def ptcache(z):
        k_ = (z.real.mid().str(30), z.imag.mid().str(30))
        if k_ not in _pc: _pc[k_] = TMe.detA(z)*L.D(z)     # same pole-free normalisation as the TM (lin_tm.py:196-198)
        return _pc[k_]
    plan = disc_plan(R1) if R_FIXED is None else disc_plan_path(upper_path(R1), lambda z: R_FIXED)
    work = [(a, b, r, 0) for (_, a, b, r) in plan]
    nd_new = 0
    while work:
        a, b, r, lev = work.pop(0)
        key = (repr(a), repr(b))
        if key in done:
            total += arb(done[key]['arg']); continue
        pa, pb = acb(a.real, a.imag), acb(b.real, b.imag); kc = (pa + pb)/2
        t0 = time.time()
        try:
            tm = TMe.D_tm(kc, r)
            inc = segment_arg(tm, pa, pb, kc)
        except ArithmeticError:
            inc = None
        if inc is None:
            if lev >= REFINE:
                raise ArithmeticError('INCONCLUSIVE: segment %s -> %s not certified at refinement %d' % (a, b, lev))
            m_ = (a + b)/2
            work[0:0] = [(a, m_, r/2, lev + 1), (m_, b, r/2, lev + 1)]
            log('refine %s -> %s to r=%g' % (a, b, r/2)); continue
        inc = sharpen(inc, lambda: ptcache(pa), lambda: ptcache(pb))
        total += inc
        nd_new += 1
        rec = dict(a=key[0], b=key[1], r=r, lev=lev, R=tm.R.str(5), arg=inc.str(40), sec=round(time.time() - t0, 1))
        with open(ledger, 'a') as f: f.write(json.dumps(rec) + '\n')
        log('disc %s r=%g  R=%s  arg %s  (%.0fs)' % (kc.str(5), r, tm.R.str(3), inc.str(8), time.time() - t0))
        if MAXDISCS and nd_new >= MAXDISCS and work:     # round-2 fix: pause only AFTER the ledger write
            log('PAUSED after %d new discs (LB_MAXDISCS); rerun the same command to resume from the ledger' % nd_new)
            return None
    w = total/arb.pi()
    n = round(float(w.mid()))
    if not (w - n).abs_upper() < arb('0.25'):
        raise ArithmeticError('Delta/pi not near an integer: %s' % w)
    parity = (1 if (D0.real > 0) == (DR.real > 0) else -1)
    log('N(R1=%g) = %d   (Delta/pi = %s; parity check (-1)^N = sign D(0) sign D(R1): %s)' % (R1, n, w.str(10), (-1)**n == parity))
    return n

def kstar(log):
    from flint import arb, acb
    import sympy as sp
    from core import sym
    import lin_stages_v2 as LS2
    L = LS2.Linear2(sp.Rational(1, 3), log=log, nsub=NSUB, tol=TOL)
    lo, hi = arb('2.7'), arb('2.9')
    flo = L.D(acb(lo)).real; fhi = L.D(acb(hi)).real
    assert (flo > 0 and fhi < 0) or (flo < 0 and fhi > 0), 'no certified sign change'
    for it in range(60):
        mid = arb(((lo + hi)/2).mid())
        fm = L.D(acb(mid)).real
        if fm.contains(0): break
        if (fm > 0) == (flo > 0): lo, flo = mid, fm
        else: hi = mid
    enc = (lo + hi)/2 + sym((hi - lo)/2)
    log('kappa* in %s   gamma = 1/kappa* in %s' % (enc.str(15), (1/enc).str(15)))
    return enc

def plan_cost(R1, sec_per_tm=None):
    p = disc_plan(R1)
    sec = 0.0
    for (_, a, b, r) in p:
        z = abs((a + b)/2)
        key = min(SEC_PER_TM, key=lambda k_: abs(k_ - z)) if sec_per_tm is None else None
        sec += (SEC_PER_TM[key] if sec_per_tm is None else sec_per_tm(z))
    return len(p), sec

if __name__ == '__main__':
    T0 = time.time()
    def log(*a): print('[%7.1fs]' % (time.time() - T0), *a, flush=True)
    mode = sys.argv[1]
    if mode == 'plan':
        R1 = float(sys.argv[2]); n, s = plan_cost(R1)
        log('R1=%g: %d discs on the upper path, projected %.1f h single process' % (R1, n, s/3600))
    elif mode == 'count':
        count(float(sys.argv[2]), sys.argv[3], log)
    elif mode == 'kstar':
        kstar(log)
