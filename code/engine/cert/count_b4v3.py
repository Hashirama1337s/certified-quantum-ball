"""B4 counting driver v3 (speedup seat, NEW file).  count_b4.py (pinned v2) is imported, never edited.

What changes against count_b4 (PROPOSAL-B4-SPEEDUP.md):
  * TM engine: lin_tm_v3.LinearTM3 (characteristic frame 'A1') replaces lin_tm_v2.LinearTM2;
  * ctx.prec = 128 after the Linear2 setup (ball arithmetic is rigorous at any precision; same printed
    digits as 200 bits, b4speed/rt3_g_p128.log);
  * optional exact step cache (b4speed/cacheL.CachedLinear2: mids bit-identical, radii rounded up);
  * configurable contour and radius rule (defaults = the speedup proposal):
      LB_R1 (40)             radius of the half-disc (Far Wall 40 <= |kappa| <= 2000: fw_cert2 40-55 + fw_cert 55-2000)
      LB_NARC (8)            chords of the circumscribed arc polygon (vertices at R1/cos(pi/(4 N)))
      LB_RULE ('2:0.25,5:0.4,90:0.5')   r(|kappa|) for Re kappa < 5
      LB_RRIGHT (1.0)        r for Re kappa >= 5
      LB_FRAME ('A1'), LB_PREC (128), LB_CACHE13 (b4speed/L13_nsub2.pkl; '' = rebuild)
The counting logic itself (disc plan, segment test, refinement REFINE=3, arg sharpening, JSONL ledger,
0.25 rule, parity check, kappa* bisection) is count_b4's, called unchanged; the two engine classes are
substituted in the modules count_b4 imports them from (lin_stages_v2.Linear2, lin_tm_v2.LinearTM2) for
the lifetime of this process only.

usage (from <LAB>, under memcap):
  count_b4v3.py plan  [R1]
  count_b4v3.py count [R1] <ledger.jsonl>        (LB_SHAPE=box for the dry-run box, as count_b4)
  count_b4v3.py kstar
  count_b4v3.py smoke <k> <parts>                 (lin_b4_smoke3.py controls with the v3 engine; its LB_* env)
"""
import sys, os, time, runpy
# Referee B1 (seats/REFEREE-B4FW.md): python-flint series built without an explicit prec truncate at
# ctx.cap (default 10).  Raise the cap BEFORE any other import.  Containment test (b4speed/r2/capchk_*.log):
# the count stack built at cap 32 is identical to the as-run objects, so this changes no number; it only
# removes the dependence on the default.  Code that sets ctx.cap itself restores it to this value.
from flint import ctx as _ctx
_ctx.cap = 32
sys.path.insert(0, '<LAB>/engine/cert'); sys.path.insert(0, '<LAB>/engine')
sys.path.insert(0, '<LAB>/b4speed')
import count_b4 as CB
import lin_stages_v2 as LS2
import lin_tm_v2 as LT2
import lin_tm_v3 as LT3
import cacheL as CL
import sympy as sp
from flint import ctx

R1_DEF = float(os.environ.get('LB_R1', '40'))
CB.N_ARC = int(os.environ.get('LB_NARC', '8'))
CB.R_RULE = [tuple(float(x) for x in p.split(':')) for p in os.environ.get('LB_RULE', '2:0.25,5:0.4,90:0.5').split(',')]
CB.R_RIGHT = float(os.environ.get('LB_RRIGHT', '1.0'))
FRAME = os.environ.get('LB_FRAME', 'A1')
PREC = int(os.environ.get('LB_PREC', '128'))
CACHE13 = os.environ.get('LB_CACHE13', '<LAB>/b4speed/L13_nsub2.pkl')
CB.SEC_PER_TM = {5: 55.0, 20: 55.0, 40: 55.0}      # measured, frame A1, prec 128 (b4speed/rt3_g_p128.log)
CACHE13_SHA = 'c239cbe5347059080ea3c4144b8c5a65133cb4c4812a280096c8f5116dfbb782'   # b4speed/L13_nsub2.pkl (round-2 pin)
import count_b4v3_seal as SEAL                      # referee B4FW2 fixes (F1-F3, F5); count_b4 stays unedited
CB.upper_path = SEAL.upper_path_v3                  # F1: inflated vertices (used by plan_cost / count / smokes)
CB.disc_plan = lambda R1, shape=None: SEAL.disc_plan_v3(R1, shape)

_Linear2 = LS2.Linear2
_LinearTM2 = LT2.LinearTM2
INFO = {}


def make_linear(k, log=print, nsub=1, tol=1e-30, **kw):
    k = sp.Rational(k)
    if k == sp.Rational(1, 3) and nsub == 2 and CACHE13 and os.path.exists(CACHE13):
        h = SEAL.sha(CACHE13)
        assert h == CACHE13_SHA, 'cache sha256 %s differs from the pin %s' % (h, CACHE13_SHA)
        L = CL.CachedLinear2(k, CACHE13, log=log, nsub=nsub, tol=tol, **kw)
        INFO['linear'] = 'CachedLinear2(%s, sha256 verified %s)' % (os.path.basename(CACHE13), h[:12])
    else:
        L = _Linear2(k, log=log, nsub=nsub, tol=tol, **kw)
        INFO['linear'] = 'Linear2 (fresh build)'
    ctx.prec = PREC
    log('[v3] engine: %s, prec set to %d after setup, TM frame %s' % (INFO['linear'], ctx.prec, FRAME))
    return L


def make_tm(L, m=16, rem_method='c1', J=4):
    return LT3.LinearTM3(L, m=m, J=J, frame=FRAME)


def install():
    LS2.Linear2 = make_linear
    LT2.LinearTM2 = make_tm


if __name__ == '__main__':
    T0 = time.time()
    def log(*a): print('[%7.1fs]' % (time.time() - T0), *a, flush=True)
    install()
    mode = sys.argv[1]
    log('[v3] count_b4v3 %s: ctx.cap %d, N_ARC %d, R_RULE %s, R_RIGHT %g, R_EDGE_RE %g, shape %s, frame %s, prec %d, m %d J %d nsub %d tol %g REFINE %d, LB_CACHE13 %s, LB_R_FIXED %r, LB_MAXDISCS %r' % (
        mode, ctx.cap, CB.N_ARC, CB.R_RULE, CB.R_RIGHT, CB.R_EDGE_RE, CB.SHAPE, FRAME, PREC, CB.M_DEG, CB.J_DEG, CB.NSUB, CB.TOL, CB.REFINE,
        CACHE13, os.environ.get('LB_R_FIXED'), os.environ.get('LB_MAXDISCS')))
    if mode == 'plan':
        R1 = float(sys.argv[2]) if len(sys.argv) > 2 else R1_DEF
        n, s = CB.plan_cost(R1)
        log('R1=%g: %d discs on the upper path, projected %.2f h single process, %.2f h on two' % (R1, n, s/3600, s/7200))
    elif mode == 'count':
        if len(sys.argv) > 3:
            R1, ledger = float(sys.argv[2]), sys.argv[3]
        else:
            R1, ledger = R1_DEF, sys.argv[2]
        SEAL.count_v3(R1, ledger, log, dict(frame=FRAME, prec=PREC, cap=ctx.cap, cache=os.path.basename(CACHE13), cache_sha=CACHE13_SHA))
    elif mode == 'pre':
        # contour + plan checks and the F3 preconditions only (no discs)
        R1 = float(sys.argv[2]) if len(sys.argv) > 2 else R1_DEF
        SEAL.check_arc(R1, log); plan = SEAL.disc_plan_v3(R1); SEAL.check_plan(plan, R1, log)
        L = LS2.Linear2(sp.Rational(1, 3), log=log, nsub=CB.NSUB, tol=CB.TOL)
        T = LT2.LinearTM2(L, m=CB.M_DEG, rem_method='c1', J=CB.J_DEG)
        SEAL.certify_preconditions(L, T, plan, log)
    elif mode == 'kstar':
        CB.kstar(log)
    elif mode == 'smoke':
        sys.argv = ['lin_b4_smoke3.py'] + sys.argv[2:]
        runpy.run_path('<LAB>/engine/cert/lin_b4_smoke3.py', run_name='__main__')
