"""build-4 QB2b: drive the PINNED build-3 instrument (imported, never edited) on the BG2-based inward tube.
In-process changes only (documented in RUN.md):
  P1  engine: lin_stages_v2.Linear2 on build-4/lin_tube_bg2_k1_3.pkl (no count_b4v3 install, no E5 cache), nsub = 2,
      tol = count_b4.TOL; its stage-L steps cached once with b4speed/cacheL (exact int tuples, no flint pickles).
  P2  `open` in the build-3 modules redirects the E5 tube path to the BG2 tube; HERE -> build-4.
  P3  RegSingAug.evaluate: the majorant radius nu' is OPTIMISED in (|t|, nu) (the build-3/engine choice
      nu' = min((nu+|t|)/2, 2|t|) forces q = |t|/nu' >= 1/2); same rigorous bound for every admissible nu'.
  P4  make_scaled / evaluate stash (for C1 at x = -0.01 on the sonic series).
"""
import sys, os, builtins, pickle, math, json, time
from flint import ctx as _ctx
_ctx.cap = 32
sys.path.insert(0, '<LAB>/engine/cert'); sys.path.insert(0, '<LAB>/engine')
sys.path.insert(0, '<LAB>/b4speed')
sys.path.insert(0, '<LAB>/quantum/seats/build-3')
import sympy as sp
from flint import arb, acb, arb_mat, ctx
from core import sym, csym
import count_b4 as CB
import lin_stages_v2 as LS2
import cacheL as CL
import qb2stages as ST
import qb2core as C
import qb2pass as PS
import qb2final as FN
import qb2lib as Q

B4 = '<LAB>/quantum/seats/build-4/'
OLD_TUBE = '<LAB>/engine/cert/lin_tube2_k1_3.pkl'
NEW_TUBE = B4 + 'lin_tube_bg2_k1_3.pkl'
CACHE = B4 + 'Lbg2t_nsub2.pkl'                 # homogeneous stage built WITH qb2b_tight (T1, T2)
if os.environ.get('QB2B_RUN1'):                  # run 1 (plain homogeneous stage): re-use its cache, no tight patches
    CACHE = B4 + 'Lbg2_nsub2.pkl'
elif os.environ.get('QB2B_NOPATCH'):            # pass processes: load the tight-built cache, no patches needed
    assert os.path.exists(CACHE), 'tight cache missing'
else:
    import qb2b_tight                             # noqa: E402  (tight jets + mean-value step tapes)


def _open(fn, *a, **kw):
    if isinstance(fn, str) and fn.replace(chr(92), '/') == OLD_TUBE:
        fn = NEW_TUBE
    return builtins.open(fn, *a, **kw)


for _m in (ST, PS, FN):
    _m.open = _open
PS.HERE = B4; FN.HERE = B4


def engine_b(log=print):
    k = sp.Rational(1, 3)
    if not os.path.exists(CACHE):
        L = LS2.Linear2(k, tube_file=NEW_TUBE, log=log, nsub=CB.NSUB, tol=CB.TOL)
        data = dict(k=str(k), nsub=CB.NSUB, tol=CB.TOL, KL=L.KL, nsplit=L.nsplit, steps=CL.ser_steps(L.steps))
        pickle.dump(data, open(CACHE, 'wb'), protocol=pickle.HIGHEST_PROTOCOL)
        log('built Linear2 on the BG2 tube and cached %d steps (nsplit %d)' % (len(L.steps), L.nsplit))
    else:
        L = CL.CachedLinear2(k, CACHE, tube_file=NEW_TUBE, log=log, nsub=CB.NSUB, tol=CB.TOL)
    ctx.prec = 128
    return L


ST.engine = engine_b

# P3: optimised majorant radius
_STASH = dict(extra_t=None, vals=[], lams=[])
_orig_ms = C.make_scaled


def _make_scaled(*a, **kw):
    r = _orig_ms(*a, **kw)
    _STASH['lams'].append(r[2])
    return r


C.make_scaled = _make_scaled


def _eval_core(self, y, t):
    K = self.K; n_ = y[0].nrows()
    at = abs(t).upper() if not isinstance(t, arb) else t.abs_upper()
    best = None
    for f in [0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95]:
        nup = at + (self.nu - at)*arb(f)
        rat = nup/self.nu
        mup = self.S*rat/(1 - rat)
        n0 = K
        ok = True
        while True:
            c = self.cn(n0 + 1)
            if c is not None and c*mup <= 1: break
            n0 += 1
            if n0 > K: ok = False; break
        if not ok:
            continue
        B = arb(0)
        for n in range(K + 1):
            B = B.max(max(abs(y[n][i, 0]).upper() for i in range(n_))*nup**n)
        qq = at/nup
        tail = B*qq**(K + 1)/(1 - qq)
        if best is None or tail < best[0]:
            best = (tail, n0, mup, f)
    if best is None:
        raise ArithmeticError('RegSingAug (optimised): no admissible nu\'')
    tail = best[0]

    def ev(tt):
        out = []
        for i in range(n_):
            s = 0*tt
            for n in range(K, -1, -1): s = s*tt + y[n][i, 0]
            out.append(s + (sym(tail) if isinstance(tt, arb) else csym(tail)))
        return out
    return ev(t), dict(tail=tail, n0=best[1], mup=best[2], nu_frac=best[3])


def _evaluate(self, y, t):
    if _STASH['extra_t'] is not None:
        v_, i_ = _eval_core(self, y, _STASH['extra_t'])
        _STASH['vals'].append((v_, i_))
    return _eval_core(self, y, t)


C.RegSingAug.evaluate = _evaluate
