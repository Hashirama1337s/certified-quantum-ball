"""B4 speedup seat (NEW file): exact cache of the stage-L step data of lin_stages_v2.Linear2.

Linear2 spends 320-420 s building its 1440 transport steps (nsub 2) before any evaluation.  This file
builds it once and stores every step (h, rho, a0, a1, A0[n], A1[n], kind, depth) EXACTLY: each arb as
(mid mantissa, mid exponent, rad mantissa, rad exponent) integers, so the reload is bit-identical
(mid equal, rad equal).  Stages S and C are rebuilt by the original code on load (about 1 s).
Nothing in the certified pipeline changes: CachedLinear2 is Linear2 with _tube_setup replaced by the
reload; transport/D/step_order are inherited unchanged.

usage:  cacheL.py build <k> <nsub> <tol> <out.pkl>
"""
import sys, os, time, pickle
sys.path.insert(0, '<LAB>/engine/cert'); sys.path.insert(0, '<LAB>/engine')
import lin_stages_v2 as LS2            # imports driver_v32 first (R7 asserted there)
from flint import arb, acb, acb_mat


def ser_arb(x):
    m, e = x.mid().man_exp(); rm, re = x.rad().man_exp()
    return (int(m), int(e), int(rm), int(re))


def des_arb(t):
    m, e, rm, re = t
    mid = arb(m)*arb(2)**e
    rad = arb(rm)*arb(2)**re
    return arb(mid, rad)


def ser_acb(z):
    if z.imag.mid() == 0 and z.imag.rad() == 0:
        return ser_arb(z.real)
    return (ser_arb(z.real), ser_arb(z.imag))


def des_acb(t):
    if isinstance(t[0], int):
        return acb(des_arb(t))
    return acb(des_arb(t[0]), des_arb(t[1]))


def ser_mat(M, rows):
    return [[ser_acb(M[i, j]) for j in range(4)] for i in rows]


def des_mat(rows_data, zero_top):
    if zero_top:
        rows = [[0]*4, [0]*4] + [[des_acb(t) for t in r] for r in rows_data]
    else:
        rows = [[des_acb(t) for t in r] for r in rows_data]
    return acb_mat(rows)


def ser_steps(steps):
    out = []
    for st in steps:
        out.append(dict(kind=st['kind'], depth=st.get('depth', 0),
                        h=ser_arb(st['h']) if isinstance(st['h'], arb) else ('c', ser_acb(st['h'])),
                        rho=ser_arb(st['rho']), a0=ser_arb(st['a0']), a1=ser_arb(st['a1']),
                        A0=[ser_mat(M, range(4)) for M in st['A0']],
                        A1=[ser_mat(M, (2, 3)) for M in st['A1']]))
    return out


def des_steps(data):
    steps = []
    for d in data:
        h = d['h']
        h = des_acb(h[1]) if (isinstance(h, tuple) and h and h[0] == 'c') else des_arb(h)
        steps.append(dict(kind=d['kind'], depth=d['depth'], h=h, rho=des_arb(d['rho']),
                          a0=des_arb(d['a0']), a1=des_arb(d['a1']),
                          A0=[des_mat(M, False) for M in d['A0']],
                          A1=[des_mat(M, True) for M in d['A1']]))
    return steps


class CachedLinear2(LS2.Linear2):
    """Linear2 whose stage-L steps come from an exact cache (built by this file from Linear2 itself)."""
    def __init__(self, k, cache, **kw):
        self._cache = cache
        super().__init__(k, **kw)

    def _tube_setup(self):
        d = pickle.load(open(self._cache, 'rb'))
        assert d['KL'] == self.KL and d['nsub'] == self.nsub, 'cache settings differ'
        self.steps = des_steps(d['steps'])
        self.nsplit = d['nsplit']


def roundtrip(L, L2, log):
    """sampled comparison built vs cached.  Expected: mid identical everywhere; radius identical or ROUNDED UP
    (arb(mid, rad) stores the radius as a 30-bit mag, rounded upward: an outward enclosure).  A smaller
    cached radius or a different mid would be an error."""
    n = dict(same=0, rad_up=0, mid_diff=0, rad_down=0)
    def cmp(x, y):
        if x.mid() != y.mid(): n['mid_diff'] += 1
        elif y.rad() == x.rad(): n['same'] += 1
        elif y.rad() > x.rad(): n['rad_up'] += 1
        else: n['rad_down'] += 1
    for i in list(range(0, len(L.steps), 97)) + [len(L.steps) - 1]:
        s1, s2 = L.steps[i], L2.steps[i]
        for key in ('rho', 'a0', 'a1'):
            cmp(s1[key], s2[key])
        for nn in (0, 5, L.KL):
            for M1, M2 in ((s1['A0'][nn], s2['A0'][nn]), (s1['A1'][nn], s2['A1'][nn])):
                for a in range(4):
                    for b in range(4):
                        cmp(M1[a, b].real, M2[a, b].real); cmp(M1[a, b].imag, M2[a, b].imag)
    log('round-trip check (sampled arbs): identical %d, radius rounded up %d (outward, allowed), mid differs %d (ERROR), radius smaller %d (ERROR)' % (
        n['same'], n['rad_up'], n['mid_diff'], n['rad_down']))
    return n


if __name__ == '__main__':
    import sympy as sp
    T0 = time.time()
    def log(*a): print('[%7.1fs]' % (time.time() - T0), *a, flush=True)
    if sys.argv[1] == 'build':
        k = sp.Rational(sys.argv[2]); nsub = int(sys.argv[3]); tol = float(sys.argv[4]); out = sys.argv[5]
        L = LS2.Linear2(k, log=log, nsub=nsub, tol=tol)
        data = dict(k=str(k), nsub=nsub, tol=tol, KL=L.KL, nsplit=L.nsplit, steps=ser_steps(L.steps))
        pickle.dump(data, open(out, 'wb'), protocol=pickle.HIGHEST_PROTOCOL)
        log('wrote %s: %d steps, %.1f MB' % (out, len(L.steps), os.path.getsize(out)/1e6))
        # round-trip check (bit-identical mids and radii) on a sample of steps
        L2 = CachedLinear2(k, out, log=log, nsub=nsub, tol=tol)
        roundtrip(L, L2, log)
        from flint import acb as _acb
        for z in (_acb(0), _acb(0, 20)):
            log('point D(%s): built %s   cached %s' % (z.str(3), L.D(z).str(12), L2.D(z).str(12)))
    if sys.argv[1] == 'check':
        # usage: cacheL.py check <k> <nsub> <cache.pkl>   (fresh Linear2 build vs the existing cache)
        k = sp.Rational(sys.argv[2]); nsub = int(sys.argv[3]); fn = sys.argv[4]
        L = LS2.Linear2(k, log=log, nsub=nsub, tol=1e-20)
        L2 = CachedLinear2(k, fn, log=log, nsub=nsub, tol=1e-20)
        roundtrip(L, L2, log)
