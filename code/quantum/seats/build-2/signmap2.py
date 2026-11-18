"""build-2 QB1b PQ2b: covering sign map of b_nn on [-12, 12] by closed cells with exact dyadic-rational endpoints.

Regions (exact breakpoints) and cell evaluators (tmcell.py, Taylor models in eps = x - x_m):
  w       [-12, -11/2]      Cartesian chart, w = z^2 series of the centre manifold with Cauchy tail (wcell)
  centre  [-11/2, -9/2]     centre manifold cm, theta(x) solved as a TM, direct s-chart
  tube    [-9/2, -1/16]     tube (E5) x-steps, ODE jets, s-chart
  sonic   [-1/16, 13/256]   sonic manifold sm (E5 tube), theta(x) solved as a TM, s-chart
  outw    [13/256, 12]      outward (BG2) x-steps, ODE jets, s-chart
A cell is certified when the TM range of b_nn (or b_nn e^{-4x} in the w region; same sign) excludes 0.  Otherwise
(or on any ArithmeticError: disc left, no contraction, ...) it is bisected, down to width 1e-6; narrower failures
are recorded as unresolved.
usage (memcap 8): py -3.11 signmap2.py <region> <lo> <hi> <n0> <P> <out.json>     (lo, hi as 'p/q' or integers)
"""
import sys, os, time, json
from fractions import Fraction as F
from flint import arb, ctx
ctx.cap = 64
import tm
import tmcell
import bg

MINW = F(1, 10**6)


def A(fr):
    return arb(fr.numerator)/arb(fr.denominator)      # exact for dyadic denominators


def main():
    region, lo, hi, n0, P, outn = sys.argv[1], F(sys.argv[2]), F(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), sys.argv[6]
    T0 = time.time()
    log = lambda *a: print('[%7.1fs]' % (time.time() - T0), *a, flush=True)
    old = ctx.prec
    try:
        ctx.prec = 256
        tm.P = P
        Bg = bg.Background()
        cells = []; unres = []; nev = 0
        work = [(lo + (hi - lo)*i/n0, lo + (hi - lo)*(i + 1)/n0) for i in range(n0)]
        for a_, b_ in work:
            assert a_.denominator & (a_.denominator - 1) == 0, 'non-dyadic endpoint'
        while work:
            a, b = work.pop(0)
            nev += 1
            s = 0; v = None; err = None
            try:
                if region == 'w':
                    o = tmcell.wcell(Bg, A(a), A(b)); v = tm.rng(o['hat_b_nn'])
                else:
                    g = tmcell.guess_theta(Bg, region, (A(a) + A(b))/2) if region in ('centre', 'sonic') else None
                    o = tmcell.cell(Bg, region, A(a), A(b), th_guess=g); v = tm.rng(o['b_nn'])
                s = 1 if v > 0 else (-1 if v < 0 else 0)
            except (ArithmeticError, AssertionError) as e:
                err = str(e)[:60]
            if s != 0:
                cells.append(dict(a=str(a), b=str(b), sign=s, enc=[v.mid().str(12, radius=False), v.rad().str(3, radius=False)]))
            elif b - a > MINW:
                m = (a + b)/2
                work[0:0] = [(a, m), (m, b)]
            else:
                unres.append(dict(a=str(a), b=str(b), err=err, enc=(v.str(5) if v is not None else None)))
            if nev % 25 == 0:
                log('%d evals, %d cells, %d unresolved, at x=%s' % (nev, len(cells), len(unres), float(a)))
        out = dict(region=region, lo=str(lo), hi=str(hi), P=P, n0=n0, evals=nev, cells=cells, unresolved=unres,
                   time_s=round(time.time() - T0, 1))
        if region == 'w':
            out['wpiece'] = tmcell._WP.get('info')
        json.dump(out, open(outn, 'w'), indent=0)
        log('done: %d evals, %d certified cells, %d unresolved -> %s' % (nev, len(cells), len(unres), outn))
    finally:
        ctx.prec = old


if __name__ == '__main__':
    main()
