# RUN.md - build-2 QB1b (Taylor-model extension of the build-1 instrument)

All commands from `<LAB>/quantum/seats/build-2`, every Python run under memcap 8, no environment flags.
The build-1 files (ajet, aser, charts, bg, ctl, qb1) are imported unedited. `tm.install()` only rebinds the
constructor names `ajet.arb`, `charts.arb` (to ArbW, which is flint.arb for numbers and passes a TM through) and
`aser.A`, inside this process. Nothing is sealed, settled or pinned.

```
M="py -3.11 <MEMCAP> 8 py -3.11"
# PQ2b sign map: 9 region runs, run in parallel (each < 1 h)
$M signmap2.py w      -12    -11/2   13  5 sm_w.json
$M signmap2.py centre -11/2  -9/2    16  5 sm_centre.json
$M signmap2.py tube   -9/2   -1/16   284 5 sm_tube.json
$M signmap2.py tube   -1/16  -13/256 3   5 sm_tube2.json   # added: sonic cells left of -0.0605 leave the sm disc
$M signmap2.py sonic  -1/16  13/256  29  5 sm_sonic.json
$M signmap2.py outw   13/256 3       151 5 sm_out1.json
$M signmap2.py outw   3      9       384 7 sm_out2.json
$M signmap2.py outw   9      21/2    384 7 sm_out3.json
$M signmap2.py outw   21/2   12      768 7 sm_out4.json
$M tasks2.py seam       # -> pq8.json
$M tasks2.py null       # -> pq6b.json
$M tasks2.py controls   # -> controls_tm.json
$M dnew2.py             # -> pq7b.json
$M finalize.py          # -> result_qb1b.json, controls_qb1b.json (merge, coverage check, zero)
```
(args: region lo hi n0 P out; cells start at width (hi-lo)/n0, bisected to >= 1e-6, endpoints dyadic rationals)

## Method
- `tm.py`: Taylor models in the cell variable eps (|eps| <= H), degree P, with a rigorous remainder.
  - Products: degrees > P are bounded by sum |d_k| H^k; cross terms by |P1| r2 + |P2| r1 + r1 r2.
  - 1/x, sqrt, log, exp, sin, cos: Taylor coefficients at c0 plus a Lagrange remainder over the TM range.
  - Comparisons are decided on the range c0 +/- (sum |c_k| H^k + r).
- `schart.py`: a direct chart in (s, r), with s = KHA time and base s0 = x0, so delta x = ln(1+dr) + ds.
  - Metric: g_ss = -A N^2 r^2, g_rr = A. Inputs are the x-series of G = A N^2, A, om and V.
  - The cell centre never enters the jets. The build-1 polar and Cartesian charts carry e^{x0}, whose TM remainder
    is amplified by the cancellations: probes showed widths of 1e4 to 1e20.
- `tmcell.py`: every cell has x = x_m + eps exactly.
  - Tube and outward cells: state TM from the step polynomial plus the order-0 record tail, then ODE x-jets in TM
    arithmetic.
  - Centre and sonic manifold cells: theta(x_m + eps) solved as a TM by a contraction (residual/(1-L),
    L = sup |1 - x'/d0| on a widened interval, x' > 0 checked), then theta-jets and reversion.
  - Deep centre, x <= -11/2: w = z^2 series. The exact w-Taylor coefficients come from the stored cm coefficients;
    the tail is Cauchy, |g_n| <= Mg/Rw^n. Mg is sup|g| over the complex box |theta| <= 0.02 (inside the cm tail disc),
    and Rw = 0.02^2 e^{2 xc}/4. A Rouche condition |e^{2X} - 1| < 3/4 is checked on that box. The evaluation uses
    the build-1 Cartesian chart in w-mode, where there is no 1/zeta factor.
- Zero (`finalize.py`): on the hull of the unresolved tube cells, b_nn' (deg-5 TM, s-chart) is checked to be
  > 0 piecewise, and the neighbouring cells must be negative on the left and positive on the right. Bisection to
  <= 1e-8 uses degenerate TMs (H = 0) on the same code path.
- PQ6b: after integration, the omega profile (Wt series, all orders) is multiplied by (1 + 1e-12) and F2, At, Vt
  are kept. Same 28 points and charts as build-1 (`qb1.route_points` list).
