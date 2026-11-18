# RUN.md - build-1 QB1 interval instrument (python-flint arb)

All commands from `<LAB>/quantum/seats/build-1`, every Python run under memcap 8. No environment flags.
Nothing here is sealed, settled or pinned; engine / sideways / gravity files are imported read-only.

```
M="py -3.11 <MEMCAP> 8 py -3.11"
$M qb1.py controls        # -> controls_qb1.json  (C0 Hayward, C3 RN, C2 Minkowski, C6 Friedmann k=1/2, sign-flip null)
$M qb1.py main            # -> result_qb1.json    (sonic, centre c0, 27 route rows, sign map + zero), main.log
$M dnew.py                # -> dnew_qb1.json      (D_new(kappa) via engine/cert/count_b4v3.py install(), kappa = 0, 1, 2)
$M null_bg.py in          # -> null_in.pkl        (off-solution V0* inward leg, lin_tube_v2 control flow)
$M null_bg.py out         # -> null_out.pkl       (off-solution V0* outward leg, outward_rec.outward_rec, BG2 params)
$M null_eval.py           # -> null_qb1.json
$M make_controls.py       # -> adds instrument sha256 + dnew/null summaries to controls_qb1.json
```

## Instrument files
- `ajet.py` arb multivariate Taylor jets (NV 3 or 4, DEG 4 or 5), curvature, Zahn T(1)/T(2), V(xi), Bach, Route A.
  Port of smoke-1 jetcurv.py.
- `aser.py` univariate arb series, reversion, `Piece` (polynomial + certified analytic tail, re-expanded at a point or
  interval; tail coefficients |e_n| <= M/R^n summed explicitly).
- `charts.py` two adapters, same profile inputs (F2 = A N^2 z^2, At = (A-1)/z^2, Wt = om/z^2, Vt = V/z):
  polar-areal (t, r) at r0 = 1, and Cartesian (t, X, Y, Z) at t0 = -1 (regular at the centre, used for c0 at r = 0).
  Identity used: s = -ln(-t), x = ln r + s turn the KHA metric into -A N^2 z^2 dt^2 + A dr^2 + r^2 dOmega^2 (z = e^x),
  so the slices and n, e_r are the KHA ones.
- `bg.py` background: tube lin_tube2_k1_3.pkl (E5; x-steps re-jetted from c +- 1.01 setrad as s1_cert; sm; cm),
  BG2 outward_C.pkl (x-steps, sharp sm). Point values: certified state ball from the record (order-0 tail), then
  x-jets from the ODE y' = f(y) (method b). Interval values: hull of the step boxes over the interval, then ODE jets.
  F2 is built from ln(N e^x)' = A - 1 - (1-k) om (KHA 18b), not from N^2 e^{2x} (catastrophic series cancellation);
  both forms are checked to overlap at every evaluation (`A_identity_check`).
- `ctl.py` controls; `qb1.py` driver; `dnew.py`; `null_bg.py`, `null_eval.py`.

## Precision rules
ctx.cap = 64 at module tops (dnew.py: 32 = the sealed count stack's value); ctx.prec = 256 set in main and restored in
`finally`. Constants such as k = 1/3 are built at call time (an import-time arb(1)/3 at 53 bits was found and fixed
before the final run). No flint object is pickled (records are strings).

## Outcome notes (details in the JSON files)
- Run order actually used: controls -> main (final, 140 s) -> dnew (29 s) -> null_bg in (510 s) / out (800 s), run in
  parallel -> null_eval (3 s) -> make_controls.
- Killed by task id (own runs): two earlier `main` runs (the first had k = 1/3 built at 53 bits at import; the second
  was an interval-cover sign map that was running away, see below).
- Sign map: interval covering was REFUSED. Naive interval evaluation of b_nn on an x-cell wraps by 1e4..1e20 per unit
  relative width (`interval_cell_attempts` in result_qb1.json; probes t_cov.py, t_tay.py also tried 2nd-4th order
  Taylor remainders). The delivered sign map is certified signs at 691 grid POINTS plus a certified zero bracket
  (bisection, width 2.4e-9) and a certified local uniqueness interval (mean-value derivative bound, half-width ~1e-5).
- D_new: the sealed count instrument matches at x_m = XBAR = -4.6 (tube end); it has no x_m parameter, so x_m = -2 is
  only available from the float lin_float (NON-COUNTING), which reproduces the certified value at -4.6 to 10 digits.
- Null: the off-solution V0* background still solves the ODE locally, so it satisfies Einstein locally: the residual
  contains 0 everywhere and Routes A and B overlap everywhere. The instrument does read the background: null b
  values are disjoint from the certified ones at x = 0 and on the outward leg 0.25..3 (b_uu shift ~1e-15).
- Probe scripts t_*.py, p_*.py are non-counting scratch.
