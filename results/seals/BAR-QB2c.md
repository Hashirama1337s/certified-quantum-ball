# BAR-QB2c: the conservation width with tight jets, the theta-theta residual over the used domain, and the replication row

Written 2026-11-18 by the lead after BAR-QB2b (373f306e9f) settled FAILED by item: PQ2b-1 only (e on the real source
encloses 0, width 6.66e-10 > 1e-10). PQ2b-2..5 HOLD. No threshold is loosened. Disclosed, calibrated_on: build-4
result_qb2b.json (e width 6.66e-10, cause = ball-recursion dependency in the homogeneous step matrices, order-0 radius
1.3e-13 vs start ball 5e-17; Psi2 = 3.4616358037141657 +/- 3.0e-12; thth on the real source at x = -0.01 only) and
seats/referee-3/RECONCILE-QB2.md (blind Psi converts exactly with 30 pi x_s^2 / A(x_V); agreement 2.4e-12 relative).

**Object.** As BAR-QB2 / QB2b. Inward background = build-4 lin_tube_bg2_k1_3.pkl (pinned under 373f306e9f).
Homogeneous stage rebuilt with build-4 qb2b_tight.py (T1 mean-value step jets, T2 mean-value step tapes), engine files
untouched. Everything else as the build-4 run-1 record.

instruments: seats/build-5/

## Predictions
- **PQ2c-1 (conservation width, the QB2 bar unchanged).** With the tight homogeneous stage: e on the real source encloses
  0 with width <= 1e-10; C8's e encloses 0 (no width bar); N-a and N-c still refuse (e excludes 0).
- **PQ2c-2 (theta-theta over the used domain).** The linearised E_thth residual (build-4 thth_sym.py expressions) on the
  real source encloses 0 at x = -6, -4.6, -2, -0.25, -0.01; the N-b variant (b_T x (1 + 1e-6)) excludes 0 at every point.
- **PQ2c-3 (headline replay).** Psi2 overlaps [3.4616358037141657 +/- 3.0e-12], excludes 0, relative width <= 1e-3.
- **PQ2c-4 (replication row).** 30 pi x_s^2 Psi_blind / A(x_V), with Psi_blind read from
  seats/blind-qb2/blind_qb2_result.json and x_s, A(x_V) from the certified enclosures, agrees with Psi2 to relative
  <= 1e-9 (seen value 2.4e-12, disclosed; the row fixes the conversion in the ledger).

## Kill conditions
- Any enclosure disjoint from the QB2b ones (Psi2, Abar_2, x_V, A(x_V)): FAILED and STOP.
- Tight jet assertion (tight jets must overlap the ball recursion) fires: FAILED and STOP.
- A single run over 2 hours or a total over 4 hours: INCONCLUSIVE for the unfinished items.

## Verdict
HOLD if PQ2c-1 .. PQ2c-4 hold.

## Seal fields
- formula_ref: DESIGN-QB2-v2.md:37-118
- calibrated_on: 373f306e9f build-4 result_qb2b.json and RUN.md; referee-3 RECONCILE-QB2.md - all values above seen before sealing (disclosed)
- depends_on: 373f306e9f, f2bcfadb3e, 4d3b12ef39, ed61d09268
- null_not_symmetry: N-a predicts a specific nonzero mismatch; N-c perturbs the background off shell; N-b breaks the theta-theta equation by scaling one stress component
- target_form: arb enclosures; overlap tests; relative agreement of two independent implementations
- parameters: e width 1e-10; thth points -6, -4.6, -2, -0.25, -0.01; N-b factor 1 + 1e-6; replication relative 1e-9
- tolerances: interval arithmetic (arb); replication float relative 1e-9
- ceiling: a fully closed certified first-order quantum push on the peak compactness of a true critical solution
- positive control: C1 manufactured phi G_2 (rerun), C8 manufactured conserved source
- toy_model: C1 / C8 manufactured cases
- files: BAR-QB2c.md
- command: the build seat's own (seats/build-5/RUN.md), under memcap 8, watch-run for anything over 10 minutes
