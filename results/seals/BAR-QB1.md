# BAR-QB1: "The Quantum Source" - the certified one-loop Bach source on the Evans-Coleman critical solution

Written 2026-11-17 by the lead (owner: "Yes, go for it"). Design: DESIGN-QB1.md (design-1). Disclosed, calibrated_on: the
float smoke seats/smoke-1/SMOKE-QB1.md (NOT SEALED) - b_nn sign pattern (negative core, one zero near x = -0.84, positive
shell), sonic values, centre coefficient c0 ~ -173.1004, D_new(2) ~ -0.1857 - was seen before sealing. Thresholds are
set from those floats plus margin; the new content is that every statement becomes a certified enclosure.

**Object.** On the certified Evans-Coleman k = 1/3 solution (E5 d0af75803f; sharper sonic box BG2 ed61d09268), the
Bach tensor B_ab (so the state-independent one-loop term u B_ab / (240 pi^2) of Zahn's Hadamard expansion, for every
curvature coupling, since R = 0 on a traceless background) in dimensionless form b_ab = r^4 B_(ab), Eulerian frame
(n, r, T) and fluid frame (b_uu). Interval arithmetic; Route A (on-shell, Ricci from the fluid stress tensor + Weyl) and
Route B (metric jets to order 4, off-shell) both evaluated.

instruments: seats/build-1/

## Predictions
- **PQ1 (not Bach-flat, certified).** b_nn and b_uu enclosures exclude 0 at the sonic point; b_uu(sonic) overlaps
  0.446673779670 +/- 1e-9 with enclosure width <= 1e-15.
- **PQ2 (sign map).** On x in [-12, 12]: b_nn < 0 certified on [-12, -1]; b_nn > 0 certified on [-0.75, 12]; exactly
  one zero of b_nn in [-1, -0.75], enclosed in an interval of width <= 1e-6 (certified sign change + certified
  monotonicity or derivative bound on the enclosing interval).
- **PQ3 (centre law).** In the blown-up centre chart, b_nn e^{-4x} -> c0 with a certified enclosure of c0 overlapping
  -173.1004 +/- 1e-3 and excluding 0 (the isotropic v^4 law; sign negative at the centre).
- **PQ4 (two routes agree).** At >= 18 points spread over [-12, 12] (centre chart, inward, sonic series, light cone,
  outward), the Route A and Route B enclosures overlap for every component; the Einstein-residual enclosure contains 0.
- **PQ5 (controls).** Through the same interval code path: Hayward xi = 0 gives (1/80, 17/240, -7/480) (enclosures
  contain the exact values); Reissner-Nordstrom (M = 1, Q = 1/2) Bach = 60 V; Minkowski patch B = 0 (enclosure contains
  0); k = 1/2 Friedmann r^4 V_uu(xi = 0) contains 2 z^4 / 243 at z = 1/2, 1, 3/2 and its Bach contains 0.
- **PQ6 (null, teeth).** An off-solution background (V0 moved to the E5 box edge + 1e-15, i.e. outside the certified box,
  integrated the same way) gives an Einstein-residual enclosure excluding 0 at some point of the sample and Route A /
  Route B enclosures that are disjoint somewhere.
- **PQ7 (bridge).** D_new(2) certified nonzero (enclosure excludes 0), negative, at matching x_m = -2.

## Kill conditions
- PQ1 fails (b enclosures contain 0 at the sonic point with width <= 1e-15): the Bach source vanishes there - FAILED,
  report; the field-independent sign story is then weakened.
- PQ5 control fails: FAILED and STOP (the code path is wrong).
- Any run over 2 hours: stop, INCONCLUSIVE.

## Not claimed
x > 12 (needs a tail lemma; next bar); the quantum state's existence (all statements read "in any Hadamard state
respecting self-similarity"); the net trapping sign (QB2/QB3 integral).

## Verdict
HOLD if PQ1-PQ7 hold.

## Seal fields
- formula_ref: DESIGN-QB1.md:26-136
- calibrated_on: seats/smoke-1/SMOKE-QB1.md float values seen before sealing (disclosed above)
- depends_on: d0af75803f, ed61d09268, b1b0ee1cde
- null_not_symmetry: PQ6 evaluates the same formulas on a metric that does not solve the field equations, so the on-shell and off-shell routes must disagree; that tests the content of Route A
- target_form: arb enclosures; certified signs on intervals; one certified zero interval; overlap tests between routes
- parameters: sonic tolerance 1e-9 width 1e-15; zero interval width 1e-6; c0 tolerance 1e-3; null offset 1e-15 outside the E5 box
- tolerances: interval arithmetic (arb)
- ceiling: the first certified one-loop stress-tensor coefficient on a true critical-collapse solution, the same shape for every free massless field
- positive control: PQ5 closed forms
- toy_model: Hayward, Reissner-Nordstrom, Minkowski, k = 1/2 Friedmann (PQ5)
- files: BAR-QB1.md
- command: the build seat's own (seats/build-1/RUN.md), under memcap 8
