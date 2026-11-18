# BAR-QB1b: "The Quantum Source", second attempt - covering sign map, a real null, the bridge at its own matching point

Written 2026-11-17 by the lead after BAR-QB1 (54caa4012c) settled FAILED by item (PQ1, PQ3, PQ4, PQ5 HOLD; PQ2, PQ6, PQ7 FAILED
as sealed). Disclosed, calibrated_on: seats/build-1/result_qb1.json - the zero bracket [-0.82611434221, -0.82611433983],
b_nn' there 1.16548010 +/- 3e-8, grid signs at 691 points, D_new(2) = -2.41262709 +/- 5.3e-9 at x_m = -4.6, and the
failed null (an off-solution V0 still solves the ODE, so it satisfies Einstein locally) were all seen before sealing.
The held items PQ1/PQ3/PQ4/PQ5 are not re-tested here; they stand on 54caa4012c.

**Object.** As BAR-QB1: b_ab = r^4 B_(ab) on the certified Evans-Coleman k = 1/3 solution, Routes A and B, same
interval instrument (seats/build-1, pinned under 54caa4012c) extended with Taylor-model arithmetic in (x - x_m) carried
through the ODE jets.

instruments: seats/build-2/

## Predictions
- **PQ2b (covering sign map).** Closed x-cells covering [-12, 12] with no gap (cell endpoints exact rationals, union
  checked), each with a certified sign of b_nn over the whole cell: negative on every cell inside [-12, -0.8262],
  positive on every cell inside [-0.8261, 12]; the cell containing the zero carries a certified derivative lower bound
  b_nn' > 0, so b_nn has exactly one zero on [-12, 12], enclosed in an interval of width <= 1e-8 overlapping the QB1
  bracket.
- **PQ6b (real null, off-shell).** Replace the certified omega profile by omega (1 + 1e-12) AFTER integration (a metric
  that does not solve the field equations), keeping everything else: the Einstein-residual enclosure excludes 0 at at
  least one of the 28 QB1 sample points, and Route A / Route B enclosures are disjoint in at least one component there.
- **PQ7b (bridge).** D_new(2) re-evaluated by the certified spherical instrument at its own matching point x_m = -4.6:
  enclosure overlaps -2.41262709 +/- 1e-7 and excludes 0 (replay of a value seen; the content - nonzero at kappa = 2 -
  is independent of the matching point because the zero set does not depend on it).
- **PQ8 (centre to grid seam).** The centre-chart enclosure of b_nn and the tube enclosure overlap at a common x in
  [-4.6, -4.4] (same x exactly).

## Kill conditions
- PQ2b finds a cell whose sign cannot be certified after refinement to cell width 1e-6: FAILED for PQ2b (Taylor
  models not enough), report.
- PQ6b Einstein residual contains 0 at all points: FAILED (instrument does not see off-shell content) - STOP.
- Any run over 2 hours: stop, INCONCLUSIVE.

## Not claimed
x > 12 (tail lemma, next bar); state existence (all statements "in any Hadamard state respecting self-similarity"); net
trapping sign (QB2/QB3).

## Verdict
HOLD if PQ2b, PQ6b, PQ7b, PQ8 hold. Together with the held items of 54caa4012c this certifies the QB1 source on [-12, 12].

## Seal fields
- formula_ref: DESIGN-QB1.md:26-136
- calibrated_on: 54caa4012c seats/build-1/result_qb1.json - zero bracket, grid signs, D_new(2) at -4.6, failed null mechanism (disclosed above)
- depends_on: 54caa4012c, d0af75803f, ed61d09268, b1b0ee1cde
- null_not_symmetry: PQ6b uses a metric that does not solve the field equations, so the on-shell Route A must disagree with the off-shell Route B; it tests the content of Route A
- target_form: arb / Taylor-model enclosures; certified signs on covering cells; one certified zero interval; route overlap tests
- parameters: zero width 1e-8; cell refinement floor 1e-6; null factor 1 + 1e-12 on omega; bridge tolerance 1e-7; seam window [-4.6, -4.4]
- tolerances: interval arithmetic (arb)
- ceiling: the first certified one-loop stress-tensor source on a true critical-collapse solution, with its full sign map
- positive control: held QB1 controls (PQ5 of 54caa4012c) re-run through the extended code path
- toy_model: Hayward, Reissner-Nordstrom, Minkowski, k = 1/2 Friedmann (carried)
- files: BAR-QB1b.md
- command: the build seat's own (seats/build-2/RUN.md), under memcap 8
