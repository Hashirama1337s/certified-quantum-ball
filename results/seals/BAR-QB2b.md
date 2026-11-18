# BAR-QB2b: closing the QB2 gaps (same headline, same thresholds, sharper inputs)

Written 2026-11-18 by the lead after BAR-QB2 (f2bcfadb3e) settled FAILED as sealed by item. Headline PQ2-8 held:
Psi2 = 3.4616358037141 +/- 1.15e-7 (seen; disclosed). This bar does NOT loosen any QB2 threshold. It re-tests the failed
items with sharper inputs and exact arguments. Disclosed, calibrated_on: seats/build-3 result_qb2.json and RUN.md
(e width 2e-5 driven by the inward E5 tube data amplified ~1e4; centre coefficients 3.5e-16 = cm record precision;
design N-a formula error |J_s|/|c_V| assumes J_c = 0 while J_c(x_m) = -0.0058227).

**Object.** As BAR-QB2. Inward background: an inward tube rebuilt from the BG2 sonic box (ed61d09268, V0 width 6.94e-29)
with the lin_tube_v2 method, or BG2 sonic-manifold coefficients intersected with the stored tube's (build seat's
choice; overlap with the old tube required).

instruments: seats/build-4/

## Predictions
- **PQ2b-1 (conservation, as sealed in QB2).** With the sharper inward background: matching mismatch e on the real
  source encloses 0 with width <= 1e-10; |D_t| <= 1e-10 and D_t encloses 0 on every inward step.
- **PQ2b-2 (centre orders, exact).** A symbolic (sympy) proof, from regularity of the background at the centre and the
  conformal-weight form of b_ab, that the theta-coefficients of b below order 4 (below 5 for b_nr) vanish identically,
  plus certified enclosures of the order-4 coefficients overlapping the QB1 centre value (c0 = -173.10044439229957).
  Replaces the numeric 1e-30 test (which the stored record cannot reach) by an exact statement.
- **PQ2b-3 (domain).** A written lemma (with the four-row system) that Psi2 is determined by data on x <= 0 only
  (centre regularity + sonic analyticity + gauge Nbar_p(0) = 0 + matching at x_m = -4.6), so the QB2 C1 points
  x = 0.5, 3 test a transport Psi2 does not use; C1 is re-run at x = -6, -4.6, -2, -0.25 and also at x = -0.01.
- **PQ2b-4 (controls completed).** C8 theta-theta residual encloses 0; C8' refuses; N-a: e excludes 0 and encloses the
  corrected prediction (J_c - J_s)/c_V; N-c refuses; C2 refuses at kappa*.
- **PQ2b-5 (headline replay).** Psi2 overlaps [3.4616358037141 +/- 1.15e-7], excludes 0, relative width <= 1e-3;
  route B overlaps.

## Kill conditions
- Any new enclosure disjoint from the QB2 / BG2 / E5 enclosures: FAILED and STOP.
- PQ2b-5 disjoint from QB2's Psi2: FAILED and STOP.
- Any single run over 2 hours, total over 5 hours: INCONCLUSIVE.

## Verdict
HOLD if PQ2b-1 .. PQ2b-5 hold.

## Seal fields
- formula_ref: DESIGN-QB2-v2.md:37-118
- calibrated_on: f2bcfadb3e seats/build-3 result_qb2.json, RUN.md - Psi2 value, e width cause, centre precision, N-a formula error seen before sealing (disclosed above)
- depends_on: f2bcfadb3e, 4d3b12ef39, ed61d09268, d0af75803f, b1b0ee1cde
- null_not_symmetry: N-c perturbs the background off shell; N-a predicts a specific nonzero mismatch; C8' scales a conserved source into a non-conserved one; C2 removes uniqueness
- target_form: arb / Taylor-model enclosures; exact sympy identities; overlap tests
- parameters: e width 1e-10; D_t bound 1e-10; Psi2 relative width 1e-3; C1 points -6, -4.6, -2, -0.25, -0.01
- tolerances: interval arithmetic (arb)
- ceiling: a fully closed certified first-order quantum push on the peak compactness of a true critical solution
- positive control: C1 manufactured phi G_2; C8 manufactured conserved source
- toy_model: C1 / C8 manufactured cases
- files: BAR-QB2b.md
- command: the build seat's own (seats/build-4/RUN.md), under memcap 8
