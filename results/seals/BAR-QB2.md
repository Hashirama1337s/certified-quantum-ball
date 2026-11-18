# BAR-QB2: "Which Way the Quantum Pushes" - the certified first-order quantum response of the critical ball of light

Written 2026-11-18 by the lead. Design: DESIGN-QB2-v2.md (referee-2 fixes applied), lemmas PROOF-QB2-LEMMAS-v2.md.
Disclosed, calibrated_on: the float smoke values of QB1 (b profiles, x_V ~ -0.2529, A(x_V) ~ 1.9143 float) and the
certified QB1b sign map (4d3b12ef39) are known. **The sign of Psi2 is NOT predicted** - it is the open question; no
float estimate of Psi2 has been computed or seen.

**Object.** The s e^{2s} (u e^{2u}) coefficient h2 of the first-order response of the certified Evans-Coleman k = 1/3
solution to the state-independent one-loop source u B_ab / (240 pi^2) (times the positive field multiplier N_c, spin <=
1), spherical sector, solved with the four-row HKA source (referee-2 K5), matched with the D_HKA(2) matrix at x_m = -4.6
in the gauge Nbar_p(0) = 0. Headline observable Psi2 = Abar_2(x_V) / A(x_V), the change of peak compactness 2m/r at the
local maximum on U = [-0.30, -0.20] (invariant under every infinitesimal gauge transformation). Psi2 > 0 means toward
trapping.

instruments: seats/build-3/

## Predictions
- **PQ2-1 (source).** Certified Taylor models of b_nn, b_nr, b_rr, b_TT on every transport step and outward cell used,
  traceless check passes, and they overlap the QB1b enclosures where both exist (node N0).
- **PQ2-2 (source rows).** The four HKA source rows rebuilt by the build match srcprop.log l.13-14 (rows 0-1) and
  r2_forced_k13.log l.21-22 (rows 2-3) exactly (node N1).
- **PQ2-3 (conservation).** On the real source, D_t and the theta-theta residual enclose 0 on [-12, 12] (N2), and the
  matching mismatch e encloses 0 with width <= 1e-10 (N6).
- **PQ2-4 (centre).** The first two theta-coefficients of all four rows of Q^-1 S vanish on the certified centre manifold
  (N4a, enclosures within [-1e-30, 1e-30]).
- **PQ2-5 (controls).** C1: manufactured phi(x) G_2 returned exactly (enclosures contain it) at x = -6, -4.6, -2, -0.25,
  0.5, 3; C2: at kappa* uniqueness is refused; C6: route B overlaps route A for Psi2; C8: manufactured conserved Sigma
  gives e and the theta-theta residual enclosing 0, and its 1e-6-scaled variant refuses.
- **PQ2-6 (nulls).** N-a: e excludes 0 and encloses the printed predicted |e|; N-c: omega (1 + 1e-6) refuses.
- **PQ2-7 (the location).** V has exactly one zero x_V in U with V' > 0 on U (N7).
- **PQ2-8 (headline).** A certified enclosure of Psi2 that EXCLUDES 0 (sign reported, not predicted), with relative
  width <= 1e-3.

## Kill conditions
- PQ2-2 mismatch, any PQ2-5 control failing, or PQ2-6 null passing: FAILED and STOP (wrong equations or no teeth).
- PQ2-3 conservation fails on the real source with controls passing: FAILED and report (would question QB1's source).
- PQ2-8 enclosure contains 0: FAILED for PQ2-8 (core/shell cancellation below resolution) - report as such.
- Any single run over 2 hours or total over 6 hours: stop, INCONCLUSIVE. Evening-only for runs over 1 hour.

## Not claimed
The direction the black-hole threshold moves (CONJECTURE: depends on x*, needs the IR-naturalness assumption incl.
|Psi1_state / Psi2| << ln(t0 / l_P), linear extrapolation, near-Planck curvature); the e^{2s} level (state-dependent);
state existence (all statements "in any Hadamard state respecting self-similarity").

## Verdict
HOLD if PQ2-1 .. PQ2-8 hold.

## Seal fields
- formula_ref: DESIGN-QB2-v2.md:37-118
- calibrated_on: 4d3b12ef39 QB1b certified source and QB1 smoke floats known; no Psi2 value computed before sealing (disclosed above)
- depends_on: 4d3b12ef39, 54caa4012c, d0af75803f, ed61d09268, b1b0ee1cde
- null_not_symmetry: N-c perturbs the background off shell so the conservation check must refuse; C2 at kappa* removes the uniqueness the solution relies on; N-a predicts a specific nonzero mismatch
- target_form: arb / Taylor-model enclosures; exact row comparison; certified Psi2 enclosure with sign
- parameters: e width 1e-10; centre coefficient bound 1e-30; Psi2 relative width 1e-3; N-c factor 1 + 1e-6; C8 scale 1e-6; U = [-0.30, -0.20]; x_m = -4.6
- tolerances: interval arithmetic (arb)
- ceiling: the first certified sign and size of the one-loop quantum push on the compactness of a true critical-collapse solution, the same for every free massless field of spin <= 1
- positive control: C1 manufactured exact solution, C8 manufactured conserved source
- toy_model: C1 / C8 manufactured cases; C5 k = 1/3 Friedmann zero floor
- files: BAR-QB2.md
- command: the build seat's own (seats/build-3/RUN.md), under memcap 8, watch-run for anything over 10 minutes
