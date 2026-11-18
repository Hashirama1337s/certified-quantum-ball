# BAR-SF-BG2: "Sharper Cone" - a tighter certified background at the light cone (sub-boss SB2c, second attempt)

Written 2026-11-17 by the lead, after BAR-SF-BG (7a68ad7998) settled FAILED (lever was wrong). Disclosed, calibrated_on: the
build-5 Step 0 diagnosis (build-5/DIAG-BG.md, result_bg.json, NOT SEALED, non-rigorous thin-input runs): the stored
centre-to-sonic tube is not the limiter; the limiter is the outward background record near the light cone
(build-1/outward_rec.pkl near root step 387, relative width 2e-16 to 9e-16; x_lc known to +/- 4.27e-12) and the light-cone
Frobenius vector (relative width 1.7e-14). With both thinned to relative width 1e-20 (not rigorous) the gate radius
was 4.2e-3; at 1e-18 it was 8.2e-2. Thresholds below are set from that ladder, not from a new run.

**Object.** The same certified Evans-Coleman solution (E5, seal d0af75803f). Rebuild, rigorously and in new files only:
a refined E5 sonic enclosure, the outward leg (driver_v2.outward_leg_v2 logic) from it, the light-cone root x_lc and the
state at x_lc, and the light-cone Frobenius vector, at whatever precision / order the build seat chooses.

instruments: build-5/

## Predictions
- **PC1 (refined sonic box).** Certified V0 enclosure of width <= 1e-24, inside the E5 box
  [0.1124394013880983362 +/- 3.31e-20].
- **PC2 (sharper cone).** Rebuilt outward record: relative half-width <= 1e-20 for every stored background quantity on
  every step with distance to the cone d in [0.012, 0.063]; the state at x_lc and the light-cone Frobenius vector also
  <= 1e-20 relative; x_lc enclosure width <= 1e-20; every new enclosure overlaps the old one (same solution).
- **PC3 (the control becomes resolvable).** With the rebuilt pieces in place of the old ones (everything else as in
  build-4 l2b_cert.py), the 16-point precision gate on the sealed control disc |kappa - (-1.8 + 6.2i)| <= 0.5 passes:
  relative radius <= 1e-2 at all 16 vertices, for at least one matching distance. Probe only; no count under this bar.
- **PC4 (positive control).** m/r at t = 0 from the rebuilt outward leg overlaps E5 PE3 [0.0595 +/- 8.24e-5]; x_lc
  overlaps 1.357941265765732 within the old enclosure; the old as-is gate (1.0) is reproduced when the old pieces are
  swapped back in.
- **PC5 (null, teeth).** Swapping in a rebuilt record started from V0 shifted by 1e-22 (inside the E5 box, outside the
  PC1 box) must make the refinement refuse, or its outward pieces must be disjoint from PC2's; it may not certify the same
  enclosure.

## Kill conditions
- PC1 or PC2 widths not reached: FAILED for that item; PC3 cannot then be claimed.
- Any new enclosure disjoint from the published E5 enclosures: FAILED and STOP, report to the owner first.
- Any single run over 2 hours, or 4 hours in total, outside the evening rule: stop, INCONCLUSIVE.

## Verdict
HOLD if PC1-PC5 hold. On HOLD the rebuilt cone pieces replace the old ones for the SB2 control re-count and SB3/SB4.

## Seal fields
- formula_ref: ../engine/cert/run_v32.py:96-97
- calibrated_on: 7a68ad7998 build-5 Step 0 (DIAG-BG.md, result_bg.json) - width ladder 1e-16/1e-18/1e-20/1e-24 seen before sealing (disclosed above)
- depends_on: d0af75803f, b1b0ee1cde, fc89881caf, 7a68ad7998
- null_not_symmetry: PC5 builds from a point that is not the solution; a correct certificate cannot reproduce the true enclosure from it
- target_form: arb enclosures with width bounds; overlap with old enclosures; gate relative radii at 16 vertices
- parameters: V0 width 1e-24; cone relative half-width 1e-20 on d in [0.012, 0.063]; gate 1e-2 at 16 points; null shift 1e-22
- tolerances: interval arithmetic (arb)
- ceiling: a background sharp enough to certify decaying-side controls, unblocking the dipole and quadrupole counts
- positive control: PC4 anchors and the swap-back reproduction
- toy_model: none - anchors are previously certified values on the same object
- files: BAR-SF-BG2.md
- command: the build seat's own (build-5/RUN-BG2.md), under memcap 8
