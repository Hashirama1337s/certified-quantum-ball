# AMENDMENT-QB2-v3: text fixes to DESIGN-QB2-v2.md and PROOF-QB2-LEMMAS-v2.md (hashed files are not edited)

Written 2026-11-18 by the lead from referee-4 (seats/referee-4/REFEREE-QB2-v2.md) and referee-3 (RECONCILE-QB2.md C4a).
No certified number changes. Each fix names the line it replaces.

## F1. The N-a expected value (PROOF-QB2-LEMMAS-v2 l.110-116; DESIGN-QB2-v2 l.147)
Replace "J_c ≡ 0, |e| = |J_s|/|c_V|" with
  e = (J_c(x_m) - J_s(x_m)) / c_V(x_m),
  J_c(x_m) = δ ∫_{-∞}^{x_m} exp(-∫_y^{x_m} A) (2 - A) S_0 dy   (J_c -> 0 at the centre),
  J_s(x_m) = -δ ∫_{x_m}^0 exp(∫_{x_m}^y A) (2 - A) S_0 dy.
For a non-conserved source the centre side carries a constraint error too. Certified, per unit δ:
- the build gives -302.35960932879925242 ± 2.06e-10, enclosing (J_c - J_s)/c_V = -302.35960932879926987 ± 2.79e-11
  (seats/build-4/result_qb2b.json, PQ2b-4 N_a);
- referee-4 recomputed it independently as [-302.35960932880 ± 3.82e-12] (r4_nums.log), with the J_c share
  [-0.03018779672 ± 1.86e-12].
Also on PROOF l.112, delete "·K": K is dropped at l.53-54.

## F2. The factor K = N_c/(30π)
The response to Σ_V = K ħ s B (per unit ε = ℓ_P²/t0²) is h = K ε e^{2s}(h1 + s h2). Accordingly:
- DESIGN-QB2-v2 l.97 and l.171: δm_U = K ε e^{2s}(sΨ2 + Ψ1) + O(ε²).
- DESIGN-QB2-v2 l.222 (node N8): "(times K ε s e^{2s}, K = N_c/(30π))", replacing "times N_c ε".
- Paper draft l.131: h = K ε e^{2s}(h1 + s·h2), applied in the draft.

## F3. Labels
- PROOF-QB2-LEMMAS-v2 l.183-186: the source orders are now SYM. Cite both proofs:
  - seats/build-4/centre_exact.log l.14-19;
  - the independent re-derivation in seats/referee-4/r4_bach_12_on.log and r4_bach_14_on.log.
- The orders, the 1/3 isotropy and b_nr_5 = (4/3) b_nn_4 are kinematic: they follow from smoothness together with
  tr B = 0 and ∇·B = 0. Only the W0-dependence is content:
  b_nn_4 = -(2/3) W0 (8 W0 - 3), c0 = b_nn_4 e^{-4 xc}.
  That is why the N-c null cannot refuse an order claim (doctrine #20). The off-shell POST HOC null in r4_bach_12_off.log
  refuses the value.
- DOMAIN-LEMMA l.4 and l.18: ρ(2) ≤ -1.9085 is PROVED from λ' (PROOF-B4-v4 l.81, CERTIFIED) and d, via PROOF-B4.md
  l.117. It is not itself certified.
- DOMAIN-LEMMA l.9-10: replace the claim about changing x > 0 data with the operational claim "the instrument consumes
  only x ≤ 0 data".

## Not checked by referee-4
- HKA rows 2-3 rest on referee-2 K5.
- n1_rows = K5.
- κ_pole = d rests on r2_pole.log.
