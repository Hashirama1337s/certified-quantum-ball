# REFEREE-QB2-v2 (referee-4, 2026-11-18)

I opened no BAR file before I wrote my scripts. I then read BAR-QB2b/QB2c (QB2c.md l.3-7) to find out what was already disclosed. My scripts are in this folder, all run under memcap 4, with peak memory ≤ 0.08 GB:
- r4_centre.py: the Einstein-fluid centre series, derived from the metric. It uses nothing from engine/ or build-3/.
- r4_bach.py: the Bach tensor computed from its definition, B = ∇∇C + ½RC. I ran it three ways: KMAX 12 on-shell, KMAX 14 on-shell, and KMAX 12 off-shell (POST HOC).
- r4_rows.py: the frame-to-coordinate source rows.
- r4_nums.py: arb rechecks of the numbers.

## In plain words
The four lemmas are sound, and the centre formulas are exactly right. I derived them from scratch: my own metric series and my own Bach tensor give the same four results. One piece of the v2 proof text is wrong: the expected size of the deliberately broken null N-a. It assumes the centre side carries no constraint error, and the certified run shows that it does. The K = N_c/(30π) fix is in the proof page and in the paper's Theorem 3. It is not yet in the design page or in one line of the paper. These are text fixes. No certified number moves.

## 1. PROOF-QB2-LEMMAS-v2

**Sec. 1, frame to coordinates and HKA rows 0, 1, σ_I: CORRECT.**
- r4_rows.log reproduces the following with difference 0:
  - σ_ss, σ_sx, σ_xx (l.55);
  - row 0 = Aβ_nn and row 1 = -A(β_nn-β_rr)/2 (l.63-64);
  - σ_I = A(Nβ_nr - β_nn) (l.67).
- My metric from ∂_s = √A r(N n - e_r) gives g_ss/r² = A(1-N²), which matches engine/derive.py l.29-32.
- Rows 2-3 were NOT CHECKED by me. They rest on referee-2's K5.

**Lemma (a): CORRECT WITH FIXES.**
- Items (i)-(iv) apply REFEREE-QB2's Fixes 1-2 as written (l.76-78, l.92-97).
- **WRONG: the N-a expected size (l.111-116).** The text says "J_c ≡ 0" and |e| = |J_s|/|c_V|.
  - With a non-conserved source, J_c is set by the centre condition, not by zero.
  - result_qb2b.json PQ2b-4 N_a:
    - J_c_xm = -0.0058226690727984285785;
    - the design formula gives 302.32942153207741039;
    - the run gives -302.35960932879925242 per unit δ;
    - "design_formula_matches": false.
  - r4_nums.log: (J_c - J_s)/c_V = [-302.35960932880 ± 3.82e-12]. The J_c share alone is [-0.03018779672 ± 1.86e-12].
  - BAR-QB2b l.7 disclosed this error. The v2 page still carries it, and so does DESIGN-QB2-v2 l.147.
  - **Fix:**
    - e = (J_c(x_m) - J_s(x_m))/c_V(x_m);
    - J_c(x_m) = δ∫_{-∞}^{x_m} exp(-∫_y^{x_m}A)(2-A)S_0 dy (J_c → 0 at the centre).
- **Minor:** l.112 writes "Aβ_nn·K", although l.53-54 drop K. Delete the "·K".

**Lemma (b): CORRECT. I checked the constants by citation only.**
- ρ(2) ≤ -1.9085 follows from Re ρ(κ) ≤ (d - Re κ)/λ' (PROOF-B4.md l.117), with λ' ∈ [1.0986, 1.0995] (PROOF-B4-v4.md l.81) and d ≤ -0.0984.
- I did not re-derive κ_pole = d. That identity rests on r2_pole.log.

**Lemma (c): CORRECT WITH FIXES.**
- **Fix:** l.183-186 still label the source orders "PAPER". build-4 centre_exact.log l.14-19 now has them SYM, and I confirm them independently (sec. 3 below). Update the label and cite both.
- POST HOC, from the exact centre values: e0 = (4/3)a0 and r0 = a0/3.
  - Applied to the l.172-178 table, the leading terms are σ_I = a0/3 and row 3 = -a0/12. Both are nonzero, so the orders are sharp.

**Lemma (d): CORRECT.** D_HKA(2) was reproduced on build-4's own transport: [-19.942295267987582258 ± 1.07e-11] (result_qb2b.json l.174-178).

**The K = N_c/(30π) fix (RECONCILE-QB2 C4a): applied only in part.**
- Applied:
  - PROOF-v2 l.53-54;
  - paper l.142 and l.340;
  - result_qb2b.json N8 "units".
- **Not applied**, and there is no v3 note or amendment on disk:
  - DESIGN-QB2-v2 l.97 and l.171 (δm_U = ε e^{2s}(…), with no K);
  - DESIGN-QB2-v2 l.222 ("times N_c ε", wrong by 1/(30π));
  - paper l.131, which writes the response to the full Σ_V as h = ε e^{2s}(h1 + s h2). With the K at l.142, that line double-counts K. It should read h = Kε e^{2s}(…).
- The scale is consistent. 8π·uB/(240π²) = uB/(30π) (DESIGN-QB1 l.~58-61), and the dimensions match: ħB ~ ℓ_P²L⁻⁴ = L⁻² = [G_ab].

## 2. DOMAIN-LEMMA (Ψ2 uses only x ≤ 0): CORRECT WITH FIXES
The logic holds:
- the bounded-at-0⁻ condition is equivalent to analyticity, since Re ρ(2) < -1.9;
- the centre condition, the matching at x_m = -4.6 and the readout at x_V ∈ U all lie in x < 0;
- det ≠ 0 on the build's own tube (above).

Fixes (wording only):
- (i) l.4 and l.18 call ρ(2) ≤ -1.9085 "CERTIFIED, PROOF-B4-v4 l.81". That line certifies λ' only. The bound itself is PROVED from λ' and d (PROOF-B4.md l.117).
- (ii) For the analytic EC background, the statement "changing data on x > 0 leaves Ψ2 unchanged" (l.9-10) is vacuous, because analytic continuation already fixes the x > 0 data. State the operational claim instead: "the instrument consumes only x ≤ 0 data".
- (iii) The claim n1_rows = K5 (l.15) was NOT CHECKED by me.

## 3. Centre-order claims (centre_exact.py / .log): CORRECT, independently re-derived
My derivation works in Schwarzschild-like coordinates with central proper time.
- In these coordinates v = z + O(z⁵), because A2/2 = c2.
- Also e^{xc} = α0 (the KHA central lapse) and W0 = 4πρ_cτ² = m0/2.
- My background series agrees with the engine's: v1 = -1/2, v3 = (8m0-21)/120 = 2W0/15 - 7/40, matching centre_exact.log l.11 (r4_centre.log).

Results, from r4_bach_12_on.log; KMAX 14 is identical:
- b_nn_4 = -2·W0·(8W0-3)/3, which is the claim. The v⁵ coefficient vanishes, so the next term is O(v⁶).
- b_nr_5 = -8·W0·(8W0-3)/9, which is the claim. The b_nr coefficient of v⁴ vanishes.
- b_rr_4/b_nn_4 = b_T_4/b_nn_4 = 1/3. The trace of B is 0.
- The free symbol v5 of the z⁶ metric terms does not enter any leading coefficient.
- c0, from the closed form on the E5 boxes: [-173.10044439229957 ± 4.55e-15] (r4_nums.log). This overlaps QB1's -173.1004443922995663248427 ± 2.99e-13.

**Instrument check.** The Weyl-trace residue appears only at r^KMAX (r¹² at KMAX 12, r¹⁴ at KMAX 14). It is truncation edge, not a bug.

**POST HOC null, outside the authors' checks.** I broke the Einstein equations off shell (A2 += 1/10, c4 += 1/7; r4_bach_12_off.log):
- the value changes, to b_nn_4 = -(2240W0² - 420W0 - 1749)/420, so the claim reads False and the null refuses;
- but the orders, the 1/3 isotropy and b_nr_5/b_nn_4 = 4/3 all survive.

**Scope fix (the claim is fine, its billing is not).** The orders, the isotropy and the ratio 4/3 are kinematic. They follow from smoothness plus the identities tr B = 0 and ∇·B = 0, which give b_nr_5 = (4/3)b_nn_4 by hand. Only the W0-dependence is content. That is also why result_qb2b.json's null_Nc rows for N4a contain 0 (l.221-240): that null cannot refuse an order claim (doctrine #20). PQ2b-2 should say so.

## What blocks release
1. The N-a formula in PROOF-v2 l.113-114 and DESIGN-QB2-v2 l.147. It must go into a new amendment file, not into an edit.
2. K is missing in DESIGN-QB2-v2 l.97, l.171 and l.222, and paper l.131 contradicts l.142.
3. PQ2b-1 has not yet held: the e width is 6.66e-10 against a bar of 1e-10. It is pending as BAR-QB2c.

**Single next strike:** write one amendment note (PROOF/DESIGN v3) that fixes the N-a expected value to (J_c - J_s)/c_V and adds K at the four lines above. Then run QB2c.
