# CLAIMS: every claim with its status, seal, result file and exact number

**In plain words.** This page lists everything the work says, one claim per row, and for each one how sure we are
(proved by computer, proved on paper, quoted from a paper, only replicated, a conjecture, or not claimed), which sealed
test decided it, which file in this repository holds the number, and the number itself. If a row says NOT CLAIMED or
CONJECTURE, the release does not assert it. Below the table, a second table names the file behind every printed
number.

Status labels: **CERTIFIED** = computer-assisted enclosure (outward ball arithmetic) from a sealed run; **PROVED** =
written or exact symbolic proof checked by at least one internal hostile referee seat; **CITED** = taken from an opened
source (locator in paper.md, References); **REPLICATION** = an independent computation without an enclosure;
**CONJECTURE** = not derived; **NOT CLAIMED**. Seal ids are the lab ledger's claim hashes (`results/ledger-quantum.json`).
Every claim holds **in any Hadamard state respecting self-similarity (Zahn eq. 4); existence of such a state on
Evans-Coleman is not proved.** Result files are given by their repository path (the instruments wrote them there;
copies of the five headline files are in `results/`). Balls are printed `[mid ± rad]`: each is either the stored
string verbatim or a shortened form written by `tools/fill_balls.py`, which rounds outward (fewer midpoint digits,
radius enlarged by the shift and rounded up) and asserts that the printed ball contains the stored one.

## Claims

| # | Claim (paper section) | Status | Seal | Result file (field) | Exact number |
|---|---|---|---|---|---|
| C1 | On R ≡ 0, Zahn's state-independent tensor is V = B/60 for every coupling ξ (Prop. 0) | CITED consequence (Zahn l.347 + R = 0); checked symbolically | - | paper.md §1, §6 | V = B/60; ⟨T⟩ ⊃ N_c ħ (s/(240π²)) B |
| C2 | Field multipliers N_c positive: 1 per scalar, 6 per Dirac fermion, 12 per Maxwell field (§3) | CITED (BO l.128, eq. 3.4a; Duff eq. 30-31) | - | - | 1 : 6 : 12 |
| C3 | Sonic point is not Bach-flat (Thm 1.2) | CERTIFIED | 54caa4012c (PQ1 HOLD) | results/result_qb1.json (sonic) | b_uu = [0.4466737796695366321 ± 1.3e-17]; b_nn = [0.4157372840770100272 ± 1.1e-17] |
| C4 | Centre law B_(nn) = c0/t⁴, isotropy b_rr/b_nn → 1/3, B_(nr) = 0 (Thm 1.3) | CERTIFIED | 54caa4012c (PQ3 HOLD) | results/result_qb1.json (centre.c0) | c0 = [−173.1004443922995663248427 ± 2.99e-13] |
| C5 | Exact closed form c0 = −(2/3)W0(8W0 − 3)e^{−4xc}; centre orders of b exact (§5, §7.3) | PROVED (sympy) + CERTIFIED value | 373f306e9f (PQ2b-2 HOLD) | results/result_qb2b.json (PQ2b-2) | c0 = [−173.1004443922995663248426 ± 8.79e-16]; referee-4 re-derived it |
| C6 | Routes A and B agree: 28 comparisons at 19 distinct x (tests the algebra, not independent physics) (§4.1) | CERTIFIED | 54caa4012c (PQ4 HOLD) | results/result_qb1.json (routes) | 28 overlaps; Einstein residual contains 0 everywhere |
| C7 | Controls: Hayward 1/80, 17/240, −7/480; RN; Minkowski; k = 1/2 Friedmann ξ-dependence; sign-flip null refuses (§6) | CERTIFIED | 54caa4012c (PQ5 HOLD) | code/quantum/seats/build-1/controls_qb1.json | as in paper §6 |
| C8 | Full sign map of b_nn on [−12, 12]: negative below x0, positive above (Thm 1.1) | CERTIFIED | 4d3b12ef39 (PQ2b HOLD) | results/result_qb1b.json (PQ2b) | 2261 cells: 370 negative, 1891 positive |
| C9 | Exactly one zero x0 of b_nn on [−12, 12] (Thm 1.1) | CERTIFIED | 4d3b12ef39 (PQ2b HOLD) | results/result_qb1b.json (PQ2b.zero) | x0 ∈ [−0.82611434347927570343, −0.82611433975398540497]; b_nn′ ≥ 1.1654794 |
| C10 | Off-shell content null refuses (§6) | CERTIFIED | 4d3b12ef39 (PQ6b HOLD) | code/quantum/seats/build-2/pq6b.json | residual excludes 0 at 24 of 28 points; routes disjoint at 8 of 28 |
| C11 | κ = 2 is not a zero of D_new (bridge to QB2) (§4.4) | CERTIFIED | 4d3b12ef39 (PQ7b HOLD) | code/quantum/seats/build-2/pq7b.json | D_new(2) = [−2.4126270857981028564 ± 1.13e-9]; D_HKA(2) = [−19.942295267987575083 ± 1.12e-8] |
| C12 | Centre series and tube agree at x = −9/2 (§4.3) | CERTIFIED | 4d3b12ef39 (PQ8 HOLD) | code/quantum/seats/build-2/pq8.json | centre b_nn = [−2.6300024457005e-6 ± 1.4e-11] |
| C13 | Forced response e^{2s}(h1 + s·h2): h2 unique modulo G_2, no s²e^{2s} term (Lemma 2) | PROVED | lemmas refereed (referee-2 v1, referee-4 v2); amendment 4cf7213207 | results/design/PROOF-QB2-LEMMAS-v2.md; results/seals/AMENDMENT-QB2-v3.md | ρ(2) ≤ −1.9085 (PROVED from λ′ and d); D_HKA(2) = [−19.942295267987582258 ± 1.07e-11] |
| C14 | Ψ2 depends only on x ≤ 0 data (domain lemma) | PROVED | 373f306e9f (PQ2b-3 HOLD) | code/quantum/seats/build-4/DOMAIN-LEMMA.md | - |
| C15 | x_V unique in U = [−0.30, −0.20], V′ > 0 on U (Thm 3) | CERTIFIED | f2bcfadb3e (PQ2-7 HOLD); replayed 373f306e9f | results/result_qb2b.json (PQ2b-5.N7) | x_V = [−0.25278710208307181878 ± 1.14e-20]; A(x_V) = [1.9143464252607319765 ± 2.91e-20] |
| C16 | **Headline: Ψ2 = Ā_2(x_V)/A(x_V), positive** (Thm 3) | CERTIFIED | f2bcfadb3e (PQ2-8 HOLD), 373f306e9f (PQ2b-5 HOLD), 7aab7efe5a (PQ2c-3 HOLD) | results/result_qb2b.json (PQ2b-5.N8.Psi2); results/result_qb2c.json (PQ2c-3) | Ψ2 = [3.4616358037141656915 ± 3.02e-12] per unit K·ε·s·e^{2s} in the Koike-Hara-Adachi clock, K = N_c/(30π), ε = ℓ_P²/t0²; Ψ2 = δ(2m/r) at x_V. First sealed run: [3.4616358037141459203 ± 1.15e-7] (results/result_qb2.json PQ2-8) |
| C17 | Route B (variation of constants + adjoint) overlaps route A (Thm 3) | CERTIFIED | 373f306e9f (PQ2b-5) | results/result_qb2b.json (PQ2b-5.C6) | [3.4616358037141656727 ± 1.03e-11] |
| C18 | Photons (N_c = 12): K·Ψ2; per unit N_c: Ψ2/(30π) | CERTIFIED (arithmetic on C16) | as C16 | as C16 | 0.440749159476; 0.036729096623 (arb: [0.440749159476 ± 4e-13], [0.0367290966230 ± 3.4e-14]) |
| C19 | KHA clock: x_s = e^{−xc} | CERTIFIED (closed form) | 7aab7efe5a (PQ2c-4) | results/result_qb2c.json (PQ2c-4.x_s) | [0.80866895631844282988 ± 6.43e-19] |
| C20 | Central-proper-time coefficient per unit N_c (u·e^{2u}) | REPLICATION-grade conversion (float) | - | code/quantum/seats/referee-3/RECONCILE-QB2.md | 0.0561653802878 |
| C21 | Constraint propagation: D_t ≡ 0 symbolically; on every inward step certified D_t series contain 0 | PROVED (sympy) + CERTIFIED | 373f306e9f (PQ2b-1, sub-bar passes) | results/result_qb2b.json (PQ2b-1.Dt) | max abs D_t ≤ 3.3606426921349536e-13 |
| C22 | θθ residual on the real source encloses 0 at x = −6, −4.6, −2, −0.25, −0.01; N-b refuses at −6, −2, −0.25, −0.01, and at −4.6 only on the centre-side representation chosen POST HOC (the sealed sonic-side N-b there, [5.81641278197e-9 ± 1.67e-8], does not exclude 0) | CERTIFIED | 7aab7efe5a (PQ2c-2 HOLD) | results/result_qb2c.json (PQ2c-2) | e.g. [0 ± 4.81e-11] at x = −6 |
| C23 | Conservation mismatch e at matching encloses 0 | CERTIFIED | 373f306e9f, 7aab7efe5a | results/result_qb2c.json (PQ2c-1.N6) | e = [2.2796218038245719129e-14 ± 3.33e-10] |
| C24 | Conservation width e ≤ 1e-10 | **NOT CLAIMED** (FAILED as sealed in QB2 PQ2-3 and QB2b PQ2b-1; QB2c PQ2c-1 FAILED, not run as sealed: a probe showed the tight stage cannot move the width, and the QB2b value is replayed) | f2bcfadb3e, 373f306e9f, 7aab7efe5a | results/result_qb2c.json (PQ2c-1) | width 6.66e-10; QB2: 1.98e-5. Floor, as a probe estimate (not a certificate): pinned tube step balls (KT = 36) |
| C25 | Controls and nulls of QB2 (C1, C2, C6, C8, C8′, N-a corrected, N-c) as predicted | CERTIFIED | 373f306e9f (PQ2b-3, PQ2b-4 HOLD); in QB2 PQ2-5, PQ2-6 FAILED as sealed (partial) | results/result_qb2b.json (PQ2b-3, PQ2b-4) | N-a: [−302.35960932879925242 ± 2.06e-10] per unit δ vs prediction [−302.35960932879926987 ± 2.79e-11]; N-c e = [−0.00011949776035180548206 ± 3.56e-10] |
| C26 | Blind replication: 30π x_s² Ψ_blind / A(x_V) agrees with Ψ2 | REPLICATION (double precision, no enclosure); conversion POST HOC | 7aab7efe5a (PQ2c-4 HOLD) | results/result_qb2c.json (PQ2c-4.rows); code/quantum/seats/blind-qb2/blind_qb2_result.json | relative differences −6.9153490713739357298e-12 ... 8.5734474244780018923e-13 (six runs); 2.0992124976238618617e-10 vs the 10-digit blind headline |
| C27 | Bach source rebuilt by the QB2 instrument agrees with QB1 | CERTIFIED | f2bcfadb3e (PQ2-1 HOLD) | code/quantum/seats/build-3/n0check.json | all overlap (sonic, centre, 8 tube points) |
| C28 | The sign of ħ's effect on the black-hole threshold: quantum effects would favour black-hole formation (lower the threshold) | **CONJECTURE** (Conjecture 4; needs IR naturalness (i)-(iii), linear extrapolation, sgn δ(2m/r) at x* = sgn at x_V) | - | paper.md §7.4 | - |
| C29 | Existence of a self-similar Hadamard state on Evans-Coleman; the state-dependent part Ψ1; spin > 1; beyond first order in ħ; x > 12; the Cauchy horizon | NOT CLAIMED | - | - | - |
| C30 | Background: E5 existence and M count (published) and the BG2 record (not published, shipped here) | CERTIFIED (inherited) | d0af75803f, b1b0ee1cde, ed61d09268 | code/sideways/build-5/result_bg2.json | BG2 V0 width 6.94e-29 inside E5 (certified box, `b1_C.log`; as serialised into the outward record, half-width 3.80e-29, width 7.58e-29, `b2_C.log`) |
| C31 | Gauge rate κ_g = 2 − A0 + (2/3)ω0, removed by the deflation (paper §2) | CERTIFIED (inherited, M control C1c) | b1b0ee1cde | verdict_M.json l.3 in the published record [MJ-EC] (results/m/runs7/verdict_M.json there); quoted in results/design/DESIGN-QB1.md l.228 | κ_g = [0.3556992037 ± 7.41e-11] |

## Numbers and their sources

| Printed | Direction | Source file (field) |
|---|---|---|
| b_nn, b_nr, b_rr, b_T, b_uu at the sonic point (paper Thm 1.2) | ball | results/result_qb1.json (sonic) |
| c0 [−173.1004443922995663248427 ± 2.99e-13]; B_(rr) = c0/3 [−57.70014813076652210828089 ± 4.11e-13] (paper prints these rounded) | ball | results/result_qb1.json (centre) |
| centre approach [−173.1004442734853028947048 ± 0.00134] at x = −12, [−173.1003931711007216411606 ± 3.33e-6] at x = −9 (paper prints these rounded) | ball | results/result_qb1.json (sign_map.rows, centre.approach) |
| 691 grid points (278 negative, 413 positive) | count | results/result_qb1.json (sign_map) |
| 2261 cells, 370 / 1891, widths 9.5e-7 to 0.5, 2049 unresolved sub-cells, zero ends [−3e-9 ± 8.87e-10], [2e-9 ± 7.33e-10] | count / ball | results/result_qb1b.json (PQ2b) |
| seam centre b_nn [−2.6300024457005e-6 ± 1.4e-11], tube [−2.6300024457005e-6 ± 2.3e-9] | ball | code/quantum/seats/build-2/pq8.json |
| 25 digits sonic b_uu, c0 −173.1004443922995663248, residuals ≤ 2.1e-32 (independent code) | exact | code/quantum/seats/referee-1/REFEREE-QB1.md, ref_ec_sonic2.log, ref_ec_centre.log |
| QB2c deep step-matrix ratio tight/plain median 0.896, max 0.995 | probe | code/quantum/seats/build-5/deep_ratio.json; result_qb2c.json (PQ2c-1.deep_ratio) |
| BG2 V0 box [0.1124394013880983362237630276 ± 3.79e-29], width 6.94e-29 | ball | code/sideways/build-5/phase2/b1_C.log; result_bg2.json |
| x_s = 0.80866895631844 (blind float); A(x_V) = 1.91434642526073 (blind float) | probe | results/result_qb2c.json (PQ2c-4) |
| Prior-art floor: arXiv 2611.23200, INSPIRE citers 2611.21356, Zenodo record 24338153 | exact | results/prior-art/PRIOR-ART-REFRESH.md |

Tokens allowed without a result-file source (definitions and sealed constants): k = 1/3; 2 (the quantum rate); 1/80,
17/240, −7/480 (Zahn's Hayward values); 1e-10, 1e-30, 1e-3, 1e-6, 1e-12, 1e-9 (sealed bars); U = [−0.30, −0.20];
x_m = −4.6; x = −6, −4.6, −2, −0.25, −0.01 (sealed points); 1 : 6 : 12 (multipliers); 30π, 240π² (normalisations).
