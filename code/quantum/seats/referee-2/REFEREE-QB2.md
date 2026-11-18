# REFEREE-QB2 (referee-2, 2026-11-18)

I did not open any BAR-*.md file. My scripts and logs are in this folder:
- r2_forced.py: my own exact sympy derivation, in HKA bar variables at general k. Logs: r2_forced_k13.log, and r2_forced_k15.log (POST HOC).
- r2_pole.py and r2_pole.log (POST HOC).

## In plain words
The four lemmas hold up in structure, and the observable is better than the design claims: it is invariant under every gauge, not only the one gauge mode. One real instrument bug needs fixing before any build. In the form the M instrument transports (the HKA pencil), the quantum source also enters the two fluid rows. The design and the lemmas say it enters only rows 0-1. A build that follows the design would solve the wrong equations, and its conservation test would then refuse and put the blame on QB1. The manufactured-solution control C1 also falls outside the centre lemma. The threshold and ceiling wording claims more than the derivation supports.

## 1(a) Constraint propagation: CORRECT WITH FIXES
- **Matches the page.**
  - My constraint row, normalised to have κ-coefficient 1, equals HKA's identity row (κ - A, ...) exactly (K1, K5).
  - Comp/D_t is free of σ (K3).
  - σ_θθ drops out of D_t but not out of D_x (K3).
- **Fix 1: the forcing coefficient.** In HKA normalisation the coefficient is **-A**, not -A²: "lam == -A : True; lam == -A^2 : False" (K3, also at k = 1/5).
  - The ratio cited from lin_sym C3 is mom/I = Acons (gen_lin_fields.log l.8).
  - So "ϑ·Comp" with Comp = -A²·D_t gives -A³. Write it as J' = c J - A·D_t.
- **Fix 2: the rate depends on which evolution system is used.**
  - Along the derive rows (the system the proof's step 2 uses), the rate is c = (8NVω + 4V²ω - 3V² + 4ω + 3)/(3(V²-1)). This is not -A, and it equals -A + I·F (K2, K5).
  - c_I = -A holds only along the HKA pencil, and only with its full source.
  - Comp itself is identical on both routes (K5 True). Both rates are regular on x ≤ 0 and tend to -1 at the centre, so (ii)-(iv) stand.
- **(iii) is right.** For a tensor with only the θθ component, the divergence is (s, x) = (2, -2)·E_θθ e^{2s-2x}, summing to 0 (K4). So E_ab = 0 on the (s, x) block forces D_t = 0, and D_x ∝ E^θ_θ.
- **(iv) is right.** Zahn l.348 checked (OPENED): "it is automatically conserved".

## 1(b) Forced sonic step: CORRECT WITH FIXES
- **Constants checked.**
  - λ' ∈ [1.0986, 1.0995]: PROOF-B4-v4.md l.81.
  - The asserts are at fw_cert.py l.385, and g = [0.08962 +/- 1.23e-6] is at runs6/fw_cert_Ka55_K2000.log l.3.
  - Recomputed: |b_hi| ≥ 0.0984552 and ρ(2) ≤ -1.9085039.
  - The steps d = b_--(0), ℓM_x'r = λ'·ℓM_s r and λ_+ - λ_- all check by hand.
- **Fix: Δ2 is affine by structure, not by a test.** The seal's check is an interval containment of the second difference (count_b4v3_seal.py l.154), not "exact".
  - Affinity is structural: κ enters only the ℓ̂P row.
- **New (POST HOC, PROVED): κ_pole = d exactly.**
  - At the sonic point, c_f·r = -8Vω·D/(9(V²-1)²), where D ∝ det M_x (r2_pole.log). So the identity row annihilates r there.
  - Hence Δ2(κ) = const·(ℓE_u r - κ ℓM_s r), with const ≠ 0 wherever c_f(0) ≠ 0 (ω > 0).
  - This explains the 13-digit "coincidence" in sec. 3, item 5. It also lets Δ2(2) ≠ 0 follow from d < 0 alone, so the lemma no longer needs the seal-identification step.
  - PROOF-B4.md l.64's float Δ2 = -0.264162 - 2.667782κ has its root at -0.0990193. That is off from d by 1e-5, so that float line is inaccurate (legacy).
- **The lemma is fine for a general analytic S.** But see the BUG below: ℓ̂S(0) ≠ 0 for the real source.

## BUG (blocks the build): the HKA source has nonzero fluid rows
- The design relies on "S only in rows 0-1" in three places: PROOF-QB2 sec. 1, DESIGN sec. 1.1 / 1.2 step 1, and failure mode 4.
- The M instrument transports the HKA pencil, Q⁻¹(P0 + κP1) (lin_mc.py l.1-10, lin_stages.py l.68).
- In that pencil, on J = 0, S = Q(F σ_I + s4 σ). My result (K5):
  - row 2: [V/N, (N+6V)/(3N), -(2N²V - N - 3V)/(3N)];
  - row 3: [1/(3N), -(NV-2)/(3N), -(2N²+NV-1)/(3N)] (columns σ_ss, σ_sx, σ_xx);
  - rows 0-1 match srcprop exactly.
- srcprop.log l.16-17 is about the derive rows, where the fluid rows carry the source implicitly through A_x and A_s.
- Node N1 compares only against srcprop, so it would miss this. C1 and C6 share the same source module, so they would also pass with the wrong source.
- Only N6/N-a would fire, and they would point at QB1 instead.
- **Fix:**
  - use all four rows (log above) in the augmented column;
  - delete the remedy in failure mode 4;
  - add the rows-2,3 comparison to N1.

## 1(c) Centre step: CORRECT WITH FIXES
- **Fix 1: the premise is false.** S_c = (S_H, S_N, 0, 0) does not hold.
- **The conclusion probably survives.** My leading-order estimate from the K5 rows (POST HOC, not a proof): Q⁻¹S has orders (v², v², v², v³). Isotropy b_rr → b_nn/3 leaves S_3 = O(v) before Q⁻¹; M_x⁻¹ ~ 1/N then gives one more power of v.
- **Fix 2: σ_I is O(v²), not "O(v³) at worst".** From the frame algebra, σ_I ∝ r²A(-Σ_nn + NΣ_nr). This is harmless, because the lemma only needs o(1).
- **Fix 3:** extend the REDUCED item in sec. 7.2 to all four rows.

## 1(d) Dimension count: CORRECT (wording fix)
- 2 + 2 - 1 = 3 and the affine argument are sound.
- verdict_M.json l.4-6 checked: PM1 N_R40 = 1, κ* = [2.810552549 +/- 3.83e-10], c_V and Ā_ss'(0) nonzero.
- pq7b.json l.13-35 checked: D_new(2) = [-2.4126270857981028564 +/- 1.13e-9].
- **Fix:** step 4 and DESIGN sec. 1.1 say the source "is the certified QB1b profile". PQ2b is a sign map of b_nn (2261 cells with signs). QB2 must still produce certified TMs of b_nn, b_nr, b_rr and b_θθ.

## 2. Split: CORRECT
- **No s² e^{2s} term.** That would need a pole of the resolvent at κ = 2. G_κ exists at every κ, so it is gauge, not resonance.
- **The e^{2s} level is non-conserved.** ∇(uV) = (∇u)V (Zahn l.348, OPENED).
- **The ambiguity of h1 lies along h2 and G_2.** Checked; this includes BO l.167 (OPENED).

## 3. Observable
- **Gauge invariance: CORRECT, and stronger than claimed.**
  - 2m/r = 1 - 1/A depends only on x, and A' = 0 at x_V (checked by hand). So ∇(2m/r) vanishes as a spacetime vector on the whole line x = x_V.
  - Hence δ(2m/r)(x_V) is invariant under every infinitesimal diffeomorphism. Ā_ss'(x_V) = 0 is the special case.
  - **Fix:** the Danskin step for the global max over ℝ needs node T1, since h2 may grow for x > 12. Define Ψ2 through the local max near x_V instead.
- **"Ψ2 > 0 means toward trapping": CORRECT.** I re-derived (∇r)² = 1/A from the inverse metric, which gives δ(2m/r) = Ā_p/A.
- **Threshold sign: CORRECT WITH FIXES.**
  - The algebra holds: sΨ2 + Ψ1 = Ψ2 ln(L*/(-t)), and -t* ~ ℓ_P(ln)^{1/2}. But L* is defined by Ψ1/Ψ2.
  - **Fix 1:** "IR naturalness" also needs |Ψ1^state/Ψ2| ≪ ln(t0/ℓ_P). Failure mode 2 itself warns that Ψ2 could be small through core/shell cancellation.
  - **Fix 2:** state two more limits. At T* the perturbation is O(1), which is the TW linear extrapolation, a model step. The curvature there is ~1/(ℓ_P² ln), only a log below Planck.
  - **Fix 3:** in sec. 10, "the direction in which Planck's constant moves the threshold" is not derived. The threshold shift depends on f at x*, not at x_V. Label it a conjecture.
- **TW reduction: CORRECT.**
  - Re-derived by hand: eqs. 6.10-6.11 give 6.12 as B = κ_c f0/(f_q(q - κ_c)) and 6.15 (tw2.txt l.4158-4212, OPENED).
  - Note: when κ_c > q and f_c < 0, existence of the AH already forces "lowered". The readout is therefore really "does f2 < 0 occur at x*", which supports the design's gauge caveat.
- **Citation fixes.**
  - LOG.md l.40 is E4 (FAILED); the E5 "one V zero" record is l.45.
  - l.17 is an uncertified first attempt.

## 4. Controls and nulls: CORRECT WITH FIXES
- **C1 is outside Lemma (c).**
  - At the centre, G_2 → (0, 1, 2, 0), because N̄_ss' → -1 and ω̄_ss' → 2.
  - So x·G_2 grows like ln v, and the source G_2(centre) ∉ Im 𝓜_c(0) = span{(-1,1,-3,0), e_V}. The O(v²) centre series cannot return it.
  - Fix: use φ(x)G_2 with φ = e^{2x}/(1+e^{2x}) - 1/2. Then the source is Qφ'G_2 = O(v²), φ(0) = 0, and every row is exercised.
- **Missing content control.** Nothing exercises the σ → S map in every row with a known answer. Add a manufactured conserved Σ: pick σ_sx and σ_xx, solve D_t = 0 for σ_ss and D_x = 0 for σ_TT. N6 must then enclose 0.
- **N-c cannot refuse.** Its 1e-12 perturbation of ω sits under the 1e-10 width (doctrine #20). Use 1e-6.
- **N-a should print its expected |e|.**
- **C5, C7 and C2 test conventions or wiring.** C3 and C4 test algebra only, as the design says. The content-bearing items are N6/N2 on the real source, N-a, and C6.

## Blocks the build (nothing is up for release)
1. The four-row HKA source.
2. C1 moved inside Lemma (c), and the manufactured conserved-Σ control added.
3. Fixes to the -A coefficient and rate wording, and the PQ2b scope.

**Next strike:** feed the K5 four-row source into N1, and certify the centre orders of all four rows of Q⁻¹S (the REDUCED item). Both are exact and take minutes.
