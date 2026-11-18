# DESIGN-QB2-v2: the forced O(hbar) response of the critical ball of light, and the functional that reads its push

thesis: main quest L2 (real gravity), LIGHTBALL Raid 3 follow-on: certified quantum back-reaction on the certified Evans-Coleman solution (SCOPE.md, approved 2026-11-17)
target: N8
reuses: E5 enclosure (lin_tube2_k1_3.pkl, BAR-LB-E5 HOLD d0af75803f, LOG.md l.45), BG2 outward record, M count instrument (count_b4v3 + CachedLinear2 L13_nsub2, seal b1b0ee1cde), QB1/QB1b Taylor-model cell code (seats/build-2/tmcell.py, schart.py), PROOF-B4 Facts 1.1-1.4 + Lemma 1.5 + Thm 1.6, PROOF-QB2-LEMMAS-v2 (a)-(d)
produces: certified forced-response instrument (augmented 5x5 transport with the four-row HKA source + forced sonic/centre series), certified TMs of b_nn b_nr b_rr b_TT along the transport, h2 profile on [-12, 12], the gauge-invariant local-max compactness functional Psi2, the conservation test for QB1 b-profiles

Design seat design-2, 2026-11-18. **v2 of DESIGN-QB2.md (v1 unedited).** Each change is marked
**[v2 2026-11-18, REFEREE-QB2 <item>]**, citing seats/referee-2/REFEREE-QB2.md.

Paper, reading and small checks, all non-counting, all in `seats/design-2/`:
- exact sympy: srcprop, xv, twlog, and in v2 centre_orders and dt_struct;
- the referee's own derivation: seats/referee-2/r2_forced_k13.log ("K1-K5") and r2_pole.log.

**Nothing on this page is certified by a QB2 run.**

## 0. In plain words

- **The question.** How does the critical ball of light answer the geometric quantum stress, which grows like s·e^{2s}?
- **The shape (PROVED, PROOF-QB2-LEMMAS-v2).** The answer is e^{2s}(h1 + s·h2). h2 is fixed by geometry, modulo the gauge
  mode G_2. h1 needs the quantum state and is ambiguous along h2.
- **The readout.**
  - (∇r)² = 1 - 2m/r = 1/A exactly, and A has a strict local maximum exactly where V = 0, at x_V ≈ -0.253.
  - Ψ2 := Ā_2(x_V)/A(x_V) is the s·e^{2s} coefficient of the first-order change of the **local** peak compactness.
  - Ψ2 is invariant under every infinitesimal gauge [v2 2026-11-18, REFEREE-QB2 sec. 3].
  - Ψ2 > 0 means toward trapping.
- **[v2 2026-11-18, REFEREE-QB2 BUG] The source column has four rows.** In the HKA pencil the M instrument transports,
  the quantum source enters the fluid rows too. v1 used only rows 0-1. A build following v1 would have solved the wrong
  equations, and its conservation test would have blamed QB1.
- **Honesty (sec. 6).** Ψ2 decides the s-secular part of the first-order push. Whether it also decides the threshold needs
  three things:
  - IR naturalness, now including |Ψ1^state/Ψ2| ≪ ln(t0/ℓ_P);
  - a linear extrapolation to O(1) perturbations;
  - curvature that is only a log below Planck.
  The threshold direction itself is a **CONJECTURE**.

## 1. What is computed

### 1.1 The forced problem and its source
- **Equations.** On x ∈ (-∞, 12], the HKA pencil with source: Q h2' = (P0 + 2P1)h2 + S(σ_V), I_2(h2) = σ_I(σ_V).
- **[v2 2026-11-18, REFEREE-QB2 BUG] Four-row HKA source** (per unit K ε s e^{2s}, K = N_c/(30π); r2_forced_k13.log
  l.19-22; PROOF-QB2-LEMMAS-v2 sec. 1), with β_ab = e^{-2x}b_ab:
  - σ_ss = A(N²β_nn - 2Nβ_nr + β_rr), σ_sx = A(Nβ_nr - β_rr), σ_xx = Aβ_rr, σ_TT = β_T;
  - S_0 = (1/N², 2/N², 1/N²)·σ = Aβ_nn;
  - S_1 = (-1/(2N²), -1/N², (N²-1)/(2N²))·σ = -(1/2)A(β_nn - β_rr);
  - S_2 = (V/N, (N+6V)/(3N), -(2N²V - N - 3V)/(3N))·σ;
  - S_3 = (1/(3N), -(NV-2)/(3N), -(2N²+NV-1)/(3N))·σ;
  - σ_I = A(-β_nn + Nβ_nr).
  - The augmented column is Q⁻¹S = (S_0, S_1, M_x⁻¹(S_2, S_3)).
  - S_2 and S_3 carry factors 1/N, never det M_x (K5 l.27). Only M_x⁻¹ meets the sonic 0/0, and only inside the pencil
    form the M instrument already handles. At x = 0 use Q Y' = ... directly, as in Lemma (b).
- **[v2 2026-11-18, REFEREE-QB2 1(d)] Source scope.**
  - The b_ab are the QB1 frame profiles r⁴B_(ab) (DESIGN-QB1 sec. 1.2).
  - QB1b's PQ2b is a **sign map of b_nn** only (result_qb1b.json). It is not a source.
  - **QB2 must produce certified TMs of b_nn, b_nr, b_rr and b_TT itself**, as Taylor models on every transport step,
    reusing tmcell.py / schart.py (node N0).
  - Tracelessness b_nn = b_rr + 2b_T is checked on every step.

### 1.2 Pipeline
1. **Augmented homogeneous system** [v2 2026-11-18, REFEREE-QB2 BUG]. Ŷ = (Y, z), z' = 0, with matrix
   [[Q⁻¹(P0+2P1), Q⁻¹S],[0,0]] and **all four rows of Q⁻¹S**. The M instrument's certified transport runs on it as a 5×5
   system.
2. **Sonic leg.** The forced Frobenius series (Lemma b).
   - Constants: |ρ(2) - n| ≥ n + 1.9085; Δ2(2) ≠ 0 from d < 0 alone via κ_pole = d; ℓM_s r ≠ 0.
   - ℓ̂S(0) ≠ 0 for the real source; it enters only the right-hand side.
3. **Centre leg.** The forced regular series in θ (Lemma c). The orders of Q⁻¹S are (v², v², v², v³) and σ_I = O(v²)
   (SYM, centre_orders.log), certified on the cm manifold in node N4a.
4. **Matching** [v2 2026-11-18: gauge stated consistently]. Gauge N̄_p(0) = 0. Sonic family y_s^p + b·s_2^HKA with
   (Ā,N̄)(y_s^p)(0) = (0,0). Then y_c^p - y_s^p = b·s_2^HKA + c·ρ1 + d·ρ2 + e·e_V at x_m = -4.6.
   - The matrix is D_HKA(2) = [-19.942295267987575083 ± 1.12e-8] (pq7b.json).
   - **e must enclose 0** (conservation; node N6).
5. **Readout.** Certify x_V in U = [-0.30, -0.20] as the unique zero of V there, with V' > 0 on U, hence A'' < 0 on U.
   Then Ψ2 = Ā_2(x_V)/A(x_V) on the sonic leg.
6. **Route B.** Variation of constants with the κ = 2 fundamental system and the adjoint (C6).

### 1.3 Interface to QB3
- h2 on [-12, 12] (the far side needs T1).
- h1 = h1^state + μh2 + cG_2, carried symbolically.

## 2. The observable

### 2.1 Trapping function
(∇r)² = 1/A exactly (PAPER; re-derived by the referee from the inverse metric, REFEREE-QB2 sec. 3). So
δ(2m/r) = Ā_p/A, and the apparent horizon is at (∇r)² = 0 (BO eq. 3.14, bo.txt l.235-241; TW eq. 6.8, tw2.txt l.4096;
both OPENED).

### 2.2 Headline functional Ψ2
- **Where A peaks.** Ā' = -(8/3)ωNV/(1-V²) (SYM, xv.log l.1), so A is stationary exactly at the zeros of V.
  - E5 certifies one V zero on the inward leg [v2 2026-11-18, REFEREE-QB2 sec. 3 citation fix: LOG.md l.45, not l.40,
    which is E4 FAILED].
  - v1's appeal to LOG.md l.17 for V > 0 outward is withdrawn: that line is an uncertified first attempt.
  - Float: x_V ∈ (-0.25456, -0.25128), A(x_V) ≈ 1.91434 (xv.log l.3).
- **[v2 2026-11-18, REFEREE-QB2 sec. 3 Fix] Definition through the local maximum.** Let U = [-0.30, -0.20], with V
  certified to have exactly one zero in U and V' > 0 on U (node N7). Then A has a unique, non-degenerate maximum on U, at
  x_V. Define:
  - **m_U(s) := max_{x∈U} (2m/r)(s, x);**
  - **δm_U = ε e^{2s}(sΨ2 + Ψ1) + O(ε²), with Ψ2 := Ā_2(x_V)/A(x_V)** (Danskin's theorem on the compact set U).
  - No global maximum over ℝ is used, so no far-side bound is needed. v1's Danskin step over ℝ would have needed T1,
    since h2 may grow for x > 12.
- **[v2 2026-11-18, REFEREE-QB2 sec. 3] Invariance under every gauge.** 2m/r = 1 - 1/A depends only on x on the
  background, and A' = 0 on the whole line x = x_V. So ∇(2m/r) = 0 there as a spacetime covector, and δ(2m/r)(x_V) is
  invariant under every infinitesimal diffeomorphism (L_ξ f = ξ·∇f = 0). Ā_ss'(x_V) = 0 (the G_2 case) is the special case.
- **Reading.** Ψ2 > 0: the log-secular quantum push raises the local peak compactness (toward trapping). Ψ2 < 0: it lowers it.
- **Field independence and conditions.** The sign is the same for every free massless field of spin ≤ 1, "in any Hadamard
  state respecting self-similarity (Zahn eq. 4); existence on Evans-Coleman not proved" (RESULT-QB1 rules 2-3).
- **Source sentence.** BO: the argument "relies heavily on the assumption that quantum matter tends to oppose black-hole
  formation" (bo.txt l.277-278, OPENED). BO's abstract: inconclusive for β < 0.5 (bo.txt l.4). For EC, β = 0.3558
  (verdict_M.json l.5).

### 2.3 Tomasevic-Wu route (secondary)
- **The extended algebra** (SYM, twlog.log; the referee re-derived 6.12 and 6.15 from tw2.txt l.4158-4212, OPENED):
  Δp = ε e^{(q-κ_c)T*}[q(f1 + T*f2) + f2]/(κ_c f_c).
- **Referee note.** For κ_c > q with f_c < 0, existence of the AH already forces "lowered". The readout is therefore
  "does f2 < 0 occur at x*".
- **The gauge caveat stands.** f2 → f2 - cĀ_ss'/A away from x_V.
- **[v2 2026-11-18, REFEREE-QB2 sec. 3 Fix 3]** The threshold shift depends on f at x*, not at x_V. Equating sgn Ψ2 with
  the direction the threshold moves is a **CONJECTURE**, not derived.

## 3. Controls

| # | control | expected | what it tests |
|---|---|---|---|
| C1 | **[v2 2026-11-18, REFEREE-QB2 sec. 4]** manufactured φ(x)·G_2, φ = e^{2x}/(1+e^{2x}) - 1/2: source S_m = φ'·Q G_2 (four rows), σ_m = 0 | solver returns φ·G_2 exactly in the gauge N̄_p(0) = 0 (φ(0) = 0), at x = -6, -2, -4.6, -0.25, 0.5, 3 | inside Lemma (c): Q⁻¹S_m = φ'G_2 = O(v²), and G_2 → (0,1,2,0) at the centre; v1's x·G_2 grew like ln v and fell outside the lemma |
| C2 | resonance null at κ* = [2.810552549 ± 3.83e-10] | uniqueness refused (determinant encloses 0) | that the D(2) ≠ 0 step has teeth |
| C3 | BO structure (bo.txt l.163-170, 229-230, 267, OPENED) | κ_q = 2; s·e^{2s} forced; BO case (ii) | algebra and structure only |
| C4 | TW closed-form Roberts (eq. 5.32, tw2.txt l.3163; PDF check before use) | the extended algebra with f2 = 0 reproduces TW's signs | algebra only: 2D, minimally coupled scalar, Roberts, no log term. Not a control of Ψ2 |
| C5 | k = 1/3 Friedmann zero floor | S ≡ 0, so h2 ≡ 0 mod G_2 | wiring |
| C6 | route B (variation of constants + adjoint) vs route A | overlap of Ψ2 | double implementation |
| C7 | linearity | Ψ2[-S] = -Ψ2[S] | wiring |
| C8 | **[v2 2026-11-18, REFEREE-QB2 sec. 4] manufactured conserved Σ** | N6's e encloses 0 (width ≤ 1e-10) and the θθ residual encloses 0 | exercises the σ → S map in **every row** with a known answer |

**C8 construction.**
- Choose explicit fake profiles β_nr = P³-type and β_rr = P²-type, for example β_rr = e^{2x}/(1+e^{2x})² and
  β_nr = e^{3x}/(1+e^{2x})².
- Solve **D_t = 0 for β_nn**. D_t contains -σ_ss'/(AN²) with a nonzero coefficient, and no σ_TT (dt_struct.log). This is
  a linear first-order ODE: take its regular centre solution, transported with the same TM machinery.
- Solve **D_x = 0 for σ_TT algebraically** (coefficient -2, dt_struct.log).
- Then the manufactured Σ is conserved by construction.
- Variant C8': scale the solved β_nn by (1 + 1e-6). Then N6 must refuse, with |e| as in N-a.

In the sources opened, there is no published closed form for the V-forced 4D response.

## 4. Content nulls

- **N-a, non-conserved source (primary)** [v2 2026-11-18, REFEREE-QB2 sec. 4: prints its expected |e|].
  - Multiply HKA row 0 of the source column by (1 + δ), δ = 1e-6.
  - **Expected:** |e| = |J_s(x_m)|/|c_V(x_m)|, with J_s(x_m) = -δ∫_{x_m}^0 exp(∫_{x_m}^y A)(2 - A)S_0 dy (PROOF-QB2-LEMMAS-v2
    sec. 2, Consequence) and c_V(x_m) = 0.19288155165634 (verdict_M.json l.6).
  - The instrument prints both. **Pass** = the printed e encloses the predicted value and excludes 0. The real source must
    give e ∋ 0 with width ≤ 1e-10.
- **N-b, radial component.** The θθ residual must enclose 0 for the real source, and exclude 0 with b_T × (1 + 1e-6).
- **N-c, wrong background** [v2 2026-11-18, REFEREE-QB2 sec. 4].
  - The ω-perturbation is now **1e-6**, not 1e-12. v1's value sat under the 1e-10 width and could not refuse (doctrine #20).
  - Rebuild the b-profiles from ω × (1 + 1e-6) as in PQ6b (seats/build-2/RUN.md l.49-50). The test must refuse.
  - An enclosure too wide to refuse is reported as "null not sharp", never as a pass.

## 5. The far side x > 12
- **Ψ2 needs no tail** [v2 2026-11-18: now exact by the local-max definition].
- The global profile and the TW readout need T1:
  - a bound |b_ab| ≤ B_far on [12, ∞), since b_ab → constant (DESIGN-QB1 sec. 2), with the asymptotics m/r → 0.05945
    (LOG.md l.7);
  - the propagator bound ‖Q⁻¹(P0 + 2P1)‖ ≤ L_far there;
  - Grönwall.
- Caveat: at fixed r with t → 0⁻, s → ∞, and the split into s- and e^{2s}-coefficients degenerates (zahn.txt l.636-637).
  Claims for x > 12 are bounds inside ε s e^{2s} ≪ 1.

## 6. Honesty: is the sign decidable from the state-independent part alone?

**No, not in general.**

1. **First order, fixed s, perturbative window.** δm_U = ε e^{2s}(sΨ2 + Ψ1). The s-coefficient is decided by Ψ2 alone,
   under Zahn's state assumption.
2. **At the threshold.**
   - sΨ2 + Ψ1 = Ψ2·ln(L*/(-t)), where L* := t0·exp(Ψ1(t0)/Ψ2) is set by the state together with the renormalised C²
     coupling a (BO eq. 3.10, bo.txt l.220, OPENED).
   - -t* ~ ℓ_P(ln)^{1/2}.
3. **The minimal extra assumption ("IR naturalness")** [v2 2026-11-18, REFEREE-QB2 sec. 3 Fix 1]. It has three parts:
   - (i) a(μ = 1/t0) = O(1);
   - (ii) Zahn's self-similar Hadamard state is normalised at the initial-data scale t0;
   - **(iii) |Ψ1^state/Ψ2| ≪ ln(t0/ℓ_P).** This is needed in its own right: Ψ2 could be small through the core/shell
     cancellation in b_nn (RESULT-QB1 rule 6; failure mode 2), and then even an O(1) state term competes.
   - Whose model: Zahn (zahn.txt l.683-684, "in a first approximation the state-dependent part could be ignored",
     OPENED) plus BO's effective action eq. 3.10-3.11.
4. **Two more limits** [v2 2026-11-18, REFEREE-QB2 sec. 3 Fix 2].
   - **(L1) Linear extrapolation.** At T* the perturbation of (∇r)² is O(1). Using the linear response there is the TW
     extrapolation (tw2.txt l.4100), a model step that the first-order theory does not justify.
   - **(L2) Near-Planck curvature.** At -t* ~ ℓ_P(ln)^{1/2}, the curvature is ~ 1/(ℓ_P² ln), only a log below the Planck
     scale. Higher loops and quantum gravity are suppressed only logarithmically there.
5. **If a is natural at the Planck scale,** then L* ~ ℓ_P, and Ψ2 does not decide the sign.
6. **In TW's 2D model** there is no log term at all (zahn.txt l.682).
7. **CONJECTURE** [v2 2026-11-18, REFEREE-QB2 sec. 3 Fix 3]. Under (3) and (4), sgn Ψ2 is the direction in which ħ moves
   the threshold. This is not derived: the shift depends on f at x*, not at x_V.

## 7. Nodes

| id | claim | source | status | expected | depends | referee | lean |
|---|---|---|---|---|---|---|---|
| D1 | D_new(2) != 0 and D_HKA(2) != 0, so S_2 ∩ R_2 = span(G_2) | runs7/verdict_M.json l.4-6; seats/build-2/pq7b.json l.13-35 | CERTIFIED | -19.9422952680 | - | REFEREE-QB2 1(d) CORRECT | no |
| Q1 | QB1 b-profiles: sign map, sonic b_uu, centre c0 (sign map only; not a source) | seats/build-2/result_qb1b.json; REFEREE-QB1 C3 | CHECK (lead verdict on BAR-QB1b pending; BAR files not opened) | -173.10044439230 | - | REFEREE-QB1 C3 CORRECT WITH FIXES | no |
| L1 | Lemma (a): J' = -A J - A D_t along the HKA pencil (c_D = -A + I F along the derive rows); conserved source keeps J = 0 | PROOF-QB2-LEMMAS-v2 sec. 2; srcprop.log l.4-11; r2_forced_k13.log l.5-12, 24-26 | PROVED (to referee; SYM, two independent codes) | n/a - identity | - | REFEREE-QB2 1(a) CORRECT WITH FIXES (fixed in v2) | no |
| L2 | Lemma (b): forced sonic step, abs(rho(2)-n) >= n+1.9085, kappa_pole = d so Delta2(2) != 0 from d < 0 | PROOF-QB2-LEMMAS-v2 sec. 3; PROOF-B4-v4.md l.81-82; fw_cert.py l.385; seats/referee-2/r2_pole.log | PROVED (to referee) | -1.90989 | L1 | REFEREE-QB2 1(b) CORRECT WITH FIXES (fixed in v2) | no |
| L3 | Lemma (c): centre step for a source with Q^-1 S = O(v^2), sigma_I = o(1) | PROOF-QB2-LEMMAS-v2 sec. 4 | PROVED (to referee) | n/a - structural | L1 | REFEREE-QB2 1(c) CORRECT WITH FIXES (fixed in v2) | no |
| L3s | centre orders of the four rows of Q^-1 S: (v^2, v^2, v^2, v^3), sigma_I = O(v^2), isotropy not needed (given b orders) | seats/design-2/centre_orders.log | PROVED (to referee; SYM) | n/a - orders | - | NEW: referee | no |
| L4 | Lemma (d): h2 exists, unique mod G_2 | PROOF-QB2-LEMMAS-v2 sec. 5 | PROVED (to referee) | n/a - structural | D1, L2, L3 | REFEREE-QB2 1(d) CORRECT | no |
| L5 | secular structure (s e^{2s} solvable; e^{2s} state-dependent, not separately conserved, ambiguous along h2) | PROOF-QB2-LEMMAS sec. 6 (v1, unchanged) | PROVED (to referee) | n/a - structural | L1, L4 | REFEREE-QB2 sec. 2 CORRECT | no |
| O1 | (grad r)^2 = 1/A; Abar' = -(8/3) om N V/(1-V^2); local max of A on U at x_V; invariant under every gauge at x_V | sec. 2.1-2.2; xv.log l.1-3 | PROVED (to referee; SYM for Abar') | -0.2529 | - | REFEREE-QB2 sec. 3 CORRECT | no |
| N0 | certified TMs of b_nn, b_nr, b_rr, b_TT on every transport step and outward cell, traceless check | tmcell.py / schart.py reuse; build seat | GAP | n/a - instrument | Q1 | - | no |
| N1 | source rows: the four HKA rows (K5) and sigma_I from b_ab, re-derived by the build; rows 0-1 vs srcprop.log l.13-14 AND rows 2-3 vs r2_forced_k13.log l.21-22 | r2_forced_k13.log l.19-25; srcprop.log l.13-15 | GAP | n/a - exact map | N0 | - | no |
| N2 | conservation of certified b: D_t and theta-theta residual enclose 0 on [-12, 12] | Lemma (a)(i),(iii) | GAP | [-1e-10, 1e-10] | N1, L1 | - | no |
| N3 | forced sonic series with the four-row source (TM tails) | Lemma (b) | GAP | n/a - instrument | L2, N1 | - | no |
| N4a | centre orders of Q^-1 S certified on the cm manifold: first two theta-coefficients of all four rows vanish | L3s; cm manifold (K 60) | GAP | [-1e-30, 1e-30] | L3s, N0 | - | no |
| N4 | forced centre series in theta, y_c^p = O(v^2), rho1, rho2 | Lemma (c) | GAP | n/a - instrument | L3, N4a | - | no |
| N5 | augmented 5x5 certified transport with all four rows of Q^-1 S | lin_stages_v2 / CachedLinear2 reuse | GAP | n/a - instrument | N3, N4 | - | no |
| N6 | matching with the D_HKA(2) matrix at x_m = -4.6, gauge Nbar_p(0) = 0; e encloses 0 | Lemma (d) practical form | GAP | [-1e-10, 1e-10] | D1, N2, N5 | - | no |
| N7 | V has exactly one zero x_V in U = [-0.30, -0.20], V' > 0 on U | E5 tube enclosure; xv.log l.3 | GAP | [-0.25456, -0.25128] | O1 | - | no |
| C1 | manufactured phi(x) G_2 reproduced (gauge Nbar_p(0) = 0) | sec. 3 C1 | GAP | n/a - exact profile phi G_2 | N5, N6 | - | no |
| C2 | resonance null at kappa*: uniqueness refused | sec. 3 C2 | GAP | n/a - must refuse | N5 | - | no |
| C6 | route B overlaps route A for Psi2 | sec. 1.2 step 6 | GAP | n/a - overlap | N5 | - | no |
| C8 | manufactured conserved Sigma: e and theta-theta residual enclose 0; scaled variant refuses | sec. 3 C8; dt_struct.log | GAP | n/a - must pass / refuse | N1, N5, N6 | - | no |
| NA | content null N-a: e excludes 0 and encloses the printed predicted abs(e) | sec. 4 | GAP | n/a - must refuse | N6 | - | no |
| NC | content null N-c at omega x (1 + 1e-6): refuses | sec. 4 | GAP | n/a - must refuse | N0, N6 | - | no |
| N8 | headline: certified sign and size of Psi2 = Abar_2(x_V)/A(x_V) (times N_c ε s e^{2s}), local-max definition on U | this design | GAP | n/a - sign not predicted | N6, N7, C1, C2, C6, C8, NA, NC, L5 | - | no |
| T1 | far-side tail for the global profile and the TW readout | sec. 5 | GAP (optional) | n/a - bound | N5 | - | no |
| TW | extended TW threshold algebra (log term); threshold direction = CONJECTURE | twlog.log; REFEREE-QB2 sec. 3 | CHECK (SYM done; threshold link conjectural) | n/a - algebra | O1 | REFEREE-QB2 sec. 3 CORRECT WITH FIXES | no |

## 8. Failure modes

1. **The source is not conserved to the instrument's width (N2/NA).** Causes: the frame-to-coordinate map, a missing
   fluid-row source, or a QB1 component error.
   - C8 (every row, known conserved answer) separates a wrong σ → S map from a bad b-profile.
   - C1 separates a solver fault.
2. **Ψ2 encloses 0.** The core/shell cancellation in b_nn. Reported as "no sign", never rounded. This also weakens
   IR-naturalness item (iii).
3. **Conditioning at the matching / width inflation** on the 5×5 transport. Remedies: move x_m, or split steps. C1
   detects it against an exact answer.
4. **[v2 2026-11-18, REFEREE-QB2 BUG: v1's remedy deleted]** v1 said "keep S in rows 0-1 outside Q⁻¹". That was wrong:
   the HKA source has nonzero fluid rows.
   - New failure mode: the build uses the derive-row source with the HKA pencil, or the reverse. This mixes systems that
     agree only on J = 0, and they differ by F·I.
   - Guard: N1 compares all four rows against r2_forced_k13.log l.19-22. C8 exercises every row.
5. **Sonic 0/0 in M_x⁻¹ acting on S_2, S_3.** Near x = 0 use the Frobenius form (Lemma b), not Q⁻¹S. Also check that the
   transported Q⁻¹S stays bounded on the switch shell |x| = 0.05.
6. **Physics-scope failure.** Without (3)-(4) of sec. 6, Ψ2 does not decide the threshold. Every write-up carries sec. 6
   and the CONJECTURE label.
7. **TW readout disagrees with Ψ2.** Logged as a tension.

## 9. Cost driver
- **The four certified b-component TMs on every step** (N0) dominate. They are now QB2's own job, since PQ2b is a sign map.
  - QB1b's sign maps took about 2 h over 2261 cells for b_nn alone.
  - Four components on about 720 tube steps plus about 1500 outward cells: projection 3-6 h, an evening run.
- The 5×5 transport is similar in cost to D_new(2) (14.9 s per κ) times the step overhead.
- The cm θ-series and N4a take minutes.

## 10. Ceiling (no odds)
- The first certified sign of the geometric (log-secular) quantum push on the peak compactness of a true critical-collapse
  solution. It is the same for every free massless field of spin ≤ 1, in any Hadamard state respecting self-similarity.
- [v2 2026-11-18, REFEREE-QB2 sec. 3 Fix 3] That this is also the direction in which ħ moves the black-hole threshold is a
  **CONJECTURE**. It holds only under IR naturalness and the TW linear extrapolation.

**Magic-trick line:** "The critical ball of light is most compact on the one sphere where its light stands still; a proof
reads, on that sphere, which way the quantum push moves it."

## 11. Curiosities (UNTESTED unless stated)
- c_I = -A along the HKA pencil, and the forcing weight is also -A. Do both hold for every k? The referee's k = 1/5 run
  agrees for the forcing (POST HOC).
- **κ_pole = d is now PROVED** (REFEREE-QB2). The remaining open question is d + λ' = 1 (PROOF-B4.md l.60).
- **S_2 and S_3 carry only 1/N, never det M_x** (K5 l.27). The source never introduces a new sonic singularity of its own.
