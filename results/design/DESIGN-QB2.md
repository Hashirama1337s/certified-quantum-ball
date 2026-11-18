# DESIGN-QB2: the forced O(hbar) response of the critical ball of light, and the functional that reads its push

thesis: main quest L2 (real gravity), LIGHTBALL Raid 3 follow-on: certified quantum back-reaction on the certified Evans-Coleman solution (SCOPE.md, approved 2026-11-17)
target: N8
reuses: E5 enclosure (lin_tube2_k1_3.pkl, BAR-LB-E5 HOLD d0af75803f), BG2 outward record, M count instrument (count_b4v3 + CachedLinear2 L13_nsub2, seal b1b0ee1cde), QB1/QB1b Taylor-model cells and b_ab profiles (seats/build-2), PROOF-B4 Facts 1.1-1.4 + Lemma 1.5 + Thm 1.6, PROOF-QB2-LEMMAS (a)-(d)
produces: certified forced-response instrument (augmented 5x5 transport + forced sonic/centre series), h2 profile on [-12, 12], the gauge-invariant compactness functional Psi2, the conservation test for QB1 b-profiles

Design + proof seat design-2, 2026-11-18. Paper, reading and three small exact or float checks (non-counting; all in
`seats/design-2/`: `srcprop.py/.log`, `xv.py/.log`, `twlog.py/.log`). **Nothing on this page is certified by a QB2 run.**
Labels: OPENED / ABSTRACT / SNIPPET for outside sources (with locator), PAPER, SYM, CERTIFIED (on file), MOCK.

## 0. In plain words

- **The question.** Quantum fields on the critical ball of light add a tiny extra stress. Its geometric part grows like
  s·e^{2s} as the collapse runs. How does the ball answer, to first order in ħ?
- **The answer's shape (PROVED in PROOF-QB2-LEMMAS).**
  - The ball's response is e^{2s}(h1 + s·h2).
  - The s·e^{2s} part h2 is fixed by geometry alone. It is unique up to a relabelling of time (the gauge mode G_2).
  - The e^{2s} part h1 needs the quantum state, which nobody has constructed on Evans-Coleman. h1 is also ambiguous by
    exactly a multiple of h2.
- **The number that reads "toward or away from a black hole".**
  - (∇r)² = 1 - 2m/r equals exactly 1/A in the KHA variables (PAPER, sec. 2.1), and A peaks exactly where the fluid
    velocity V crosses zero (PAPER + SYM: Ā' = -(8/3)ωNV/(1-V²)).
  - So the maximum compactness of each slice, 2m/r at x_V ≈ -0.253 (float, xv.log), is a clean scalar. The gauge mode
    does not move it, because Ā_ss'(x_V) = 0 exactly.
  - QB2's headline is the sign of Ψ2 := Ā_2(x_V)/A(x_V), the s·e^{2s} coefficient of the first-order change in the maximum
    compactness. Ψ2 > 0 means the geometric quantum push makes the ball more compact (toward trapping); Ψ2 < 0 means
    less compact (away).
- **Honesty, up front (sec. 6).** Ψ2 alone decides the sign of the *full* O(ħ) push at the threshold only under one extra
  assumption: the state and the renormalised Weyl² coupling are of natural size at the scale of the initial data, not at
  the Planck scale. Without it, the sign is not decidable from the state-independent part.

## 1. What is computed, certified

### 1.1 The forced problem (PROOF-QB2-LEMMAS sec. 1, 6)
- **Equations.** L(2)h2 = S(σ_V) and I_2(h2) = σ_I(σ_V), on x ∈ (-∞, 12], plus a tail on x > 12 (T1, optional).
- **Source rows (SYM, srcprop.log l.13-15; units: per unit ε s e^{2s}, ε = ℓ_P²/t0², times the field multiplier N_c > 0).**
  - Hamiltonian row: S_H = (N_c/30π)·A e^{-2x} b_nn.
  - Slicing row: S_N = -(N_c/60π)·A e^{-2x}(b_nn - b_rr).
  - Momentum row: σ_I is built from b_nn, b_nr, b_rr with the derive-normalised coefficients of srcprop.log l.15,
    divided by the HKA/derive ratio ϑ (lin_sym.py C3).
  - The fluid rows carry no source (srcprop.log l.16-17).
  - The b_ab are the QB1 profiles r⁴B_(ab) in the Eulerian frame (DESIGN-QB1 sec. 1.2), certified on [-12, 12] by QB1b
    (seats/build-2/result_qb1b.json: PQ2b "coverage complete, no gap").
- **Frame to coordinates.** The PAPER steps are in PROOF-QB2-LEMMAS sec. 1.
  - Σ(∂_s+∂_x, ∂_s+∂_x) = A r²Σ_(nn) per N², and Σ_xx = A r² Σ_(rr).
  - The build seat re-derives the full map symbolically from the KHA metric and checks it against srcprop.log l.13-15.
    This check is node N1.

### 1.2 Pipeline (smarter-first: reuse the M instrument unchanged)
1. **Augmented homogeneous system.** Ŷ = (Y, z), z' = 0, with matrix [[Q⁻¹(P0+2P1), Q⁻¹S],[0,0]].
   - The certified linear transport of the M instrument (lin_stages_v2 / CachedLinear2) then runs on a 5×5 system, with no
     new transport code. The only new input is the source column as a Taylor model on each step.
   - Build-2's tmcell.py already yields b_ab as degree-P Taylor models in ε = x - x_m per cell (RUN.md l.37-45).
2. **Sonic leg.** Forced Frobenius series at x = 0 (Lemma b).
   - The constants are certified and exact: |ρ(2) - n| ≥ n + 1.9085, Δ2(2) ≠ 0, ℓM_s r ≠ 0.
   - Run as a TM recursion with the M instrument's sonic_tm and the source series added to R_n, giving y_s^p and the
     homogeneous s_2, s_2^HKA.
   - Transport from x = 0 inward to x_m = -4.6 and outward to 12.
3. **Centre leg.** Forced regular series in θ on the certified centre manifold (Lemma c): y_c^p = O(v²) plus ρ1, ρ2.
   Transport outward to x_m.
4. **Matching.** Solve y_c^p - y_s^p = a s_2 + c ρ1 + d ρ2 + e e_V at x_m by Cramer's rule with the D(2) matrix.
   - |D_new(2)| ≥ 2.4126 (pq7b.json l.35) gives the conditioning.
   - **e must enclose 0.** This is the conservation test (N6).
   - Gauge: h2 is taken with Ā_p(0) = Ā^p_s(0), i.e. no s_2^HKA. Any other choice changes h2 by c·G_2, which Ψ2 cannot see.
5. **Readout.**
   - Certified enclosure of x_V (the unique zero of V on ℝ) and of V'(x_V) > 0.
   - Then Ψ2 = Ā_2(x_V)/A(x_V) on the sonic leg (x_V ∈ (-0.26, -0.25) lies between x_m and 0).
   - The x_V enclosure width enters only through |dĀ_2/dx|·width (TM).
6. **Route B (independent).** Variation of constants with the homogeneous fundamental system at κ = 2 from the M instrument.
   - Ψ2 = Σ_k ∫ w_k(x) S(x) dx, with Wronskian-weighted kernels built from {s_2, s_2^HKA, ρ1, ρ2} and the adjoint left
     null vector (sec. 3, C6).
   - It must overlap route A.

### 1.3 Interface to QB3
- The global profile h2 on [-12, 12] supports the Tomasevic-Wu-type threshold readout (sec. 2.3).
- The unsolvable h1 is carried symbolically as h1 = h1^{state} + μ·h2 + c·G_2, with μ the renormalisation-scale freedom
  (PROOF-QB2-LEMMAS sec. 6 item 3).

## 2. The observable: which functional decides "toward or away from black-hole formation"

### 2.1 The trapping function is 1/A (PAPER, exact)
- On ds² = r²{A[(1-N²)ds² - 2ds dx + dx²] + dΩ²} (engine/derive.py l.31-37) with r = e^{x-s}:
  g^{ss} = g^{sx} = -1/(r²AN²), g^{xx} = -(1-N²)/(r²AN²), ∂_s r = -r, ∂_x r = r.
- Hence **(∇r)² = 1/A**, i.e. A = 1/(1 - 2m/r), so A is KHA's a².
- The apparent horizon is (∇r)² = 0, the function used by both Brady-Ottewill (eq. 3.14, bo.txt l.235-241, OPENED) and
  Tomasevic-Wu (eq. 6.8, tw2.txt l.4096, OPENED).
- To first order, δ(∇r)² = -Ā_p/A and δ(2m/r) = Ā_p/A.

### 2.2 The headline functional Ψ2 (PAPER + SYM)
- **Where A peaks.** On shell, Ā' = 1 - A + 2ω(1+V²/3)/(1-V²) = **-(8/3)ωNV/(1-V²)** (SYM, xv.log l.1).
  - ω, N > 0, so A is stationary exactly at the zeros of V.
  - V has exactly one zero on x ≤ 0 (E5 PE2, LOG.md l.40 "1 V zero"), and V > 0 on the outward leg (LOG.md l.17).
  - So A increases up to x_V and decreases after it: **x_V is the unique global maximum of A**, with A''(x_V) < 0 iff V'(x_V) > 0.
  - Float: x_V ∈ (-0.25456, -0.25128), A(x_V) ≈ 1.91434 (xv.log l.3), matching "sup A = 1.914" (PROOF-B4.md l.210).
- **Gauge invariance.** G_2 = (Ā_ss', N̄_ss' + 2, ω̄_ss', V_ss') (PROOF-B4.md l.85), and Ā_ss'(x_V) = 0 exactly.
  - So Ā_2(x_V) is the same for every representative h2 + cG_2.
  - It is also independent of the slicing-origin freedom f(x) (PROOF-QB2-LEMMAS sec. 6).
- **Envelope.** max_x(2m/r) is a scalar on each slice s, and to first order its change is δ(2m/r)(x_V) (Danskin, unique
  non-degenerate maximiser). Therefore
  **δ max_x(2m/r) = ε e^{2s}[ s·Ψ2 + Ψ1 ] + O(ε²),  Ψ2 := Ā_2(x_V)/A(x_V),  Ψ1 := Ā_1(x_V)/A(x_V).**
  - **Ψ2 is state-independent; Ψ1 is not.**
  - **Reading:** Ψ2 > 0 means the log-secular quantum push raises the peak compactness of the critical ball (toward
    trapping); Ψ2 < 0 means it lowers it (away).
- **Field independence.** Ψ2 ∝ N_c > 0 times the response to B/(240π²). So its sign is the same for every free massless
  field of spin ≤ 1 (scalar via Zahn's theorem; spin 1/2 and Maxwell via BO's conformal transformation law, RESULT-QB1
  rule 2). This holds "in any Hadamard state respecting self-similarity (Zahn eq. 4); existence on Evans-Coleman not
  proved" (RESULT-QB1 rule 3).
- **Source sentence (why this is the open physics).** Brady-Ottewill end their argument on an assumption: it "relies
  heavily on the assumption that quantum matter tends to oppose black-hole formation" (bo.txt l.277-278, OPENED). Their
  abstract: "When β < 0.5 our argument is inconclusive" (bo.txt l.4); for EC, β = 0.3558 (verdict_M.json l.5).
  Ψ2 replaces that assumption, for the state-independent part, with a computed sign.

### 2.3 Tomasevic-Wu route (QB3 readout, secondary; a model step)
- **TW's rule.** "sgn(Δp) = sgn(f_q/f_c)" at the earliest apparent horizon (tw2.txt l.4208, OPENED; the general-rule
  remark is footnote 45, l.4323).
- **TW's setting.** (∇r)² ≈ f0 + (p - p*)e^{κ_c T}f_c + ħ e^{qT}f_q (eq. 6.8, l.4096), with **no** T e^{qT} term. Their 2D
  anomaly model has no geometric log (zahn.txt l.682, OPENED).
- **Extended here (SYM, twlog.log).** With f_q → f1 + T·f2 (our f2 = -Ā_2/A, f1 = -Ā_1/A), the extremality condition gives
  Δp = ε e^{(q-κ_c)T*}[q(f1 + T*f2) + f2]/(κ_c f_c).
  - With f2 = 0 this reproduces TW (6.12) and (6.15) exactly (twlog.log T1).
  - Leading log: sgn Δp = sgn(f2/f_c).
  - Since κ_c = 2.81 > q = 2, the AH branch at large T* exists only where (1 - q/κ_c)T f2 < 0, i.e. where f2 < 0 (twlog.log T3).
- **Gauge caveat (PAPER).**
  - TW's pointwise test "f2(x) < 0 somewhere" is **not** gauge-invariant in our setting: f2 → f2 - cĀ_ss'/A, and
    Ā_ss' takes both signs. So for any profile, some c makes f2 negative somewhere.
  - TW's extrapolation lets (∇r)² perturbations reach O(1) (tw2.txt l.4100), which in areal gauge is an O(1) metric change.
  - So the TW readout is reported only (i) in the stated gauge (HKA, N̄_p(0) = 0) and (ii) next to Ψ2. If the two disagree
    in sign, that is a tension to log, not a result.
  - Ψ2, at the one point where the gauge term vanishes, is the certified headline.

## 3. Controls

| # | control | expected | what it tests |
|---|---|---|---|
| C1 | manufactured solution: source S_m = Q·G_2, σ_m = 0 (exact: L(2)(x·G_2) = Q G_2, I(xG_2) = 0) | solver returns x·G_2(x) (HKA gauge, N̄_p(0) = 0 holds since x = 0 there) at x = -1, -2, -4.6, 0.5, 3 | forced sonic series, centre series, transport and matching on an exactly known answer |
| C2 | resonance null: same pipeline at κ ∈ κ* = [2.810552549 ± 3.83e-10] (verdict_M.json l.5) | matching determinant encloses 0; uniqueness **refused** | that Lemma (d)'s D(2) ≠ 0 step has teeth |
| C3 | Brady-Ottewill structure (bo.txt l.163-170, l.229-230, OPENED) | (i) rate κ_q = 2 from the u-scaling; (ii) the s e^{2s} term is forced (no s² term since D(2) ≠ 0); (iii) BO case (ii) applies (κ_c = 2.81 > 2, bo.txt l.267) | structure; PAPER (PROOF-QB2-LEMMAS sec. 6) |
| C4 | Tomasevic-Wu closed-form Roberts | TW's 2D Roberts response F_q = -r_q (eq. 5.32, tw2.txt l.3160-3168, OPENED; garbled extraction, PDF check before use) and footnote 45 (f_q > 0 at f_c > 0 → lowered) are reproduced by the extended algebra with f2 = 0 | **algebra only.** 2D Polyakov s-wave, minimally coupled scalar, Roberts background, TW's ζ = -3/4 state, no log term: **not a control of Ψ2.** 2D vs 4D: "dimensional reduction and quantization do not commute in general" (tw.txt l.409, OPENED) |
| C5 | zero floor: k = 1/3 flat Friedmann in the KHA chart (Bach ≡ 0, DESIGN-QB1 C5) | S ≡ 0, so h2 ≡ 0 mod G_2 | the source module |
| C6 | route B (variation of constants with the adjoint) vs route A (direct forced transport + matching) | overlap of Ψ2 enclosures | double implementation of the headline |
| C7 | linearity | Ψ2[2S] = 2Ψ2[S]; Ψ2[-S] = -Ψ2[S] at the point level | readout sign path |

No published closed form exists for the **V-forced** (4D, log-secular) response on any CSS background in the sources
opened: Zahn gives V on Roberts (zahn.txt l.612-621) but no response; TW's responses (Garfinkle, Roberts) have no log; BO
give structure only. C1 is the exact-answer control in its place.

## 4. Content nulls

- **N-a, non-conserved source (primary).** Multiply only the Hamiltonian-row source S_H by (1 + 10⁻⁶), keeping S_N and σ_I.
  - By PROOF-QB2-LEMMAS Lemma (a)(i), Comp ≠ 0. The matching coefficient e (or J transported from the sonic side to x_m)
    must **exclude 0**.
  - The real source must give e ∋ 0 with width ≤ 10⁻¹⁰.
  - This tests that the instrument reads conservation, and that QB1's certified b-profiles satisfy (∇·V)_t = 0 on EC.
    That is a non-trivial check on QB1.
- **N-b, radial component.** The θθ Einstein residual (Lemma (a)(iii)) must enclose 0 for the real source, and exclude 0
  when b_rr is scaled by (1 + 10⁻⁶).
- **N-c, wrong background.** Feed the PQ6b off-shell b-profiles (omega × (1 + 10⁻¹²), seats/build-2/RUN.md l.49-50).
  - They are not the Bach tensor of a solution, so (∇·B)_t ≠ 0 and N-a's test must refuse. An enclosure that is too wide
    to refuse is reported as "null not sharp", not as a pass.

## 5. The far side x > 12

- **Ψ2 needs no tail.** x_V ≈ -0.25. The BVP lives on x ≤ 0 (sonic point plus centre), and x_V lies inside it.
- **The global profile and the TW readout need two bounds (node T1, optional for the headline).**
  1. **Source tail.** |b_ab(x)| ≤ B_far on [12, ∞). B is finite at fixed r > 0 as t → 0⁻ (smooth there) and has
     dimension length⁻⁴, so b_ab = r⁴B_ab → constant, and e^{-2x}b_ab → 0 (DESIGN-QB1 sec. 2 table). It is certified from
     the far-side asymptotics (N → 0, ω → ω∞, A → 1/(1 - 2·0.05945), LOG.md l.7) by a tail lemma of the S1b kind.
  2. **Propagator tail.** ‖Q⁻¹(P0 + 2P1)‖ ≤ L_far on [12, ∞), with no sonic point there (D_sonic < 0 continues; to be
     shown). Grönwall then gives ‖h2(x)‖ ≤ (‖h2(12)‖ + B_far'/(L_far + 2))e^{L_far(x-12)}.
- **Physical caveat (PAPER).** At fixed r as t → 0⁻, s → ∞, so ε s e^{2s} diverges. The split into s-coefficient and
  e^{2s}-coefficient degenerates there: Zahn needs r e^{2(u-r)}-type terms on Roberts for the same reason (zahn.txt
  l.636-637, OPENED). Claims for x > 12 are therefore bounds inside the window ε s e^{2s} ≪ 1, never values at t = 0.

## 6. Honesty: is the sign decidable from the state-independent part alone?

**No, not in general.** Plainly:

1. **At first order and fixed s, in the regime where perturbation theory holds,** the change in peak compactness is
   ε e^{2s}(sΨ2 + Ψ1). The s-coefficient is decided by Ψ2 alone (PROVED structure + QB2 computation). That is a clean,
   gauge-invariant, field-independent statement, conditional on Zahn's state assumption.
2. **What decides the push at the threshold** is the sign of sΨ2 + Ψ1 at the trapping time T*.
   - Ψ1 depends on the state, the renormalised Weyl² coupling a (BO eq. 3.10, bo.txt l.220, OPENED), and the origin of s.
     All three shift Ψ1 by multiples of Ψ2 (PROOF-QB2-LEMMAS sec. 6 item 3).
   - The invariant form is sΨ2 + Ψ1 = Ψ2·ln(L*/(-t)). L* is a length fixed by the state together with a.
   - At the threshold, -t* ~ ℓ_P·(ln)^{1/2} (from ε T e^{2T}|f2| ~ f0, TW-style).
   - **So the log-secular term decides the sign iff L* ≫ ℓ_P(ln)^{1/2}.**
3. **The minimal extra assumption that decides it.** "IR naturalness":
   - the state-dependent T(y) and the renormalised C² coupling a are O(1) when normalised at the initial-data scale t0, so
     L* ~ t0 up to e^{O(1)};
   - then T* ≈ ln(t0/ℓ_P) and sgn(push) = sgn(Ψ2), with corrections O(1/ln(t0/ℓ_P));
   - **whose model:** Zahn's Hadamard state respecting self-similarity, normalised at the initial-data scale. This is Zahn's
     own suggested approximation: "in a first approximation the state-dependent part could be ignored" (zahn.txt l.683-684,
     OPENED). Add BO's effective action (eq. 3.10-3.11) with a(μ = 1/t0) not anomalously large.
4. **If instead a is natural at the Planck scale** (a(1/ℓ_P) = O(1)), then L* ~ ℓ_P, the log at the threshold is O(1), and
   QB2's state-independent part does **not** decide the sign. The state would be needed.
5. **In TW's 2D model** there is no log term at all (zahn.txt l.682), so the whole push is state-dependent there. Their
   result and Ψ2 answer different questions.

## 7. Nodes

| id | claim | source | status | expected | depends | referee | lean |
|---|---|---|---|---|---|---|---|
| D1 | D_new(2) != 0, so S_2 ∩ R_2 = span(G_2) | runs7/verdict_M.json l.4-6; seats/build-2/pq7b.json l.13-35 | CERTIFIED | -2.41262708580 | - | REFEREE-QB1 C4 CORRECT WITH FIXES | no |
| Q1 | QB1 source b_ab certified on [-12,12], centre c0, sonic b_uu | seats/build-2/result_qb1b.json (PQ2b, PQ8); REFEREE-QB1 C3 | CHECK (lead verdict on BAR-QB1b pending; this seat did not open BAR files) | -173.10044439230 | - | REFEREE-QB1 C3 CORRECT WITH FIXES | no |
| L1 | Lemma (a): J' = -A J + ϑ·Comp, Comp = -A^2 (div Sigma)_t; conserved source keeps J = 0 | PROOF-QB2-LEMMAS sec. 2; seats/design-2/srcprop.log l.4-11 | PROVED (to referee; SYM P1-P3) | n/a - identity | - | NEW: referee | no |
| L2 | Lemma (b): forced sonic Frobenius at kappa=2, abs(rho(2)-n) >= n+1.9085, Delta2(2) != 0, affine 2-plane | PROOF-QB2-LEMMAS sec. 3; PROOF-B4-v4.md l.81-82; fw_cert.py l.385; verdict_M.json l.6 | PROVED (to referee) | -1.90989 | L1 | NEW: referee | no |
| L3 | Lemma (c): centre step, O(v^2) source, affine 2-plane, J = 0 | PROOF-QB2-LEMMAS sec. 4; PROOF-B4.md l.70-81 | PROVED (to referee) | n/a - structural | L1 | NEW: referee | no |
| L4 | Lemma (d): forced h2 exists, unique mod G_2 (2+2-1 = 3) | PROOF-QB2-LEMMAS sec. 5 | PROVED (to referee) | n/a - structural | D1, L2, L3 | NEW: referee | no |
| L5 | secular structure: s e^{2s} state-independent and solvable; e^{2s} state-dependent, not separately conserved, ambiguous along h2 | PROOF-QB2-LEMMAS sec. 6; zahn.txt l.136, l.152, l.348 (OPENED) | PROVED (to referee) | n/a - structural | L1, L4 | NEW: referee | no |
| O1 | (grad r)^2 = 1/A; Abar' = -(8/3) om N V/(1-V^2); x_V unique max of A; Abar_ss'(x_V) = 0 | sec. 2.1-2.2; seats/design-2/xv.log l.1-3 | PROVED (to referee; SYM for the Abar' identity) | -0.2529 | - | NEW: referee | no |
| N1 | source rows S_H, S_N, sigma_I from b_ab (frame to KHA), re-derived and matched to srcprop.log l.13-15 | srcprop.log l.13-15; build seat script | GAP | n/a - exact map | Q1 | - | no |
| N2 | conservation of certified b: Comp(sigma_V) and theta-theta residual enclose 0 on [-12,12] | Lemma (a)(i),(iii); build seat | GAP | [-1e-10, 1e-10] | N1, L1 | - | no |
| N3 | forced sonic series, TM tails, y_s^p and s_2 at kappa=2 | Lemma (b); count_b4v3 sonic_tm + source | GAP | n/a - instrument | L2, N1 | - | no |
| N4 | forced centre series in theta, y_c^p = O(v^2), rho1, rho2 | Lemma (c); cm manifold | GAP | n/a - instrument | L3, N1 | - | no |
| N5 | augmented 5x5 certified transport on tube [-4.6,-0.05] and outward [0.05,12] | lin_stages_v2 / CachedLinear2 reuse | GAP | n/a - instrument | N3, N4 | - | no |
| N6 | matching via D(2) matrix at x_m = -4.6; constraint coefficient e encloses 0 | Lemma (d) practical form | GAP | [-1e-10, 1e-10] | D1, N2, N5 | - | no |
| N7 | certified enclosure of x_V, V'(x_V) > 0, V > 0 on (x_V, 12] | E5 tube + BG2 outward record; xv.log l.3 | GAP | [-0.25456, -0.25128] | O1 | - | no |
| C1 | manufactured solution x·G_2 reproduced | sec. 3 C1 | GAP | n/a - exact profile x·G_2 | N5, N6 | - | no |
| C2 | resonance null at kappa*: uniqueness refused | sec. 3 C2; verdict_M.json l.5 | GAP | n/a - must refuse | N5 | - | no |
| C6 | route B (variation of constants, adjoint) overlaps route A for Psi2 | sec. 1.2 step 6 | GAP | n/a - overlap | N5 | - | no |
| NA | content null N-a: non-conserved source refused (e excludes 0) | sec. 4 | GAP | n/a - must refuse | N6 | - | no |
| N8 | headline: certified sign and size of Psi2 = Abar_2(x_V)/A(x_V) (times N_c ε s e^{2s}) | this design | GAP | n/a - sign not predicted | N6, N7, C1, C2, C6, NA, L5 | - | no |
| T1 | far-side tail: abs(b_ab) <= B_far and propagator bound on [12, inf) for the global profile and TW readout | sec. 5 | GAP (optional) | n/a - bound | N5 | - | no |
| TW | extended TW threshold algebra (log term), reduces to TW 6.12/6.15 | seats/design-2/twlog.log T1-T3 | CHECK (SYM done; PDF check of tw2 eq. 5.32 before C4 numbers) | n/a - algebra | O1 | - | no |

## 8. Failure modes

1. **The source is not conserved to the instrument's width (N2/NA).** The QB1 b-profiles might fail (∇·B)_t = 0 beyond the
   width, from a convention slip in the frame-to-coordinate map or a QB1 component error.
   - Then e excludes 0 on the real source and the forced problem refuses.
   - C1 separates a solver fault from a source fault; N-a shows the test can refuse.
2. **Ψ2 encloses 0.** h2 is fine but Ā_2(x_V) is too small or too wide to sign. Causes: cancellation between the negative
   core and positive shell of b_nn (RESULT-QB1 rule 6), and the x_V enclosure.
   - Remedy: sharper TMs near x_V, route B as a second enclosure. Reported as "no sign", never rounded to a sign.
3. **Conditioning at the matching.** |D_new(2)| ≥ 2.41 is healthy, but the centre forced series and the 5×5 transport may
   inflate widths (the M instrument was tuned for 4×4 at complex κ).
   - Remedy: move x_m (the zero set is x_m-independent wherever c_V(x_m) ≠ 0; REFEREE-QB1 C4), or split the transport.
     C1 detects width blow-up against an exact answer.
4. **Sonic-point removable singularities.** The source column Q⁻¹S has a det M_x denominator in the fluid block. It is zero
   there, but S enters only rows 0-1, and Q is the identity on them. Any build that passes S through Q⁻¹ symbolically meets
   0/0. Remedy: keep the source outside Q⁻¹ (rows 0-1 only). C1 exercises the sonic leg.
5. **Physics-scope failure (not numerical).** Even a certified Ψ2 is the s-secular coefficient, under Zahn's state
   assumption. Without IR naturalness (sec. 6) it does not decide the threshold push. Every write-up must carry sec. 6.
6. **TW readout disagrees with Ψ2.** This is gauge or extrapolation dependence (sec. 2.3). Log it as a tension, and keep Ψ2
   as the headline.

## 9. Cost driver

- **The source column as Taylor models on every transport step,** not the transport itself.
  - D_new(2) on the M instrument took 14.9 s (pq7b.json).
  - The QB1b b-maps took about 2 h in total over 2261 cells (result_qb1b.json, sum of time_s).
  - QB2 needs b_ab jets along the E5 tube steps (651 x-steps plus 69 τ-steps) and about 1500 outward cells. Projection:
    1-3 h for one pass, which fits an evening run.
- The forced centre series needs b_ab as θ-series from the cm coefficients (K 60) once. That costs minutes.
- Route B costs about the same as route A (one more transport of the adjoint).

## 10. Ceiling (no odds)

- The first certified sign of the geometric (log-secular) quantum push on a true critical-collapse solution. It is the same
  for every free massless field of spin ≤ 1, in any Hadamard state respecting self-similarity.
- Under IR naturalness, this is also the direction in which Planck's constant moves the black-hole threshold for radiation.

**Magic-trick line:** "The critical ball of light is most compact on the one sphere where its light stands still; a proof
reads, on that sphere, which way the quantum push moves it."

## 11. Curiosities to surface (UNTESTED, labelled)

- **c_I = -A exactly**, and the forcing enters the constraint as -A²(∇·Σ)_t (SYM). Is the constraint-propagation rate equal
  to minus the metric function for every k, i.e. a general KHA identity? One sympy run at symbolic k would tell.
- **A peaks exactly where V = 0.** The KHA identity Ā' = -(8/3)ωNV/(1-V²) (k = 1/3; the general-k factor is 2(1+k)) means
  peak compactness sits on the sphere where the fluid is at rest relative to the Eulerian observers. Compare the Weyl-zero
  locus A - 1 = 2ω/3 (DESIGN-QB1 sec. 11). At x_V the peak condition reads A = 1 + 2ω, which differs.
- **The certified Δ2 root equals the float d to 13 digits** (PROOF-QB2-LEMMAS sec. 3). Combined with "d + λ' = 1"
  (PROOF-B4.md l.60), are the sonic constants exactly κ_pole = d = 1 - λ'? Then ρ(κ) = (1 - λ' - κ)/λ' in closed form.
- **Gauge sensitivity of the TW pointwise rule** in areal gauge (sec. 2.3): their numerical results used other gauges, so a
  vein, not a criticism, until checked against tw2 sec. 6 in full.
