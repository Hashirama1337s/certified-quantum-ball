# PROOF-QB2-LEMMAS-v2: the forced O(hbar) linear response at rate 2 is solvable and unique modulo gauge

Proof seat design-2, 2026-11-18. v2 of PROOF-QB2-LEMMAS.md; v1 is left unedited. Every change from v1 is marked
**[v2 2026-11-18, REFEREE-QB2 <item>]**, citing seats/referee-2/REFEREE-QB2.md. These are the four lemmas that
RESULT-QB1.md rule 7 (l.16) and REFEREE-QB1 C4 fix 2 (seats/referee-1/REFEREE-QB1.md l.69-73) say QB2 owes. The
u e^{2u} secular structure is in sec. 6.

Scripts and logs are in `seats/design-2/`; all exact sympy, all non-counting, all run under memcap 8:
- v1: srcprop.py/.log, xv.py/.log, twlog.py/.log;
- v2: centre_orders.py/.log, dt_struct.py/.log.

The referee's independent derivation is seats/referee-2/r2_forced.py with r2_forced_k13.log ("K1-K5"), and r2_pole.py/.log.

**Labels.**
- **PROVED**: complete on this page; still needs a referee.
- **SYM**: exact sympy identity, log named.
- **CERTIFIED**: an interval enclosure already on file (file:line).
- **CITED (OPENED / ABSTRACT / SNIPPET)**: an outside source, with locator.
- **REDUCED**: a finite list of ball evaluations, not yet run.
- **CONJECTURE**: not derived.

## 0. In plain words

- **The response.** The quantum source grows like s·e^{2s}, and the ball answers with e^{2s}(h1(x) + s·h2(x)). Only h2 is
  fixed by geometry alone.
- **What is proved.** The h2 problem has exactly one solution modulo the gauge solution G_2:
  - (a) a conserved source keeps the constraint solved;
  - (b) the sonic point passes a 2-parameter family;
  - (c) so does the regular centre;
  - (d) the two families meet modulo G_2, because 2 is not a growth rate.
- **What is not solvable.** h1 is state-dependent, and its state-independent part is not conserved on its own.
- **[v2 2026-11-18, REFEREE-QB2 BUG] The source reaches all four HKA rows.** In the pencil the M instrument transports
  (HKA, Q⁻¹(P0 + κP1)), the source enters the fluid rows 2 and 3 as well, not only rows 0-1. v1 sec. 1 said otherwise,
  and that statement was true only for the derive rows.
- **Exact by-products.**
  - The rate is c_I = -A along the HKA pencil.
  - The forcing enters the constraint as -A·(∇·Σ)_t [v2 2026-11-18, REFEREE-QB2 1(a) Fix 1: v1 said -A²].
  - κ_pole = d exactly [v2 2026-11-18, REFEREE-QB2 1(b), POST HOC].

## 1. Setting

- **Background and unknowns.** HKA bar variables Y = (Ā_p, N̄_p, ω̄_p, V_p)(x). Pencil Q Y' = (P0 + κP1)Y, Q = diag(1,1,M_x),
  identity row I(Y) = c(x,κ)·Y (PROOF-B4.md l.24-31; engine/cert/lin_sym.py l.43-57).
- **Rows used (engine/derive.py l.101, l.147).**
  - Hamiltonian, slicing and momentum are solved from E_ss, E_sx, E_xx.
  - Two fluid rows come from div T = 0.
  - E_θθ is unused; it follows from the Bianchi identity.
- **Forced problem.** Σ is added to 8πT in Einstein's equations; the fluid is unchanged.
  - Σ = 8π⟨T⟩_ren, with state-independent part Σ_V = (N_c ħ/(30π)) s B[g] (B = 60V on R = 0; RESULT-QB1 rule 1).
  - N_c > 0 is the field multiplier.
- **κ form.** Σ = ε e^{κs}σ(x) and h = ε e^{κs}Y(x).
- **Frame to coordinates (PAPER; signs checked against REFEREE-QB2 1(c) Fix 2).** Use ∂_s = -t∂_t - r∂_r, ∂_x = r∂_r,
  ∂_t = αn, ∂_r = a e_r, α = N a e^x, tα = -N a r. Write β_ab := e^{-2x} b_ab (b_ab = r⁴B_(ab)), and drop the common factor
  K = N_c/(30π). Then
  σ_ss = A(N²β_nn - 2Nβ_nr + β_rr),  σ_sx = A(Nβ_nr - β_rr),  σ_xx = Aβ_rr,
  and σ_TT = β_T, since Σ_θθ = r²Σ_(θθ) and b_T = r⁴B_(θθ). σ_TT appears only in (∇·Σ)_x.
- **[v2 2026-11-18, REFEREE-QB2 BUG] Two equivalent forced systems.**
  - **(D) The derive rows.** Four evolution rows with source C·σ (srcprop.log l.13-17). Here the fluid rows have no
    explicit σ, because A_x and A_s carry it implicitly.
  - **(H) The HKA pencil.** M_4 - M_HKA = F·I, with F ≠ 0 only in rows 2 and 3 (r2_forced_k13.log l.17-18).
  - On J = 0, (D) equals (H) with source Q⁻¹S = s4·σ + F·σ_I. The HKA source rows S (columns σ_ss, σ_sx, σ_xx), per
    r2_forced_k13.log l.19-22:
    - row 0: (1/N², 2/N², 1/N²); this equals srcprop.log l.13 divided by A, giving A r²Σ_(nn) = Aβ_nn per K;
    - row 1: (-1/(2N²), -1/N², (N²-1)/(2N²)); this equals srcprop.log l.14, giving -(1/2)A(β_nn - β_rr);
    - **row 2: (V/N, (N+6V)/(3N), -(2N²V - N - 3V)/(3N))**;
    - **row 3: (1/(3N), -(NV-2)/(3N), -(2N²+NV-1)/(3N))**.
  - Momentum: σ_I = (-1/N², (N²-2)/N², (N²-1)/N²)·σ = **A(-β_nn + Nβ_nr)**. The referee printed the coefficients at
    K5 l.25; the frame form is my reduction, PAPER.
  - **The build transports (H), so it must use all four rows.**
- κ = 2 unless stated. W_κ = {I ≡ 0}, of dimension 3 (PROOF-B4.md l.47).

## 2. Lemma (a): constraint propagation with a conserved source

**Lemma (a)** [v2 2026-11-18, REFEREE-QB2 1(a) Fix 1 and Fix 2: coefficient and rate restated per system].
Let J := I(Y) - σ_I, in HKA normalisation (κ-coefficient 1).
- (i-H) Along the HKA pencil with its full four-row source: **J' = -A·J - A·D_t(σ)**.
- (i-D) Along the derive evolution rows: **J' = c_D·J - A·D_t(σ)**, with c_D = (8NVω + 4V²ω - 3V² + 4ω + 3)/(3(V²-1)) =
  -A + I·F (r2_forced_k13.log l.6, l.26).
  - Here D_t := e^{2(x-s)}e^{-κs}[(∇·Σ)_s + (∇·Σ)_x], the component of div Σ along ∂_s + ∂_x = -t∂_t|_r.
  - The forcing term is the same on both routes (r2_forced_k13.log l.24).
- (ii) Both rates are regular on x ≤ 0 and tend to -1 at the centre. Suppose D_t = 0. Then J(x0) = 0 at one point gives
  J ≡ 0, and the forced solutions with J ≡ 0 form an affine 3-space with direction W_κ.
- (iii) The unused θθ equation holds iff (∇·Σ)_x = 0 as well.
- (iv) For the QB2 s-coefficient problem both components vanish, because ∇·V[g] = 0.

**Proof.**
1. **Homogeneous rates.**
   - HKA pencil: c_I = (8NVω + 2V²ω - 3V² + 6ω + 3)/(3(V²-1)) (PROOF-B4.md l.41) = -Acons (xv.log l.2). Confirmed
     independently: r2_forced K1, l.5.
   - Derive rows: c_D ≠ -A (K2, l.6). The difference is I·F, which vanishes on I = 0 (K5, l.26).
   - My v1 srcprop P1 (srcprop.log l.4) is the derive-variable version, with its own normalisation.
2. **Forcing coefficient.**
   - SYM P2 (srcprop.log l.7-10) gives Comp = l·D_t with l = -A², in the normalisation of the derive momentum row.
   - That row is Acons × the HKA identity row (the lin_sym C3 ratio, gen_lin_fields.log l.8, per REFEREE-QB2 1(a)).
   - In HKA normalisation, therefore, the coefficient is -A²/A = **-A**. The referee confirms this independently:
     "lam == -A : True; lam == -A^2 : False" (r2_forced_k13.log l.11), also at k = 1/5 (POST HOC).
   - **v1's "ϑ·Comp, Comp = -A²·D_t" double-counted the ratio. It is withdrawn.**
3. **σ_θθ.** It does not enter D_t but does enter D_x (r2_forced_k13.log l.12; dt_struct.log: D_x carries -2σ_TT
   algebraically, and D_t carries no σ_TT). Consistent: σ_θθ enters none of the three rows used.
4. **(ii)** is a linear ODE with a regular coefficient on x ≤ 0, the sonic point included.
5. **(iii)** For a tensor with only a θθ component, the divergence has (s, x) components (2, -2)·E_θθ e^{2s-2x}, summing to 0
   (r2_forced K4, per REFEREE-QB2 1(a)). So with the (s,x) block of E zero, Bianchi gives D_t = 0 and D_x ∝ E^θ_θ.
   - This upgrades v1 step 4 from PAPER to SYM, via the referee's K4.
6. **(iv)** The identity is polynomial in κ with s-independent coefficients. So the s·e^{2s} coefficient obeys the κ = 2
   identity with σ2 alone.
   - ∇·V[g] = 0 identically: V "is automatically conserved" (zahn.txt l.348, OPENED; re-opened by the referee).
   - The remainder (∇s)·V[g] sits at the e^{2s} level (same line). ∎

**Consequence (content null).** A source with D_t ≢ 0 gives J' = -AJ - A·D_t. The sonic solution, which has J(0) = 0,
reaches the centre with J ≠ 0, growing like e^{-x} there, so the matching refuses.
- [v2 2026-11-18, REFEREE-QB2 4] **Expected size for the DESIGN null N-a.** Multiply HKA row 0 of the source by (1 + δ).
  - The change in I·Q⁻¹S is δ·(κ - A)·S_0 = δ(2 - A)·Aβ_nn·K.
  - Hence J_s(x_m) = -δ ∫_{x_m}^0 exp(∫_{x_m}^y A) (2 - A(y)) S_0(y) dy, while J_c ≡ 0.
  - So |e| = |J_s(x_m)|/|c_V(x_m)|, with c_V(x_m) = [0.19288155165634 ± 1.41e-15] (verdict_M.json l.6).
  - The instrument prints this predicted |e| from the same TMs, and passes only if the printed e encloses it and
    excludes 0.

## 3. Lemma (b): the forced sonic Frobenius step at κ = 2

**Lemma (b).** Suppose S and σ_I are analytic at x = 0 (all four rows of S; no assumption that ℓ̂S(0) = 0)
[v2 2026-11-18, REFEREE-QB2 1(b) last line]. Then the forced solutions that are analytic at 0 and have J ≡ 0 form an
affine 2-plane S_2^F = y_s^p + S_2, parametrised by (Ā_p, N̄_p)(0).

**Constants** (certified or exact; recomputed by the referee: |b_hi| ≥ 0.0984552 and ρ(2) ≤ -1.9085039).
1. λ' ∈ [1.0986, 1.0995] (CERTIFIED, PROOF-B4-v4.md l.81; asserts at fw_cert.py l.385).
2. d ≤ -0.0984 < 0.
   - g = [0.08962 ± 1.23e-6] (runs6/fw_cert_Ka55_K2000.log l.3), with b_-- ≤ b_hi < 0 on [-s0, 0].
   - d = b_--(0), PROVED (v1 argument, checked by the referee).
3. ℓM_s r ≠ 0, PROVED. It follows from λ_+ - λ_- = 2cN(1-V²)/(1-c²V²) > 0 at the sonic point. Then ℓM_x'r = λ'ℓM_s r ≠ 0.
4. |ρ(2) - n| ≥ n + 1.9085 for all n ≥ 0.
5. **Δ2(2) ≠ 0** [v2 2026-11-18, REFEREE-QB2 1(b) Fix + POST HOC identity].
   - **Affinity is structural:** κ enters only the ℓ̂P row (P1 sits in the fluid block). v1 wrote "exact second
     difference", which is wrong: the seal checks that an interval encloses 0 (count_b4v3_seal.py l.154).
   - **κ_pole = d exactly, PROVED** (REFEREE-QB2 1(b); r2_pole.log). At the sonic point, M_x r = 0, and the identity row's
     fluid part satisfies c_f·r = -8Vω·D/(9(V²-1)²), with D ∝ det M_x = 0.
     - So I annihilates r̂ there, and Δ2(κ) = const·(ℓE_u r - κℓM_s r) = const·ℓM_s r·(d - κ).
     - The constant is nonzero wherever c_f(0) ≠ 0 (ω > 0).
     - Hence **Δ2(2) ≠ 0 follows from d < 0 alone** (constant 2). The seal-identification step of v1 is no longer needed.
   - This explains v1's 13-digit "coincidence" (certified detA root = [-0.0990290820550 ± 2.91e-14], verdict_M.json l.6).
   - **Legacy note:** PROOF-B4.md l.64's float line Δ2 = -0.264162 - 2.667782κ puts its root at -0.0990193, off d by about
     1e-5. That line is inaccurate (legacy float) and is no longer used.

**Proof.**
1. **Order 0.** ℓ̂(P(0)y_0 + S(0)) = 0 and c(0)·y_0 = σ_I(0). S(0) includes rows 2 and 3, so ℓ̂S(0) ≠ 0 in general; it
   enters only the right-hand side. The 2×2 matrix is Δ2(2) ≠ 0.
2. **Order n ≥ 1.** The coefficient of β_n stays (ℓM_x'r)(ρ(2) - n) ≠ 0.
3. **Convergence.** By majorant, or via the augmented system Ŷ = (Y, z), whose sonic exponents are {0,0,0,0, ρ(2)}
   (Coddington-Levinson ch. 4, as in PROOF-B4.md l.62).
4. **Constraint.** J is analytic, J' = -AJ (Lemma (a)(i-H), conserved source), and J(0) = 0, so J ≡ 0.
   Uniqueness holds modulo S_2. ∎

The singular branch |x|^{ρ(2)} is unbounded (Re ρ(2) < -1.9), as in PROOF-B4.md l.117.

## 4. Lemma (c): the centre step

[v2 2026-11-18, REFEREE-QB2 1(c) Fix 1-3: v1's premise S_c = (S_H, S_N, 0, 0) was false and is withdrawn. Redone with the
true orders.]

**Lemma (c).** In the blown-up chart v = p e^x, suppose Q⁻¹S(v) and σ_I(v) are analytic at v = 0, with Q⁻¹S = O(v²) in
every row and σ_I = o(1). The real source has (O(v²), O(v²), O(v²), O(v³)) and σ_I = O(v²) (table below). The C1 control
source φ'G_2 is O(v²) in every row. Then the bounded forced solutions with Ā_p → 0 form an affine 2-plane
R_2^F = y_c^p + R_2, with y_c^p = O(v²) and J ≡ 0.

**The orders (SYM, exact; seats/design-2/centre_orders.log). This is REFEREE-QB2's next strike, done.**
- **Inputs:**
  - the four HKA rows above;
  - Q = diag(1,1,M_x);
  - the chart N = 1/P, ω = WP², V = UP with U → -1/2 and W0 free (lin_sym C6);
  - the source orders of a regular centre: β_nn = P²(a0 + ...), β_nr = P³(e0 + ...), β_rr = P²(r0 + ...).
- **Result, generic coefficients:**

| row of Q⁻¹S | order | leading coefficient (isotropic r0 = a0/3) | leading coefficient (r0 free) |
|---|---|---|---|
| 0 | v² | a0 | a0 |
| 1 | v² | -a0/3 | (r0 - a0)/2 |
| 2 | v² | a0/3 | a0 - 2r0 |
| 3 | v³ | -(5a0 - 3e0)/12 | (-3a0 + e0 + 4r0)/4 |
| σ_I | v² | e0 - a0 | e0 - a0 |

- Before Q⁻¹, row 3 of S is O(v) (coefficient a0/9 isotropic), and M_x⁻¹ ~ 1/N adds one power. The orders (v², v², v², v³)
  hold **without** isotropy (last column), confirming the referee's POST HOC estimate.
- σ_I is O(v²), not "O(v³) at worst" as v1 said (REFEREE-QB2 1(c) Fix 2). The lemma needs only o(1).
- **Premises still PAPER:** the source orders β_nn = O(v²), β_nr = O(v³), β_rr = O(v²).
  - B_(ab) is smooth at a regular centre; B_(nr) is odd in r; b = r⁴B = O(v⁴), b_nr = O(v⁵) (DESIGN-QB1 sec. 2).
  - QB1 numerics agree: b_nn e^{-4x} → c0 and b_rr/b_nn → 1/3 (REFEREE-QB1 l.47-52).
  - P = v + O(v²) is analytic, so orders in P are orders in v.
  - The certified check is REDUCED: evaluate the first θ-coefficients of the four rows on the cm manifold (DESIGN-QB2-v2 node N4a).

**Proof.**
1. v dY/dv = 𝓜_c(v,2)Y + Q⁻¹S, with 𝓜_c(0) diagonalisable: eigenvalues -1, 0, 0, -2, kernel span(e_N̄, e_ω̄).
2. Order 0: Q⁻¹S(0) = 0, so y_0 ∈ ker. Order 1: Q⁻¹S has no v¹ term, so y_1^p = 0.
3. For n ≥ 1, (n - 𝓜_c(0)) is invertible. The particular series starts at v² and converges by majorant.
4. The other directions grow like v⁻¹ and v⁻², and the source cannot cancel them.
5. Constraint: J = J(x1)exp(-∫A), so |J| ~ Ce^{-x}. On R_2^F, J → 0 (c_N, c_ω, c_V → 0; Ā_p → 0; σ_I = O(v²) → 0), hence
   C = 0. ∎

## 5. Lemma (d): the dimension count

**Lemma (d).** S_2^F ∩ R_2^F ≠ ∅, and it is unique modulo c·G_2.

**Inputs.**
- D_new(2) ≠ 0, two ways:
  - logically: verdict_M.json l.4-5 (PM1, PM2), with the preconditions at l.6;
  - directly: pq7b.json l.13-35: D_new(2) = [-2.4126270857981028564 ± 1.13e-9] and
    D_HKA(2) = [-19.942295267987575083 ± 1.12e-8].
  - The referee re-checked both.
- dim S_2 = dim R_2 = 2 and dim W_2 = 3.

**Proof.** 2 + 2 - 1 = 3 gives S_2 + R_2 = W_2. The affine argument is as in v1. Two elements of the intersection differ
by an element of span(G_2). ∎

**Continuation and source scope** [v2 2026-11-18, REFEREE-QB2 1(d) Fix].
- The sonic point is crossed analytically, and the coefficients are analytic on (0, 12] (no second sonic point; E5,
  LOG.md l.45).
- **The certified source on [-12, 12] is not on file.**
  - QB1b's PQ2b (seats/build-2/result_qb1b.json) is a **sign map of b_nn** only: 2261 cells with signs.
  - **QB2 must itself produce certified Taylor models of b_nn, b_nr, b_rr and b_θθ** on every transport step (DESIGN-QB2-v2
    node N0). v1's wording "the source is the certified QB1b profile" is withdrawn.
- x > 12 is DESIGN node T1. It is not needed for Ψ2.

**Practical form** [v2 2026-11-18: gauge stated consistently with control C1].
- **Gauge.** Fix the gauge by N̄_p(0) = 0 (HKA normalisation). This is possible since N̄ of G_2 at 0 is 2 - κ_g ≠ 0.
- **Sonic family.** y_s^p + b·s_2^HKA, with y_s^p chosen so that (Ā, N̄)(0) = (0, 0).
- **Matching.** Solve at x_m = -4.6 with the D_HKA(2) matrix [s_2^HKA, ρ1, ρ2, e_V]:
  y_c^p - y_s^p = b·s_2^HKA + c·ρ1 + d·ρ2 + e·e_V.
  - |D_HKA(2)| ≥ 19.94 (pq7b.json).
  - e·c_V(x_m) = J_c - J_s = 0 for a conserved source: the conservation check.

## 6. The u e^{2u} secular structure

Unchanged from v1 sec. 6. The referee rated it CORRECT (REFEREE-QB2 sec. 2):
- s·e^{2s}: solvable as posed, unique modulo G_2. There is no s²·e^{2s} term: that would need a pole of the resolvent at
  κ = 2, and G_κ exists at every κ (it is gauge, not a resonance).
- e^{2s}: state-dependent, with a non-conserved state-independent part, since ∇(uV) = (∇u)V (zahn.txt l.348).
- The ambiguity of h1 lies along span(h2, G_2), as in BO l.167.
- BO's ansatz eq. 3.12 lacks the s·e^{2s} term; Zahn l.152 says so too.

The table and items 1-4 of v1 sec. 6 apply verbatim. For "L(2)" read the HKA pencil with the four-row source.

## 7. Still owed

1. A certified version of the centre source orders on the cm manifold (node N4a; the SYM structure is done).
2. The convergence constants in (b) and (c) are majorant-level. The build replaces them with TM tails.
3. A referee pass on this v2 page.
