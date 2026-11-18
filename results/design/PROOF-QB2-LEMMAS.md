# PROOF-QB2-LEMMAS: the forced O(hbar) linear response at rate 2 is solvable and unique modulo gauge

Proof seat design-2, 2026-11-18. Paper, reading and two small exact sympy checks (non-counting; scripts and logs in
`seats/design-2/`). These are the four lemmas that RESULT-QB1.md rule 7 (l.16) and REFEREE-QB1 C4 fix 2
(seats/referee-1/REFEREE-QB1.md l.68-73) say QB2 owes. The u e^{2u} secular structure is in sec. 6.

**Labels.**
- **PROVED**: complete on this page; still needs a referee.
- **SYM**: exact sympy identity, run here (log named).
- **CERTIFIED**: an interval enclosure already on file in a sealed/HOLD run (file:line given).
- **CITED (OPENED / ABSTRACT / SNIPPET)**: outside fact, with locator and how much of the source was read.
- **REDUCED**: reduced to a finite list of ball evaluations, not yet run.

## 0. In plain words

- **What is solved.** The quantum source grows like s e^{2s} (s = self-similar time). The metric and fluid answer it with a
  response e^{2s}(h1(x) + s·h2(x)). Only h2 is fixed by geometry alone. h1 needs the quantum state.
- **What is proved here.** The h2 problem has exactly one solution up to the gauge solution G_2, which is a relabelling of
  time, not physics. Four steps:
  - (a) the constraint stays solved if the source is conserved;
  - (b) the sonic point lets through exactly a 2-parameter family;
  - (c) so does the regular centre;
  - (d) the two families meet in exactly one point modulo G_2, because 2 is not a growth rate of the critical solution.
- **What is not solvable.** The h1 problem has a state-dependent source, and its state-independent part is not
  conserved on its own. It has no meaning without the state.
- **Two exact by-products (SYM, sec. 1).**
  - The constraint propagation rate is exactly c_I = -A.
  - The forcing enters the constraint only through the time component of div(Sigma), with weight -A².

## 1. Setting

- **Background and unknowns.** HKA bar variables Y = (Ā_p, N̄_p, ω̄_p, V_p)(x). Pencil Q Y' = (P0 + κP1)Y, Q = diag(1,1,M_x),
  identity row I(Y) = c(x,κ)·Y (PROOF-B4.md l.24-31; engine/cert/lin_sym.py l.43-57).
- **The 5 PDE rows (engine/derive.py l.101, l.147).**
  - Hamiltonian (A_x), slicing (N_x) and momentum (A_s) are solved from the Einstein components E_ss, E_sx, E_xx.
  - Two fluid rows come from div T = 0.
  - E_θθ is not used; it follows from the contracted Bianchi identity.
- **Forced problem.** Add a symmetric source Σ to 8πT in the Einstein equations and leave the fluid rows alone.
  - The quantum field is a separate conserved system, so div T_fluid = 0 still holds.
  - Here Σ = 8π⟨T⟩_ren. Its state-independent part is Σ_V = (N_c ħ/(30π)) s B[g], with B[g] = 60 V[g] on R = 0
    (DESIGN-QB1.md sec. 1.2; RESULT-QB1 rule 1). N_c > 0 is the field multiplier.
- **The κ form.** Write Σ = ε e^{κs} σ(x) and h = ε e^{κs} Y(x). The rows become
  Q Y' = (P0 + κP1)Y + S(σ),  I(Y) = σ_I(σ).
  - S has nonzero entries only in rows 0 and 1; the fluid rows carry no source (SYM, srcprop.log l.16-17).
- **Source rows (SYM, seats/design-2/srcprop.log l.13-15; frame translation PAPER).** In derive normalisation:
  - Hamiltonian: Ā_p' gets + (1/N²)(σ_ss + 2σ_sx + σ_xx) = A r² Σ_(nn).
    Here ∂_s + ∂_x = -t∂_t|_r, and N² = α²/(A e^{2x}).
    This agrees with DESIGN-QB1 sec. 1.2: S_H = (N_c/(30π)) A e^{-2x} b_nn per unit ε s e^{2s}.
  - Slicing: N̄_p' gets -(1/(2N²))(σ_ss + 2σ_sx + σ_xx - N²σ_xx) = -(1/2) A r²(Σ_(nn) - Σ_(rr)).
    Control: on the background this reproduces the matter term -(2/3)ω of KHA (18b) for T_(nn) - T_(rr) = ρ - p at V = 0.
  - Momentum: the coefficients are printed in srcprop.log l.15.
- **Notation.** κ = 2 throughout unless stated. W_κ = {solutions with I ≡ 0}, dim 3 (PROOF-B4.md l.47).

## 2. Lemma (a): constraint propagation with a conserved source

**Lemma (a).** Let Y solve the four evolution rows with source S(σ), and set J := I(Y) - σ_I.
- (i) The following holds identically:
  J' = -A·J + ϑ(x)·Comp(σ),  Comp(σ) = -A² · e^{2(x-s)} e^{-κs} [ (∇·Σ)_s + (∇·Σ)_x ].
  - ϑ ≠ 0 is the fixed ratio between the HKA identity row and the derive.py momentum row (lin_sym.py C3, l.86-88).
  - (∇·Σ)_s + (∇·Σ)_x is the component of div Σ along ∂_s + ∂_x = -t∂_t|_r, the polar time component.
- (ii) If (∇·Σ)_t = 0, then J(x) = J(x0)·exp(-∫_{x0}^x A). So J(x0) = 0 at one point gives J ≡ 0, and
  {forced solutions with J ≡ 0} is an affine 3-space with direction W_κ.
- (iii) The unused θθ Einstein equation holds iff, in addition, (∇·Σ)_r = 0.
- (iv) For the QB2 top-order problem (the s-coefficient, κ = 2) both components vanish, because ∇·V[g] = 0 identically.

**Proof.**
1. **Homogeneous identity.** In derive variables, d/dx(P_2·hp) = ct·(P_2·hp) along the evolution rows, for every κ.
   SYM P1 (srcprop.log l.4) re-derives lin_sym.py C4 (l.89-95).
   - In HKA variables the rate is c_I = (8NVω + 2V²ω - 3V² + 6ω + 3)/(3(V² - 1)) (PROOF-B4.md l.41).
   - The same numerator gives A = Acons = -(8NVω + 2V²ω - 3V² + 6ω + 3)/(3(V² - 1)) (seats/design-2/xv.log l.2).
   - Hence **c_I = -A exactly** (PAPER, by comparing two printed exact expressions).
2. **Forced identity.** Write m = P_2·hp + C_2·σ and hp' = -L_q⁻¹(L_p hp + L_c σ). Then
   m' - ct·m = -(P_2 L_q⁻¹ L_c)σ + (C_2)'σ + C_2σ' - ct·C_2σ =: Comp(σ).
   - This is an exact first-order expression in σ = (σ_ss, σ_sx, σ_xx).
   - SYM P2 (srcprop.log l.7, l.10): Comp(σ) = l_s·D_s + l_x·D_x exactly. D_ν is the covariant divergence of
     Σ_μν = e^{κs}σ_μν(x) on the on-shell background (σ_θθ included), with the weight e^{2(x-s)}e^{-κs} stripped.
   - l_s = l_x = -(8NVω + 2V²ω - 3V² + 6ω + 3)²/(9(V² - 1)²) = **-A²** ≠ 0 (srcprop.log l.8-9). They are κ-free (l.11).
   - σ_θθ enters D_s and D_x but cancels in D_s + D_x, as it must: σ_θθ does not enter the three rows used.
   - Moving to HKA rows multiplies by ϑ. This gives (i).
3. **(ii)** is the linear ODE J' = -AJ with a regular coefficient on x ≤ 0, at the sonic point included (A is analytic
   there, E5). The affine-space statement follows from Cor. 1.2 (PROOF-B4.md l.44-47).
4. **(iii)** By the contracted Bianchi identity and div T_fluid = 0, ∇^μ(E - Σ)_μν = -(∇·Σ)_ν. Here E = G - 8πT.
   - When the three used rows and J = 0 hold, E_ss = E_sx = E_xx = 0.
   - The ν = x component then reads E^θ_θ·(2Γ^θ_θx-type coefficient) = -(∇·Σ)_x. The coefficient is ≠ 0 because r_x = r.
   - So E_θθ = 0 iff (∇·Σ)_r = 0, with x ↔ r at fixed t because Σ_xx and Σ_sx are already fixed.
   - This last step is PAPER, not SYM. The θθ residual is listed as an instrument check in DESIGN-QB2 node N2.
5. **(iv) Physical-space version.** The identity in (i) is polynomial in κ with s-independent coefficients. It is
   therefore the Laplace image of an operator identity in ∂_s. Apply it to Σ = s e^{2s}σ2(x) + e^{2s}σ1(x):
   - The s e^{2s} coefficient obeys the κ = 2 identity with σ2 alone, because P(∂_s)[s e^{2s}f] = s e^{2s}P(2)f + e^{2s}P'(2)f.
   - For Σ_V the s-coefficient of ∇·(sV[g]) is ∇·V[g]. This vanishes identically: V is the metric variation of a local
     action (CITED OPENED, zahn.txt l.348: V "is automatically conserved").
   - The remainder (∇s)·V[g] sits at the e^{2s} level only (zahn.txt l.348, same sentence). ∎

**Consequence (content null, used in DESIGN-QB2).** A source that is not conserved in its t component produces
J' = -AJ + ϑ Comp with Comp ≢ 0. The forced sonic solution (J(0) = 0) then reaches the centre with J ≠ 0, which is
unbounded like e^{-x} there. The matching must therefore refuse. Doing this on the real Σ_V tests that QB1's certified
b-profiles are conserved, which is a non-trivial check on QB1.

## 3. Lemma (b): the forced sonic Frobenius step at κ = 2

**Lemma (b).** Let S and σ_I be analytic at x = 0 with radius ≥ ρ_S. Then the set of forced solutions that are
analytic at x = 0 and have J ≡ 0 is an affine 2-plane S_2^F = y_s^p + S_2. Its direction S_2 is the homogeneous
sonic-analytic space of PROOF-B4.md l.63. It is parametrised by (Ā_p, N̄_p)(0) ∈ ℂ², and the series converge on
|x| < min(ρ_S, ρ_bg).

**Constants (all certified or exact, none float).**
1. **λ' := λ_-'(0) ∈ [1.0986, 1.0995].** CERTIFIED: PROOF-B4-v4.md l.81, from engine/cert/fw_cert.py on the E5 enclosure.
   The run asserts dλ_-/dx > 0 on the whole sonic piece [-s0, 0] (fw_cert.py l.385; log runs6/fw_cert_Ka55_K2000.log l.3).
2. **d ≤ -0.0984 < 0.** fw_cert.py l.385 asserts b_-- ≤ b_hi < 0 on [-s0, 0], which contains x = 0. The log prints
   g = |b_hi|_lower/λ'_hi = [0.08962 ± 1.23e-6] (fw_cert.py l.386; runs6/fw_cert_Ka55_K2000.log l.3), so
   |b_hi| ≥ 0.089618 × 1.0986 > 0.0984.
   - **d = b_--(0) (PROVED, one line).** At x = 0, λ_- = 0, so b_--(0) = l_- M_s⁻¹ E_u r_- (PROOF-B4.md l.159).
   - l_- Λ0 = 0 forces l_- M_s⁻¹ ∝ ℓ (the left null vector of M_x(0)), and l_- r_- = 1 gives l_- M_s⁻¹ = ℓ/(ℓ M_s r).
   - Hence b_--(0) = ℓE_u r/(ℓM_s r), which is the d of Fact 1.3 (PROOF-B4.md l.54).
3. **ℓ M_s r ≠ 0 (PROVED, exact).**
   - λ_+ - λ_- = 2cN(1 - V²)/(1 - c²V²) with c = 1/√3 (from C5, PROOF-B4.md l.42; algebra
     (V+c)(1-cV) - (V-c)(1+cV) = 2c(1-V²)). At the sonic point (λ_- = 0, N > 0, |V| < 1), λ_+(0) > 0.
   - So 0 is a simple eigenvalue of Λ0 = M_s⁻¹M_x, and for a simple eigenvalue the left-right pairing ℓM_s r ≠ 0.
     M_s is invertible, since det M_s = γ(1 - (γ-1)V²)/(1 - V²) > 0.
   - Then ℓM_x'r = λ'·ℓM_s r ≠ 0, by item 1.
4. **ρ(2) = (d - 2)/λ' ≤ -(2 + 0.0984)/1.0995 < -1.9085.** So for every n ≥ 0,
   **|ρ(2) - n| ≥ n + 1.9085.**
   - This needs only the sign of d and the upper end of λ'. No float constant from Fact 1.3 (PROOF-B4.md l.56-61) is used.
   - The float value is ρ(2) = -1.90989.
   - Curiosity (UNTESTED, from PROOF-B4.md l.60): if d + λ' = 1 exactly, then ρ(2) = -1 - 1/λ' ∈ [-1.91025, -1.90951].
5. **Δ2(2) ≠ 0.** Δ2 is affine in κ (PROOF-B4.md l.64). The M seal checks affinity (exact second difference) and that the
   root is real and negative: detA root = [-0.0990290820550 ± 2.91e-14] (CERTIFIED, runs7/verdict_M.json l.6 PM3;
   engine/cert/count_b4v3_seal.py l.12-13, l.149-155).
   - Identifying the seal's "detA" with PROOF-B4's Δ2 is by reading those two files (same definition: the ω̄ and V
     columns of [ℓ̂P(0); c(0)], Cramer form PROOF-B4.md l.63-68).
   - An affine function with its only root at κ < 0 does not vanish at κ = 2.
   - Side note (UNTESTED): the certified root equals the float d = -0.0990290820550211564 (PROOF-B4.md l.59) to 13 digits.

**Proof.**
1. **Order 0.** Left-multiply Q Y' = (P0 + 2P1)Y + S by ℓ̂ = (0,0,ℓ). Since ℓ̂ Q(0) = 0, analyticity forces
   ℓ̂(P(0)y_0 + S(0)) = 0.
   - Together with J(0) = 0, i.e. c(0)·y_0 = σ_I(0), this is a 2×2 linear system for (ω̄_p, V_p)(0).
   - Its matrix is Δ2(2) ≠ 0 (constant 5), for any given (Ā_p, N̄_p)(0). So y_0 = y_0^p + (homogeneous y_0), and the
     homogeneous data form the 2-dimensional set of PROOF-B4.md l.63.
2. **Order n ≥ 1.** The recursion is PROOF-B4.md l.52-54 with σ_n added to R_n. The source never touches the
   coefficient of the free direction β_n r̂, which stays (ℓM_x'r)(ρ(2) - n). By constants 3 and 4 this is nonzero for
   every n ≥ 1 (indeed for every n ≥ 0). So β_n, and hence y_n, is uniquely determined.
3. **Convergence.** Majorant argument.
   - Let C_bg, ρ_bg bound the Taylor coefficients of Q, P0, P1 at 0 (E5 sonic manifold, radius about 1 in x;
     DESIGN-QB1 sec. 2), and C_S, ρ_S those of S and σ_I.
   - Solving step 2 gives |y_{n+1}| ≤ K·(Σ_j C_bg ρ_bg^{-j}|y_{n-j}| + C_S ρ_S^{-n})/(n + 1.9085).
     K = max(1/|Q-block pivots|, |r̂|/|ℓM_x'r|) is finite by constant 3.
   - The comparison series solves a scalar first-order regular-singular equation with exponent < 0 and an analytic
     right side. It converges on |x| < min(ρ_bg, ρ_S).
   - Equivalent form: augment Ŷ = (Y, z), z' = 0, and the source becomes a column S·z. The augmented system is
     homogeneous, its sonic exponents are {0,0,0,0, ρ(2)}, and none of them differ by a positive integer from ρ(2).
     Convergence is then the cited homogeneous theorem (Coddington-Levinson 1955 ch. 4, as in PROOF-B4.md l.62).
4. **Constraint.** J is analytic at 0, satisfies J' = -AJ for x ≠ 0 (Lemma (a), conserved source), and J(0) = 0.
   By continuity J ≡ 0. So y_s^p + S_2 ⊂ {J ≡ 0}.
   - Uniqueness: two such forced solutions differ by a homogeneous sonic-analytic solution with I ≡ 0, i.e. by an element
     of S_2. ∎

**Why the singular branch is excluded at κ = 2.** Re ρ(2) < -1.9 < 0, so |x|^{ρ(2)} is unbounded at x = 0. A forced
solution that is bounded (or in H^1) near the sonic point is therefore in S_2^F, the same argument as PROOF-B4.md l.117.

## 4. Lemma (c): the centre step for an O(v²) source

**Lemma (c).** In the blown-up chart v = p e^x (PROOF-B4.md l.70-73), suppose:
- the source column S_c(v) := Q⁻¹S = (S_H, S_N, 0, 0)(v) and σ_I(v) are analytic at v = 0 (radius ρ_c > 0);
- S_c(v) = O(v²) and σ_I(v) = o(1).

Then the forced solutions that are bounded as v → 0 with Ā_p → 0 form an affine 2-plane R_2^F = y_c^p + R_2, with
y_c^p = O(v²). Every element has J ≡ 0.

**The source orders (PAPER, from QB1).**
- B is a smooth tensor at a regular centre with B_(nn) → c0 t⁻⁴ (DESIGN-QB1 sec. 2 table). So b_ab = r⁴B_ab = O(e^{4x}) = O(v⁴),
  and S_H = (N_c/30π) A e^{-2x} b_nn = O(v²). Q acts as the identity on rows 0 and 1.
- The certified value is c0 = -173.1004443922995663248 (REFEREE-QB1 l.51). It was obtained by Richardson extrapolation in
  e^{2x}, and the deep-centre series in w = z² is used in seats/build-2/RUN.md l.42-45.
- S_N uses Σ_(nn) - Σ_(rr), which is also O(v²). B_(nr) → 0 like r (DESIGN-QB1 sec. 2), so σ_I = O(v³) at worst.
- Analyticity in v: b_ab is built from the centre manifold (P, W, U, X)(θ), which is analytic in θ (cm, K 60). Also
  v = θ e^{X(θ)} is an analytic change of variable with v'(0) = 1. Both statements are PAPER from the chart definition in
  seats/build-1/bg.py l.5-6.

**Proof.**
1. The forced system is v dY/dv = 𝓜_c(v,2)Y + S_c(v). Here 𝓜_c(0) = [[-1,0,0,0],[1,0,0,0],[-3,0,0,0],[0,0,0,-2]]
   (SYM, lin_sym.py C6; PROOF-B4.md l.73). It is diagonalisable:
   - eigenvalue -1 with eigenvector (1, -1, 3, 0);
   - eigenvalue 0 with eigenvectors e_N̄, e_ω̄;
   - eigenvalue -2 with eigenvector e_V.
2. Seek Y = Σ y_n v^n. At order 0: 𝓜_c(0)y_0 = S_c(0) = 0, so y_0 ∈ ker = span(e_N̄, e_ω̄). Hence Ā_p(0) = 0 and V_p(0) = 0.
3. At order n ≥ 1: (n - 𝓜_c(0))y_n = Σ_{j≥1} 𝓜_j y_{n-j} + s_n. The eigenvalues n, n, n+1, n+2 are all ≥ 1 (Fact 1.4,
   PROOF-B4.md l.77), so this is uniquely solvable.
   - With y_0 = 0, s_0 = s_1 = 0 gives y_1 = 0, and the particular series starts at v².
   - Convergence: ‖(n - 𝓜_c(0))⁻¹‖ ≤ κ_E/n, where κ_E is the condition number of the eigenvector matrix (a fixed
     rational number). The majorant is as in Lemma (b), step 3. It converges on |v| < min(ρ_c, centre series radius 0.51;
     PROOF-B4.md l.144).
4. **Everything else is unbounded.** The two remaining directions grow like v⁻¹ (point mass in Ā_p) and v⁻² (V_p), as in
   Fact 1.4 (l.80). The source is O(v²) and so cannot cancel them.
5. **Constraint.** By Lemma (a), J = J(x1)exp(-∫_{x1}^x A), and A → 1 at the centre (A = 1 + O(v²), Acons at P = 0), so
   |J| ~ C e^{-x}. On R_2^F, J = I(Y) - σ_I → 0:
   - c_N, c_ω, c_V → 0 (PROOF-B4.md l.46);
   - Ā_p → 0 and the first entry κ - A is bounded;
   - σ_I → 0.
   Hence C = 0 and J ≡ 0. This is Cor. 1.2(b) (l.46) with the source added. ∎

**Remark.** An O(v) source would also do: only S_c(0) ∈ Im 𝓜_c(0) is needed at order 0. O(v²) leaves a margin of two
orders. A source with S_c(0) ∉ Im 𝓜_c(0) would force a log or a point mass at the centre (the Ayal-Piran type centre
problem; DESIGN-QB1 sec. 2 table). For the local V source this cannot happen.

## 5. Lemma (d): the dimension count, solvable and unique modulo G_2

**Lemma (d).** With S, σ_I as in (b) and (c) and conserved (Lemma (a)(iv)):
- S_2^F ∩ R_2^F ≠ ∅;
- any two elements of it differ by c·G_2.

So the top-order response h2 exists on x ≤ 0, is unique modulo the gauge solution, and continues analytically to ℝ
wherever the coefficients and the source are regular.

**Inputs.**
- **D_new(2) ≠ 0**, which gives S_2 ∩ R_2 = span(G_2) (Lemma 1.5(c) and Thm 1.6(b), PROOF-B4.md l.91-93, l.105-108).
  It is CERTIFIED two ways:
  - (i) logically: runs7/verdict_M.json l.4 (PM1: N_R40 = 1, a single zero in {Re κ ≥ 0, |κ| < 40}) and l.5
    (κ* ∈ [2.810552549 ± 3.83e-10] ∌ 2), with the preconditions c_V(x_m) ≠ 0 and Ā_ss'(0) ≠ 0 at l.6;
  - (ii) directly: D_new(2) = [-2.4126270857981028564 ± 1.13e-9], |D_new(2)| ≥ 2.41262708466
    (seats/build-2/pq7b.json l.13-35).
  - Wording: REFEREE-QB1 C4 fix 1 (l.68) says "certified" here means computed with the sealed count instrument, and
    PQ7 as first sealed FAILED. The logical route (i) does not depend on PQ7/PQ7b.
- dim S_2 = 2 (Lemma (b)); dim R_2 = 2 (Fact 1.4, PROOF-B4.md l.79); both lie in W_2, which has dim 3 (Cor. 1.2).

**Proof.**
1. dim(S_2 + R_2) = 2 + 2 - dim(S_2 ∩ R_2) = 2 + 2 - 1 = 3 = dim W_2. Hence S_2 + R_2 = W_2.
2. Take y_s ∈ S_2^F and y_c ∈ R_2^F (Lemmas (b), (c)). Both satisfy the forced rows and J ≡ 0, so y_s - y_c ∈ W_2.
   Write y_s - y_c = a - b with a ∈ S_2 and b ∈ R_2. Then y* := y_s - a = y_c - b lies in S_2^F ∩ R_2^F.
3. If y*, y** both lie in the intersection, their difference lies in S_2 ∩ R_2 = span(G_2).
4. **Continuation to x > 0.** The sonic point is crossed analytically (Lemma (b)). On (0, 12] the coefficients are
   analytic: the certified outward leg has no second sonic point (PROOF-B4.md l.86; BG2 record).
   - On [-12, 12] the source is the certified QB1b profile (seats/build-2/result_qb1b.json, coverage
     "complete, no gap").
   - For x > 12 a tail bound is needed. That is DESIGN-QB2 node T1, and it is not needed for the headline functional. ∎

**Practical form (used by the build).** The matching at x_m uses exactly the D(2) matrix [s_2, ρ1, ρ2, e_V](x_m)
(PROOF-B4.md l.100):
  y_c^p(x_m) - y_s^p(x_m) = a·s_2 + c·ρ1 + d·ρ2 + e·e_V   (gauge fixed by dropping s_2^HKA, i.e. Ā_p(0) = Ā^p(0)).
- Cramer's rule with |D_new(2)| ≥ 2.41 gives (a, c, d, e).
- **e is a constraint residual.** Lemma (a) and Thm 1.6(b) (l.106) force e·c_V(x_m) = J-difference = 0. A certified
  e ∋ 0 of tiny width is the conservation check, and a certified e ∌ 0 is a refusal. This is DESIGN-QB2 node N6.

## 6. The u e^{2u} secular structure: what is forced at s e^{2s} and what at e^{2s}

**Ansatz.** h = ε e^{2s}(h1(x) + s·h2(x)). This is forced by the source structure Σ = ε e^{2s}(sσ_V + σ_T), Zahn eq. 6
(CITED OPENED, zahn.txt l.136). Zahn notes that the plain e^{qu} ansatz "has to be slightly modified to also include a
supplementary linear growth in u" (zahn.txt l.152).

Write L(κ) := Q∂_x - (P0 + κP1) on the κ form. Then L acting on s e^{2s}h2 gives e^{2s}(s·L(2)h2 + ∂_κL·h2), with
∂_κL = -P1. Also ∂_κI = e_Ā, from the κ - A entry.

| order | equation | source | status |
|---|---|---|---|
| **s e^{2s}** | L(2)h2 = S(σ_V), I_2(h2) = σ_I(σ_V) | σ_V = (N_c/(30π))·(rows of B[g]) per unit ε s e^{2s}. **State-independent and scheme-independent**: V is "state-independent and independent of the renormalization scheme" (zahn.txt l.133-134, OPENED); with R = 0, B/60 for every ξ (RESULT-QB1 rule 1) | **solvable as posed**: Lemmas (a)-(d) give h2 unique mod G_2 |
| **e^{2s}** | L(2)h1 = S(σ_T) + P1·h2 - S(f·σ_V) + 2a-term; I_2(h1) + Ā_{2}·1 = σ_I(σ_T) + ... | σ_T = the state-dependent T(y) of Zahn eq. 6; f(x) from s = u + f(x) if Zahn's u differs from KHA s; the renormalised C² coupling a (BO eq. 3.10-3.11, bo.txt l.220-226, OPENED) adds 2aB, i.e. a multiple of σ_V | **not solvable as posed**: T(y) is not determined by geometry (zahn.txt l.134) and no state is constructed on EC (zahn.txt l.126, l.131) |

**Exact statements (PROVED, from the table).**
1. **No s² e^{2s} term.** A source with s¹ forces s² only if L(2) has a kernel that obstructs (resonance). By Lemma (d),
   the kernel at κ = 2 is span(G_2), and the forced problem is solvable without a secular s² term. If 2 were a mode,
   s² e^{2s} would appear.
2. **The e^{2s} source is not conserved by itself, and its state-independent part is not either.**
   - ∇·(sV[g]) = (∇s)·V[g] ≠ 0 lands at the e^{2s} level. Only the total, including the state term, is conserved: Zahn
     says the state term's divergence cancels it (zahn.txt l.348; on Roberts, zahn.txt l.636-637).
   - So "the state-independent part of h1" (for example the ∂_κ-trick piece -L(2)⁻¹∂_κL h2) has a non-conserved source.
     By Lemma (a)(i) it violates the constraint. It is **not** a well-posed object.
   - Corollary: QB2 can certify h2 and nothing else at this order.
3. **The ambiguity of h1 is exactly along h2 (and G_2).**
   - A shift of the origin of s by s0 maps e^{2s}(h1 + s h2) to e^{2s0}e^{2s}((h1 + s0 h2) + s h2).
   - The finite renormalisation T(y) → T(y) + bV(y) (on EC, T(1) = 0 and T(2) = -240V; DESIGN-QB1 sec. 1.2 Consequence 3)
     adds a multiple of h2 to h1.
   - The C² coupling a does the same. It is BO's "arbitrary renormalization scale hidden in the logarithm" (bo.txt l.167, OPENED).
   - So even with the state known, h1 is defined only modulo span(h2, G_2), while h2 is canonical modulo G_2. This is the
     QB2 form of DESIGN-QB1's "ambiguous only along V itself".
4. **Brady-Ottewill's ansatz.** BO eq. 3.12-3.13 (bo.txt l.229, OPENED) use g_q e^{κ_q u} with no u e^{κ_q u} term. Their
   κ_q = 2 matching (bo.txt l.230) is correct. But with their own u e^{4u}B term (eq. 3.9, l.163-170) the response must
   carry s e^{2s}, by item 1 and Zahn l.152. BO's case analysis (l.255-270) uses only the exponents and is unaffected.

## 7. What is still owed (the gaps these lemmas leave)

1. Lemma (a)(iii) step 4 (the θθ equation via the r component) is PAPER only. A SYM check, or the instrument's θθ residual,
   closes it.
2. Lemma (c) needs a written bound S_c = O(v²) from the certified centre series: b_nn e^{-4x} → c0 and the θ-analyticity.
   It is a 1-line REDUCED item: evaluate the first two θ-coefficients of S_H and S_N and show they are 0 (structural).
3. The convergence constants in (b) and (c) are majorant-level. The certified build replaces them with Taylor-model tail
   bounds (DESIGN-QB2 nodes N3 and N4).
4. A referee seat on this page.
