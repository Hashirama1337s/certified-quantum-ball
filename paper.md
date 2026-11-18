# The quantum fingerprint on the critical ball of light: a certified one-loop Bach source on the Evans-Coleman solution, and the sign of its secular push on the peak compactness

**Moki & Julio** (AI-assisted)

Version 1.0, 2026. Every number below is traced to a result file named in `CLAIMS.md`. Balls are printed
as `[mid ± rad]`. Each printed ball is either the stored string of the certifying program, verbatim, or a shortened form
written by `tools/fill_balls.py`, which rounds outward: it keeps fewer midpoint digits, enlarges the radius by the
midpoint shift and rounds it up, and asserts that the printed ball contains the stored one. Interval ends are printed as
stored.

---

> **The magic trick.** Every free massless quantum field of spin up to 1 leaves the same fingerprint on the critical ball
> of light: the Bach tensor. A computer proof maps that fingerprint from x = −12 (next to the centre) out to x = 12. It is negative in
> the core, positive in the shell, and it has exactly one sign change on that range. Then the proof reads it on the
> sphere near x ≈ −0.25 where the gas of light stands still and 2m/r peaks locally. There it finds that the slowly growing quantum push
> raises (toward trapping) the ball's peak compactness.

> **In plain words**
>
> When a ball of light (a gas of radiation) collapses under its own gravity, a tiny change in the starting push decides
> its fate: it becomes a black hole or it flies apart. Exactly on the dividing line sits a special collapse that keeps
> its shape while it shrinks. This is the Evans-Coleman critical solution. In earlier work we proved that it exists and
> that it has exactly one way to fall.
>
> Einstein's theory knows nothing of Planck's constant ħ. Quantum fields living on the collapsing ball add a small extra
> stress. Most of that stress depends on which quantum state the fields are in, which nobody knows here. One piece does
> not: it is fixed by the geometry alone, and it grows like s·e^{2s} as the collapse proceeds (s is the logarithm of
> the time left). On this background that piece is the Bach tensor, a curvature quantity, times a positive number that
> depends only on which fields are present. This form follows in one line from a theorem of Zahn (2025) and the fact that
> radiation has zero trace, and, for spin 1/2 and photons, from the conformal transformation law of Brady and Ottewill (1998).
>
> **What we computed.** We evaluated the Bach tensor on the certified critical solution, with a guaranteed error bar on
> every number.
>
> - Its energy component is negative near the centre, positive further out, and changes sign exactly once on the range
>   we mapped. The change happens at x = ln(r/(−t)) ≈ −0.8261143, the sphere r ≈ 0.4377·|t| in the time of Koike, Hara and Adachi
>   (r ≈ 0.354·|t| in the proper time of the centre).
> - We then solved for how the ball responds to this push, to first order in ħ.
> - We read the response at the sphere near x ≈ −0.25 where 2m/r peaks locally. That reading does not depend on any choice of
>   coordinates. Its sign is positive: the slowly growing quantum push raises (toward trapping) the peak compactness.
>
> **What we did not compute.** We did not compute whether quantum effects move the black-hole threshold. That needs the
> unknown state-dependent part. We state the link as a conjecture under named conditions (Section 7.4).
>
> **The idea in one picture.** Think of a fingerprint pressed onto a shrinking balloon. Every kind of quantum field
> presses with the same finger; only how hard it presses differs. **Where the picture breaks.** The "finger" also
> depends on a quantum state. We assume that state is itself self-similar, and nobody has proved that such a state
> exists on this background.

---

## Abstract

We study the one-loop (first order in ħ) quantum stress on the Evans-Coleman critical solution. This is the continuously
self-similar threshold solution for a p = ρ/3 perfect fluid in spherical symmetry, whose existence and single growing
mode we certified earlier [MJ-EC].

For a free massless scalar of any curvature coupling, in any Hadamard state respecting self-similarity, Zahn proved
that the renormalised stress tensor has the form e^{2u}[T(y) + (u/4π²)V(y)], with V local, state-independent and
conserved [Zahn]. On this background R ≡ 0, and V = B/60 then follows in one line from Zahn's remark that, at conformal
coupling, V reduces to a multiple of the Bach tensor B (on R ≡ 0 the reduction holds for every ξ). For conformal spin-1/2 and Maxwell fields, the same structure with a positive
field-dependent multiplier follows from the exact conformal transformation law of Brady and Ottewill [BO]. We do not
claim either reduction as new. The new content is the evaluation of B on Evans-Coleman.

**QB1/QB1b (sealed, HOLD).** With ball arithmetic and Taylor models we certify the profile b_nn(x) = r⁴B_(nn) of the
Bach energy component over x ∈ [−12, 12], in 2261 cells:

- the full sign map: b_nn < 0 for x below the zero and b_nn > 0 above it;
- exactly one sign change, inside x ∈ [−0.82611434347927570343, −0.82611433975398540497] (width 3.7e-9), where b_nn′ ≥ 1.165;
- the sonic-point fluid-frame value b_uu = [0.4466737796695366321 ± 1.3e-17];
- the centre coefficient c0 = lim b_nn e^{−4x} = [−173.1004443922995663 ± 2.99e-13], with exact isotropy b_rr/b_nn → 1/3.

An off-shell content null has teeth, the sign of B is the same for every such field at every point, and an independent
code reproduces the sonic value to 25 digits and c0 to 22.

**QB2, QB2b, QB2c (sealed; each FAILED by item, with the headline item HOLD in all three).** We prove that the forced linear response at rate 2 is e^{2s}(h1 + s·h2), with h2
fixed by geometry and unique modulo gauge. We certify

- the gauge-invariant first-order change of the local peak compactness, Ψ2 ∈ [3.4616358037141656915 ± 3.02e-12] (per unit K·ε·s·e^{2s} in the Koike-Hara-Adachi clock, K = N_c/(30π), ε = ℓ_P²/t0²; Ψ2 = δ(2m/r) at x_V);
- its sign: positive.

**The threshold.** That sgn Ψ2 is also the direction in which ħ moves the black-hole threshold is a conjecture. It holds
only under stated infrared-naturalness conditions and a linear extrapolation (Section 7.4).

Every claim is conditional on the existence of a Hadamard state respecting self-similarity on Evans-Coleman, which is
not proved.

---

## 1. Results

Throughout, "every field" means every free massless field of spin ≤ 1: a scalar of any coupling ξ (by Zahn's
theorem), and conformal spin-1/2 and Maxwell fields (by the Brady-Ottewill conformal transformation law). Each claim
below holds **in any Hadamard state respecting self-similarity (Zahn eq. 4); existence of such a state on Evans-Coleman
is not proved.**

**Proposition 0 (cited consequence, not a finding).** Let g be any spacetime with R ≡ 0. Then Zahn's state-independent
tensor V_μν (eq. 19 of [Zahn]) equals B_μν/60 for every ξ, where

B_μν = ∇^ρ∇^σ C_μρνσ + ½ R^ρσ C_μρνσ

is the Bach tensor (Wald conventions, signature −+++).

- *Proof.* Every term of T(1) contains R, so T(1) ≡ 0, and T(2) − T(1)/3 = −4B holds identically [Zahn l.250-256, l.347].
- *Checks.* This was checked symbolically on a generic two-variable metric, on Reissner-Nordström, and on a pp-wave
  (Section 6).
- *Application.* The Evans-Coleman fluid is traceless (p = ρ/3), so R ≡ 0 there, and the leading secular term is
  ⟨T_μν⟩ ⊃ N_c ħ (s/(240π²)) B_μν.
- *Multiplier.* N_c = 1 per scalar, 6 per Dirac fermion, 12 per Maxwell field. These are all positive [BO l.128, eq. 3.4a;
  Duff eq. 30-31].

**Theorem 1 (QB1/QB1b, computer-assisted, HOLD).** On the certified Evans-Coleman solution [MJ-EC], let
b_ab(x) = r⁴B_(ab). Here (ab) are components in the Eulerian orthonormal frame (n, e_r), u is the fluid 4-velocity,
x = ln(r/(−t)) is the self-similar variable in the time normalisation of Koike, Hara and Adachi, and t < 0.

1. **Sign map.** b_nn is certified negative on [−12, x0) and positive on (x0, 12], with
   x0 ∈ [−0.82611434347927570343, −0.82611433975398540497] (the dyadic interval [−443516761, −443516759]/2²⁹).
   The zero x0 is the only zero of b_nn in [−12, 12]: b_nn′ ≥ 1.1654794 on the 9.5e-7-wide cell that contains it.
2. **Sonic point (x = 0).**
   - b_nn = [0.4157372840770100272 ± 1.1e-17];
   - b_nr = [0.1190421959584470894 ± 5.8e-18];
   - b_rr = [−0.1171195084663042197 ± 2e-17];
   - b_T = r⁴B_(θθ) = [0.2664283962716571235 ± 7.3e-18] = (b_nn − b_rr)/2 (traceless);
   - b_uu = r⁴B_ab u^a u^b = [0.4466737796695366321 ± 1.3e-17].
3. **Centre.**
   - B_(nn)(r = 0, t) = c0/t⁴, with c0 = [−173.1004443922995663 ± 2.99e-13];
   - B_(rr) = B_(θθ) = c0/(3t⁴) (certified [−57.70014813076652211 ± 4.2e-13]);
   - B_(nr) = 0 at the centre.
4. **Field independence of the sign.** At every point, the leading secular quantum source has the same sign for every
   field, because the multipliers N_c are positive. Theorem 1 says **where** that common sign is negative (the core) and
   where it is positive (the shell). It does not by itself give a single "direction of the quantum push".

**Lemma 2 (QB2 structure, PROVED, internal referee: CORRECT WITH FIXES, fixes applied in v2; v2 refereed by referee-4: (a) and (c) CORRECT WITH FIXES, (b) and (d) CORRECT; text fixes in AMENDMENT-QB2-v3, sealed 4cf7213207, no certified number changed).**
Write the first-order response to the source Σ_V = (N_c ħ/(30π)) s B = K ħ s B as h = K ε e^{2s}(h1(x) + s·h2(x)), K = N_c/(30π).

- The s·e^{2s} problem for h2 has a solution that is analytic through the sonic point and regular at the centre. It is
  unique modulo the gauge mode G_2.
- There is no s²e^{2s} term.
- h1 depends on the quantum state, and it is ambiguous along h2 and G_2.

**Theorem 3 (QB2/QB2b/QB2c, computer-assisted; headline item HOLD in all three bars: PQ2-8, PQ2b-5, PQ2c-3).** Let U = [−0.30, −0.20].

- The fluid velocity V has exactly one zero x_V ∈ [−0.25278710208307181878 ± 1.14e-20] in U, with V′ > 0 on U.
- Hence A, with (∇r)² = 1/A, has a unique non-degenerate maximum on U at x_V, where A(x_V) ∈ [1.9143464252607319765 ± 2.91e-20].
- The local peak compactness m_U(s) = max_U 2m/r changes, to first order, by δm_U = K ε e^{2s}(sΨ2 + Ψ1) + O(ε²), K = N_c/(30π), with
  **Ψ2 = Ā_2(x_V)/A(x_V) ∈ [3.4616358037141656915 ± 3.02e-12]** (per unit K·ε·s·e^{2s} in the Koike-Hara-Adachi clock, K = N_c/(30π), ε = ℓ_P²/t0²; Ψ2 = δ(2m/r) at x_V).
- **Ψ2 is positive.**
- Ψ2 is invariant under every infinitesimal change of coordinates (Section 7.2). Its sign is the same for every field.
- For the photon field alone (N_c = 12), the coefficient is K·Ψ2 = 0.440749159476 (per unit ε·s·e^{2s}). Per unit N_c the coefficient is Ψ2/(30π) = 0.036729096623; in central proper time u = −ln(−t_c) the u·e^{2u} coefficient per unit N_c is 0.0561653802878 (referee-3, RECONCILE-QB2). The KHA clock is t_KHA = x_s·t_c with t_c the central proper time and x_s = e^{−xc} = [0.80866895631844282988 ± 6.43e-19] (closed form, QB2c PQ2c-4); the sonic point sits at x = 0.
- Not claimed: a conservation width ≤ 1e-10. The matching check e encloses 0 at width 6.66e-10; the sealed 1e-10 width bar FAILED in QB2 and QB2b, and in QB2c it was not run as sealed (a probe showed the tight stage cannot move it; Section 7.3).

**Conjecture 4 (not derived).** Assume the infrared-naturalness conditions (i)-(iii) of Section 7.4, the linear
extrapolation of Tomasevic and Wu [TW1, l.4100], and that the sign of δ(2m/r) at the extremum x* of the trapping
function equals its sign at x_V (Section 7.4, item 6). Then sgn Ψ2 is the direction in which ħ moves the black-hole threshold
for radiation: quantum effects would favour black-hole formation (lower the threshold).

**To our knowledge**, after a prior-art search (17 November 2026; arXiv to 2611.10320, INSPIRE citers to 2611.07007,
Zenodo to 18 November 2026; refreshed 18 November 2026 to arXiv 2611.23200, INSPIRE citers 2611.21356 and Zenodo record 24338153), no earlier work:

- evaluates the one-loop state-independent source on the Evans-Coleman background, in any named stress-tensor model;
- gives its forced rate-2 response or a sign derived from it;
- certifies any semiclassical quantity on a true (not toy) critical-collapse solution.

The structure of the quantum correction is due to Brady and Ottewill (1998) and Zahn (2025). Closed-form one-loop
results on toy self-similar spacetimes are due to Tomasevic and Wu (2025-26). Section 9 gives the credit in detail.

## 2. Background and conventions

- **Background.** We use the critical solution certified in [MJ-EC]:
  - sonic velocity V0 ∈ [0.1124394013880983362 ± 3.31e-20];
  - growth rate κ* ∈ [2.810552549 ± 3.83e-10];
  - gauge rate κ_g ∈ [0.3556992037 ± 7.41e-11].
- **Variables and metric.** We use the self-similar variables of Koike, Hara and Adachi (KHA): N, ω = 4πr²Aρ, V, A.
  The metric is ds² = r²{A[(1 − N²)ds² − 2ds dx + dx²] + dΩ²}, with s = −ln(−t) and x = ln r + s.
- **Inputs.** The inward enclosure is a Lohner tube with certified sonic and centre series. The outward record runs to
  x = 12, and its sonic-manifold box has V0 width 6.94e-29 (the certified box, `b1_C.log`). The same box
  as serialised into the outward record has half-width 3.80e-29, so width 7.58e-29 (`b2_C.log`); both describe one box.
- **Conventions.** Curvature conventions are Wald's, with signature −+++, the same as Zahn [Zahn l.227].
  - Perturbations grow like e^{κs}, the convention of [BO], [Zahn] and [TW1].
  - The quantum rate is κ_q = 2 [BO l.230; Zahn eq. 6].
  - For radiation κ* = 2.81 > 2. This is Brady-Ottewill's case (ii), where "we can say less" [BO l.267].
- **The Bach profiles** are b_ab(x) = r⁴B_(ab). They are functions of x alone, because B has dimension length⁻⁴ and
  conformal weight −2.

## 3. The source in every field model (QB1 design)

- **Scalar, any ξ.** Zahn's eq. 6 for a CSS spacetime is ⟨T_μν⟩(u, y) = e^{2u}[T_μν(y) + (u/4π²)V_μν(y)], with
  V = (1/720 − (6ξ−1)²/288)T(1) − (1/240)T(2) [Zahn eq. 19, l.329-343].
  - The extracted text loses the π². We fix the normalisation 1/(4π²) in three ways: from the parametrix normalisation,
    by exact reproduction of Zahn's Hayward values (eq. 33), and by agreement with Brady-Ottewill's log coefficient
    (below).
  - On R ≡ 0, V = B/60 for all ξ (Proposition 0).
- **Spin 1/2 and Maxwell.** Brady and Ottewill's exact conformal transformation law [BO eq. 3.1, 3.8-3.9, l.100-106,
  l.163-170] gives a log term 8α·u·B for each conformally invariant field.
  - α is the Weyl² anomaly coefficient. One conformal scalar has 8α = 1/(240π²) [BO l.128], which is exactly
    (1/4π²)(1/60).
  - The field ratios are 1 : 6 : 12 (scalar : Dirac : Maxwell). This agrees with Duff's b coefficient [Duff eq. 30-31],
    whose contributions "are all positive".
  - This route needs conformal fields in a state invariant under the rescaled time translation. It uses the exact law,
    not the Page state approximation criticised in [TW1, App. D, l.6315-6320].
- **Sign conventions.** Brady and Ottewill's sign conventions are not in the text we extracted. The absolute sign
  agreement between [BO] and [Zahn] is therefore UNTESTED. This does not affect the claim that the multipliers are
  positive: inside [BO] every field carries 8α_s B with α_s > 0. The sign of B itself comes from Zahn's conventions and
  our own code (Section 6).
- **Renormalisation.** On R ≡ 0 the finite-renormalisation freedom aT(1) + bT(2) collapses onto multiples of V. So does
  a shift of the origin of u [BO l.167]. The s·e^{2s} coefficient is therefore unambiguous. Brady and Ottewill already
  note this for conformal fields; we record it as known.

## 4. Certified evaluation (QB1, QB1b)

### 4.1 Two routes, two charts
- **Route A (on-shell).** V = (1/120)□R_μν + (1/60)R^ρσC_μρνσ + (1/240)g_μν R_ρσR^ρσ − (1/60)R_μ^ρR_ρν with
  R_μν = 8πT_μν. It needs jets of (N, ω, V) up to order 2, obtained from the ODE y′ = f(y).
- **Route B (off-shell).** Zahn's full eq. 19, with ξ symbolic, from metric jets up to order 4. No field equation is
  used.
- **Charts.** A polar-areal chart (t, r), and a Cartesian chart (t, X, Y, Z) that is regular at the centre.
- **Comparison.** The two routes were compared in **28 comparisons at 19 distinct x** (9 values of x appear in both
  charts). All overlap, and the Einstein residual contains 0 at every point. Both routes use the same ODE-derived jets,
  so their agreement tests the algebra, not independent physics.

### 4.2 Sonic point and centre
- **Sonic point.** The values come from exact Taylor coefficients of the certified sonic manifold. The widths are at
  most 3.9e-17, and they overlap the inward tube's sonic manifold.
- **Centre.** The values come from the stored centre-manifold coefficients, in the Cartesian chart at r = 0.
  - The w = z² series has no odd part, as regularity requires.
  - The centre limit is approached from x = −12, where b_nn e^{−4x} = [−173.10044427 ± 0.0014], and from
    x = −9, where it is [−173.1003932 ± 3.4e-6].

### 4.3 The sign map (QB1b)
- **Grid points (QB1).** A first instrument certified signs at 691 grid points (278 negative, 413 positive), with one
  change between x = −0.83 and −0.82.
  - Covering whole cells with plain interval evaluation failed: the evaluation of b_nn on an x-cell overestimates by
    factors of 1e4 to 1e20.
  - The sealed item asked for a covering map, so it is recorded FAILED (Section 8).
- **Cells (QB1b).** The re-test evaluates b_nn as a Taylor model in ε = x − x_m on each cell, with degree P = 5 or 7 and
  a rigorous remainder.
  - It uses a direct chart (s, r) in which the cell centre never enters the jets.
  - In the deep centre, x ≤ −11/2, it uses the w = z² series with a Cauchy tail. There the Rouché condition
    |e^{2X} − 1| < 3/4 is checked on a complex box.
  - Nine regions cover [−12, 12] with no gap:
    - 2261 certified cells (370 negative, 1891 positive);
    - widths from 9.5e-7 to 0.5;
    - total compute about 1.9 hours of single-core time, in parallel runs.
  - Of the 2049 unresolved sub-cells, 2048 lie inside other certified cells. The remaining one is the hull of the zero,
    where b_nn′ > 0 is certified piecewise, with lower bound 1.1654794.
  - Bisection with degenerate Taylor models, on the same code path, gives the zero enclosure of width 3.7e-9. Both ends
    have certified signs: [−3e-9 ± 8.9e-10] and [2e-9 ± 7.4e-10].
- **Seam.** At x = −9/2 the centre series and the tube agree in every component, in both the direct chart and the
  first instrument's polar chart. For example, centre b_nn = [−2.6300024457005e-6 ± 1.4e-11] and tube
  b_nn = [−2.6300024457005e-6 ± 2.3e-9].

### 4.4 The bridge to QB2
- The κ = 2 response is unique modulo gauge iff 2 is not a zero of the deflated characteristic function D_new.
- This already follows from the certified count of [MJ-EC]. That count has one zero in the half-disc, at κ*, and 2 is
  not in κ*'s enclosure.
- QB1b certifies the value directly at the matching point x_m = −4.6:
  - D_new(2) = [−2.4126270857981028564 ± 1.13e-9];
  - D_HKA(2) = [−19.942295267987575083 ± 1.12e-8], after exact gauge deflation.

## 5. Structure of the forced response (QB2 lemmas)

- **Setting.** Σ = 8π⟨T⟩_ren is added to Einstein's equations; the fluid itself is unchanged. The state-independent part
  is Σ_V = (N_c ħ/(30π)) s B.
- **Frame to coordinates.** Write β_ab = e^{−2x}b_ab and drop K = N_c/(30π). Then
  - σ_ss = A(N²β_nn − 2Nβ_nr + β_rr);
  - σ_sx = A(Nβ_nr − β_rr);
  - σ_xx = Aβ_rr;
  - σ_TT = β_T.
- **The four-row source.** In the pencil the count instrument transports (KHA form, Q Y′ = (P0 + κP1)Y), the source
  enters all four rows, the two fluid rows included. The momentum constraint carries σ_I = A(−β_nn + Nβ_nr).

The four lemmas (PROOF-QB2-LEMMAS-v2, refereed in v1; v2 refereed by referee-4: (a) and (c) CORRECT WITH FIXES, (b) and (d) CORRECT; text fixes in AMENDMENT-QB2-v3, sealed 4cf7213207, no certified number changed):

- **(a) Constraint propagation.** Along the KHA pencil, J := I(Y) − σ_I obeys J′ = −A·J − A·D_t(σ), where D_t is the
  divergence of Σ along −t∂_t at fixed r.
  - For the s-coefficient problem D_t = 0, because V is "automatically conserved" [Zahn l.348]. So the constraint stays
    solved.
  - The θθ equation holds iff the radial divergence also vanishes.
- **(b) Sonic step.** Every analytic forced solution with J ≡ 0 lies in an affine 2-plane. The constants:
  - the singular exponent satisfies ρ(2) ≤ −1.9085, so no resonance occurs;
  - Δ2(2) ≠ 0 follows from d < 0 alone, through the identity κ_pole = d (proved; found by the referee, POST HOC).
- **(c) Centre step.** The source rows have orders (v², v², v², v³) and σ_I = O(v²), with or without isotropy
  (exact symbolic computation). So the bounded forced solutions with Ā_p → 0 form an affine 2-plane starting at v².
- **(d) Dimension count.** 2 + 2 − 1 = 3 together with D_HKA(2) ≠ 0 gives a forced solution that is unique modulo
  G_2.
- **Secular structure.** There is no s²e^{2s} term, because G_κ exists at every κ (it is gauge, not a resonance).
  - The e^{2s} level is state-dependent, and its state-independent part is not conserved on its own
    (∇(uV) = (∇u)V [Zahn l.348]).
  - Brady and Ottewill's ansatz (eq. 3.12) lacks the s·e^{2s} term; Zahn notes this too [Zahn l.152-153].

## 6. Controls and nulls

All controls use the same metric-in, tensor-out module as the main run.

| control | expected | result |
|---|---|---|
| C0 Hayward, ξ = 0 [Zahn eq. 33] | V_ττ = 1/80, V_uu = 17/240, V_θθ = −7/480 | PASS (point and Taylor-model versions) |
| sign-flip null (flip the R^ρσR_μρνσ term) | C0 must fail | REFUSES (−1/32 in place of −7/480) |
| C3 Reissner-Nordström M = 1, Q = 1/2 (R = 0, Bach ≠ 0) | closed forms, the same for ξ = 0 and ξ = 1/6; Bach = 60V; route A = V | PASS |
| C2 Minkowski self-similar patch | V ≡ 0 | PASS (zero floor) |
| C6 flat Friedmann k = 1/2 (R ≠ 0) | V_ab u^a u^b = 2(6ξ−1)²/(243τ⁴): nonzero at ξ = 0, zero at ξ = 1/6 | PASS (ξ-dependence detected where it must appear) |
| PQ6b off-shell null: ω × (1 + 1e-12) | Einstein residual must exclude 0; routes A and B must disagree | REFUSES as required: residual excludes 0 at 24 of 28 points, routes disjoint at 8 of 28; certified background contains 0 everywhere |
| QB2 C1: manufactured φ(x)·G_2 | returned exactly | PASS. Matching coefficients b = [−4.7963240638274756700e-20 ± 8.83e-14] ∋ 0, c = [0.49999999999999999879 ± 2.19e-12] ∋ 1/2, d = [0.99999999999999999675 ± 6.15e-12] ∋ 1, e ∋ 0; the solver and the exact φG_2 overlap componentwise at x = −6, −4.6, −2, −0.25 and −0.01 (QB2b PQ2b-3). The first QB2 run printed the wrong expected signs c = −1/2, d = −1 in a sanity line; corrected before the record run (build-3 RUN.md) |
| QB2 C2: resonance null at κ* | uniqueness refused | REFUSES as required: D_HKA(κ*) = [0 ± 4.55e-7] and D_new(κ*) = [0 ± 1.52e-7] contain 0; at κ = 2, D_HKA(2) = [−19.942295267987582258 ± 1.07e-11] excludes 0 (QB2b PQ2b-4) |
| QB2 C6: route B (variation of constants + adjoint) | overlaps route A for Ψ2 | PASS: route B Ψ2 = [3.4616358037141656727 ± 1.03e-11] overlaps route A [3.4616358037141656915 ± 3.02e-12] (QB2b PQ2b-5; replayed in QB2c) |
| QB2 C8 / C8′: manufactured conserved Σ / the same with β_nn × (1 + 1e-6) | e ∋ 0 / refuses | C8: e = [−1.8804402479588588903e-15 ± 4.68e-9] ∋ 0 and its θθ residual [0 ± 3.06e-10] ∋ 0 at x = −0.01; C8′: e = [0.0012111376937769522267 ± 4.68e-9] excludes 0 and encloses the prediction [0.0012111376937788326689 ± 1.11e-15] (QB2b PQ2b-4) |
| QB2 N-a: non-conserved source, row 0 × (1 + 1e-6) | e excludes 0 and encloses the predicted value | REFUSES, and encloses the corrected prediction e = (J_c(x_m) − J_s(x_m))/c_V(x_m) (AMENDMENT-QB2-v3 F1): per unit δ the run gives [−302.35960932879925242 ± 2.06e-10] against the prediction [−302.35960932879926987 ± 2.79e-11]. The sealed design formula |J_s|/|c_V| (302.32942153207741039) assumed J_c = 0 and was wrong by 1e-4 relative, so this null was scored FAILED in QB2 and HOLD (corrected) in QB2b |
| QB2 N-c: wrong background, ω × (1 + 1e-6) | refuses | REFUSES: e = [−0.00011949776035180548206 ± 3.56e-10] excludes 0 (QB2b PQ2b-4) |
| QB2c N-b: θθ null, b_T × (1 + 1e-6) | excludes 0 at x = −6, −4.6, −2, −0.25, −0.01 | REFUSES at −6, −2, −0.25, −0.01; at −4.6 only on the centre-side representation, chosen POST HOC (the sealed sonic-side one does not exclude 0) (QB2c PQ2c-2) |
| QB2c replication row (PQ2c-4): blind seat, own formulation, double precision, no enclosure | 30π x_s² Ψ_blind / A(x_V) agrees with Ψ2 to relative ≤ 1e-9 | PASS: relative differences from −6.9153490713739357298e-12 to 8.5734474244780018923e-13 over six runs; 2.0992124976238618617e-10 against the blind 10-digit headline 0.1075199950. The conversion was found after both numbers existed (POST HOC, referee-3); a replication, not a certificate |

- **Independent recomputation (internal referee).** A separate code, written from scratch in Wald conventions with its
  own series reversion, used only the certified manifold coefficients as input. It reproduced:
  - sonic b_uu = 0.4466737796695366321323105, identical to 25 digits;
  - c0 = −173.1004443922995663248;
  - Einstein residuals ≤ 2.1e-32, and isotropy b_rr/b_nn → 1/3.
- **Post hoc support.** The same referee code verified Proposition 0 on a non-static pp-wave. There V_uu =
  −(1/240)∇⁴H matches a hand prediction. This check was POST HOC.

## 7. The observable and its reading (QB2)

### 7.1 Trapping function
- (∇r)² = 1 − 2m/r = 1/A exactly, so δ(2m/r) = Ā_p/A.
- The apparent horizon is where (∇r)² = 0 [BO eq. 3.14, l.235-241; TW1 eq. 6.8, l.4096].
- On the background, A′ = −(8/3)ωNV/(1 − V²), so A is stationary exactly where the fluid velocity V vanishes.

### 7.2 Ψ2 and its invariance
- **Definition.** Ψ2 is defined through the local maximum on U = [−0.30, −0.20], via Danskin's theorem on a compact
  set. No global maximum over all x is needed, so no bound beyond x = 12 is needed either.
- **Invariance.** On the background, 2m/r = 1 − 1/A depends only on x, and A′ = 0 on the whole line x = x_V. So
  ∇(2m/r) vanishes there as a spacetime covector, and δ(2m/r)(x_V) is invariant under every infinitesimal
  diffeomorphism. The gauge mode G_2 is a special case.

### 7.3 Result
- x_V ∈ [−0.25278710208307181878 ± 1.14e-20] and A(x_V) ∈ [1.9143464252607319765 ± 2.91e-20] (QB2b N7, N8).
- The conservation check at matching gives e = [2.2796218038245719129e-14 ± 3.33e-10] (QB2b, replayed unchanged in QB2c). It encloses 0, as it must, but its width 6.66e-10 is above the sealed bar 1e-10, so that item FAILED as sealed in QB2 (width 1.98e-5) and QB2b. In QB2c it was not run as sealed: a probe showed that the tight homogeneous stage cannot move the width, and the QB2b value is replayed (FAILED). We do not claim a conservation width ≤ 1e-10. Our estimate of the floor (a probe estimate, not a certificate) is the step balls of the pinned inward background tube (Taylor order KT = 36) times the state sensitivity of the transported HKA system: a tighter homogeneous stage leaves the order-0 step-matrix radii almost unchanged (ratio tight/plain: median 0.896, at most 0.995 on the deep substeps; build-5 `deep_ratio.json`). An earlier note that blamed interval dependency was corrected by that probe. The nulls keep their teeth at this width: N-a gives e = [−0.00030235960930600302070 ± 3.33e-10] and N-c gives [−0.00011949776035180548206 ± 3.56e-10], both excluding 0. The source conservation along the
  transport: D_t ≡ 0 and (div Σ)_x ≡ 0 hold identically for the Bach tensor of every metric of this form (exact sympy, `bachsym.log`). On every inward step the certified D_t series coefficients contain 0, with |D_t| ≤ 3.3606426921349536e-13 (QB2b PQ2b-1; that sub-bar passes). The linearised θθ residual on the real source encloses 0 at all five sealed points (QB2c PQ2c-2 HOLD): [0 ± 4.81e-11] at x = −6, [0 ± 8.00e-12] at x = −4.6 (centre-side representation), [3.19839640883e-18 ± 1.14e-10] at x = −2, [4.33680868994e-19 ± 1.60e-10] at x = −0.25 and [0 ± 1.10e-9] at x = −0.01. The N-b null (b_T × (1 + 1e-6)) excludes 0 at x = −6, −2, −0.25 and −0.01. At x = −4.6 the sealed, transported sonic-side representation (|Y| ≈ 600) does not exclude 0 (N-b [5.81641278197e-9 ± 1.67e-8]); the centre-side representation of the same solution, chosen after that (POST HOC; build-5 RUN.md), excludes it ([5.81641272299e-9 ± 8.00e-12]). The item was scored HOLD on that basis. The outward range (0, 12] is not covered (domain cut, below).
- The Bach source rebuilt by the QB2 instrument agrees with the QB1 instrument: all overlap (`code/quantum/seats/build-3/n0check.json`): at the sonic point b_nn [0.415737284 ± 2.32e-10] against QB1 [0.4157372840770100 ± 3.82e-17] (ratio [1.0000000 ± 3.73e-10]), and likewise b_nr, b_rr, b_T; the centre c0 [−173.100444392300 ± 5.67e-13] against [−173.100444392300 ± 7.33e-13]; all four components at 8 tube points from x = −0.2489 to −4.5.
- **Ψ2 ∈ [3.4616358037141656915 ± 3.02e-12] (per unit K·ε·s·e^{2s} in the Koike-Hara-Adachi clock, K = N_c/(30π), ε = ℓ_P²/t0²; Ψ2 = δ(2m/r) at x_V). Sign: positive.**
- **Reading.** The s-secular (log-growing) part of the first-order quantum push raises (toward trapping) the local peak
  compactness of the critical ball of light. The sign is the same for every free massless field of spin ≤ 1, in any
  Hadamard state respecting self-similarity.
- x_V lies in the region where b_nn > 0 (Theorem 1, since x_V > x0). Ψ2 nonetheless integrates the source over the whole
  inward solution, core and shell. So its sign is not implied by the sign of b_nn at x_V, and it has to be computed.
- **Disclosed deviations.** (i) Domain cut: the instrument consumes only x ≤ 0 data (centre regularity, sonic analyticity, the gauge N̄_p(0) = 0, matching at x_m = −4.6, readout at x_V), a written lemma sealed as QB2b PQ2b-3 (`code/quantum/seats/build-4/DOMAIN-LEMMA.md`). So no h2 profile for x > 0 was built, and the QB2 C1 points x = 0.5, 3 were not run. (ii) The sealed numeric centre test (coefficients within 1e-30) could not be reached by the stored record (3.5e-16) and FAILED in QB2; it was replaced in QB2b by an exact symbolic proof of the orders and the closed form c0 = −(2/3)W0(8W0 − 3)e^{−4xc}, certified [−173.1004443922995663248426 ± 8.79e-16] (PQ2b-2 HOLD; re-derived independently by referee-4). The orders and the 1/3 isotropy are kinematic (smoothness plus tr B = 0 and ∇·B = 0); only the W0-dependence is content. (iii) The QB2b inward tube is rebuilt from the BG2 sonic box (V0 width 6.94e-29), same 720-step partition, overlapping the E5 tube at every step (`tube_compare.json`).

### 7.4 Is the sign decidable from the geometric part alone? Not in general
1. **At fixed s.** At first order and fixed s, inside the window ε s e^{2s} ≪ 1, δm_U = K ε e^{2s}(sΨ2 + Ψ1), K = N_c/(30π). The
   s-coefficient is Ψ2 alone, under Zahn's state assumption.
2. **At the threshold.** sΨ2 + Ψ1 = Ψ2 ln(L*/(−t)), where L* is set by the state together with the renormalised C²
   coupling [BO eq. 3.10, l.220]. The trapping time satisfies −t* ~ ℓ_P (ln)^{1/2}.
3. **Infrared naturalness, our stated assumption.**
   - (i) the renormalised coupling is O(1) at the initial-data scale t0;
   - (ii) the self-similar Hadamard state is normalised at t0;
   - (iii) |Ψ1^state/Ψ2| ≪ ln(t0/ℓ_P).
   - Item (iii) matters because Ψ2 could be small through the cancellation between core and shell.
   - Zahn writes that "in a first approximation the state-dependent part could be ignored" [Zahn l.683-684].
   - Tomasevic and Wu caution that no theorem bounds the state-dependent part by the state-independent one, which can
     be of "comparable or even larger order" [TW1 l.6376].
4. **Two further limits.**
   - (L1) At the trapping time the perturbation of (∇r)² is O(1), so using the linear response there is the linear
     extrapolation of Tomasevic and Wu [TW1, l.4100], a model step.
   - (L2) The curvature there is only a logarithm below the Planck scale.
5. **If the coupling is natural at the Planck scale**, then L* ~ ℓ_P and Ψ2 does not decide the sign.
6. **The threshold shift depends on the trapping function at its extremum x*, not at x_V.** Equating sgn Ψ2 with the
   threshold direction is therefore Conjecture 4, not a derived result.
7. **Context.** In the 2D s-wave models of Tomasevic and Wu there is no log term at all [Zahn l.682].
   - Brady and Ottewill assumed that "quantum matter tends to oppose black-hole formation" [BO l.277-278].
   - Tomasevic and Wu found that quantum corrections "lower the critical point" in toy scalar models [TW2 l.275-278].
   - Wu found, in a large-D model of a discretely self-similar (DSS) scalar collapse, that the direct quantum
     contribution "can favor or oppose trapping" in a large-D DSS scalar model [Wu2 l.200].
   - Our Ψ2 is positive. Under Conjecture 4 it would favour black-hole formation (lower the threshold). This direction matches the toy-model result of Tomasevic and Wu and is opposite to the assumption of Brady and Ottewill. Both comparisons are across different models; neither is a test of Ψ2.

## 8. Sealed tests and verification

Every test below was sealed (prediction and pass bar hashed) before its run, and scored as written. Failed items are
listed as failed.

| test | seal | item | what | verdict |
|---|---|---|---|---|
| QB1 | 54caa4012c | PQ1, PQ3, PQ4, PQ5 | certified source at the sonic point and centre, 28 route comparisons, positive controls (PQ1 sonic point not Bach-flat; PQ3 centre law c0; PQ4 routes A and B at ≥ 18 points; PQ5 controls) | HOLD |
| QB1 | 54caa4012c | PQ2 | covering sign map on [−12, 12] | FAILED (grid points only: cell covering overestimated by factors of 1e4 to 1e20) |
| QB1 | 54caa4012c | PQ6 | off-solution content null (V0 + 1e-15) | FAILED (ill-posed: the shifted background still solves the ODE locally, so the Einstein residual cannot refuse) |
| QB1 | 54caa4012c | PQ7 | D_new(2) bridge | FAILED (computed at the wrong matching point x_m) |
| QB1b | 4d3b12ef39 | PQ2b | covering sign map with Taylor models, unique zero | HOLD |
| QB1b | 4d3b12ef39 | PQ6b | off-shell null ω × (1 + 1e-12) | HOLD (refuses at 24 of 28 points) |
| QB1b | 4d3b12ef39 | PQ7b | D_new(2), D_HKA(2) at x_m = −4.6 | HOLD |
| QB1b | 4d3b12ef39 | PQ8 | centre-tube seam at x = −9/2 | HOLD |
| QB2 | f2bcfadb3e | PQ2-1 | certified source Taylor models, overlap with QB1b | HOLD |
| QB2 | f2bcfadb3e | PQ2-2 | four HKA source rows, exact | HOLD |
| QB2 | f2bcfadb3e | PQ2-3 | conservation: D_t and θθ on [−12, 12]; e ∋ 0 with width ≤ 1e-10 | FAILED (e ∋ 0 but width 1.98e-5; (0, 12] not covered) |
| QB2 | f2bcfadb3e | PQ2-4 | centre coefficients within 1e-30 | FAILED (3.5e-16, the precision of the stored record) |
| QB2 | f2bcfadb3e | PQ2-5 | controls C1, C2, C6, C8, C8′ | FAILED (partial: C1 at x = 0.5, 3 not run; C8 θθ residual not computed; the rest passed) |
| QB2 | f2bcfadb3e | PQ2-6 | nulls N-a, N-c | FAILED (partial: both refuse, but the printed N-a prediction formula was wrong by 1e-4 relative) |
| QB2 | f2bcfadb3e | PQ2-7 | x_V unique in U, V′ > 0 | HOLD |
| QB2 | f2bcfadb3e | PQ2-8 | headline Ψ2, excludes 0, relative width ≤ 1e-3 | HOLD ([3.4616358037141459203 ± 1.15e-7]) |
| QB2b | 373f306e9f | PQ2b-1 | conservation with the BG2 inward tube: e width ≤ 1e-10; abs D_t ≤ 1e-10 | FAILED (e ∋ 0, width 6.66e-10; abs D_t ≤ 3.4e-13 passes) |
| QB2b | 373f306e9f | PQ2b-2 | exact centre orders + certified order-4 coefficient | HOLD |
| QB2b | 373f306e9f | PQ2b-3 | domain lemma (x ≤ 0 only) + C1 at five points | HOLD |
| QB2b | 373f306e9f | PQ2b-4 | C8 θθ, C8′, N-a (corrected), N-c, C2 | HOLD |
| QB2b | 373f306e9f | PQ2b-5 | headline replay, route B overlaps | HOLD ([3.4616358037141656915 ± 3.02e-12]) |
| QB2c | 7aab7efe5a | PQ2c-1 | conservation width ≤ 1e-10 with a tighter homogeneous stage | FAILED (not run as sealed: a probe showed the tight stage cannot move the width; the QB2b value is replayed) |
| QB2c | 7aab7efe5a | PQ2c-2 | θθ residual at five points; N-b refuses | HOLD (N-b at x = −4.6 on the centre-side representation, POST HOC) |
| QB2c | 7aab7efe5a | PQ2c-3 | headline replay | HOLD ([3.4616358037141656915 ± 3.02e-12]) |
| QB2c | 7aab7efe5a | PQ2c-4 | replication row, relative ≤ 1e-9 | HOLD (≤ 6.9e-12 over six runs; POST HOC conversion) |

Ledger for this work: QB1 FAILED by item (4 HOLD, 3 FAILED); QB1b HOLD; QB2 FAILED by item (4 HOLD, 4 FAILED); QB2b FAILED by item (4 HOLD, 1 FAILED); QB2c FAILED by item (3 HOLD, 1 FAILED). The headline Ψ2 held in all three QB2 bars. The conservation width bar failed in QB2 and QB2b and was not run as sealed in QB2c (FAILED). The QB1b settle row carries an automatic pin warning ("STALE: result_qb1b.json names instrument hash 582e24236c8f…"): that hash is `engine/cert/lin_stages.py`, recorded by the instrument as provenance, unchanged and pinned under b1b0ee1cde.

Inputs inherited from earlier sealed work [MJ-EC]:

- existence E5 (d0af75803f, HOLD);
- mode count M (b1b0ee1cde, HOLD);
- the sharper outward record BG2 (ed61d09268, HOLD). BG2 is not in any published record: its sonic box, outward record (`outward_C.pkl`), sealed bar and regeneration scripts are shipped here (`code/sideways/build-5/`, `results/seals/BAR-SF-BG2.md`).

**Internal referee seats.** "Seats" are separate internal working sessions of the lab, each with its own files. Referee
seats are internal checks, not external peer review.

| seat | object | findings | what was done |
|---|---|---|---|
| referee-1 | QB1 design and results | C1-C5 all CORRECT WITH FIXES. Independent code reproduces sonic b_uu (25 digits) and c0 (22 digits). Wording fixes: "28 comparisons at 19 distinct x"; spin 1/2 and Maxwell rest on BO's conformal law; V = B/60 is a cited consequence; no "exactly one zero" until PQ2b; no single "sign of the push" from QB1 | PQ2b, PQ6b and PQ7b sealed as QB1b (HOLD); wording rules adopted throughout this paper |
| referee-2 | QB2 lemmas and design (v1) | lemmas CORRECT WITH FIXES. BUG: the source also enters the fluid rows of the transported pencil. Observable invariant under every gauge (stronger than claimed). Threshold wording overclaimed. Control C1 outside the centre lemma; null N-c could not refuse | four-row source; C1 moved inside the lemma; new control C8; N-c at 1e-6; threshold labelled CONJECTURE; v2 written before the build |
| referee-3 | QB2 result against the blind seat | Build Ψ2 CORRECT; blind Ψ CORRECT in its own convention; the two convert exactly by 30π x_s²/A(x_V) (clock, source unit, readout), relative agreement 2.4e-12 (POST HOC). Text fixes: K = N_c/(30π) missing in Theorem 3 and design node N8; the unit must name the KHA clock | K added (AMENDMENT-QB2-v3 F2); unit statement in Theorem 3; replication row sealed as QB2c PQ2c-4 (HOLD) |
| referee-4 | QB2 lemmas v2, domain lemma, centre formulas | Lemmas sound; centre formulas re-derived from scratch (own metric series and own Bach tensor), c0 closed form exact. WRONG: the N-a expected value in the v2 text (J_c ≠ 0). Fixes: K in the design page; ρ(2) ≤ −1.9085 is PROVED, not certified; the domain claim stated operationally; the centre orders are kinematic | AMENDMENT-QB2-v3 (sealed 4cf7213207) carries every fix; no certified number changed; QB2c re-tested the open width item |
| referee-5 (two rounds) | this text (paper, README, CLAIMS, SCORECARD, CITATION, figures) | Science CORRECT: Ψ2, the photon and per-N_c values re-derived in arb; the scorecard matches the ledger; the figures regenerate byte for byte; the literature comparison is fair. Round 1: 3 blocking text fixes (printed balls not all outward; QB2c PQ2c-1 described as measured although it was not run as sealed; magic-trick scope) and 11 smaller ones. Round 2: CORRECT WITH FIXES (3 sentence-final balls missed by the printing tool; minor wording) | every printed ball is now generated by tools/fill_balls.py (113 balls, 0 not verbatim-or-outward); PQ2c-1 worded as not run as sealed; magic trick scoped to x in [-12, 12] and the local peak; all smaller fixes applied |

**Clean-room replay.** On 18 November 2026 the repository was copied into fresh folders and replayed with the documented commands only (tools/cleanroom.sh). Quick tier: all 16 steps IDENTICAL to the shipped files (canonical content), ALL CHECKS PASS, 20.1 min. BG2 tier: the sonic box (44.0 min) and the outward record (10.1 min) IDENTICAL; the gate step first differed only in two wall-clock fields (time_s, path_info.setup_s), the comparator was fixed to ignore them (PATCHES.md), and the gate step re-run from the shipped outward record is IDENTICAL, ALL CHECKS PASS. The full tier (about 7 hours: the transported QB2b record) was not replayed in a clean room; the quick tier replays the matching, controls, nulls and readout on the shipped transported record.

## 9. Fairness notes, credit and corrections

### 9.1 Credit
- **Brady and Ottewill (1998)** [BO] found the structure of the one-loop correction for continuously self-similar
  collapse, perfect fluids included:
  - the quantum rate κ_q = 2;
  - the u-linear Bach term with coefficient 8α;
  - the renormalisation freedom along the Bach tensor (l.167).
  They also stated plainly that the sign was an assumption that "can not be addressed without a complete calculation"
  (l.278), and that for κ_c > κ_q "we can say less" (l.267). Everything in Sections 3 and 5 is built on their
  structure.
- **Zahn (2025)** [Zahn] proved the Hadamard scaling theorem behind V for a scalar with any coupling. He identified V as
  state-independent, scheme-independent and conserved, and noted that, at conformal coupling, it reduces to a multiple of the
  Bach tensor (l.347). Proposition 0 is a one-line consequence of his work.
- **Tomasevic and Wu (2025-26)** [TW1, TW2] solved the one-loop problem in closed form on toy self-similar spacetimes
  (Garfinkle, Roberts). They gave the sign rule sgn(Δp) = sgn(f_q/f_c) [TW1 l.4196-4218], and named radiation-fluid
  collapse "the important next step", to be done by simulation [TW2 l.417]. Ours is a certified linear-order
  coefficient on a proved background. It is a different kind of result from their planned simulation, and it does not
  replace it.
- **Wu (2026)** [Wu1, Wu2] restated the radiation case as open [Wu1 l.3103], and found that the direct quantum
  contribution "can favor or oppose trapping" in a large-D DSS scalar model [Wu2 l.200].
- **Ayal and Piran (1997)** [AP] did an early numerical semiclassical study of scalar-field critical collapse. They found that evaporation
  acts like a decrease of the initial amplitude (abstract), and pointed out honestly that their 2D-derived tensor
  "diverges at the origin" (l.145-151). Our local tensor is finite at the centre (Theorem 1, item 3).
- **Berczi, Saffin and Zhou (2021)** [BSZ] simulated semiclassical Choptuik (scalar-field) critical collapse numerically
  and found that Choptuik scaling survives (abstract). Theirs is the closest numerical semiclassical study of critical
  collapse; it is a simulation of a scalar field, not a certified quantity on a fluid critical solution.
- **Chiba and Siino (1997)** [CS] state that semiclassical equations "do not admit continuously self-similar
  solutions". We know this work **only from its abstract**; we did not open the full text. We claim no semiclassical
  self-similar solution: the source breaks self-similarity through its s·e^{2s} growth, consistent with their
  statement. Two early readers describe their model as two-dimensional: "write down their own 2-dimensional toy
  model" [G97 l.1980-1981], and Tomasevic and Wu speak of "the same two-dimensional
  model" [TW1 l.869]. Their centre divergence therefore concerns a 2D-derived tensor, not the 4D local Bach term.
- **Centre structure is not new.** Smoothness at a regular centre forces the radial flux to vanish linearly in r and
  radial minus angular pressure to vanish there [Wu2 l.2226], and self-similarity gives the 1/t⁴ law. What is new here
  is the evaluation on Evans-Coleman: the finite centre limit of B with the exact coefficient
  c0 = −(2/3)W0(8W0 − 3)e^{−4xc} (proved symbolically, re-derived by an independent referee) and its certified value.
- **Gundlach, Hilditch and Martín-García** [GHM] summarise Brady-Ottewill's perturbative approach in their review of
  critical collapse (l.1541-1549). We used their treatment of dimensionful parameters as relevant perturbations
  (l.458-505) as a guide to how the Planck length enters.
- **Duff** [Duff] gives the standard spin weights of the Weyl² anomaly coefficient that we use for N_c.
- **Evans and Coleman (1994)** [EC] discovered the critical solution itself.

### 9.2 Our own corrections (none of this was published)
- An internal prior-art note wrote Zahn's log coefficient as u/4π. The correct factor is u/(4π²); the π² is lost in
  the text extraction.
- The same note read Brady-Ottewill's coefficient as 8β. It is 8α, the Weyl² coefficient. With β the two sources would
  differ by a factor of 3 and a sign.
- An internal roadmap said QB2 needs "2 between the gauge rate and κ*". The correct condition is "2 is not a zero of
  D_new"; the gauge rate is deflated out.
- QB1 items PQ2, PQ6 and PQ7 FAILED as sealed and were re-tested as QB1b (Section 8).
- The QB2 design v1 omitted the fluid-row source. The internal referee caught this before any build ran, and v2 fixed
  it.
- The v1 null N-c used a perturbation (1e-12) below the instrument width, so it could not refuse. v2 uses 1e-6.
- The internal referee's own first sonic run had a frame-sign bug on its side. It found the bug, fixed it, and logged
  it.

## 10. Limitations (what we do not claim)

- **The state.** Every claim holds in any Hadamard state respecting self-similarity (Zahn eq. 4). The existence of
  such a state on Evans-Coleman is not proved. Zahn proves existence for a Minkowski patch, Hayward and Roberts, and
  sketches others [Zahn l.126, l.649-663]. Tomasevic and Wu write that the proposal "would benefit from further
  clarification and development" [TW1 l.6375-6382].
- **Fields.** The scalar case rests on Zahn's theorem. Spin 1/2 and Maxwell rest on the Brady-Ottewill conformal law,
  for conformal fields only. Nothing is claimed for spin > 1.
- **The radiation fluid is not quantised.** The quantum fields are independent free fields on a classical background.
  Reading N_c = 12 as "the fluid's own photons" is an interpretation. Double counting between the fluid and its quantum
  field is not addressed.
- **Order and window.**
  - All results are first order in ħ, inside ε s e^{2s} ≪ 1.
  - The source breaks exact self-similarity, consistent with Chiba and Siino [CS] and with Zahn (l.153).
  - Higher loops and quantum gravity are not addressed.
- **Region.**
  - Everything is at t < 0. The Cauchy horizon (t > 0), where quantum fluxes are known to diverge in self-similar
    collapse, is outside this work.
  - The sign map covers x ∈ [−12, 12]. The centre limit c0 covers x → −∞.
  - There is no certified bound for x > 12. Ψ2 does not need one, but the global response profile does.
- **What QB1 alone does not say.** b_nn changes sign. So QB1 gives no single "sign of the quantum push", only the field
  independence of the sign at each point.
- **The state-dependent part** e^{2s}T(y), and so Ψ1, is not computed.
- **Threshold.** The direction in which ħ moves the black-hole threshold is Conjecture 4, not a result.
- **2D models** of the Tomasevic-Wu type have no log term [Zahn l.682]. They are a different object, not a control of Ψ2.
- **Route agreement** (28 comparisons at 19 distinct x) tests the algebra, not independent physics.
- **Trusted layer.** We trust ball arithmetic (python-flint arb) with outward rounding. Classical theorems are used
  without proof: Frobenius theory, majorant convergence, Coddington-Levinson for regular singular points, Danskin's
  theorem, the Rouché condition. No step is checked in Lean.
- **Agreement** between our two instruments and the internal referee's code was checked after the fact (POST HOC).
- **Review status.** Four internal hostile referee seats on the science (referee-1 on QB1; referee-2 on the QB2 design v1; referee-3 on the QB2 result against the blind seat; referee-4 on the QB2 lemmas v2 and the centre formulas), every fix applied; release-text referee: referee-5, two rounds, final verdict CORRECT WITH FIXES, all fixes applied. Prior-art refresh to 18 November 2026. These are internal checks, not external peer review.

## 11. Open doors

What this work leaves open. "Open according to" quotes the literature; "UNTESTED" marks an observation of ours that we
did not test.

- **The threshold itself.** The sign rule sgn(Δp) = sgn(f_q/f_c) of Tomasevic and Wu [TW1 l.4196-4218] turns the
  threshold question into one sign at the trapping extremum x*. That sign needs the state-dependent part.
  - *Status:* open according to Tomasevic and Wu, who name radiation-fluid collapse "the important next step" [TW2
    l.417]; Wu [Wu2 l.200] finds, for a large-D DSS scalar model, that the direct contribution "can favor or oppose trapping"; Brady and Ottewill [BO l.278]
    call the sign an assumption.
  - *Start here:* h2 on [−12, 12] from this work, and the extended threshold algebra with the log term (Section 7.4).
  - *Ceiling:* the direction in which ħ moves black-hole formation for radiation, under each named stress-tensor model.
- **A self-similar Hadamard state on Evans-Coleman.** Every claim here is conditional on one. Zahn proves existence on
  three backgrounds and sketches a route for others [Zahn l.126, l.649-663].
  - *Status:* open according to Zahn, who argues persistence "though not rigorously" [Zahn l.131].
  - *Ceiling:* the condition dropped from every theorem in this paper.
- **4D log term against 2D models on the same background.** Zahn notes that 2D results "differ ... in the coefficient
  of e^{2u}" because they have no log term [Zahn l.682], and Tomasevic and Wu note that "dimensional reduction and
  quantization do not commute in general" [TW2 l.409].
  - *Status:* UNTESTED observation from this work. Nobody we found has put both on the same true critical solution.
  - *Start here:* the certified background and the four-row source of Section 5. A 2D s-wave source would replace σ.
  - *Ceiling:* a direct, certified measurement of how much the 2D reduction changes the quantum push on a real critical
    collapse.

## Figures (drawn from shipped data only by `figures/make_figures.py`; the scorecard by `tools/make_scorecard.py`)

1. **`figures/bach-profile.svg`: the quantum fingerprint and its one sign change.**
   - *Data:* `results/result_qb1.json` (`sign_map.rows`: 691 certified point values of b_nn, with the scale tag
     `hat` = b_nn e^{−4x} in the centre region); `results/result_qb1b.json` (`PQ2b.cells`: 2261 cells with sign and
     region; `PQ2b.zero`).
   - *Top panel:* b_nn against x on [−3, 12], linear axis, with points joined. Mark the sonic point x = 0
     (b_nn = 0.4157) and the zero x0 = −0.8261143 with a vertical line.
   - *Middle panel:* the centre approach b_nn e^{−4x} on [−12, −3], tending to c0 = −173.1004 (dashed horizontal line).
   - *Bottom strip:* the 2261 certified cells as coloured bars, negative in one colour and positive in another.
     Annotate "370 negative · 1891 positive · one zero, width 3.7e-9".
   - *Caption must say:* the error bars are far smaller than the line; x > 12 is not certified.
2. **`figures/core-and-shell.svg`: where the quantum source points, at one instant.**
   - *Data:* x0 from `PQ2b.zero`, the sonic point x = 0, x_V from the QB2 result file ([−0.25278710208307181878 ± 1.14e-20]).
   - *Drawing:* a cross-section of the ball at t = −1 (so r = e^x), as concentric circles:
     - the core r < 0.43775 shaded "b_nn < 0";
     - the shell r > 0.43775 shaded "b_nn > 0";
     - the sonic sphere r = 1 (dashed);
     - the sphere of peak compactness r = e^{x_V} ≈ 0.776 (solid), labelled "Ψ2 read here".
   - *Caption must say:* the picture shrinks self-similarly as t → 0⁻; the sign is the same for every field of spin ≤ 1;
     conditional on the self-similar Hadamard state.
3. **`figures/psi2-number-line.svg`**: Ψ2's certified enclosure on a number line with 0 marked, and route A, route B, the first QB2 run and the six converted blind values stacked in a zoom.
   - *Data:* `results/result_qb2b.json` (PQ2b-5), `results/result_qb2.json` (PQ2-8), `results/result_qb2c.json` (PQ2c-4).
4. **`scorecard.svg`:** the HOLD/FAILED chip table of Section 8, failed items included.

## References

Labels: OPENED means the lab read the full text; locators are line numbers in the lab's text copy, named per entry
("copy:"; the copies are not shipped). ABSTRACT means the
abstract or record only. arXiv identifiers are given; journal data only where the opened text shows it.

- [AP] S. Ayal, T. Piran, "Spherical collapse of a mass-less scalar field with semi-classical corrections",
  arXiv:gr-qc/9704027, abstract (l.6), l.145-151. OPENED; copy: gr-qc_9704027.plain.txt.
- [BO] P. R. Brady, A. C. Ottewill, "Quantum corrections to critical phenomena in gravitational collapse",
  arXiv:gr-qc/9804058 (Phys. Rev. D 58, 024006), l.32, l.100-106 (eq. 3.1), l.128
  (eq. 3.4a), l.163-170 (eq. 3.8-3.9), l.167, l.220 (eq. 3.10), l.230, l.235-241 (eq. 3.14), l.260, l.267, l.277-278.
  OPENED; copy: bo.txt.
- [CS] Chiba, Siino, Mod. Phys. Lett. A 12 (1997) 709, DOI 10.1142/S0217732397000741, INSPIRE record 444860.
  ABSTRACT (full text not opened).
- [Duff] M. J. Duff, arXiv:hep-th/9308075, eq. 30-31. OPENED by an internal referee (quoted in REFEREE-QB1.md); no
  lab copy kept, so the release referee did not re-check it.
- [EC] C. R. Evans, J. S. Coleman, arXiv:gr-qc/9402041. OPENED (in our earlier work [MJ-EC]).
- [GHM] C. Gundlach, D. Hilditch, J. M. Martín-García, review of critical phenomena in gravitational collapse,
  arXiv:2507.07636 v1, l.458-505, l.1520-1549. OPENED; locators not re-checked by the
  release referee.
- [BSZ] B. Berczi, P. M. Saffin, S.-Y. Zhou, "Gravitational collapse of quantum fields and Choptuik scaling",
  arXiv:2111.11400 (2021). ABSTRACT (title and authors checked against the arXiv record).
- [MJ-EC] Moki & Julio, A certified Evans-Coleman critical solution: existence, one growing spherical mode, and its
  exponent, Zenodo, DOI 10.5281/zenodo.24303877.
- [TW1] M. Tomašević, C.-H. Wu, "Unveiling horizons in quantum critical collapse", arXiv:2509.03584 v3 (JHEP),
  l.4096 (eq. 6.8), l.4100, l.4196-4218, l.4270, l.6310-6320, l.6375-6382. OPENED; copy: tw2.txt.
- [TW2] M. Tomašević, C.-H. Wu, "Quantum critical collapse abhors a naked singularity", arXiv:2509.03587 v3, l.80,
  l.275-278, l.406, l.409, l.417. OPENED; copy: tw.txt.
- [Wu1] C.-H. Wu, "Quantum fate of the Choptuik naked singularity", arXiv:2606.17135 v2, l.3103. OPENED; copy: 2606.17135.plain.txt.
- [Wu2] C.-H. Wu, "Quantum critical collapse in the large D limit", arXiv:2610.20705 v1, l.200, l.1760, l.2226. OPENED; copy: 2610.20705.plain.txt.
- [G97] C. Gundlach, "Critical phenomena in gravitational collapse" (review), arXiv:gr-qc/9712084, l.1980-1983. OPENED; copy: gr-qc_9712084.txt (prior-art refresh).
- [Zahn] J. Zahn, "Quantum fields on self-similar spacetimes", arXiv:2509.10360 v3, l.122 (eq. 4), l.126, l.131, l.133-139 (eq. 5-6), l.152-153, l.227, l.250-256
  (eq. 11-12), l.329-343 (eq. 19), l.347, l.348, l.490-496 (eq. 33), l.649-663, l.682-684. OPENED; copy: zahn.txt.

---

*Prepared with AI assistance; what was checked, and by what, is stated above.*
