# DESIGN-QB1: the state-independent quantum source on the certified Evans-Coleman solution

Design seat design-1, 2026-11-17. Paper and reading only. Scratch files and logs are in `seats/design-1/` (`vcheck.py`,
`rncheck.py`, `weylcheck.py` plus `.log`). These are tiny sympy checks under memcap 8 and are NON-COUNTING. **Nothing on
this page is certified.** Labels: OPENED = full text read, copy in `src/`; ABSTRACT; SNIPPET; PAPER = derived here by
hand; SYM = a sympy identity checked here (exact, symbolic); MOCK = the check reproduced a published closed form exactly.

## 0. In plain words

- **The question.** QB1 asks for the one piece of the quantum stress tensor on the critical ball of light that does not
  depend on the quantum state. That piece is Zahn's local tensor V, the coefficient of u e^{2u}.
- **The finding (PAPER + SYM).** On the Evans-Coleman background the Ricci scalar is exactly zero, because radiation is
  traceless. Zahn's V then loses all dependence on the curvature coupling xi, and it equals **Bach/60** exactly.
  - So "which field model?" has a clean answer for QB1. Every free massless field has the same tensor shape here:
    - a scalar with any coupling;
    - by Brady-Ottewill's Page law, conformal spin 1/2 and Maxwell too.
  - The field content changes only one positive multiplier (the "c" anomaly coefficient).
  - **The sign of QB1 is the same for every field.**
- **What gets computed.** QB1 is therefore one certified object: the Bach tensor of the E5 solution. In practice that is
  three dimensionless profiles in the KHA variable x, built from N, omega, V, A and their first two x-derivatives.
- **Where it could fail.** Possibly the Bach tensor of Evans-Coleman is identically zero ("Bach-flat"). Then QB1 is zero
  and the leading quantum push is the state-dependent part. Node N7 decides this first.
- **Dropped control.** The brief suggested a wrong-coupling content null. That null is void on this background, since V
  cannot depend on xi when R = 0. It moves to the k = 1/2 Friedmann control, where R != 0 and xi must show up.

## 1. Exactly which tensor, in which field model

### 1.1 The general formula (OPENED, Zahn arXiv:2509.10360v3, src/zahn.txt)
- **Field.** A free massless scalar with arbitrary coupling xi, (box - xi R)phi = 0 (l.78-80). The state must respect
  self-similarity: w(u+δ, y; u'+δ, y') = e^{2δ} w(u, y; u', y') (eq. 4, l.122).
- **Structure.**
  - Eq. 5 (l.133): <T>(u+δ) = e^{2δ}<T>(u) + (δ/4π²)V.
  - For CSS, eq. 6 (l.136): <T_μν>(u, y) = e^{2u} [ T_μν(y) + (u/4π²) V_μν(y) ].
  - The π² is lost in the extraction. Three things fix it at 1/(4π²) here (and correct PRIOR-ART-QB.md, which writes "u/4pi"):
    - the parametrix normalisation 1/(8π²) with ln σ in eq. 14-16 (l.271-305) gives δ/(4π²);
    - the match with Brady-Ottewill in N2 below;
    - the Hayward MOCK, N0.
- **The tensor, eq. 19 (l.329-343):**
  V_μν = (1/720 - (6ξ-1)²/288) T(1)_μν - (1/240) T(2)_μν.
  - T(1)_μν = 4∇_μ∇_ν R - 4 g_μν □R - 4 R R_μν + g_μν R² (eq. 11, l.250).
  - T(2)_μν = 2∇_μ∇_ν R - g_μν □R - 2 □R_μν + g_μν R_ρσR^ρσ - 4 R^ρσ R_μρνσ (eq. 12, l.256).
  - Conventions: mostly plus, Wald curvature conventions (l.227).
  - The extracted text keeps only the denominators. The reading above is pinned by **MOCK N0**: `vcheck.py H` reproduces
    Zahn's Hayward values, eq. 33 (l.490-496), exactly: V_ττ = 1/80, V_uu = 17/240, V_θθ = -7/480, traceless.
    It also reproduces eq. 32 term by term (l.470-482, 7 coefficients from 2 unknowns).
- **What is state-independent (Zahn l.134, l.153).**
  - V_μν is state-independent and scheme-independent, and depends on xi. It is local in the curvature, conserved
    ("renormalization ambiguity ... automatically conserved", l.348), and traceless whenever □R = 0 (l.498-500).
  - T_μν(y) is state-dependent, and it is ambiguous by the finite renormalisations a T(1) + b T(2) and by the choice of
    origin of u (renormalisation scale).

### 1.2 Reduction on Evans-Coleman (PAPER + SYM)
- **R ≡ 0 on the background.** For k = 1/3, T^μ_μ = -ρ + 3p = 0. The Einstein equations then give R = 0 and
  R_μν = 8π T_μν = 8πρ[(4/3) u_μu_ν + (1/3) g_μν], exactly on the solution.
- **Consequence 1: xi drops out.** Every term of T(1) contains R, so T(1) ≡ 0 and **V = -(1/240) T(2) for every ξ**.
  - SYM check on Reissner-Nordström (R = 0, not conformally flat), `rncheck.py`: ∂V/∂ξ = 0 in all components.
- **Consequence 2: V is the Bach tensor.** SYM checks:
  - `vcheck.py B`, generic static metric: T(2) - T(1)/3 = -4B, with B_μν = ∇^ρ∇^σ C_μρνσ + (1/2) R^ρσ C_μρνσ.
  - Hence V(ξ = 1/6) = B/60 on any background (Zahn l.347 says "a multiple of the Bach tensor").
  - On R = 0, V = B/60 for all ξ (`rncheck.py`: V - B/60 = 0).
  - **Leading term on EC: <T_μν> ⊃ (u/(240π²)) B_μν, for any ξ.**
- **Consequence 3: one renormalisation ambiguity, not two.**
  - On EC, T(1) = 0 and T(2) = -240 V. The finite-renormalisation freedom in T(y) is therefore a multiple of the same tensor V.
  - A shift of the origin of u is also c·V (Brady-Ottewill l.167: a constant rescaling "changes the RSET by the addition of a multiple of the Bach tensor").
  - So the u e^{2u} coefficient is unambiguous, and the e^{2u} coefficient is ambiguous **only along V itself** (PAPER).
- **Route-A formula (SYM on Reissner-Nordström: exact).** With R = 0:
  V_μν = (1/120) □R_μν + (1/60) R^ρσ C_μρνσ + (1/240) g_μν R_ρσR^ρσ - (1/60) R_μ^ρ R_ρν.
  - Insert R_μν = 8πT_μν. The result needs ρ, the 4-velocity and the Weyl tensor (algebraic in the metric's second
    derivatives), plus □ of the fluid tensor. That is **jets of y = (N, ω, V) up to order 2 only**.
  - y' = f(y) is the sealed ODE (engine/system.py), so on the solution V(x) is a rational function of y. Its denominators
    are powers of N, A, (1 - V²) and det M (the sonic factor).
- **Field content (OPENED structure, garbled coefficients; Brady-Ottewill gr-qc/9804058, src/bo.txt).**
  - For conformally invariant fields, the Page law (l.100-106, eq. 3.1) gives eq. 3.8-3.9 (l.163, l.170):
    <T> = e^{4u}<T~>(y) + 8[α]u e^{4u} B(y) + e^{4u} S(y).
  - The Greek coefficient is lost in the extraction. PRIOR-ART-QB.md writes "8 beta u B". **Our check says it must be α**,
    the Weyl² ("c") coefficient. With α (h0 = 1) = 12/(2⁹·45π²) from l.128, 8α = 1/(240π²). That is **exactly** Zahn's
    1/(4π²) · (1/60). With β the two sides would differ by a factor of 3 and a sign. (Node N2; agreement in sign assumes
    BO and Zahn use the same curvature conventions; BO's conventions are not in the extracted text, UNTESTED.)
  - Multipliers, relative to one scalar, read from l.128 (12h0 + 18h_{1/2} + 72h1): scalar 1 (any ξ, by Consequence 1);
    each fermion helicity state 3/2; Maxwell (h1 = 2) 12.
  - So the QB1 tensor for any field content is N_c · B/(240π²) with N_c > 0.
  - **For "the radiation fluid's own quanta" (photons) N_c = 12.** That multiplier inherits the BO coefficient reading,
    which needs a PDF check of l.128-132. It is not needed for the sign.
- **Which tensor QB1 computes (decision).** B_μν of the E5 solution, delivered as dimensionless self-similar profiles in
  the Eulerian orthonormal frame (n = α⁻¹∂_t, e_r = a⁻¹∂_r):
  - b_nn(x) = r⁴ B_(nn), b_nr(x) = r⁴ B_(nr), b_rr(x) = r⁴ B_(rr), b_T(x) = r⁴ B_(θθ) = (b_nn - b_rr)/2 (traceless);
  - the fluid-frame invariant b_uu(x) = r⁴ B_ab u^a u^b, as the slicing-free headline scalar.
  - These are functions of x alone: B has dimension length⁻⁴ and r = e^{x-s}.
- **Interface to QB2 (PAPER).** KHA (18a) carries matter as 8π r² A T_nn: since ω = 4πr²Aρ, 2ω(1+kV²)/(1-V²) = 8π r² A T_nn.
  - The ħ source in the Hamiltonian row is therefore
    ε s e^{2s} · S_H(x), with S_H(x) = (2/π)·(N_c/60)·A(x) e^{-2x} b_nn(x), and ε = ℓ_P²/t₀².
  - The units are those of KHA's t, so e^{2s} = t₀²/t². The slicing and momentum rows take b_rr and b_nr in the same way.
  - Since s = u + f(x) for any two homothety times, the coefficient of s is slicing-independent (PAPER). The f·V shift
    lands in the ambiguous T(y).
- **Side node QB1b: the trace anomaly (state-independent only for conformal fields).**
  - Basis (C², E, □R) with BO's α, β (l.128-130; reading consistent with N2). On R = 0, the γ□R term and its scheme
    ambiguity ("48 or -72", l.132) drop out, and E = C² - 2R_μνR^μν. So:
    r⁴<T^μ_μ> = (α+β)·r⁴C² - 2β·r⁴R_μνR^μν.
  - SYM (`weylcheck.py`, generic time-dependent spherical metric): r⁴C² = 48(m/r - (4π/3)r²(T^θ_θ - T^A_A))².
  - In KHA variables (PAPER): r⁴C² = 12(A - 1 - 2ω/3)²/A² and r⁴R_μνR^μν = (16/3)ω²/A². These are closed forms in the
    E5 enclosure.
  - Scalar: 2880π² r⁴<T> = 12(A-1-2ω/3)²/A² + (16/3)ω²/A² > 0, positive everywhere.
  - Maxwell: 23040π² r⁴<T> = -104 r⁴C² + 496 r⁴Ric². The sign is not fixed: it is **negative at the sonic point** (mock below).
  - For ξ ≠ 1/6 the trace also contains a state-dependent (ξ - 1/6)□<φ²> term, so it is not a QB1 object.
  - The trace is one component of T(y). It is not a consistent source by itself.

## 2. What the certified background must supply, and where it could be singular

**Inputs.**
- E5 (HOLD d0af75803f; LOG.md l.45): V0, W0, x_c enclosures; the certified sonic series (radius ≈ 1 in x); the centre
  Frobenius series in v = p eˣ; the Lohner tube between them.
- The outward record to t = 0: build-1 outward_rec.pkl, and the sharper BG2 record (HOLD ed61d09268, LOG.md l.107;
  sideways/build-5/phase2/outward_C.pkl).
- Gear: gravity_engine C10 (css_linear.py l.89-106). It already writes the background as
  ds² = r²{A[(1-N²)ds² - 2ds dx + dx²] + dΩ²} with jets F0..F4, plus Christoffel and Ricci routines (l.108-128).
  QB1 adds □R_μν, the Weyl contraction and the frame projection.

**Two independent implementations (doctrine: every link double-implemented).**
- **Route A (on-shell, primary).** Use the route-A formula with R_μν = 8πT_μν. Jets ≤ 2, via y' = f(y), y'' = (∂f)f
  away from x = 0 and the centre. At the sonic point and the centre, use direct differentiation of the certified Taylor
  series.
- **Route B (off-shell, check).** Use Zahn's full eq. 19 with ξ symbolic, from metric jets N, A up to order 4, with no
  field equation used (the C10 style).
  - It must return T(1) ⊂ [-w, w] and R ⊂ [-w, w], where w is set by the enclosure width.
  - It must return a ξ-coefficient ⊂ [-w, w], and route-B minus route-A ⊂ [-w, w].

**Singular or delicate places.**

| place | QB1 tensor (local, Bach) | why | treatment |
|---|---|---|---|
| centre r = 0 (x → -∞) | **finite.** It is a smooth local tensor on a smooth regular centre. The Ayal-Piran divergence ("Unfortunately this tensor diverges at the origin", src/gr-qc_9704027.plain.txt l.145) belongs to their ad hoc 2D/4πr² tensor, not to a covariant local tensor. | the KHA chart degenerates (N → ∞); in r⁴ scaling b_ab = O(v⁴) cancels from O(1) pieces | evaluate in the blown-up chart (P, W, U)(v) from the certified centre series. Expected (PAPER, isotropy): B_(nn) → c0 t⁻⁴, B_(rr) = B_(θθ) → c0/3 t⁻⁴, B_(nr) → 0 like r. Certify c0 and b_ab ~ v⁴ |
| sonic point x = 0 | **analytic** (background analytic, E5) | the rational form via f(y) has det M in the denominator: a removable 0/0 | use the sonic Taylor series jets on |x| ≤ ρ_s and the rational form outside; the two must overlap on a shared shell (control) |
| light cone N = 1, x_lc = 1.357941265765732 (BAR-SF-L0 HOLD 596c7146ff, LOG.md l.84) | **nothing singular.** The light cone matters for the 2D-model state (TW regularity at the past light cone, src/tw2.txt l.1247) and for l ≥ 2 perturbations, not for a local tensor | none | the outward record covers it (BG2) |
| far side x → +∞ (t → 0⁻, r > 0) | finite. Curvature ~ r⁻² at t = 0, so b_ab → const and e^{-2x}b_ab → 0 | the outward record is not certified past x ≈ 12 (sideways/BAR-SF-S1b.md l.5) | tail lemma: bound |b_ab| on x ≥ X_far from the far-side invariant region (S1b) and the N → 0 asymptotics. QB1 claims a bound there, not values |
| Cauchy horizon, t > 0 | **outside QB1.** x ∈ ℝ covers only t < 0 | the quantum flux diverges there (PRIOR-ART-QB obstruction 5) | out of scope; say so in every claim |

## 3. Positive controls

All run through the **same** generic metric-in, V-out module as EC, except C5 and C6, which run in the KHA chart.

| # | control | expected value | status here |
|---|---|---|---|
| C0 | Zahn Hayward, ξ = 0, metric e^{-2u}(-2dτ² + 2du² + dΩ²) (eq. 23 reading; R = -e^{2u}, matching l.500) | V_ττ = 1/80, V_uu = 17/240, V_θθ = -7/480, traceless (eq. 33, l.490-496) | **MOCK PASS** (vcheck.log) |
| C1 | Zahn critical Roberts, metric e^{-2(u-r)}[2(1-e^{-2r})du² - 4du dr + dΩ²] (eq. 36 reading, l.526) | V_rr = e^{-2r}(5 - 17e^{-2r})/60, V_ur = -e^{-2r}(17 - 58e^{-2r} + 45e^{-4r})/240, V_θθ = -e^{-2r}(7 + 26e^{-2r} - 69e^{-4r})/480 (eq. 46, l.612-621) | **MOCK PASS** on these three components; V_uu still to add |
| C2 | Minkowski self-similar patch (Zahn l.358-377); in KHA, N = e^{-x}, A = 1, ω = V = 0 | V ≡ 0 | trivial; run as the zero-floor check |
| C3 | Reissner-Nordström M = 1, Q = 1/2 (static, R = 0, Bach ≠ 0): a nonzero control for the ξ-blindness and route A | V_tt = (4r²-8r+1)²/(640 r¹⁰), V_rr = -1/(120 r⁶), V_θθ = (4r²-8r+1)/(240 r⁶), for every ξ | **SYM** (rncheck.log) |
| C4 | Brady-Ottewill structure: (i) κ_q = 2 by matching the Einstein-tensor u-dependence (bo.txt l.230); (ii) the log coefficient 8α_BO = (1/4π²)(1/60) for one conformal scalar | (i) the forced rate equals 2 exactly in the QB2 ansatz; (ii) 1/(240π²) on both sides | (ii) **PAPER + SYM** (T(2) - T(1)/3 = -4B); (i) is a QB2 node |
| C5 | flat Friedmann k = 1/3 in the KHA chart (homogeneous self-similar radiation; conformally flat, R = 0) | V ≡ 0 for every ξ, Bach ≡ 0 | **SYM** in comoving coordinates (vcheck F). A KHA-chart closed form is needed (engine/frw.py keeps k general) |
| C6 | flat Friedmann k = 1/2 in the KHA chart (AMEND-LB-E-1 closed form: V0 = -2√2/5, N0 = 7√2/9, ω0 = 12/17, A0 = 25/17; LOG.md l.10) | with τ the cosmic time and a ∝ τ^{4/9}: V_ττ = 2(6ξ-1)²/(243 τ⁴), V_rr = 4(6ξ-1)²/(243 τ^{28/9}) (comoving), so V_ab u^a u^b = 2(6ξ-1)²/(243 τ⁴). Nonzero and ξ-dependent; zero at ξ = 1/6 | **SYM** (vcheck F2) |
| C7 | Tomasevic-Wu closed-form Roberts (2509.03584 Sec. 5.2, src/tw2.txt l.2831ff) | a control for the **2D s-wave model** (QB3), not for V: TW's tensor has no u e^{2u} term (Zahn l.682) | **wrong object for QB1.** Keep it for QB3, and record the disagreement as a vein |

Internal instrument checks:
- conservation ∇^μV_μν ⊂ [-w, w], two ODE identities in x;
- tracelessness;
- route A = route B;
- series and rational forms overlap on the sonic shell.

## 4. Content nulls

- **The "wrong conformal coupling" null is void on EC.** It is a prediction instead: V(ξ = 0) - V(ξ = 1/6) ⊂ [-w, w] on
  E5 (route B). That is the sealed **ξ-blindness** item.
- **ξ-null where it should break (C6).** The same pipeline on k = 1/2 Friedmann must give V(0) - V(1/6) = -(1/288)T(1)
  ≠ 0. For example V_ab u^a u^b(ξ = 0) = 2/(243τ⁴) against 0 at ξ = 1/6, detected at many σ. A pipeline that silently
  drops T(1) fails here.
- **Perturbed-background null (primary content null).** Run on E5 data with V0 moved outside its certified box, by
  1e-8 (well outside the 3.3e-20 radius) and integrated as a non-solution.
  - The Einstein equations then fail. Route B must report R and ∇^μV_μν enclosures that **exclude** 0.
  - Route A (which assumes R = 0) and route B must **disagree** detectably.
  - This tests that the instrument actually reads the background (the E5 null template: LOG.md l.40, "content nulls
    W0/xc/V0 all refuse").
- **Sign-flip null.** Flip the sign of the R^ρσR_μρνσ term in T(2). C0 and C1 must FAIL (cheap and immediate).

## 5. Honest obstruction list

1. **Chiba-Siino, Mod. Phys. Lett. A 12 (1997) 709 (ABSTRACT via INSPIRE 444860; PRIOR-ART-QB table).** "semiclassical
   equations of motion do not admit continuously self-similar solutions"; "regularity at the origin implies the
   disappearance of the echo" (scalar field, trace-anomaly method).
   - *Sidestepped:* QB1 claims no semiclassical CSS solution. It evaluates an O(ħ) source on the classical CSS
     background. The source grows like s e^{2s}, so it breaks CSS, which is exactly Chiba-Siino's statement. Zahn says
     the same: "the spacetime will in general cease to be self-similar" (l.153).
   - *Not sidestepped:* (a) validity needs ε s e^{2s} ≪ 1, a finite s-window; every claim is first order. (b) Their
     centre statement may concern a non-local tensor with a centre problem. For local V the centre is regular (sec. 2),
     but we cannot rule on their argument from an abstract. The owner may fetch DOI 10.1142/S0217732397000741.
2. **The Page-approximation objection (OPENED, TW 2509.03584 App. D, src/tw2.txt l.6315-6320).** Page's method
   "requires that (i) the physical spacetime be static, (ii) the Weyl rescaling factor be time-independent, and (iii) the
   background metric be an Einstein metric", and "cannot consistently describe this perturbative process".
   - *Effect on QB1:* none on the coefficient. QB1 uses Zahn's Hadamard scaling argument, not Page. The Page-derived BO
     coefficient agrees with it exactly (C4(ii)).
   - *Remaining bite:* the state. Zahn's statement holds only for states respecting self-similarity (eq. 4).
     - Existence is proved only for the Minkowski patch, Hayward and Roberts (l.126), with a sketch for others (l.649-663).
     - Persistence for other states is argued, "though not rigorously prove[d]" (l.131; l.681 "highly implausible").
     - TW dispute that such a state is physical and say q = D - 2 "is effectively built into the construction" (l.6382).
   - **Every QB1 claim must read: "in any Hadamard state respecting the self-similarity (Zahn eq. 4); existence of such a
     state on Evans-Coleman not proved."**
3. **2D vs 4D do not commute (OPENED).** TW: "dimensional reduction and quantization do not commute in general"
   (src/tw.txt l.409). Zahn: the 2D result "differ[s] ... in the coefficient of e^{2u}, as there is no geometric
   component growing linearly in u" (l.682).
   - QB1 is a 4D-Hadamard object with **no** 2D counterpart. QB3's s-wave model will not reproduce it, and the
     comparison is a tension to log, not a control.
   - Spin content: V is derived for scalars (Zahn). For spin 1/2 and 1 we rely on the BO/Page structure (conformal fields only).
4. **The radiation fluid is not a quantum field.** QB1 treats an independent free field on the background. "The fluid's
   own quanta" (photons, N_c = 12) is an interpretation: the vacuum polarisation of the Maxwell field on a background
   sourced by a classical photon gas. Double counting between the fluid and its quantum field is not addressed (UNTESTED).
5. **Finite-u dominance.** The u e^{2u} term dominates the e^{2u} term only once s ≫ |T(y)/V(y)| (Zahn l.683
   "generically the dominant one at large u"). Growth in s is logarithmic, so QB1 alone does not settle the QB3 sign at
   a finite trapping time. It settles the sign of the leading secular term.

## 6. The QB2 bridge: D(2) ≠ 0

- **What is needed (PAPER).** Use the ansatz h = h_ss + ε e^{2s}(h1(x) + s·h2(x)). Two forced problems result:
  - top order s e^{2s}: L(2) h2 = S_V. This is consistent by itself because ∇·V = 0;
  - next order e^{2s}: L(2) h1 + ∂_κL·h2 = S_T. This is state-dependent, and its divergence compensates ∇·(sV)
    (Zahn l.348).
  - Both need the homogeneous problem at κ = 2 to have only the gauge solution:
    - S_2 ∩ R_2 = span(G_2), which holds iff D(2) ≠ 0 (PROOF-B4.md Lemma 1.5(c) l.91-93, Thm 1.6(b) l.103-108);
    - then S_2 + R_2 = W_2 by dimension count (2 + 2 - 1 = 3);
    - so the forced solution exists and is unique modulo c·G_2. The gauge part cannot affect a geometric trapping sign.
  - Further requirements:
    - Sonic point: the singular exponent ρ(2) = (d - 2)/λ' = (-0.09903 - 2)/1.09903 ≈ -1.910 is not in {0, 1, 2, ...}
      (Fact 1.3, PROOF-B4.md l.56-64), so the forced analytic branch exists uniquely.
    - Centre: the source e^{-2x}b ~ v² vanishes at the centre, inside the regular class. n - 𝓜_c(0) is invertible for n ≥ 1
      (Fact 1.4, l.75-80), so there is no centre resonance.
- **Is it already certified? Yes, logically, from the M files.**
  - runs7/verdict_M.json l.4: "PM1": N_R40 = 1, Δ/π = [1.0000 ± 2.42e-5], parity true. D_new has exactly one zero,
    counted with multiplicity, in {Re κ ≥ 0, |κ| < 40}.
  - verdict_M.json l.5: κ* ∈ [2.810552549 ± 3.83e-10].
  - κ = 2 is interior to the half-disc, off its contour, and 2 ∉ κ*'s enclosure. So D_new(2) ≠ 0.
  - Then D_HKA(2) = -(2 - κ_g)D_new(2)/Ā_ss'(0) ≠ 0 too (Lemma 1.5(d), l.94), with κ_g = [0.3556992037 ± 7.41e-11]
    (verdict l.2) and Ā_ss'(0) = [-0.198928186801772 ± 4.57e-16] (verdict l.6).
  - The seal wording in BAR-LB-M.md l.36-37 and l.61 covers it.
  - Correction to ROADMAP QB2: the needed fact is "2 is not a zero of D_new", not "2 lies between κ_g and κ*". The gauge
    rate is deflated out and plays no role.
- **Still to add (cheap):** a direct certified point value D_new(2) with the M instrument (count_b4v3 point mode). It is
  belt and braces, and QB2 needs a lower bound |D_new(2)| for conditioning anyway. No float value for D(2) exists in the
  files: PROOF-B4 l.129 prints D(0) = -7.958 and D(1) = -0.544, and D(κ*) = 0.

## 7. Node table

| node | content | source | status | expected mock value |
|---|---|---|---|---|
| N0 | eq. 19 reading and conventions | zahn.txt l.329-343, l.470-496 | MOCK PASS | Hayward (1/80, 17/240, -7/480) exact |
| N1 | Roberts control | zahn.txt l.526, l.612-621 | MOCK PASS (3 of 4 components) | eq. 46 exact |
| N2 | BO = Zahn log coefficient | bo.txt l.128, l.148-170; vcheck B | PAPER + SYM | 8α = 1/(240π²) = (1/4π²)(1/60); T(2) - T(1)/3 = -4.000 B |
| N3 | R = 0 gives V = B/60 for all ξ | zahn.txt l.250-256; rncheck | SYM | ∂_ξV = 0, V - B/60 = 0 |
| N4 | Friedmann controls C5 and C6 in the KHA chart | engine/frw.py; AMEND-LB-E-1 | SYM in comoving coordinates; KHA chart TO BUILD | C5: 0; C6: V_ab u^a u^b = 2(6ξ-1)²/(243τ⁴) |
| N5 | E5 background jets to order 2 (route A) and order 4 (route B), centre and sonic charts | E5 d0af75803f; BG2 ed61d09268 | TO BUILD | widths ≲ 1e-15 on the tube (E5 tube width 6.8e-15 class, LOG.md l.13) |
| N6 | profiles b_nn, b_nr, b_rr, b_uu on x ∈ (-∞, X_far], plus the tail bound | this design | TO BUILD | traceless; conserved; b ~ v⁴ at the centre; e^{-2x}b → 0 at both ends; **sign not predicted** |
| N7 | B ≢ 0 on EC: certified b_uu(x0) ≠ 0 at ≥ 1 point | N6 | TO BUILD (first smoke) | nonzero expected but not implied: Weyl ≠ 0 at the sonic point (r⁴C² = 1.43769477, below) does **not** imply Bach ≠ 0 |
| N8 | centre coefficient c0 | centre series | TO BUILD | one number; no prediction |
| N9 | QB1b trace anomaly, scalar and Maxwell | closed form, sec. 1.2 | PAPER + SYM; certify | at x = 0: r⁴C² = 1.4376947722, r⁴Ric² = 0.1632723809; scalar 2880π² r⁴<T> = 1.6009671531 (> 0); Maxwell 23040π² r⁴<T> = -68.537155397 (< 0) |
| N10 | ξ-blindness on E5 (route B) | sec. 4 | TO SEAL | ⊂ [-w, w] |
| N11 | D_new(2) ≠ 0 | verdict_M.json l.4-5 + Thm 1.6 | CERTIFIED (logical consequence); direct point value TO RUN | sign of D_new(2) not on file |

The N9 mock uses the float sonic data from DESIGN.md (A0 = 1.861426722571081642101, ω0 = 0.325688889423068404479). It
is not certified.

## 8. Failure modes

1. **EC is Bach-flat (B ≡ 0).** Then QB1 = 0. The leading quantum term is the state-dependent e^{2s}T(y), and the
   field-independent sign story is gone.
   - Intervals cannot prove ≡ 0. If N7 refuses everywhere, we need a symbolic argument: B is a rational function of y
     on a 1-parameter curve.
   - This is the decisive first smoke.
2. **Interval blow-up at the removable singularities.** The sonic 0/0 and the centre O(v⁴) cancellation give refusals
   near x = 0 or x → -∞. Remedy: series jets on both charts, with an overlap shell. Failure means refusal, not a wrong
   number.
3. **Convention or normalisation slip.** Candidates: 4π² vs 4π, index placement, the frame sign of b_nr, the α/β
   misreading in BO.
   - C0, C1 and C4(ii) catch the first three. The BO multipliers (N_c) need a PDF check of bo.txt l.128-132 before any
     Maxwell number is quoted.
4. **Route B cannot see R = 0.** The order-4 jets on the tube are too wide, so the T(1) contamination is O(w·jets⁴) and
   is not small. Then route B is a consistency check only near the series charts.
5. **Far tail.** Without a certified bound on x ≥ X_far, QB2's response is undefined near t = 0. The claim must then be
   restricted to x ≤ X_far.
6. **Physical-scope failure (not numerical).** No self-similar Hadamard state is shown to exist on EC. QB1 stays a
   conditional statement, and that must be in the headline wording.

## 9. Cost driver

- **The size of the symbolic expression, not the cell count.**
  - Route A is □ of an algebraic tensor with jets ≤ 2. It is a modest rational function of y, and evaluation over
    ~550 tube cells (E5 Lohner step count, LOG.md l.13) plus the series charts takes minutes.
  - Route B has 4-jets of N and A and a full □Ric_μν. It is the expensive one: the sympy `cancel` on 4th-order jets
    grows like the C10 linearisation. Generate it once, pin it, and evaluate it on point jets (not boxes) with a mean-value
    remainder, to avoid dependency blow-up.
- Projection: under 1 h total, with no evening run needed. The heavy part is the one-time symbolic build.

## 10. Ceiling (no odds)

- The first certified value of a one-loop quantum stress-tensor coefficient on a true critical-collapse solution.
- If N7 is nonzero, also a sign that is the same for every free massless field: the direction of the leading secular
  quantum source on the critical ball of light.

**Magic-trick line:** "On the knife-edge ball of light, every kind of massless quantum field leaves the same
fingerprint, the Bach tensor, and a proof reads its sign."

## 11. Curiosities to surface (UNTESTED, labelled)

- **Weyl-zero locus.** r⁴C² ∝ (A - 1 - 2ω/3)², so the Weyl tensor of EC vanishes wherever A - 1 = 2ω/3. At the sonic
  point that quantity is +0.6443. Does EC have a conformally flat sphere other than the centre? One certified sign scan
  would decide it.
- **The renormalisation ambiguity collapses onto V itself on any R = 0 CSS background.** The QB3 sign at large s is then
  scheme-free. Is this already noted in the literature for radiation? It was not seen in the sources opened.
- **Field-dependent trace sign.** On EC the anomaly trace is positive for scalars and negative for Maxwell at the sonic
  point, while the Bach term has one sign for all fields. This is a vein if QB3 models built on the trace (BO-style)
  and on V disagree.
