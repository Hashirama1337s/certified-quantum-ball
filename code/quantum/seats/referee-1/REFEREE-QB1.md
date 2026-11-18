# REFEREE-QB1 (referee-1, hostile)

Process note: this seat sent an early hand-back by mistake. It contained only the text "PLACEHOLDER - not used". This file is the real report.

## In plain words
The main physics claim holds. On the Evans-Coleman background R = 0, so the state-independent quantum term is Bach/60 for every coupling. I checked this with my own symbolic code. The two headline numbers reproduce to 22-25 digits from my own geometry code: sonic b_uu and centre c0. What needs fixing is the wording, not the mathematics. Three things are overstated: "28 points" (19 distinct x), "certified" for the bridge (no bar has held it yet), and the spin-1/2 / Maxwell part, which rests on Brady-Ottewill's conformal argument, not on Zahn's theorem. QB1 also still has no content null that passed (PQ6 FAILED).

## C1 - CORRECT WITH FIXES
- **Own code.** `ref_geom.py` and `ref_tests.py`, written from scratch in Wald conventions. Zahn eq. 11, 12 and 19 were read at src/zahn.txt l.250, l.256 and l.329-343; conventions at l.227.
  - Generic diagonal metric depending on 2 coordinates, R != 0 (R = 85638304669382887875/23194084542881794592 at a rational point): T2 - T1/3 + 4B = 0 exactly in all 16 components, and tr B = 0.
  - So V(xi = 1/6) = B/60 on any background. On R = 0, T1 = 0 term by term, so V = B/60 for every xi.
  - Hayward (eq. 23 reading, R = -e^{2u}): V(xi = 0) = (1/80, 17/240, -7/480), which is eq. 33 (l.488-496) exactly. A sign-flip null on the R^{cd}R_{acbd} term refuses: it gives -1/32.
  - Reissner-Nordstrom: dV/dxi = 0 and V - B/60 = 0. Design C3 closed forms reproduced exactly.
- **POST HOC.** A pp-wave with H = x^4 + x y^3 + sin(u) x^2 y^2 (R = 0, not Ricci-flat, not static):
  - dV/dxi = 0 and V = B/60;
  - V_uu = -(1/240) Lap^2 H = -sin(u)/30 - 1/10, matching a hand prediction;
  - all other components are 0.
- **Friedmann.** k = 1/2: V_tautau = 2(6xi-1)^2/(243 tau^4) and V_xx = 4(6xi-1)^2/(243 tau^(28/9)), Bach = 0. k = 1/3: V = 0 for all xi. (ref_frw.log)
- **Normalisation.** Brady-Ottewill (OPENED, src/bo.txt l.128, eq. 3.4a): "= 12h0 + 18h 1/2 + 72h1 /(2^9 45 pi^2)". This gives 8 alpha(h0 = 1) = 1/(240 pi^2) = (1/4 pi^2)(1/60), as the design says.
- **Spin content.** Duff, hep-th/9308075 (OPENED, PDF text, eq. 30-31): "b = 1/(120(4pi)^2)[NS + 6NF + 12NV]", with "NF ... spin 1/2 Dirac fermions", and "the contributions to the b coefficient are all positive".
  - This agrees with BO: per Dirac fermion 18·4/12 = 6, Maxwell 72·2/12 = 12, scalar 1/(1920 pi^2) in both.
- **Convention check (Zahn vs BO).** BO's conventions are not in the extracted text (grep for convention/signature: 0 hits), so the absolute sign agreement stays UNTESTED.
  - It does not matter for the "positive multiplier" claim. Inside BO every field carries 8 alpha_s B with alpha_s > 0, so the ratios between fields are positive in any convention.
  - BO write B = C_{..;..} + (1/2) R C. Under a flip of the Riemann sign the relative sign of these two terms would flip, so this form fits a Wald-sign Riemann. The overall sign under a signature flip cannot be fixed from the text.
- **Fixes.**
  1. Zahn's theorem covers only the scalar. The spin-1/2 and Maxwell statements rest on BO's exact conformal transformation law (eq. 3.1, 3.8) for conformal fields, in a state invariant under d_u of the rescaled metric. Label them that way.
  2. Restrict "every free massless field" to spin <= 1.
  3. Novelty. V = B/60 on R = 0 follows in one line from Zahn l.347 ("reduces to (a multiple of) the Bach tensor") plus R = 0. Present it as a cited consequence, not a finding. The new content is the evaluation of B on Evans-Coleman.

## C2 - CORRECT WITH FIXES
- In the massless case the ambiguity is a T1 + b T2 (Zahn l.250-256). T1 vanishes identically when R = 0, and T2 = -240 V. A shift of the origin of u adds a multiple of V (Zahn eq. 6; BO l.167).
- **Fix.** Say what this does not cover. It holds on the background only:
  - The linearised T1 (4 DD dR - 4 g box dR) is nonzero, but it enters only at O(eps^2) in QB2, so this is harmless.
  - The remaining state-dependent part T(y) is not fixed by it.
  - Brady-Ottewill already say this for conformal fields (l.167). Call it a folklore-level observation.

## C3 - CORRECT WITH FIXES
- **Own recomputation (`ref_ec.py`).** The pieces that are mine:
  - my curvature code on -F dt^2 + G dr^2 + r^2 dOmega^2, with F = A N^2 z^2 and G = A;
  - my own mpmath series reversion;
  - the frame n = d_t/sqrt(F), e_r = d_r/sqrt(G), with u = gamma(n + V e_r).
  - Input: only the certified manifold coefficients (outward_C.pkl sm, lin_tube2 cm, as mid strings).
- **Sonic point** (ref_ec_sonic2.log):
  - b_nn = 0.4157372840770100272197622, b_nr = 0.1190421959584470894048207, b_rr = -0.1171195084663042197300078;
  - **b_uu = 0.4466737796695366321323105**, identical to build-1's 25 digits;
  - traceless to 3e-51.
- **Einstein check on my reading of the background.** The residual r^2(G - 8 pi T) is <= 2.1e-32 in all components. With the sign of V flipped, the nr component gives -0.10627. So the V sign convention is confirmed.
  - Note: my first run carried a T_(nr) sign bug on my side (u.n = -gamma). I found and fixed it, and the result is logged.
- **Centre** (ref_ec_centre.log):
  - b_nn e^{-4x} at x = -7..-12 runs from -173.0976478 to -173.1004442653;
  - Richardson in e^{2x} gives **c0 = -173.1004443922995663248**, matching build-1's -173.1004443922995663248427;
  - b_rr/b_nn -> 1/3 (0.333333333425 at x = -12), which is isotropy;
  - residual <= 5.6e-32.
- **Zero bracket.** [-0.82611434221, -0.82611433983] is a certified sign change: at least one zero, and unique only inside [-0.826124, -0.826105] by the derivative bound.
  - Fix: PQ2 FAILED as sealed, so the wording must not imply "exactly one zero on [-12, 12]" until PQ2b holds.
- **Routes.** The 28 route rows sit at **19 distinct x** (9 x values appear in both the polar and Cartesian charts).
  - Fix the wording to "28 comparisons at 19 distinct x". PQ4's threshold of >= 18 points is still met.
  - Both routes use the same ODE-derived jets, so their agreement tests the algebra, not independent physics.
- **Controls.** Pass (controls_qb1.json all true). My code reproduces the RN values at r = 3: V_tt 4.471921624413622584633101e-6 and V_rr -1.143118427069044352994970e-5.
- **Fix: the null.** The content null PQ6 FAILED (the off-solution V0 still solves the ODE). Every C3 number stands without a passing null until PQ6b holds. My independent agreement is POST HOC support, not a substitute.

## C4 - CORRECT WITH FIXES
- **The value.** D_new(2) = -2.4126270857981028564 +/- 1.13e-9 (dnew_qb1.json). Rounding to -2.41262709 gives a total error of 5.33e-9, so the quoted +/- 5.3e-9 is honest.
  - Gauge deflation D_HKA = -(kappa - kappa_g) D_new / Abar'(0) reproduces the file's D_HKA at kappa = 0, 1, 2 to relative 3e-11, 2e-11 and 7e-12.
- **Non-vanishing is already implied.** verdict_M.json PM1 has N_R40 = 1 with kappa* in [2.810552549 +/- 3.83e-10]. The zero set does not depend on x_m wherever c_V(x_m) != 0 (Thm 1.6(b) proof).
  - The ratio D(-4.6)/D(-2) is kappa-dependent (13.64, 13.31, 12.99), as expected for a non-Wronskian determinant.
- **Fixes.**
  1. "Certified" here means computed with the sealed count instrument. PQ7 FAILED as sealed (wrong x_m) and PQ7b has not run.
  2. Thm 1.6(b) concerns only the homogeneous problem. "Solvable, unique modulo G_2" for the forced top-order (V-sourced) problem needs four further PAPER steps, to be written out:
     - constraint propagation with a conserved source;
     - the forced sonic Frobenius step, using rho(2) ~ -1.91, whose constants are still [float]/REDUCED in Fact 1.3;
     - the centre step for an O(v^2) source;
     - the dimension count 2 + 2 - 1 = 3.
  3. The e^{2s} (h1) problem has a state-dependent source and is not solvable as posed.
  4. The far side x > 12 has no tail bound.

## C5 - CORRECT WITH FIXES
- The obstruction list is honest. Chiba-Siino is correctly labelled ABSTRACT, and the state-existence conditional is stated.
- **Not sufficient as written, for three reasons:**
  1. The rebuttal "QB1 uses Zahn, not Page" is true only for the scalar. For spin 1/2 and Maxwell QB1 does use BO's Page-type transformation law. It is the exact law, not the Page state approximation that TW criticise (tw2.txt l.6315-6320), but the claim must say so.
  2. The wording must read: "in any Hadamard state respecting self-similarity (Zahn eq. 4); existence on Evans-Coleman not proved; spin 1/2 and Maxwell via the conformal transformation law."
  3. The magic-trick line "a proof reads its sign" overclaims, because b_nn changes sign. Only the field-independence of the sign at each point is supported.

## Blocks release
- PQ6b (a passing content null) and PQ2b (global uniqueness of the zero) have not run.
- The C4 forced-problem lemmas are not written.
- The spin-1/2 / Maxwell attribution and the "19 distinct x" wording must be fixed.

## Single next strike
Run BAR-QB1b PQ6b, the off-shell omega(1 + 1e-12) null. It is the only missing teeth test on numbers that two independent codes now agree on.

Files (all in <LAB>/quantum/seats/referee-1/):
- code: ref_geom.py, ref_tests.py, ref_frw.py, ref_ec.py, chk_defl.py, routes.py
- logs: ref_tests_1.log, ref_tests_2.log, ref_frw.log, ref_ec_sonic.log, ref_ec_sonic2.log, ref_ec_centre.log
