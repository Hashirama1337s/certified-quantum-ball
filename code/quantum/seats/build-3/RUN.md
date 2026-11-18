# RUN.md - build-3 QB2 (certified first-order quantum response at rate 2)

All commands from `<LAB>/quantum/seats/build-3`, every Python run under memcap 8, no environment flags.
Engine, build-1/2, sideways and referee files are imported or read, never edited. `engine/cert/lin_mc.py` is NOT imported:
on import it rewrites `engine/cert/lin_fields2.pkl`. Its `hka()` is copied verbatim into `hka_copy.py`, and a diff check
is in the transcript (IDENTICAL). Nothing is sealed, settled or pinned. No BAR file was opened.

```
M="py -3.11 <MEMCAP> 8 py -3.11"
# exact sympy (non-counting derivations whose outputs the certified code evaluates)
$M bachsym.py        # -> bach_exprs.txt, bachsym.log   b_ab = Bach_(ab)[gt]; traceless; D_t == 0, (div Sigma)_s,x == 0 identically
$M bach_u.py         # -> bach_u.txt (exp(l0) -> 1/u, u = 1/N; numeric identity check)
$M n1_rows.py        # -> n1_rows.txt, n1_rows.log      N1: four HKA rows + sigma_I by undetermined coefficients
$M qb2sym2.py        # -> qb2_exprs.pkl, qb2sym2.log    per-unit-beta rows (generic + centre chart, u-orders)
$M a3.py             # D_t structure for C8 (coefficients of beta', beta)
# certified (python-flint arb)
$M n0check.py        # -> n0check.json   N0 cross-check of the Bach source against the QB1 instrument (result_qb1.json)
$M qb2stages.py      # -> stages.json, stages.log   N3 forced sonic series, N4a + N4 forced centre series
./run_pass.sh        # -> pass_ckpt.json, pass_steps.jsonl, pass.log   N5 transport, N7, route B, N2 on tube cells
                     #    (loops: $M qb2pass.py k k+240, exit 3 = resume, 1440 engine substeps; each chunk < 20 min)
$M qb2final.py       # -> final.json   N6 matching, C1 C2 C6 C8 C8' N-a N-c, gate, then N7, N8
$M make_outputs.py   # -> result_qb2.json, controls_qb2.json (instrument sha256)
```
Probe and debug scripts p_*.py and _patch*.py are non-counting scratch. `_patch*.py` applied the source edits recorded in the
transcript; the final instrument files are the hashed ones.

## Method (one paragraph per node)
- **N0 source.** On the background, g = r^2 gt with gt = -N^2 A ds^2 + A(dx - ds)^2 + dOmega^2, which depends on x only. The
  Bach tensor has conformal weight -2 in 4D, so b_ab = r^4 B_(ab)[g] = Bt_(ab)[gt] exactly. bachsym.py derives it symbolically
  in (u = 1/N, A, ..., A'''', l', ..., l''''), with l = ln N and primes = d/dx.
  - Along every chart parameter t, the x-derivatives are series operations D = (1/x') d/dt.
  - b is therefore a certified Taylor series on every transport substep: to order 40, plus a disc sup.
  - The normalisation matches QB1 with ratio exactly 1: sonic point, centre c0, and 10 tube points all overlap
    (n0check.json).
- **N1.** Undetermined coefficients. Substitute Y' = Q^-1(P0 + kP1)Y + Shat sigma into the linearised E_ss, E_sx, E_xx and the
  two fluid rows, and require each residual to be lam_r (I.Y - sigma_I.sigma). This gives 15 linear equations in 15 unknowns,
  with one solution. The rows are compared exactly with K5 and srcprop.
- **N3 / N4.** Augmented regular-singular series (12 components: Y, J, Phi, six constant drivers) with the engine's sonic and
  centre matrices M0 + 2 M1. Forcing columns come from the qb2core builders. The tail is the majorant of
  lin_stages.RegSing.evaluate. The power-of-2 diagonal rescaling (make_scaled) keeps the majorant constant small.
  - Centre: b = O(theta^4) and b_nr = O(theta^5) are used as exact structural zeros. This is the L3s PAPER premise; N4a
    certifies that those coefficients enclose 0.
- **N5.** The engine substeps (CachedLinear2, sha-pinned) supply the homogeneous block. On the same substeps, with start balls
  replicated from lin_stages_v2 nsub = 2, this build adds the forcing series and the sup of every entry on the disc
  rho_T = rho/2.
  - Sup bound: Taylor part plus a Cauchy tail from a crude 4x4-cell sup on the engine disc rho. Crude cell sups of the Bach
    columns wrap by 1e4 to 1e10, and the factor 0.5^41 removes this.
  - Transport: 12x12 Taylor series, the optimal-radius remainder (rem_opt), tolerance 1e-16 and order <= 40, and a real
    Lohner set per vector.
  - Near the sonic point, the Bach series is formed in mean-value form over the start ball (point jets plus Dual
    derivatives). This happens wherever the plain ball run gives a radius above 1e-13.
- **N6.** The 4x4 matrix [s^HKA, rho1, rho2, e_V] at x_m = -4.6 is the D_HKA(2) matrix; its determinant is printed.
  - For C8, Phi continuity at x_m first fixes the sonic homogeneous amount.
- **N7.** On every x-substep overlapping U = [-0.30, -0.20]: V' > 0 by series plus derivative tail, the end signs, and
  bisection.
- **N8.** Psi2 = Abar_2(x_V)/A(x_V) on the sonic leg: y_s^p + b s^HKA, evaluated inside the substep that contains x_V.
- **C6 route B.** Variation of constants with the adjoint series Psi' = -Psi A (bound ||Psi|| <= e^{4a|t|}). It runs from
  the centre-side solution at x_m back to x_V, accumulated during the forward pass as G = product of Psi(h) and S = sum of
  G w.

## Lead's speed levers (2026-11-18)
1. **Domain cut: USED.** Lemmas (b)-(d) need only the sonic point, x <= 0 and the centre: matching at x_m = -4.6, gauge
   Nbar(0) = 0, x_V in U inside (-4.6, 0). No outward cells or outward transport were built. Consequences, stated as gaps:
   - C1 at x = 0.5 and 3 is not checked;
   - N2 covers [-4.6, 0] cells, the sonic series and the centre chart, not (0, 12];
   - no h2 profile for x > 0.
2. **Float preview: NOT USED.** The certified pass is ~1.5 h in 240-substep chunks. Controls and Psi2 come out of the same
   certified pass, and qb2final.py gates Psi2 on N1, C1 and C8 before printing it. A separate float implementation would have
   been a third code path with no certified value.
3. **Lean precision: USED.**
   - prec 128 (not 200);
   - transport tolerance 1e-16 (not the count's 1e-20);
   - series order <= 40 (not 60);
   - disc rho/2;
   - mean-value Bach only where the plain ball radius exceeds 1e-13.

## Run record (2026-11-18)
- Order used: bachsym, bach_u, n1_rows (N1 PASS) -> qb2sym2 -> n0check (N0 cross-check PASS) -> qb2stages -> probes of qb2pass
  (k 0..40, kept as the start of the final pass; same code) -> run_pass.sh (k 40..1440, 6 chunks, 03:48 done) -> c2check ->
  qb2final -> make_outputs.
- Settled pass rate: 0.8 s per substep (plain) and 3.4-5.9 s (mean-value Bach); about 75 min wall for k 40..1440.
- Own PIDs stopped: the first n1_rows.py run (Stop-Process on its 4 PIDs). The fluid rows were too slow under expand;
  it was rewritten to take coefficients by differentiation.
- **qb2final.py, first run (final_run1.json): gate CLOSED by my own C1 sanity line.** I had written the expected matching
  coefficients as c = -1/2, d = -1. The correct values are c = +1/2, d = +1: phi G_2 = y_c^p - (1/2) rho1 - rho2 at the
  centre, so y_c^p - y_s^p = +(1/2) rho1 + rho2. The design's C1 criterion, the point equality at x = -6, -4.6, -2, -0.25,
  had already passed in run 1. The line was corrected (comment in qb2final.py), and the second run (final.json) is the
  record.
- Display trap: `arb.str(n)` adds the decimal-truncation error to the printed radius. All widths quoted here use `.rad()`.
