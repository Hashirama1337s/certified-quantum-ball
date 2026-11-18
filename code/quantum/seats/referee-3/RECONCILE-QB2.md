# RECONCILE-QB2 (referee-3): build-3 Psi2 against blind-qb2 Psi

Seat referee-3, 2026-11-18. Wrote only in seats/referee-3/. Scripts: conv.py, reconcile_qb2.py (log: reconcile_qb2.log),
both run under memcap 8 (peak 0.02 GB). Inputs are read straight from build-3/result_qb2.json (PQ2-8.N8) and
blind-qb2/blind_qb2_result.json.

## In plain words

The two numbers are the same physical result, written in three different conventions. The two seats use different clocks,
source units and readouts. Converted exactly, the blind value 0.10751999497716666 becomes 3.461635803705687 in the build's
units, against the build's certified 3.4616358037141459203 +/- 1.15e-7. That is a relative difference of 2.4e-12, about
four orders of magnitude inside the build's own enclosure. Neither code is wrong. The text around the number has two
fixes, and those block release until made.

## The exact conversion (each factor read from code, not fitted)

**Psi2(build) = 30 pi * x_s^2 * Psi(blind) / A(x_V)**

There are three factors.

1. **Clock: factor x_s^2.**
   - Build: metric g = e^{2(x-s)}[A(1-N^2)ds^2 - 2A ds dx + A dx^2 + dOmega^2] (n1_rows.py l.33-34). So r = e^{x-s}, with x = 0
     at the sonic point. The KHA time therefore obeys r_s = -t_K.
   - Blind: t is the central proper time, alpha(0) = 1, and the sonic point sits at r/(-t) = x_s = 0.80866895631844.
   - Hence t_K = x_s t_c, s = u - ln x_s, and e^{2s} = e^{2u}/x_s^2.
   - Independent check, from the files alone:
     - x_V(build) - ln(x_peak/x_s)(blind) = -4.5e-15;
     - A(build) - 1/(1-2m/r)(blind) = 1.6e-14.
   - **The s shift does not mix the levels that matter.**
     - s e^{2s} = (u - ln x_s) e^{2u}/x_s^2. The u e^{2u} coefficient is exactly (build coefficient)/x_s^2.
     - The -ln x_s piece goes to the state-dependent e^{2u} level only.
     - The same holds for the source: s B = u B - (ln x_s) B.
2. **Source unit: factor 30 pi.**
   - Build: residual G - 8 pi T - Sigma, with Sigma = eps e^{kap s} sigma (n1_rows.py docstring l.5-7). sigma is built
     from beta = e^{-2x} b and b = r^4 B_(ab) (qb2core.py l.6), with no 1/(30 pi). The design divides out
     K = N_c/(30 pi) (PROOF-QB2-LEMMAS-v2 sec. 1).
   - So the build's per-unit-K source is Sigma = eps s B.
   - Blind: rhs = Bhat/(30 pi) per unit K (qb2_pert.py l.74-80). That is, Sigma = eps u B/(30 pi) per unit K.
   - Both K's describe the same physics:
     - one conformal scalar has <T> = u B/(240 pi^2) (DESIGN-QB1 l.61, l.81; blind "K = 1 is Zahn's u V/(4 pi^2)");
     - so K_blind = N_c and K_build = N_c/(30 pi).
   - The 8 pi is inside both: 8 pi/(240 pi^2) = 1/(30 pi).
3. **Readout: factor 1/A.**
   - Build: A = Ab(1 + eps a e^{kap s}) (n1_rows.py l.31), so Abar = deltaA/A. Then Psi2 = Abar_2/A = deltaA/A^2 =
     delta(2m/r) (qb2final.py l.155).
   - Blind: Psi = delta(2m/r)/(1-2m/r) = A delta(2m/r) (qb2_pert.py l.110-116).
   - The blind's brief called this "A_2/A". That label collides with the build's "Abar_2/A".
   - Cross-check: Abar2_xV/(30 pi x_s^2) = 6.6267701263947/(30 pi x_s^2) = 0.107519994977429 = blind Psi (N60 run
     0.10751999497716666).

**Agreement, all six blind runs** (reconcile_qb2.log):

| quantity | value |
|---|---|
| relative differences | -6.92e-12, -3.39e-12, -2.44e-12, -9.03e-13, 8.56e-13, -6.08e-13 |
| build relative width | 3.32e-8 |
| implied factor/pi | 30.0000000000733 |

The lead's two guesses fail on this scale:
- the "5 pi^2" guess gives 3.46976872434, which is off by 2.3e-3 relative;
- no clean factor showed in the lead's quick look because the 1/A factor was missing.

**Physical numbers, per unit N_c eps:**
- delta(2m/r) at x_V, in KHA s e^{2s}: Psi2/(30 pi) = 0.036729096623;
- the same in central-proper-time u e^{2u}: 0.0561653802878;
- for photons (N_c = 12), K * Psi2 = 0.440749159476.

## Verdicts

**C1. Build Psi2 = 3.4616358037141 +/- 1.15e-7 (KHA units, per K eps s e^{2s}): CORRECT.**
- The blind replication agrees to 2.4e-12 relative after an exact conversion.
- All three factors were read from the code and checked against independent file quantities.
- Sign: positive in both.

**C2. Blind Psi = +0.1075199950 per unit K eps: CORRECT**, in its own stated convention (central proper time, source
u B/(240 pi^2) inside 8 pi <T>, readout A delta(2m/r)).

**C3. The lead's reading "same quantity, size disagrees": WRONG.**
- It assumed the same readout, and the readouts differ by A(x_V) = 1.91434642526073.
- It assumed the same source unit, and the units differ by 30 pi.

**C4. Text that states the unit: CORRECT WITH FIXES.** Both fixes block release.
- **(a)** paper.md Theorem 3 (l.142) and DESIGN-QB2-v2 sec. 2.2 write delta m_U = eps e^{2s}(s Psi2 + Psi1). This is
  missing K.
  - Fix: write delta m_U = K eps e^{2s}(s Psi2 + Psi1), K = N_c/(30 pi).
  - Or define Psi2 per unit N_c and quote 0.036729096623.
  - DESIGN-QB2-v2 sec. 7 node N8 says "(times N_c eps s e^{2s})". That is wrong by 1/(30 pi). Fix it in a v3 note, not by
    editing.
- **(b)** {{QB2_PSI2_UNIT}} must state all three conventions:
  - the KHA clock (x = ln(r/(-t)), sonic at x = 0, t_KHA = x_s * central proper time, x_s = 0.80866895631844);
  - the s e^{2s} coefficient of delta(2m/r);
  - K = N_c/(30 pi).
  - Add the converted blind value with this formula as the replication row.
  - Without these, a reader comparing with any proper-time literature value is off by x_s^2 = 0.654.
- {{QB2_PSI2_PHOTONS}} = 0.440749159476 (K * Psi2 at N_c = 12), but only if (a) is fixed so that "coefficient" means
  K * Psi2.

**C5. The blind seat's own note "Psi_lead = c^2 Psi_here": CORRECT.** It is the clock factor alone. It cannot close the gap
without the 30 pi and A factors.

## What blocks release

The numbers do not block release. Two things do:
- the missing-K statement in Theorem 3 / design sec. 2.2 / node N8;
- an unspecified unit placeholder.

A secondary issue: the blind result is double precision with no enclosure. It is a replication, not a certificate. Say so
in the scorecard row.

## Single next strike

Seal a one-line replication row:
- Psi2 = 30 pi x_s^2 Psi_blind / A(x_V);
- pass if the relative difference is <= 1e-9;
- record it as POST HOC, since this reconciliation was done after both numbers existed.

Then fill the paper's unit placeholder with the three conventions above.
