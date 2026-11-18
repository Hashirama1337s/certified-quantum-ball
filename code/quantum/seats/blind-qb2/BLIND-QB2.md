# BLIND-QB2: blind replication of the quantum-forced response of the Evans-Coleman solution

Seat: blind-qb2. I wrote only in this folder.

## Brief hygiene and independence

- The brief gave control values: kappa* ~ 2.8106 and gamma ~ 0.3558. It also gave the identity V_ab = B_ab/60 when R = 0, and the source normalisation. It did not give a target value, sign or functional form for Psi. I checked V = B/60 myself (see below).
- I did not open any forbidden path: build-3/, release/, any BAR-, RESULT- or TARGETS file, referee-2/, design-2/, DESIGN-QB2*, PROOF-QB2*, or lab memory. A directory listing of quantum/ showed some of those file names, but I opened none of them.
- I did read some lab background code that the brief allows: engine/system.py, engine/eqs_gen.py, the head of engine/bg.py and the head of engine/cert/lin_float.py. lin_float.py prints the background constants V0 = 0.11243940138809833622, W0 = 8.901323277755643519117 and PSCALE = 0.808668956318442829477.
  - I used W0 only to choose where to scan my own shooting parameter. My first scan had found only a band of nodal sonic solutions.
  - My bisection and collocation then found the solution on their own. I use these constants afterwards as a cross-check of the background, not as input.

## Formulation (mine)

- **Coordinates.** I use self-similar coordinates (u, x, theta, phi), with u = -ln(-t) and x = r/(-t). Here t is the polar time, normalised so that alpha(0) = 1, so t is the proper time of the central observer.
- **Metric.** g = e^{-2u}[-a du^2 + 2b du dx + c dx^2 + x^2 dOmega^2].
- **Background.** It is polar-areal written in (u, x): a = alpha^2 - A x^2, b = -A x, c = A. The fluid is p = rho/3, with 8 pi rho = e^{2u} w(x) and radial 3-velocity v relative to the polar normal observers.
- **Perturbation gauge (not polar-areal).** delta g_ux = 0 and delta g_thth = 0 in the (u, x) chart. The unknowns are da, dc, dw and dQ, where dQ is the perturbation of the fluid u^x. The residual gauge freedom is one constant: xi^u = e^{lam u} f(x) and xi^x = x xi^u, where f obeys a first-order ODE whose coefficient is -alpha^2, which never vanishes. I fix it with the gauge condition da(0) = 0.
- **Equations.** All equations come from sympy applied directly to the metric ansatz (qb2_derive.py): the Einstein tensor, nabla_a T^ab, and their linearisation with the ansatz eps e^{lam u}(...). Nothing was hand-derived.
- **Why the s e^{2s} level stands alone.** I write the response as eps[u e^{2u} h2 + e^{2u} h1]. Collecting the u e^{2u} coefficient gives L_{lam=2} h2 = source. Here L_lam is the operator for the ansatz e^{lam u}. The source V_ab = e^{2u} Vhat_ab is conserved, so the fluid equations at this level are homogeneous. The state-dependent part only enters at the e^{2u} level.
- **Quantum source.** I built it in two independent ways:
  - (i) From Zahn 2509.10360 (src/zahn.txt), eq. (19) with (11)-(12), at xi = 1/6 and R = 0: V_ab = (1/120) Box R_ab - (1/240) g_ab R_cd R^cd + (1/60) R_acbd R^cd. I reconstructed the coefficients from the garbled text and cross-checked all eight of them against Zahn eq. (32), the xi = 0 form (l.470-482). R_ab = 8 pi T_ab on shell.
  - (ii) The Bach tensor B_ab = nabla^c nabla^d C_acbd + (1/2) R^cd C_acbd, computed from the metric only, using spectral 4th derivatives (qb2_bach.py).
  - **Result:** Bhat/Vhat = 60.0000 pointwise, within 2e-6 at N = 60. So V = B/60 holds with the standard Bach sign.
- **Normalisation used.** The equation is dG_ab - 8 pi dT^fluid_ab = 8 pi K u B_ab/(240 pi^2), which is K u Bhat_ab/(30 pi) e^{2u}. K = 1 is Zahn's u V/(4 pi^2).
- **Numerics.** Chebyshev collocation in zeta = x^2 on [0, X], with X = 0.95 or 0.90, which is beyond the sonic point x_s = 0.8087.
  - Polynomials in zeta give a regular centre with the right parity. A polynomial basis can only represent the analytic branch at the sonic point.
  - The background is solved by Gauss-Newton, starting from my own centre shooting.
  - The perturbation uses all 6 equations (4 Einstein components and 2 fluid equations) plus the gauge row, solved by least squares. The system is overdetermined, so consistency is a built-in check: residuals are about 1e-10 against a source of size about 4.
- **Observable.** 1 - 2m/r = g^{ab} d_a r d_b r, with r = x e^{-u} left unperturbed in my gauge. Psi = -A * delta(1 - 2m/r) at x_peak, which equals A_2/A in the brief's notation.
- **Gauge invariance.** 2m/r is a spacetime scalar, and the background value of 2m/r depends only on x. At x_peak, d_x A = 0; in my equations A' is proportional to v, so A' vanishes exactly where v = 0. So the gradient of the background 2m/r vanishes there, and the Lie derivative along any gauge vector is zero (Stewart-Walker).
  - Numerical check: adding the gauge mode (da(0) = 1 instead of 0) changes Psi at x_peak by at most 8e-12.
  - At x = 0.3 the same shift changes Psi from 0.0227 to 0.2797, and at x = 0.75 from 0.0969 to 0.0232. So the test has teeth.

## Controls

1. **Background against the lab's stored constants (computed elsewhere, KHA chart).** I matched them as follows.

   | Quantity | Mine | Lab |
   |---|---|---|
   | x_s | 0.80866895631844 | PSCALE 0.808668956318443 |
   | v_s | 0.1124394013881 | V0 0.1124394013881 |
   | w0/2 | 8.90132327775 | W0 8.90132327775 |

2. **Growing mode.** kappa* = 2.8105525488 (N = 40, 50, 60 agree to 2e-11), so gamma = 1/kappa* = 0.3558019224. This matches the brief's check values 2.8106 and 0.3558. It is a positive control on the whole linear operator, including the sonic-point and centre handling.
   - smin(lam = 2) is about 1e-8, at the same level as away from any mode, against about 1e-17 at kappa*. So lam = 2 is not a mode and the forced solution is unique modulo gauge.
3. **Manufactured source.** I chose analytic regular fields: da = 0.7 x^2 e^{-x^2}, dc = 0.4 x^2/(1+x^2), dw = 2 + cos(1.3 x^2), dQ = x(0.5 - 0.2 x^2). I built the source from the symbolic operator applied to their exact derivatives, then solved.
   - Fields were recovered to within 2e-11.
   - Psi_exact = 0.059103772656349, Psi_solver = 0.059103772656343 (N = 40) and 0.059103772655859 (N = 60).
4. **Regular centre emerges unforced:** dc(0) is about 1e-13.
5. **Two forms of the source agree:** the Bach form and the fluid-on-shell form of V give the same Psi to about 1e-14.

## Result

**Psi = +0.1075199950 per unit K eps** (10 significant digits). The SIGN is positive.

| run | Psi |
|---|---|
| N=40, X=0.95 | 0.10751999497668573 |
| N=50, X=0.95 | 0.10751999497706545 |
| N=60, X=0.95 | 0.10751999497716666 |
| N=40, X=0.90 | 0.10751999497733232 |
| N=50, X=0.90 | 0.10751999497752145 |
| N=60, X=0.90 | 0.107519994977364 |

- The spread is 8.4e-13 and the gauge-shift changes are at most 8e-12. I quote an error of 1e-11 and trust 10 digits.
- The runs are double precision only. The two independent axes were resolution N and domain X.
- At the peak (x_peak = 0.62803916560721, background 2m/r = 0.47762850714765): delta(2m/r) = eps K u e^{2u} * 0.10752 * (1 - 2m/r) + (e^{2u}-level, state-dependent terms). So the quantum term raises 2m/r at its peak.

**Convention dependence (important for the lead).**
- u = -ln(-t), with t the central proper time. Rescaling t by t_lead = c t_here gives Psi_lead = c^2 Psi_here.
- Example: if the lead's chart puts the sonic point at x = 1, then c = x_s and Psi = 0.0703122148.
- The source is K u B_ab/(240 pi^2), with B in the convention above and mostly-plus signature with Wald curvature signs. K = 1 equals Zahn's u V_ab/(4 pi^2), because V = B/60 here.

## Surprises

- In the polar-areal chart, centre shooting finds a whole band of w0 (about 0.1 to 0.78) whose solutions reach the sonic line through a node. The Evans-Coleman solution sits elsewhere, at w0 = 17.8026465555, where the shooting function changes regime.
- The 2m/r peak (v = 0) lies inside the sonic point.
- No digits were identified, so identify.py was not used.

## Files

These are all in <LAB>/quantum/seats/blind-qb2/.

- **Derivation:** qb2_derive.py and qb2_bach.py, producing qb2_exprs.pkl and qb2_bach.pkl (sympy srepr strings).
- **Solver:** qb2_common.py, qb2_bg.py (shooting), qb2_coll.py (collocation), qb2_bgrun.py (background), qb2_pert.py (operator, source, observable), qb2_main.py (kappa and forced runs), qb2_control.py (manufactured source).
- **Results:** blind_qb2_result.json, main_all_X0.950_40-50-60.json, main_forced_X0.900_40-50-60.json, bg_X0.9*.json, control_manufactured_X0.950.json.
- **Exploration scripts** (not used for results): qb2_shoot_scan*.py, qb2_look.py, qb2_cont.py, qb2_probe_rhs.py.
- **Checkpoints:** none, so a resume test is not applicable. The pickles hold only lists and strings, no flint objects.
