"""Certified gravity engine: the reusable code behind the lab's certified gravity results.

In plain words: each main-quest project (K3, BHSLICE, ...) leaves behind the piece of code that made its certificate
work. Those pieces live here, so the next project imports them instead of re-deriving or copying them.
CAPABILITIES.md maps what the engine can certify today; HARVEST.md ranks what is still to be extracted.

How a module gets in (doctrine #121):
- only after its certificate holds: the project's sealed test is settled HOLD and the referee round is done;
- the lead extracts it from the project folder; project code is never moved here before its verdict;
- the module docstring names its DESIGN file, the seal id of the certificate, and the CERTIFIED nodes it backs;
- the project's DESIGN lists it in `produces: engine:gravity_engine/<module>.py`, and later projects that use it list it
  in `reuses: engine:gravity_engine/<module>.py`. `lab.py dryrun --release` checks that the file exists and backs a
  CERTIFIED node; `lab.py status` counts modules and reuse.

Versions: ENGINE_VERSION below is the package version (bumped whenever a module is added or changed). Each module may
carry its own version; lichnerowicz_axisym.ENGINE_VERSION stays 2 (its file is unchanged since v2, so runners that pin
its sha256, e.g. BHSLICE p2c.py, are unaffected by package versions 3 and 4).

Modules:
- lichnerowicz_axisym.py (2026-11-14): certified axisymmetric Lichnerowicz / Hamiltonian-constraint sandwiches; first problem: Bowen-York spinning puncture. Backs BHSLICE DESIGN-P2 nodes R1-R5, P1, K1, K2 (seal d62bfc0e0b; refereed REFEREE-P2a; blind-replicated VERDICT-REP2). Tests: gravity_engine/tests (reproduce P2a byte for byte). Module version 2 (2026-11-14): fit-only fit_ops / cache_ops options, see README CHANGELOG.
- heat_kernel_green.py (2026-11-14, package v3): certified oscillation (Green) constants on closed Riemannian 4-manifolds: Lemma G (Li-Yau + Doeblin) and Lemma G' + Z' (Andrews-Clutterbuck spectral tail). From K3 stage2/s2.py and wp5/wp5.py; backs S2b (seal 1e719251da) and WP5 (seal 055009634d); refereed REFEREE-S2, REFEREE-S2CODE, REFEREE-WP5. Tests: tests/test_heat_kernel_green.py (123 stored values byte for byte).
- kahler_route2.py (2026-11-14, package v3): certified Route-2 continuity assembly (Theorem 2) for a near-Ricci-flat Kahler surface, uniform in the gluing parameter; Kummer K3 diameter bounds L5 / L5'. From wp5/wp5.py and s2.py; backs WP5 V2-V4 and S2b S2-S3 (the published K3 certificate, DOI 10.5281/zenodo.24216647). Tests: tests/test_k3_route2.py (60 intervals + 4 leading values byte for byte).
- arbkit.py: shared arb helpers for the two K3 modules (verbatim from wp5.py l.32-131, run at 160 bits).
- bowen_york_sqrt.py (2026-11-14, package v4): certified sqrt-map sandwich for the Bowen-York spinning puncture (second, independent certifier of the same PDE as lichnerowicz_axisym). Rebuilt from BHSLICE rep4/rep4.py with the REFEREE-REP4 fixes. Entry rule: backed by the REP4 sealed claim ad62cca292, verdict HOLD (a blind-REPLICATION verdict: REP4 reproduced P2c's delta(3000), delta(10000)); refereed REFEREE-REP4 (method sound, fixes applied); comparison lemma BHSLICE LEMMA-REP4.md (written for this module, not yet refereed). Backs REP4's certified nodes delta(3000), delta(10000) and K1. Tests: tests/test_bowen_york_sqrt.py (REP4 resolution A reproduced with identical exact rationals). Module version 2 (2026-11-14, package stays v4): adds the tangent certificate of DESIGN-P2d (certify_tangent, TL cell/tail bounds; tests/test_bowen_york_tangent.py) and REFEREE-LEMMA FIX F (residual v-degree max(8K, 1)); the tangent certificate backs no sealed node yet (first use: BHSLICE P2d, claim 9677804a51).
- argument_count.py (2026-11-17, package v5): certified zero count of a real-symmetric analytic function in a right half-disc (argument principle; disc plan with an exact-Fraction arc check, half-plane segment test, arg sharpening, 0.25 rule, parity check, resumable JSONL ledger). Verbatim from LIGHTBALL count_b4.py / count_b4v3_seal.py / core.py (published certified-evans-coleman, DOI 10.5281/zenodo.24303877); backs the M count N(40) = 1 (seal b1b0ee1cde HOLD; refereed REFEREE-B4FW, B4FW2, B4TAIL). Tests: tests/test_argument_count.py (published 99-disc ledger replayed; exact-polynomial controls, null, refusals).
- frobenius_levels.py (2026-11-17, package v6): level-shift lemma (symbolic) + certified exponent-collision counter rho = (d - kappa)/lambda at a singular point, checked at the level actually integrated. Verbatim from LIGHTBALL sideways/build-2/x1_cert.py; backs X1 (seal cb9af3b657 HOLD; referee-1, referee-2 REFEREE-SF-X1). Tests: tests/test_frobenius_levels.py (seven X1 levels replayed).
- css_linear.py (2026-11-17, package v6): symbolic first-order Einstein-Euler (p = k rho) perturbations of a CSS spherical background, all l, polar and axial. Verbatim LIGHTBALL sideways/build-1/lin_core.py; backs X1 (same seal and referees; known limits: l = 1 sonic identity with background row 0 only; rank condition (C41, C44) != 0). Tests: tests/test_css_linear.py (pure-gauge exactness + perturbed-background null).
- greybody.py (2026-11-17, package v7): certified Schwarzschild greybody factors (A_in, A_out, Gamma; Regge-Wheeler any s, Zerilli s = 2), point mode and Taylor models in w. Verbatim HAWKING build-1/hkcore.py; backs HK-1b (seal 2040a9985f HOLD; PROOF-HK.md; referee-1 twice + release referee; blind Teukolsky/MST replication; published in certified-page-time, DOI 10.5281/zenodo.24327159). Importing sets ctx.cap = 64. Tests: tests/test_greybody.py.
- greybody_bounds.py (2026-11-17, package v7): one-sided emission bounds (sech^2 band bound, Planck tails, l > 60 lemma). Verbatim HAWKING build-1/hkbounds.py except its import of the core. Same seal and tests.
"""
ENGINE_VERSION = 7
