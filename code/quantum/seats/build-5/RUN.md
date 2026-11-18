# RUN.md - build-5 BAR-QB2c (7aab7efe5a)

All commands run from `<LAB>/quantum/seats/build-5`, every Python run under memcap 8. build-4 and build-3 are
imported and read, never edited (build-4 runs with `QB2B_RUN1=1`, i.e. its run-1 configuration). BAR-QB2c.md was read on
the lead's instruction. Nothing is sealed, settled or pinned.

```
M="py -3.11 <MEMCAP> 8 py -3.11"
$M probe_tight.py        # 2-min probe: plain vs tight (build-4 qb2b_tight T1+T2) cost and radius per tube step -> probe_tight.log
$M deep_ratio.py 480     # PQ2c-1 evidence: tight vs plain step-matrix radii on all substeps of tube steps 480..719 -> deep_ratio.json
$M thth_c.py             # PQ2c-2: E_thth residual, real source + N-b, x = -6, -4.6, -2, -0.25, -0.01 -> thth_c.json
$M final_c.py            # PQ2c-3 replay (qb2final on the QB2b run-1 record) + PQ2c-4 replication row -> final_c.json, conv_c.json
$M make_outputs_c.py     # -> result_qb2c.json, controls_qb2c.json
```

## PQ2c-1: the lever decision (smarter first)
- **Probe (probe_tight.log).** Tight costs 4.5-7.5 s per tube step, against 0.3-0.5 s plain. The full rebuild is about
  75-90 min, plus about 85 min for the pass: roughly 2h50 of the 4 h budget.
- **But the order-0 step-matrix radius is the same at depth.**
  - Tube step 600: 9.0e-16 tight vs 1.0e-15 plain.
  - Tube step 690: 8.0e-14 vs 8.1e-14.
  - Tube step 715: 3.1e-13 vs 3.1e-13.
  - These track the pinned record's own step balls: setrad 2.5e-18 at tube step 600 and 4.7e-17 at 715.
- **Conclusion.** The deep radius is the record ball times the genuine state sensitivity of the HKA matrix (N ~ 1e2 deep),
  not interval dependency.
  - My build-4 attribution ("ball-recursion dependency") was wrong; this is a correction.
  - T1+T2 cannot move e from 6.66e-10 to ≤ 1e-10 on the sealed object (the pinned tube).
  - The full tight rebuild plus pass was therefore NOT run. deep_ratio.py quantifies the ratio over all deep substeps
    (deep_ratio.json).
  - The e, C8, N-a and N-c reported are the QB2b run-1 record, replayed (final_c.json).
- **Lever that would work (NOT run; it changes the sealed object).** Rebuild the inward tube with a higher Taylor order
  in `lin_tube_v2` (module constant KT = 36, remainder q^{KT+1} = 0.25^37 ≈ 1e-22 per step). That shrinks setrad and hence
  the step-matrix radii.

## PQ2c-2: method
- Y at each point is the full real-source h2: the centre side at -6; the pass check points at -2 and -0.25; the sonic
  series at -0.01; and both representations at -4.6.
- Y' and Y'' come from the ODE at the point: Y' = MY + F, Y'' = M_x Y + MY' + F_x, with M the engine's HKA step tape on
  the local ODE x-series.
- The x = -4.6 sonic-side representation (|Y| ≈ 600, after 1440 transport steps) encloses 0 but is too wide for N-b.
  The centre-side representation of the same solution is sharp, and it is used for N-b at -4.6. Both are in
  thth_c.json.

## PQ2c-4: the clock factor, certified
- At the centre the lapse is sqrt(A N^2 e^{2x}) → N e^x → e^{xc}. This uses A → 1, and P = v + O(v^4) with
  v = e^{x-xc} (build-4 centre_exact).
- So t_central = e^{xc} t_KHA, and the sonic point (x = 0, r = -t_KHA) sits at r/(-t_central) = e^{-xc}.
- Hence x_s = e^{-xc}, taken from the E5 box of the cm record. It agrees with the blind float x_sonic to 3.5e-15.
- The conversion is Ψ2 = 30π x_s² Ψ_blind / A(x_V) (referee-3 RECONCILE-QB2.md).

## Budget
- Phase-3 compute: about 40 min in total (probes 4 min, deep_ratio ~20 min, thth 2 runs, final/conv 2 min).
- Only my own PIDs.

## deep_ratio.json result (PQ2c-1 evidence)
- Over all 480 substeps of tube steps 480-719 (k 960-1439), the tight/plain order-0 step-matrix radius ratio is: min
  0.46, median 0.90, max 0.995.
- By band (median): k 960-1100: 0.59; k 1100-1300: 0.90; k 1300-1440: 0.99. The deepest band is where the radii are
  largest (1e-14 to 3e-13).
- The T1 overlap assertion never fired.
- Estimate (not a certificate): even if every per-step homogeneous contribution shrank by the best ratio (0.46), the e
  width would stay ≥ about 3e-10. The source columns and the start widths do not change. The needed factor is 6.7.
- So the sealed PQ2c-1 prediction cannot hold on the pinned tube. The full tight rebuild plus pass (about 2h50) was not
  spent.
