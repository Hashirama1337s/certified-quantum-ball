# RUN.md - build-4 QB2b (phase 2 of QB2; build-3 imported, never edited; no BAR file opened)

All commands from `<LAB>/quantum/seats/build-4`, every Python run under memcap 8.
Nothing is sealed, settled or pinned.
`engine/cert/lin_tube2_k1_3.pkl` is verified unchanged after the tube build: sha f6bff9e2…, mtime Nov 16.

```
M="py -3.11 <MEMCAP> 8 py -3.11"
$M tube_bg2.py          # item 1: inward tube from the BG2 V0 box -> lin_tube_bg2_k1_3.pkl, tube_bg2.log (515 s)
$M tube_compare.py      # BG2 tube vs E5 tube on every step -> tube_compare.json
$M centre_exact.py      # item 2: exact centre orders + certified c0 -> centre_exact.json, centre_exact.log
$M thth_sym.py          # item 4: linearised E_thth (to Y''), D_x, sigma_TT(D_x = 0) -> thth_exprs.txt
$M stages_b.py          # engine Linear2 on the BG2 tube (built once, cached Lbg2_nsub2.pkl) + N3/N4 + C1 x=-0.01 -> stages.json
$M make_pass_b.py       # generates pass_b.py from build-3 qb2pass.py (7 asserted text substitutions)
./run_pass_b.sh         # N5 transport + N7 + route B + N2 on all 1440 inward substeps -> pass_ckpt.json, pass_steps.jsonl
$M final_b.py           # N6, C1 (5 points), C2, C6, C8, C8', N-a, N-c, gate, N7, N8 -> final.json
$M thth_eval.py         # item 4: C8 / real / N-b theta-theta residual at x = -0.01 -> thth.json
$M make_outputs_b.py    # -> result_qb2b.json (PQ2b-1..5), controls_qb2b.json (sha256)
```
Helpers: `qb2b.py` (in-process patches, below) and `stages_b_helpers.py`. `_p*.py` files are applied source edits
(transcript).

## What changed from build-3, and why (all in build-4 files)
1. **Inward tube from the BG2 sonic box.** The pinned `lin_tube_v2.build` is called unchanged, except:
   - V0 is the BG2 box (outward_C.pkl, ±3.79e-29, inside the E5 box);
   - its pickle is redirected to build-4;
   - its module constants `XBAR = arb('-4.6')` and `X_SWITCH = arb('-0.05')` are re-created at prec 300 with the same
     values. They had been built at import precision, 53 bits, radius about 3e-18 and 3.6e-16. A first build without
     this fix (stopped, own PIDs) showed the x-leg width jumping to 1e-17 at the switch.
   - W0 and xc stay the E5 boxes.
   - Result: same 720-step partition, every step overlaps the E5 tube, max setrad 5.2e-17 (was 1.0e-13), sonic-manifold
     coefficient radii 1e9 smaller (tube_compare.json).
2. **Homogeneous stage** rebuilt by `lin_stages_v2.Linear2` on the new tube: nsub 2, tol 1e-20, depth 0 on every step.
   It is cached with `b4speed/cacheL` ser_steps (exact ints).
3. **Series tails** (qb2b.P3). The RegSing majorant radius ν' is optimised in (|t|, ν) instead of the engine's
   min((ν+|t|)/2, 2|t|), which forces q ≥ 1/2. This is the same bound, valid for every admissible ν'. Centre forced
   tails fell from 5.8e-11 to 5e-25.
4. **Transport settings** (pass_b.py): series order 50 (was 40), remainder tolerance 1e-24 (was 1e-16), mean-value Bach
   trigger 1e-16 (was 1e-13).
   - Probes, kept in this file's history: order 60 everywhere gives 7.7 s per substep.
   - The mean-value trigger at 1e-12 gives near-sonic D_t bounds of 7e-11 and width growth of 3e-15 per step. Both were
     rejected; the final pass ran from k = 0 with one uniform setting.
5. **C1 at x = -0.01** is evaluated on the forced sonic series (the τ-steps have no check points), with s^HKA from the
   engine sonic RegSing at the same θ.

## Speed levers
- Lever 1: domain cut, justified in DOMAIN-LEMMA.md.
- Lever 3: lean precision; the settings above are the smallest that meet the targets in the probes.
- Lever 2 (float preview) was not needed.

## Run record (2026-11-18)
- **Run 1 is the record.** Files: final.json, pass_ckpt.json, stages.json, thth.json; copies are in run1/.
  - It used the plain homogeneous stage on the BG2 tube (cache Lbg2_nsub2.pkl).
  - qb2b.py's current `QB2B_RUN1=1` switch reproduces that configuration. qb2b_tight.py and the env switches were added
    after run 1.
  - C8 θθ: `QB2B_RUN1=1 QB2B_DATA=<build-4>/run1/ $M thth_eval.py` (uses only the run-1 cache, no rebuild).
- **Run 2 (STOPPED by budget, own PIDs).** A width budget at x_m showed the dominant source of the e width: the
  homogeneous step matrices built by the engine's ball recursion. Their order-0 radius is 1.3e-13 at depth, while the
  start ball is 5e-17.
  - The fix is qb2b_tight.py (T1: tight jets; T2: mean-value step tapes).
  - Linear2 rebuilt with it ran more than 40 min (the plain build: 5.7 min) and was stopped. A full rerun (rebuild plus
    the ~85 min pass) would have exceeded the 5 h phase-2 budget.
  - Log: stages_b_tight_STOPPED.log. NOT a result.
- Compute, phase 2: about 2h55. Single runs: the pass chunks ≤ 25 min each; the stopped rebuild 41 min.
