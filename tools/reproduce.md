# Reproduce: re-run the certificates

Requirements: Python 3.11, python-flint 0.9.0, sympy, mpmath, numpy (scipy for the blind seat only; matplotlib 3.11 for
the figures only). Every step is a single process; the largest peak in the lab was about 1.3 GB. The lab ran every step
under an 8 GB memory cap.

**Status of this file.** The commands and lab times below come from the lab's run records (`RUN.md` in each seat folder,
the shipped logs). On 18 November 2026 the repository was copied into fresh folders and replayed with the documented commands only (tools/cleanroom.sh). Quick tier: all 16 steps IDENTICAL to the shipped files (canonical content), ALL CHECKS PASS, 20.1 min. BG2 tier: the sonic box (44.0 min) and the outward record (10.1 min) IDENTICAL; the gate step first differed only in two wall-clock fields (time_s, path_info.setup_s), the comparator was fixed to ignore them (PATCHES.md), and the gate step re-run from the shipped outward record is IDENTICAL, ALL CHECKS PASS. The full tier (about 7 hours: the transported QB2b record) was not replayed in a clean room; the quick tier replays the matching, controls, nulls and readout on the shipped transported record.

## Commands

From the repository root:

```
python tools/check_pins.py                   # every pinned file against its ledger pin; every file against MANIFEST.sha256
python tools/reproduce.py --tier quick       # about 30 minutes on one core (lab times)
python tools/reproduce.py --tier full        # about 7 hours on one core (lab times)
python tools/reproduce.py --tier bg2         # about 70 minutes: the BG2 sonic box, outward record and gate
python tools/reproduce.py --list             # every step, its tier, lab time and command
sh tools/cleanroom.sh quick                  # fresh copy of the repository, then the quick tier (or: full, bg2)
```

Expected last line: `[reproduce] ALL CHECKS PASS`. The runner works in `./work` (or `--work DIR`, which must be empty
or absent); `--python "CMD"` picks another interpreter; `--only NAME ...` runs named steps (the data files they need are
prepared first); `--regen` regenerates data files instead of using the shipped ones. A JSON report is written to
`work/_out/reproduce_report.json`; each step's output is in `work/_out/<step>.log` and the launcher's substitution
lines in `work/_out/<step>.launcher.log`.

## Steps

"Canonical" = compared with the shipped file by `tools/canon.py`: elapsed times, peak memory, run-time keys
(`time_s`, `seconds`, ...) and the project folder are removed; every other character (each enclosure, count, sign and
sha256 provenance field) must agree. Data files are compared byte for byte.

### Quick tier

| Step | Claim | Command (lab path; run from its seat folder) | Lab time | Compared with |
|---|---|---|---|---|
| pins | all | `tools/check_pins.py` | seconds | `tools/pins.json`, `tools/path_patch.json`, `MANIFEST.sha256` |
| qb1_main | QB1 PQ1, PQ3, PQ4, grid sign map | `quantum/seats/build-1/qb1.py main` | 2.3 min | `result_qb1.json`, `main.log` |
| qb1b_sm_w | QB1b PQ2b, deep-centre region | `quantum/seats/build-2/signmap2.py w -12 -11/2 13 5 sm_w.json` | 1.9 min | `sm_w.json`, `sm_w.log` |
| qb1b_seam, qb1b_null, qb1b_controls | QB1b PQ8, PQ6b, controls | `quantum/seats/build-2/tasks2.py seam / null / controls` | 0.2 min each | `pq8.json`, `pq6b.json`, `controls_tm.json`, `t2_*.log` |
| qb1b_finalize | QB1b verdict file (from the shipped `sm_*.json`) | `quantum/seats/build-2/finalize.py` | 0.3 min | `result_qb1b.json` |
| qb2_bachsym, qb2_bach_u, qb2_n1_rows | QB2 exact source and rows | `quantum/seats/build-3/bachsym.py`, `bach_u.py`, `n1_rows.py` | 1.5, 0.5, 0.4 min | `bach_exprs.txt`, `bach_u.txt`, `n1_rows.txt` and logs |
| qb2_n0check | QB2 PQ2-1 cross-check with QB1 | `quantum/seats/build-3/n0check.py` | about 1 min | `n0check.json` |
| qb2b_centre | QB2b PQ2b-2 exact centre | `quantum/seats/build-4/centre_exact.py` | 3.4 min | `centre_exact.json` |
| qb2b_thth_sym | θθ linearisation (exact) | `quantum/seats/build-4/thth_sym.py` | about 2 min | `thth_exprs.txt` |
| qb2b_stages | Lbg2 step cache + N3/N4 (needed by the QB2c replay) | `quantum/seats/build-4/stages_b.py` | 5.7 min | `stages.json`, `Lbg2_nsub2.pkl` (sha256) |
| qb2c_final | QB2c PQ2c-3 headline replay, PQ2c-4 replication row | `quantum/seats/build-5/final_c.py` | 2 min | `final_c.json`, `conv_c.json` |
| qb2c_thth | QB2c PQ2c-2 θθ at five points, N-b | `quantum/seats/build-5/thth_c.py` | about 4 min | `thth_c.json` |
| qb2c_outputs | QB2c verdict file | `quantum/seats/build-5/make_outputs_c.py` | seconds | `result_qb2c.json` |

The quick tier replays the matching, controls, nulls and readout on the shipped transported record of QB2b run 1
(`code/quantum/seats/build-4/run1/`, `pass_ckpt.json`, `stages.json`); the full tier recomputes that record.

### Full tier (adds)

| Step | Claim | Command | Lab time |
|---|---|---|---|
| gen_eqs ... gen_tube13, gen_cacheL | engine data (as certified-evans-coleman) | `engine/derive.py`, `engine/cert/fields.py`, `lin_sym.py`, `lin_mc.py`, `lin_tube_v2.py 1/3`, `b4speed/cacheL.py build 1/3 2 1e-20 b4speed/L13_nsub2.pkl` | about 22 min |
| gen_qb2_exprs | QB2 per-unit-β rows | `quantum/seats/build-3/qb2sym2.py` | about 1 min |
| gen_bg2_tube | QB2b inward tube from the BG2 box | `quantum/seats/build-4/tube_bg2.py` | 8.6 min |
| gen_null_in, gen_null_out | QB1 PQ6 off-solution background | `quantum/seats/build-1/null_bg.py in / out` | 8.5, 13.3 min |
| qb1_controls, qb1_dnew, qb1_null_eval, qb1_make_controls | QB1 PQ5, PQ7, PQ6 | `quantum/seats/build-1/qb1.py controls`, `dnew.py`, `null_eval.py`, `make_controls.py` | about 2 min |
| qb1b_sm_* (8 regions) | QB1b PQ2b sign map | `quantum/seats/build-2/signmap2.py ...` (the 9 region lines of `build-2/RUN.md`) | about 1.9 h total |
| qb1b_dnew2 | QB1b PQ7b | `quantum/seats/build-2/dnew2.py` | 0.5 min |
| qb2_stages, qb2_pass, qb2_c2check, qb2_final, qb2_outputs | QB2 (f2bcfadb3e) | `quantum/seats/build-3/qb2stages.py`; `qb2pass.py k k+240` chunks (exit 3 = resume); `c2check.py`; `qb2final.py`; `make_outputs.py` | about 85 min |
| qb2b_tube, qb2b_pass, qb2b_final, qb2b_thth, qb2b_outputs | QB2b (373f306e9f) | `tube_compare.py`; `stages_b.py`, `pass_b.py k k+240` chunks, `final_b.py` and `thth_eval.py`, all with `QB2B_RUN1=1` (the run-1 record; the shipped `run_pass_b.sh` sets `QB2B_NOPATCH=1`, the configuration of the stopped run 2); `make_outputs_b.py` | about 95 min |
| qb2c_deep | QB2c PQ2c-1 evidence | `quantum/seats/build-5/deep_ratio.py 480` | about 20 min |

### BG2 tier (ed61d09268; not in any published record)

| Step | Command (from `sideways/build-5/phase2`) | Lab time | Compared with |
|---|---|---|---|
| bg2_box | `b1_box.py C 256 50 24 --nonewton` (LB_NGRID = 24 as sealed) | 49 min | `b1_C.json` |
| bg2_outward | `b2_outward.py C b1_C.json 256 50` | 15 min | `b2_C.json`; `outward_C.pkl` byte for byte |
| bg2_gate | `b3_gate.py C outward_C.pkl` | 4.3 min | `b3_C.json` |

The gate imports the lab's shared count engine `gravity_engine/argument_count.py`, shipped in `code/_pocketlab/` (the
launcher points the lab's tool folder there). The collector `make_bg2.py C C C nullshift` and the controls of
`sideways/build-5/RUN-BG2.md` (null box, null-shift outward leg, E5 control leg, old-record gate) are shipped but not run by
the runner; their recorded outputs are in `code/sideways/build-5/phase2/`. Disclosed in RUN-BG2.md: `b1_box.py` gained its
`--ladder` option after runs C and NULL had started (default ladder unchanged).

## Data files (pickles)

Loading a pickle can execute code. Shipped pickles (`tools/layout.json`, `shipped: true`) are copied into the work folder
only after their sha256 matches; the full tier regenerates them. Not shipped (regenerated by the runner and checked
against the lab file's sha256): `b4speed/L13_nsub2.pkl` (41 MB, ledger pin under b1b0ee1cde), `build-4/Lbg2_nsub2.pkl`
(44 MB), `build-1/null_in.pkl`, `build-1/null_out.pkl`. Not shipped and not used by any step: the blind seat's pickles
and arrays (its scripts regenerate them).

## What the runner does (and why)

- **Work folder = the lab tree.** The instruments find each other and their data by paths under one project folder;
  `code/` is that tree, and `code/_pocketlab/` holds the two files used from the lab's tool folder.
- **Original bytes.** The shipped copies are path-patched (PATCHES.md); the runner restores the originals in the work
  folder first, so every file has exactly the bytes the ledger pinned.
- **Lab folder literals.** `tools/run_sealed.py` compiles every module loaded from the work folder with only the lab-root
  literal (and the lab's tool-folder literal) pointed at the work folder; files on disk stay unchanged, so every recorded
  sha256 is the original. Opens of lab paths are redirected into the work folder, and an audit hook refuses any open
  outside it.
- **Environment.** `LB_*`, `FW*` and `QB2B_*` flags are cleared before each step; only the sealed ones are set
  (`LB_NGRID = 24` for BG2; `QB2B_RUN1 = 1` for every QB2b step, which reproduces the run-1 record; the
  build-5 scripts set `QB2B_RUN1` themselves).
- **Known limits of this runner (not yet exercised).** Some instruments append to their logs or record the folder they
  ran in; such logs are not compared (the JSON products are). The order in which an instrument with two `bg.py`
  candidates (engine and build-1) imports them follows `sys.path` exactly as in the lab.

## Running one step by hand

```
python tools/reproduce.py --only qb2_bach_u --work run   # sets up run/ and runs one step
cd run/quantum/seats/build-3                              # then by hand, from the seat folder (absolute paths):
python <repo>/tools/run_sealed.py --root <repo>/run <repo>/run/quantum/seats/build-3/bach_u.py
```

## Printed numbers (optional; not run by tools/reproduce.py)

`python tools/fill_balls.py --check` reads every ball printed in paper.md, README.md, CLAIMS.md, SCORECARD.md and
CITATION.cff, finds the stored ball it stands for in the shipped JSON files and logs, and reports whether the print is
the stored string verbatim or contains it (outward rounding). Expected last line: `... 0 not verbatim-or-outward`.

## Figures and scorecard (optional; not run by tools/reproduce.py)

`python figures/make_figures.py .` redraws `figures/*.svg` from the shipped JSON result files only (seconds; two runs give
byte-identical SVGs). `python tools/make_scorecard.py` rebuilds `scorecard.svg` and `SCORECARD.md`, checked against the
ledger export `results/ledger-quantum.json`.

## Smoke steps run during assembly

- `--only qb2_bach_u qb1b_seam`: both IDENTICAL (canonical), 0.1 min.
- `--only qb2c_final`: qb2b_stages regenerated `Lbg2_nsub2.pkl` byte for byte and `stages.json` canonically identical (5.0 min), then
  final_c.py reproduced `final_c.json` and `conv_c.json` canonically identical (0.2 min): ALL CHECKS PASS.
