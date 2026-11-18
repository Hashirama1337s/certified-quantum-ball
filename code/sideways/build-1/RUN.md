# RUN.md - build-1 (B7 SB0 + SB1): main-run commands for the lead

All commands from `<LAB>/sideways/build-1`, under memcap 8. No environment flags.
Pin first (sha256 in controls_sf.json): lin_core.py, polar_sys.py, eqs_l.py, eqs_l.pkl, bgnum.py, outward_rec.py,
outward_rec.pkl, s1_cert.py, l0_cert.py. Inputs read (not modified): engine/cert/lin_tube2_k1_3.pkl (f6bff9e2...),
engine/system.py, engine/eqs_gen.py, engine/cert/fields.pkl and the engine/cert modules imported by the instruments.

## Main runs (lead only)

1. BAR-SF-S1 (seal 2c30c7615f):
   `py -3.11 <MEMCAP> 8 py -3.11 s1_cert.py --main`
   writes `result_s1.json` (PS1-PS4, enclosure strings, pass flags, instrument sha256).
   Expected run time: about 1-3 min (dry run: 23.6 s total with a 40/720 tube-step slice; tube loop 1.5 s per 40 steps).
2. BAR-SF-L0 (seal 596c7146ff):
   `py -3.11 <MEMCAP> 8 py -3.11 l0_cert.py --main`
   writes `result_l0.json` (PL1-PL5, PL3b, kill-switch block, instrument sha256).
   Expected run time: under 1 min (full dry run 3.1 s).

Both are far under the 30 min cap; peak memory in every probe was at most 0.16 GB.

## Already run by the build seat (not main runs)

- `py -3.11 .../memcap.py 8 py -3.11 eqs_l.py` -> eqs_l.pkl, eqs_checks.json (69 s; 41 symbolic checks, all true)
- `py -3.11 .../memcap.py 8 py -3.11 outward_rec.py` -> outward_rec.pkl (363 s; 731 steps; m/r overlaps PE3)
- `s1_cert.py --controls` (24 s) and `l0_cert.py --controls` (3 s) -> controls_s1.json, controls_l0.json
- `s1_cert.py --dryrun`, `l0_cert.py --dryrun`: full code path, values suppressed, nothing written
- `t_gauge.py`, `t_fm1.py`, `t_l0path.py`, then `controls_sf.py` -> controls_sf.json

To rebuild everything from scratch (about 8 min): eqs_l.py, outward_rec.py, then the two main runs.
