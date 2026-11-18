# RUN-BG2.md - build-5 phase 2 (rebuild of the outward leg near the light cone)

All scripts in `<LAB>/sideways/build-5/phase2/`, run from that folder under memcap 8, no environment flags
(LB_NGRID = 24 as sealed). New files only; build-1, build-4 and engine/cert files are imported read-only.

## B1 refined V0 box (E5 Krawczyk / R8 method, parameters changed only)
    py -3.11 <MEMCAP> 8 py -3.11 b1_box.py C 256 50 24 --nonewton        -> b1_C.json
Probes that fixed the parameters (kept): b1_A_Kt36.* (prec 200, K_t 36, K_s 16: point residual [+/- 3.9e-21],
Lohner-truncation floor, stopped by PID after the point run), b1_B_Ks16.* (prec 256, K_t 50, K_s 16: point residual
[+/- 3.9e-30], Newton step [+/- 7.9e-30], but error bound E = 1.67e-24 from the sonic tail eps 1.9e-32 x growth; stopped
by PID before its rungs), p_eps.py / p_eps0.py (eps vs K_s: 1.9e-32 / 2.0e-40 / 2.7e-48 for K_s 16 / 20 / 24).

## Null (V0 shifted by 1e-22, Krawczyk content null at the B1 rung radius)
    py -3.11 <MEMCAP> 8 py -3.11 b1_box.py NULL 256 50 24 --nonewton --null 1e-27   -> b1_NULL.json
    py -3.11 <MEMCAP> 8 py -3.11 b2_outward.py nullshift b1_C.json 256 50 --v0shift 1e-22

## B2 outward leg + x_lc + state at x_lc
    py -3.11 <MEMCAP> 8 py -3.11 b2_outward.py C b1_C.json 256 50            -> outward_C.pkl, b2_C.json
Control (same code, OLD E5 box; must reproduce the old widths):
    py -3.11 <MEMCAP> 8 py -3.11 b2_outward.py e5ctl e5box.json 256 50        -> outward_e5ctl.pkl, b2_e5ctl.json

## B3 rigorous gate (probe only, no count) and the old-record control
    py -3.11 <MEMCAP> 8 py -3.11 b3_gate.py C outward_C.pkl                   -> b3_C.json
    py -3.11 <MEMCAP> 8 py -3.11 b3_gate.py old <LAB>/sideways/build-1/outward_rec.pkl  -> b3_old.json

## Collect
    py -3.11 <MEMCAP> 8 py -3.11 make_bg2.py C C C nullshift   -> ../result_bg2.json, ../controls_bg2.json

## Added during the run
    py -3.11 <MEMCAP> 8 py -3.11 b1_box.py E 256 50 24 --ladder 1e-25,1e-24   # replicate: Newton-centred
    py -3.11 <MEMCAP> 8 py -3.11 p_dover.py 0     # D at the 16 vertices, old vs new record: must overlap
Killed by PID (own runs, superseded): b1_D (Newton + ladder 1e-27..1e-25; replaced by E after C showed the K-radius floor
2.6e-26 in the W0 component, so rungs below 1e-25 cannot contract).
Disclosure: b1_box.py gained the --ladder option after runs C and NULL had started (default ladder unchanged); the pinned
hash in controls_bg2.json is the final file.  content_null (driver_v31.py:367-386, pinned) prints K/X with 20 digits, so
the radii shown in b1_NULL.json are display-limited; the outcome fields (contracts, inside, miss_over_r) are exact.
