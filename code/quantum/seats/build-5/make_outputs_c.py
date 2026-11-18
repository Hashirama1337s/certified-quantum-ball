"""build-5: result_qb2c.json (PQ2c-1..PQ2c-4) and controls_qb2c.json (+ sha256).  usage: py -3.11 make_outputs_c.py"""
import json, hashlib, os, glob
H = '<LAB>/quantum/seats/build-5/'
B4 = '<LAB>/quantum/seats/build-4/'
J = lambda f: json.load(open(H + f))
F = J('final_c.json'); CV = J('conv_c.json'); TH = J('thth_c.json'); DR = J('deep_ratio.json')
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()
rows = DR['rows']
rat0 = sorted(r['tight0']/r['plain0'] for r in rows)
deep = [r for r in rows if r['k'] >= 1300]
res = {
 'PQ2c-1': dict(
    verdict_by_build='NOT RUN AS SEALED (tight homogeneous stage + full pass): probe shows T1+T2 cannot change the deep '
                     'step-matrix radii that set the e width; the e below is the QB2b run-1 record replayed',
    probe=dict(file='probe_tight.log', tight_sec_per_tube_step='4.5-7.5', plain_sec_per_tube_step='0.3-0.5',
               full_tight_rebuild_ETA_min='~75-90 (720 tube steps)', plus_pass_ETA_min='~85'),
    deep_ratio=dict(file='deep_ratio.json', substeps=len(rows), k_range=[rows[0]['k'], rows[-1]['k']],
                    order0_radius_ratio_tight_over_plain=dict(min=rat0[0], median=rat0[len(rat0)//2], max=rat0[-1]),
                    k_ge_1300_ratio_max=max(r['tight0']/r['plain0'] for r in deep) if deep else None,
                    tight_jet_overlap_assertion='never fired (assert in build-4 qb2b_tight.py T1)'),
    cause='floor = the pinned tube record step balls (setrad 2.5e-18 at tube step 600 -> 4.7e-17 at 715) times the genuine '
          'state sensitivity of the HKA step matrix (order-0 radius 8e-14 .. 3e-13 at depth), not interval dependency',
    N6=F['N6'], C8=F['C8'], N_a=F['N-a'], N_c=F['N-c']),
 'PQ2c-2': TH,
 'PQ2c-3': dict(CV['PQ2c3'], N8=F['N8'], C6=F['C6'], note='replay of the final stage on the QB2b run-1 transported record'),
 'PQ2c-4': CV,
}
ctl = dict(C1=F.get('C1'), C2=F.get('C2'), C6=F.get('C6'), C8=F.get('C8'), C8prime=F.get("C8'"), N_a=F.get('N-a'),
           N_c=F.get('N-c'), N_b_thth=TH, gate=F.get('gate_N1_C1_C8'),
           build5_sha256={os.path.basename(p): sha(p) for p in sorted(glob.glob(H + '*.py') + glob.glob(H + '*.md') + glob.glob(H + '*.log'))},
           build4_read_sha256={f: sha(B4 + f) for f in ('qb2b.py', 'qb2b_tight.py', 'stages_b_helpers.py', 'thth_exprs.txt',
                                                        'lin_tube_bg2_k1_3.pkl', 'Lbg2_nsub2.pkl', 'final.json', 'stages.json', 'pass_ckpt.json')},
           blind_sha256=sha('<LAB>/quantum/seats/blind-qb2/blind_qb2_result.json'))
json.dump(res, open(H + 'result_qb2c.json', 'w'), indent=1, default=str)
json.dump(ctl, open(H + 'controls_qb2c.json', 'w'), indent=1, default=str)
print('wrote result_qb2c.json, controls_qb2c.json; deep ratio min/median/max', rat0[0], rat0[len(rat0)//2], rat0[-1])
