"""make_bg2.py - collects phase-2 outputs into ../result_bg2.json and ../controls_bg2.json (no computation)."""
import json, os, hashlib, sys
HERE = os.path.dirname(os.path.abspath(__file__)); UP = os.path.dirname(HERE)


def J(name):
    p = os.path.join(HERE, name)
    return json.load(open(p)) if os.path.exists(p) else None


def sha(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()


B1TAG, B2TAG, B3TAG, NULLB2 = sys.argv[1:5]
b1 = J('b1_%s.json' % B1TAG); b2 = J('b2_%s.json' % B2TAG); b3 = J('b3_%s.json' % B3TAG); b3old = J('b3_old.json')
nul = J('b1_NULL.json'); nb2 = J('b2_%s.json' % NULLB2) if NULLB2 != '-' else None
B1 = dict(b1['B1'], prec=b1['prec'], Kt=b1['Kt'], Ks=b1['Ks'], nu_s=b1['nu_s'], ladder=b1['ladder'], tried=b1.get('ladder_tried', b1.get('ladder')),
          centre=b1['newton'], time_s=b1.get('time_s'),
          earlier_runs=dict(A_Kt36=J('b1_A_Kt36.json').get('newton'), B_Ks16=J('b1_B_Ks16.json').get('newton')),
          replicate_E_newton_centred=dict(V0=J('b1_E.json')['B1']['V0'], width=J('b1_E.json')['B1']['width'], inside_E5=J('b1_E.json')['B1']['inside_E5'],
                                          overlaps_C=None))
from flint import arb as _arb, ctx as _ctx
_ctx.prec = 256
B1['replicate_E_newton_centred']['overlaps_C'] = bool(_arb(B1['V0']).overlaps(_arb(B1['replicate_E_newton_centred']['V0'])))
B2 = dict(record='phase2/outward_%s.pkl' % B2TAG, record_sha256=sha(os.path.join(HERE, 'outward_%s.pkl' % B2TAG)),
          prec=b2['prec'], Kt=b2['Kt'], steps=b2['steps'], lemma=b2['lemma'], ok_D=b2['ok_D'], ok_V=b2['ok_V'],
          near_cone_max_rel_halfwidth=b2['near_cone_max_rel'], near_cone_steps=len(b2['near_cone_steps']),
          near_cone_per_step=b2['near_cone_steps'],
          x_lc=b2['x_lc'], x_lc_width=b2['x_lc_width'], root_step=b2['root_step'],
          state_at_xlc=b2['state_at_xlc'], state_at_xlc_rel=b2['state_at_xlc_rel'],
          frobenius_Ylc_rel_per_vertex=(b3 or {}).get('Ylc_rel_per_vertex'),
          overlap_old=b2['overlap_old'], time_s=b2['time_s'])
B3 = dict(outward=(b3 or {}).get('outward'), tried=(b3 or {}).get('tried'), gate_tol='1e-2', note='probe only, no count',
          D_old_vs_new_k0=J('p_dover_k0.json'))
anchors = dict(m_over_r_t0=b2.get('m_over_r_t0'), m_over_r_old=b2.get('m_over_r_old'), m_over_r_overlaps_old=b2.get('m_over_r_overlaps_old'),
               x_lc_new=b2['x_lc'], x_lc_old=b2['overlap_old']['x_lc_old'], x_lc_inside_old=b2['overlap_old']['x_lc_inside_old'],
               old_pieces_gate=dict(tried=[(t['match_k'], t['max_rel']) for t in b3old['tried']] if b3old else None,
                                    reproduces_1_0=bool(b3old and all(t['max_rel'] == 1.0 for t in b3old['tried']))),
               e5ctl_record=J('b2_e5ctl.json') and dict(x_lc=J('b2_e5ctl.json')['x_lc'], near_cone_max_rel=J('b2_e5ctl.json')['near_cone_max_rel'],
                                                        note='same rebuild code from the OLD E5 box: reproduces the old widths'))
null = dict(shift='1e-22', krawczyk=(nul or {}).get('null'),
            rebuild_from_shifted=(nb2 and dict(x_lc=nb2['x_lc'], x_lc_overlaps_B2=None, m_over_r=nb2.get('m_over_r_t0'))))
if nb2:
    from flint import arb
    a, b = arb(nb2['x_lc']), arb(b2['x_lc'])
    null['rebuild_from_shifted']['x_lc_overlaps_B2'] = bool(a.overlaps(b))
    null['rebuild_from_shifted']['x_lc_separation'] = (a - b).str(5)
res = dict(spec='build-5 phase 2 (DESIGN-SF-BG re-targeted: outward leg near the light cone)', B1=B1, B2=B2, B3=B3, anchors=anchors, null=null)
json.dump(res, open(os.path.join(UP, 'result_bg2.json'), 'w'), indent=1)
files = ['b1_box.py', 'b2_outward.py', 'b3_gate.py', 'make_bg2.py', 'p_xlc.py', 'p_eps.py', 'p_eps0.py', 'p_dover.py']
ctl = dict(instrument_sha256={'phase2/' + f: sha(os.path.join(HERE, f)) for f in files if os.path.exists(os.path.join(HERE, f))},
           imported_read_only_sha256={p: sha(p) for p in ['<LAB>/sideways/build-1/outward_rec.py', '<LAB>/sideways/build-1/l0_cert.py',
                                                         '<LAB>/sideways/build-4/l2_lib.py', '<LAB>/sideways/build-4/l2b_ext.py',
                                                         '<LAB>/sideways/build-4/l2_cert.py', '<LAB>/engine/cert/driver_v2.py',
                                                         '<LAB>/engine/cert/driver_v31.py', '<LAB>/engine/cert/driver_v32.py',
                                                         '<LAB>/engine/cert/driver2.py', '<LAB>/engine/cert/lin_tube2_k1_3.pkl',
                                                         '<LAB>/sideways/build-1/outward_rec.pkl']},
           records_sha256={'phase2/outward_%s.pkl' % t: sha(os.path.join(HERE, 'outward_%s.pkl' % t)) for t in (B2TAG, 'e5ctl', NULLB2) if t != '-' and os.path.exists(os.path.join(HERE, 'outward_%s.pkl' % t))},
           controls=dict(old_record_gate_reproduces_1_0=anchors['old_pieces_gate'], e5ctl=anchors['e5ctl_record'],
                         D_point_values_old_new_overlap_k0=(J('p_dover_k0.json') or {}).get('all_overlap'),
                         B1_replicate_overlap=B1['replicate_E_newton_centred']['overlaps_C'], b2_overlap_old=b2['overlap_old']))
json.dump(ctl, open(os.path.join(UP, 'controls_bg2.json'), 'w'), indent=1)
print('wrote result_bg2.json, controls_bg2.json')
