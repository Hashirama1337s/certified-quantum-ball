"""build-1 QB1: add instrument sha256 (all seat files), input-file sha256, D_new and null summaries to controls_qb1.json.
No computation.  usage (memcap 8): py -3.11 make_controls.py"""
import json, os, hashlib
HERE = os.path.dirname(os.path.abspath(__file__))


def sha(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()


c = json.load(open(os.path.join(HERE, 'controls_qb1.json')))
files = ['ajet.py', 'aser.py', 'charts.py', 'bg.py', 'ctl.py', 'qb1.py', 'dnew.py', 'null_bg.py', 'null_eval.py', 'make_controls.py']
c['instrument_sha256_final'] = {f: sha(os.path.join(HERE, f)) for f in files}
c['inputs_sha256'] = {p: sha(p) for p in ['<LAB>/engine/cert/lin_tube2_k1_3.pkl',
                                          '<LAB>/sideways/build-5/phase2/outward_C.pkl',
                                          '<LAB>/engine/cert/fields.pkl',
                                          '<LAB>/engine/cert/taylor.py',
                                          '<LAB>/engine/cert/count_b4v3.py']}
for fn, key in (('dnew_qb1.json', 'dnew'), ('null_qb1.json', 'null')):
    p = os.path.join(HERE, fn)
    if os.path.exists(p):
        d = json.load(open(p))
        c[key] = d.get('summary') if key == 'null' else dict(x_m=d.get('x_m'), values=d.get('values'),
                                                                 D_new_2_nonzero_certified=d.get('D_new_2_nonzero_certified'),
                                                                 float_crosscheck=d.get('float_lin_float_D_new_NONCOUNTING'))
        c[key + '_file_sha256'] = sha(p)
r = os.path.join(HERE, 'result_qb1.json')
if os.path.exists(r):
    c['result_qb1_sha256'] = sha(r)
json.dump(c, open(os.path.join(HERE, 'controls_qb1.json'), 'w'), indent=1)
print('controls_qb1.json updated')
