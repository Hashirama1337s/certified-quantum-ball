"""build-1 QB1 item 6: run the unchanged QB1 evaluators (qb1.py / bg.py / charts.py) on the OFF-SOLUTION background
null_in.pkl + null_out.pkl (null_bg.py, V0* = E5 upper edge + 1e-15) and on the certified background at the same
points.  Reports: Einstein residual excludes 0 anywhere?  Route A and Route B disjoint anywhere?  Null b vs certified b
disjoint where?
usage (memcap 8): py -3.11 null_eval.py -> null_qb1.json
"""
import sys, os, time, json, pickle, hashlib
from flint import arb, ctx
ctx.cap = 64
import bg
import qb1
from aser import Piece
HERE = os.path.dirname(os.path.abspath(__file__))


class NullBackground(bg.Background):
    def __init__(self):
        self.tb = pickle.load(open(os.path.join(HERE, 'null_in.pkl'), 'rb'))
        self.ow = pickle.load(open(os.path.join(HERE, 'null_out.pkl'), 'rb'))
        self._tube_pieces = None

    def sm_piece(self, which):
        src = self.ow if which == 'bg2' else self.tb
        sm = src['sm']
        Y = [[bg.deser(p) for p in sm['Y'][i]] for i in range(3)]
        X = [bg.deser(p) for p in sm['X']]
        eps = bg.deser(sm['eps']).abs_upper(); epsX = bg.deser(sm['epsX']).abs_upper()
        nus = sm['nu_s'] if which == 'bg2' else src['nu_s']
        nu = arb(nus).lower()*(1 - arb('1e-9'))
        return Piece(Y + [X], [eps, eps, eps, epsX], arb(nu), 'sch')


def row(o, x0, tag):
    E = [o['E_nn'], o['E_nr'], o['E_rr'], o['E_T']]
    return dict(tag=tag, x=qb1.S(x0), b={c: qb1.S(o['b_' + c]) for c in qb1.COMP}, a={c: qb1.S(o['a_' + c]) for c in qb1.COMP},
                AB_disjoint_any=any(not o['b_' + c].overlaps(o['a_' + c]) for c in qb1.COMP),
                einstein={k: qb1.S(e) for k, e in zip(('nn', 'nr', 'rr', 'TT'), E)},
                einstein_excludes_0_any=any(not e.contains(0) for e in E), R=qb1.S(o['R']), _o=o)


def main():
    T0 = time.time()
    log = lambda *a: print('[%7.1fs]' % (time.time() - T0), *a, flush=True)
    old = ctx.prec
    try:
        ctx.prec = 256
        Nb = NullBackground(); Cb = bg.Background()
        pts = []
        for xt in ('-4.5', '-3', '-2', '-1.5', '-1', '-0.84', '-0.5', '-0.2489'):
            pts.append(('tube', arb(xt)))
        pts.append(('sonic', arb(0)))
        for xt in ('0.2489', '0.5', '1', '2', '3', '6', '9', '12'):
            pts.append(('outw', arb(xt)))
        rows = []
        for which, x in pts:
            out = {}
            for nm, B in (('null', Nb), ('cert', Cb)):
                if which == 'sonic':
                    x0, o = qb1.ev_sonic(B, 'bg2', x)
                    x1, o1 = qb1.ev_sonic(B, 'tube', x)
                    out[nm + '_sm_in'] = row(o1, x1, nm + ' sonic (inward-run sm)')
                else:
                    x0, o = qb1.ev_record_point(B, which, x, 'polar')
                out[nm] = row(o, x0, '%s %s x=%s' % (nm, which, x.str(6)))
            n, c = out['null']['_o'], out['cert']['_o']
            out['null_vs_cert_disjoint'] = {cc: bool(not n['b_' + cc].overlaps(c['b_' + cc])) for cc in qb1.COMP}
            out['null_minus_cert_b_uu'] = qb1.S(n['b_uu'] - c['b_uu'])
            for k in list(out):
                if isinstance(out[k], dict): out[k].pop('_o', None)
            rows.append(out)
            log('%s x=%s: null Eres excl0 %s, null A/B disjoint %s, null-vs-cert disjoint %s, d b_uu %s' % (
                which, x.str(6), out['null']['einstein_excludes_0_any'], out['null']['AB_disjoint_any'],
                [k for k, v in out['null_vs_cert_disjoint'].items() if v], (n['b_uu'] - c['b_uu']).str(5)))
        summ = dict(V0_star='0.1124394013880983362 + 3.31e-20 + 1e-15 (radius 1e-30)',
                    einstein_residual_excludes_0_anywhere=any(r['null']['einstein_excludes_0_any'] for r in rows),
                    routeA_routeB_disjoint_anywhere=any(r['null']['AB_disjoint_any'] for r in rows),
                    null_b_disjoint_from_certified_at=[r['null']['x'][0][:10] for r in rows if any(r['null_vs_cert_disjoint'].values())],
                    null_record_flags=dict(outward_ok_D=Nb.ow.get('ok_D'), outward_ok_V=Nb.ow.get('ok_V'), outward_lemma=Nb.ow.get('lemma')))
        res = dict(summary=summ, rows=rows, sha256={f: hashlib.sha256(open(os.path.join(HERE, f), 'rb').read()).hexdigest()
                                                    for f in ('null_bg.py', 'null_eval.py', 'null_in.pkl', 'null_out.pkl')},
                   time_s=round(time.time() - T0, 1))
        json.dump(res, open(os.path.join(HERE, 'null_qb1.json'), 'w'), indent=1)
        log('summary', summ)
    finally:
        ctx.prec = old


if __name__ == '__main__':
    main()
