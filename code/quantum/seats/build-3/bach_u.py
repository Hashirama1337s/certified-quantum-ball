"""build-3: b_ab of gt with exp(l0) -> 1/u (u = 1/N); writes bach_u.txt (srepr, one per line). Exact sympy."""
import sympy as sp
Bx = [eval(l, sp.__dict__.copy()) for l in open('<LAB>/quantum/seats/build-3/bach_exprs.txt')]
l0 = sp.Symbol('l0'); uu = sp.Symbol('u')
out = []
for nm, e in zip(('nn', 'nr', 'rr', 'T'), Bx):
    e2 = e.subs(l0, -sp.log(uu))
    e2 = sp.cancel(sp.together(sp.expand(sp.powsimp(sp.expand_log(e2, force=True), force=True))))
    assert l0 not in e2.free_symbols and not e2.has(sp.log) and not e2.has(sp.exp), (nm, e2.free_symbols)
    print(nm, 'free', sorted(map(str, e2.free_symbols)), ' denom', sp.factor(sp.denom(e2)))
    out.append(sp.srepr(e2))
    # numeric identity check at a random rational point against the original
    import random
    pt = {sp.Symbol('A%d' % i): sp.Rational(random.randint(5, 30), 7) for i in range(5)}
    pt.update({sp.Symbol('l%d' % i): sp.Rational(random.randint(-9, 9), 5) for i in range(1, 5)})
    pt[l0] = sp.Rational(3, 7)
    v1 = sp.N(e.subs(pt), 40); v2 = sp.N(e2.subs(pt).subs(uu, sp.exp(-sp.Rational(3, 7))), 40)
    print('   check', abs(v1 - v2) < 1e-30)
open('<LAB>/quantum/seats/build-3/bach_u.txt', 'w').write('\n'.join(out) + '\n')
