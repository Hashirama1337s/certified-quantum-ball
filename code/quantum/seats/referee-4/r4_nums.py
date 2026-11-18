"""referee-4: numeric rechecks (arb). c0 from the closed form on the E5 boxes; N-a arithmetic from result_qb2b.json."""
from flint import arb, ctx
ctx.prec = 200
W0 = arb('8.90132327775564352') + arb(0, arb('7.92e-18'))
xc = arb('0.212365646758945385') + arb(0, arb('7.95e-19'))
c0 = -arb(2)/3*W0*(8*W0 - 3)*(-4*xc).exp()
print('c0 closed form =', c0.str(22))
print('b_nr_5 e^{-4xc}... ratio b_nr5/b_nn4 = 4/3 -> hat_b_nr centre limit (x^5 coefficient) =', (c0*4/3).str(18))
Js = arb('58.313767936470775906') + arb(0, arb('1.65e-13'))
Jc = arb('-0.0058226690727984285785') + arb(0, arb('4.51e-18'))
cV = arb('0.19288155165634') + arb(0, arb('1.41e-15'))
print('design |J_s|/|c_V| =', (Js/cV).str(15))
print('(J_c - J_s)/c_V    =', ((Jc - Js)/cV).str(15))
print('J_c/c_V share       =', (Jc/cV).str(10))
