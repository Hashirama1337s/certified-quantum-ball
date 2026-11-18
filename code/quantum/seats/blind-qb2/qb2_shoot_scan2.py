import sys; sys.path.insert(0, '<LAB>/quantum/seats/blind-qb2')
from qb2_bg import *
import numpy as np
for w0v in np.linspace(0.09, 0.8, 30):
    xe, F, r = shoot(float(w0v))
    y = r.y[:, -1]
    d = rhs(r.t[-1], y)
    print('w0 %.5f  xend %.8f  F %.3e  v %.6f  dw %.3e dv %.3e  det %.3e' % (w0v, r.t[-1], F if F is not None else np.nan, y[3], d[2], d[3], detM(r.t[-1], y)), flush=True)
