import sys; sys.path.insert(0, '<LAB>/quantum/seats/blind-qb2')
from qb2_bg import *
import numpy as np
for w0v in np.linspace(15, 21, 25):
    xe, F, r = shoot(float(w0v))
    y = r.y[:, -1]
    print('w0 %.4f xend %.6f F %.4e v %.5f det %.3e st %d' % (w0v, r.t[-1], F if F is not None else np.nan, y[3], detM(r.t[-1], y), r.status), flush=True)
