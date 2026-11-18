import sys; sys.path.insert(0, '<LAB>/quantum/seats/blind-qb2')
from qb2_bg import *
import numpy as np
for w0v in np.geomspace(0.05, 20, 25):
    xe, F, r = shoot(float(w0v))
    print('w0 %.5g  xsonic %s  F %s  status %s xend %.4g' % (w0v, xe, F, r.status, r.t[-1]), flush=True)
