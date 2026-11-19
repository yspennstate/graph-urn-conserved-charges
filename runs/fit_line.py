"""Fits to the line runs of Section 5: log sd = b - a log log t (free exponent), and the heuristic's form
sd = C (log t + tau0)^(-1/4) with tau0 scanned and C fitted, both by least squares in log sd. Two ranges: from the
row t = 169.4 on, and t >= 10^4. Run from this folder: python fit_line.py"""
import os
import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
for lo in (169.0, 1e4):
    print("rows with t >= %g" % lo)
    for s in (1, 2, 3):
        R = np.loadtxt(os.path.join(here, "line_L10001_s%d.txt" % s), comments="#")
        R = R[R[:, 0] >= lo]
        L, y = np.log(R[:, 0]), np.log(R[:, 2])
        X = np.vstack([np.ones_like(L), np.log(L)]).T
        k = np.linalg.lstsq(X, y, rcond=None)[0]
        rms_free = np.sqrt(np.mean((y - X @ k) ** 2))
        best = None
        for tau0 in np.linspace(-4.0, 40.0, 4401):
            z = y + 0.25 * np.log(L + tau0)
            e = np.sqrt(np.mean((z - z.mean()) ** 2))
            if best is None or e < best[0]:
                best = (e, tau0)
        print("  seed %d: first row t = %.1f, %d rows | free exponent a = %.4f, rms %.5f | exponent 1/4: tau0 = %.2f, rms %.5f"
              % (s, R[0, 0], len(R), -k[1], rms_free, best[1], best[0]))
