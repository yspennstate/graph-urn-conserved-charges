"""Summary of the small-graph runs of Section 5 (code/warm_g.c, 2000 runs each): per time T, the standard
deviation over runs of each X_e and the mean and maximum over runs of max_v |S_v/t - 2|. Run from this folder."""
import os
import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
for g in ("K4", "triangle"):
    for T in (1000, 10000, 100000):
        rows = [l.split() for l in open(os.path.join(here, "%s_T%d.txt" % (g, T)))]
        X = np.array([[float(x) for x in r[1:-1]] for r in rows])
        dev = np.array([float(r[-1].split("=")[1]) for r in rows])
        print("%-8s T=%6d runs %d | sd of X_e over runs %s | mean X %.3f | max_v|S_v/t-2|: mean %.3f, max %.3f"
              % (g, T, len(X), " ".join("%.3f" % s for s in X.std(axis=0)), X.mean(), dev.mean(), dev.max()))
