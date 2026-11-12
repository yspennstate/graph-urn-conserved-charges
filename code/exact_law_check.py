"""Compare the charge of the urn on the even cycle C_L with its exact terminal law.

Runs R independent copies of the embedded jump chain (each step a uniformly chosen vertex rings) for K rings,
records the alternating charge M = sum_e (-1)^e H(N_e - 1), and compares it with a sum of L/2 independent
standard logistic variables: variance L*pi^2/6, excess kurtosis 12/(5L), and a two-sample Kolmogorov-Smirnov
distance. For L = 4 it also compares Var((1 + u)/2), u = tanh(M/L), with pi^2/32 - 1/4.

usage: python exact_law_check.py L K R seed
"""
import sys
import numpy as np

L, K, R, seed = (int(a) for a in sys.argv[1:5])
rng = np.random.default_rng(seed)
N = np.ones((R, L))
M = np.zeros(R)
sign = np.where(np.arange(L) % 2 == 0, 1.0, -1.0)
rows = np.arange(R)
for _ in range(K):
    v = rng.integers(0, L, R)            # ringing vertex; it sits between edges v-1 and v
    le, re = (v - 1) % L, v
    a, b = N[rows, le], N[rows, re]
    left = rng.random(R) * (a + b) < a
    e = np.where(left, le, re)
    n = np.where(left, a, b)
    M += sign[e] / n
    N[rows, e] += 1.0

ref = rng.logistic(size=(10**6, L // 2)).sum(axis=1)
var, ex = M.var(), ((M - M.mean())**4).mean() / M.var()**2 - 3
grid = np.sort(np.concatenate([M, ref]))
ks = np.abs(np.searchsorted(np.sort(M), grid, 'right') / R - np.searchsorted(np.sort(ref), grid, 'right') / len(ref)).max()
print(f"L={L} rings={K} runs={R} seed={seed}")
print(f"mean M      {M.mean():+.4f}  (s.e. {M.std() / R**0.5:.4f}; exact 0)")
print(f"Var M       {var:.4f}  (exact L*pi^2/6 = {L * np.pi**2 / 6:.4f}; s.e. about {var * (2 / R)**0.5:.4f})")
print(f"excess kurt {ex:+.4f}  (exact 12/(5L) = {12 / (5 * L):.4f})")
print(f"KS distance {ks:.4f}  (about {1.36 * (1 / R)**0.5:.4f} at the 5% level)")
print(f"min count   {N.min():.0f}")
if L == 4:
    p = (1 + np.tanh(M / L)) / 2
    print(f"Var p       {p.var():.5f}  (exact pi^2/32 - 1/4 = {np.pi**2 / 32 - 0.25:.5f})")
