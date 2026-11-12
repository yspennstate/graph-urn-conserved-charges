# Conserved charges in linearly reinforced graph urns

Yitzchak Shmalo

This repository contains a paper on the graph-based Pólya urn with linear reinforcement (WARM with alpha = 1),
together with the programs used for its simulations.

Paper: [main.pdf](main.pdf) (source [main.tex](main.tex)).

## Summary

Each vertex of a finite graph rings at rate one and reinforces one of its edges with probability proportional to
the edge's count. The paper proves:

* for every sigma in the kernel of the unsigned vertex-edge incidence matrix, sum_e sigma_e H(N_e - 1) is an exact
  martingale, where H is the harmonic number; up to a vanishing error it equals sum_e sigma_e log(N_e / t), which
  the mean-field flow conserves;
* on every finite graph the terminal value of the charge has the law of sum_e sigma_e log W_e, with W_e
  independent standard exponential variables;
* on the even cycle C_L the normalised counts converge almost surely to (1 + u, 1 - u, ...), and
  L artanh(u) equals the terminal value of the charge, so u = tanh(M / L) with M a sum of L/2 independent
  standard logistic variables; sqrt(L) artanh(u) converges in law to the normal law with variance pi^2/6;
* on the integers the expected quadratic variation of the charge per vertex is pi^2/6.

This gives the law of the random limit for the continuous-time, one-vertex-at-a-time analogue of the urn studied
by Chen and Lucas (2014) and by Lima (2016, Question 1); no transfer to the synchronous urn is proved. The last
section derives the local conservation law on the line, a formal envelope equation for the alternating mode, and
the logarithmic decay rate it predicts for Kleptsyn's conjecture on the critical urn on Z.

## Code

The C programs are plain C99 (`cc -O2 -o warm warm.c -lm`) and print to standard output.

* `code/warm.c L Tmax seed` runs the urn on the cycle of length L and prints, at log-spaced times, the mean,
  spread and extremes of N_e / t, the alternating amplitude and the neighbour correlation.
* `code/warm2.c L T R seed` runs R independent copies on the even cycle to the first ring after T and prints,
  for each, the charge M, the ring-sampled bracket (the sum over rings of 1/(N_{v-1} N_v) for the ringing vertex,
  whose mean is the expected predictable quadratic variation) and the alternating amplitude u read from
  L artanh(u) = sum_e (-1)^e log(N_e / t).
* `code/warm_g.c T R seed nV u0 w0 u1 w1 ...` runs R copies on a small graph given by its edge list
  (at most 16 vertices) and prints the final N_e / t and the largest deviation of a vertex sum from 2.
* `code/exact_law_check.py L K R seed` (Python 3 with NumPy) runs R copies of the jump chain on C_L for K rings
  and compares the charge with the exact terminal law (variance, excess kurtosis, Kolmogorov-Smirnov distance to
  10^6 samples of a sum of L/2 logistic variables). The runs reported in the paper are in `exact_law_runs.txt`.

## License

The text is released under CC BY 4.0 and the code under the MIT license.
