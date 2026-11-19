# Conserved charges in linearly reinforced graph urns

Yitzchak Shmalo

This repository contains a paper on the graph-based Pólya urn with linear reinforcement (WARM with alpha = 1),
together with the programs used for its simulations.

Paper: [main.pdf](main.pdf) (source [main.tex](main.tex)).

## Summary

Each vertex of a finite graph rings at rate one and reinforces one of its edges with probability proportional to
the edge's count (the (W,A)-reinforcement model, WARM, with alpha = 1). The paper gives written proofs of the
following statements:

* for every sigma in the kernel of the unsigned vertex-edge incidence matrix, sum_e sigma_e H(N_e - 1) is a
  martingale (the charge), where H is the harmonic number and all counts start at 1;
* on every finite graph the limit of the charge has the law of sum_e sigma_e log W_e, with W_e independent standard
  exponential variables;
* on the even cycle C_L, L >= 4, the ratios N_e / t converge almost surely to (1 + u, 1 - u, ...), and
  L artanh(u) equals the limit of the charge, which has the law of a sum of L/2 independent standard logistic
  variables;
* on the integers the expectation of the integral over time of 1/(N_{-1} N_0) is pi^2/6.

The martingales and the form of the limit on even cycles are due to M. Holmes and V. Kleptsyn, whose work on the
critical regime is cited as "Infinite WARM graphs II: critical regime" in Hirsch, Holmes and Kleptsyn
(arXiv:2009.07682) and in Infinite WARM graphs III (Nonlinearity, 2023). The model is due to van der Hofstad,
Holmes, Kuznetsov and Ruszel (2016); the discrete-time urn and its convergence at linearity to Benaim, Benjamini,
Chen and Lima (2015), Chen and Lucas (2014) and Lima (2016). The paper's AI declaration says which system wrote
which part.

Section 4 also gives the local conservation law on the line and a heuristic envelope equation for the alternating
mode, which suggests a decay like (log t)^(-1/4) on Z.

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
  10^6 samples of a sum of L/2 logistic variables). Three runs are in `exact_law_runs.txt`.

## License

The text is released under CC BY 4.0 and the code under the MIT license.
