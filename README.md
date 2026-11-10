# Conserved charges in linearly reinforced graph urns

Yitzchak Shmalo

This repository contains a paper on the graph-based Pólya urn with linear reinforcement (WARM with alpha = 1),
together with the C programs used for its simulations.

Paper: [main.pdf](main.pdf) (source [main.tex](main.tex)).

## Summary

Each vertex of a finite graph rings at rate one and reinforces one of its edges with probability proportional to
the edge's count. The paper proves:

* for every sigma in the kernel of the unsigned vertex-edge incidence matrix, sum_e sigma_e H(N_e - 1) is an exact
  martingale, where H is the harmonic number; up to a vanishing error it equals sum_e sigma_e log(N_e / t), which
  the mean-field flow conserves;
* on the even cycle C_L the normalised counts converge almost surely to (1 + u, 1 - u, ...), and
  L artanh(u) equals the terminal value of the charge;
* as L tends to infinity through even values, sqrt(L) artanh(u) converges in law to N(0, c), where c is the
  expected quadratic variation per vertex of the process on the integers (numerically c is about 1.64).

This gives the large-L law of the random limit for the continuous-time form of the question asked by Chen and
Lucas (2014) and by Lima (2016, Question 1). The last section derives the local conservation law on the line,
a formal envelope equation for the alternating mode, and the logarithmic decay rate it predicts for Kleptsyn's
conjecture on the critical urn on Z.

## Code

Each program is plain C99 (`cc -O2 -o warm warm.c -lm`) and prints to standard output.

* `code/warm.c L Tmax seed` runs the urn on the cycle of length L and prints, at log-spaced times, the mean,
  spread and extremes of N_e / t, the alternating amplitude and the neighbour correlation.
* `code/warm2.c L T R seed` runs R independent copies on the even cycle and prints, for each, the charge M,
  its quadratic variation and the alternating amplitude u read from L artanh(u) = sum_e (-1)^e log(N_e / t).
* `code/warm_g.c T R seed nV u0 w0 u1 w1 ...` runs R copies on a small graph given by its edge list
  (at most 16 vertices) and prints the final N_e / t and the largest deviation of a vertex sum from 2.

## License

The text is released under CC BY 4.0 and the code under the MIT license.
