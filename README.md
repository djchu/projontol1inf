## Implementation for "Semismooth Newton Algorithm for Efficient Projections onto $\ell_{1, \infty}$-norm Ball" to appear in ICML 2020 by Dejun Chu.

Our algorithm is implemented in C with a Matlab interface. To run the demo, you first need to compile the mex file with Matlab terminal:
   >> mex -output myssnewton mex-ssnewton.c
