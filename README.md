# Adaptive Rational Krylov for Sylvester Equations
This repository contains the implementation of a **block rational Krylov subspace method** for the solution of **Sylvester equations** of the form $AX - XB = UV^T$, with, $U,V$ tall and skinny, and $A$ and $B$ large and square.
The method is described in detail in [1], and has the following features:
 * Builds on **matrix-vector products and linear solves** with shifted copies of $A$ and $B$. In particular, the methods can take advantage of sparse or data-sparse matrices.
 * **Adaptive pole selection**: optimal poles are chosen during the construction of the projection subspace automatically. The pole selection is based on a novel criterion that enables a better understanding of the convergence of block Krylov methods, which is not just a straightforward extension of "scalar" Krylov methods.

The examples included in the repository can be used to test the method on two examples: the discretizations of a 2D Poisson problem, and a 2D convection-diffusion equation.

The algorithm computes rational Krylov subspaces using the rktoolbox described in [2]. The code has been tested using the version 2.9 of rktoolbox.

# References

[1]. Casulli, A. & Robol, L., [An effcient block rational Krylov solver for Sylvester equations with adaptive pole selection](https://arxiv.org/abs/2301.08103), arXiv.

[2]. Berljafa, M. & S. Elsworth, S. &  Guttel, S., [A rational Krylov toolbox for Matlab (2014)](http://guettel.com/rktoolbox/).

