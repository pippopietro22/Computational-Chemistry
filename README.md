# Computational-Chemistry
Set of basic programs written in Fortran 90 for Computational Chemistry.
This repo contains:
1) Cholensky decomposition: an implementation in fortran 90 of the Cholensky decomposition of a matrix "A" in two triangular matrix L1 and L2 (L1 = transpose(L2)).
2) Jacobi: an implementation in Fortran 90 of the Jacobi algorithm used to solve linear systems, but only if associeted to a diagonal dominant matrix.
3) Jacobi - DIIS: an implementation in fortran 90 of the Jacobi algorithm with the DIIS (direct inversion in the iterative subspace) that impoves our algorithm for any kind of matrices.
4) Conjugate Gradient: an implementation in fortran 90 of the Conjugate Gradient algorithm, used to solve linear systems. There is another file (Conjugate Gradient - Preconditioned) which implements the same algorithm but with a difference: there's a preconditioning of the vector used to find the solution to the linear system. This improves the efficency of the algorithm.
5) Davidson:  an implementation in fortran 90 of the Davidson algorithm, used to find initial few eigenvector and eigenvalues of a Huge matrix (dimension >> 1000).

All this programs, with the possible exception of the Cholensky Decomposition, must be compiled with a reference to the blas, lapack libraries: "gfortran [name_program.f90] -o [name] -lblas -llapack" 

BONUS: There is a program that works with second quantization: fils a binary vector that rapresents the spin-orbitalis using operators of creation and destruction.
