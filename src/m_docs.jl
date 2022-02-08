# Portions extracted from SLICOT-Reference distribution
# Copyright (c) 2002-2020 NICONET e.V.


"""

Compute the complex square root YR + i*YI of a complex number
XR + i*XI  in real arithmetic.  The returned result is so that
YR >= 0.0  and  SIGN(YI) = SIGN(XI).



See the SLICOT documentation for details.
"""
function ma01ad! end

"""

Compute the general product of K real scalars without over-
or underflow.



See the SLICOT documentation for details.
"""
function ma01bd! end

"""

Compute the general product of K complex scalars trying to
avoid over- and underflow.



See the SLICOT documentation for details.
"""
function ma01bz! end

"""

Compute, without over- or underflow, the sign of the sum of two
real numbers represented using integer powers of a base (usually,
the machine base). Any base can be used, but it should the same
for both numbers. The result is an integer with value 1, 0, or -1,
depending on the sum being found as positive, zero, or negative,
respectively.



See the SLICOT documentation for details.
"""
function ma01cd! end

"""

Transpose all or part of a two-dimensional matrix A into
another matrix B.



See the SLICOT documentation for details.
"""
function ma02ad! end

"""

Reverse the order of rows and/or columns of a given matrix A
by pre-multiplying and/or post-multiplying it, respectively, with
a permutation matrix P, where P is a square matrix of appropriate
order, with ones down the secondary diagonal.



See the SLICOT documentation for details.
"""
function ma02bd! end

"""

Reverse the order of rows and/or columns of a given matrix A
by pre-multiplying and/or post-multiplying it, respectively, with
a permutation matrix P, where P is a square matrix of appropriate
order, with ones down the secondary diagonal.



See the SLICOT documentation for details.
"""
function ma02bz! end

"""

Compute the pertranspose of a central band of a square matrix.



See the SLICOT documentation for details.
"""
function ma02cd! end

"""

Compute the pertranspose of a central band of a square matrix.



See the SLICOT documentation for details.
"""
function ma02cz! end

"""

Pack/unpack the upper or lower triangle of a symmetric matrix.
The packed matrix is stored column-wise in the one-dimensional
array AP.



See the SLICOT documentation for details.
"""
function ma02dd! end

"""

Store by symmetry the upper or lower triangle of a symmetric
matrix, given the other triangle.



See the SLICOT documentation for details.
"""
function ma02ed! end

"""

Store by skew-symmetry the upper or lower triangle of a
skew-symmetric matrix, given the other triangle. The diagonal
entries are set to zero.



See the SLICOT documentation for details.
"""
function ma02es! end

"""

Store by (skew-)symmetry the upper or lower triangle of a
(skew-)symmetric/Hermitian complex matrix, given the other
triangle.



See the SLICOT documentation for details.
"""
function ma02ez! end

"""

Compute the coefficients c and s (c^2 + s^2 = 1) for a modified
hyperbolic plane rotation, such that,

    `y1 := 1/c * x1 - s/c * x2 = sqrt(x1^2 - x2^2),`
    `y2 :=  -s * y1 +  c  * x2 = 0,`

given two real numbers x1 and x2, satisfying either x1 = x2 = 0,
or abs(x2) < abs(x1).



See the SLICOT documentation for details.
"""
function ma02fd! end

"""

Perform a series of column interchanges on the matrix A.
One column interchange is initiated for each of columns K1 through
K2 of A. This is useful for solving linear systems X*A = B, when
the matrix A has already been factored by LAPACK Library routine
DGETRF.



See the SLICOT documentation for details.
"""
function ma02gd! end

"""

Perform a series of column interchanges on the matrix A.
One column interchange is initiated for each of columns K1 through
K2 of A. This is useful for solving linear systems X*A = B, when
the matrix A has already been factored by LAPACK Library routine
DGETRF.



See the SLICOT documentation for details.
"""
function ma02gz! end

"""

Check if A = DIAG*I, where I is an M-by-N matrix with ones on
the diagonal and zeros elsewhere.



See the SLICOT documentation for details.
"""
function ma02hd! end

"""

Check if A = DIAG*I, where I is an M-by-N matrix with ones on
the diagonal and zeros elsewhere, A is a complex matrix and DIAG
is a complex scalar.



See the SLICOT documentation for details.
"""
function ma02hz! end

"""

Compute the value of the one norm, or the Frobenius norm, or
the infinity norm, or the element of largest absolute value
of a real skew-Hamiltonian matrix

              `[  A   G  ]          T         T`
        `X  =  [       T ],   G = -G,   Q = -Q,`
              `[  Q   A  ]`

or of a real Hamiltonian matrix

              `[  A   G  ]          T         T`
        `X  =  [       T ],   G =  G,   Q =  Q,`
              `[  Q  -A  ]`

where A, G and Q are real n-by-n matrices.

Note that for this kind of matrices the infinity norm is equal
to the one norm.



See the SLICOT documentation for details.
"""
function ma02id! end

"""

Compute the value of the one norm, or the Frobenius norm, or
the infinity norm, or the element of largest absolute value
of a complex skew-Hamiltonian matrix

              `[  A   G  ]          H         H`
        `X  =  [       H ],   G = -G,   Q = -Q,`
              `[  Q   A  ]`

or of a complex Hamiltonian matrix

              `[  A   G  ]          H         H`
        `X  =  [       H ],   G =  G,   Q =  Q,`
              `[  Q  -A  ]`

where A, G and Q are complex n-by-n matrices.

Note that for this kind of matrices the infinity norm is equal
to the one norm.



See the SLICOT documentation for details.
"""
function ma02iz! end

"""

Compute || Q^T Q - I ||_F for a matrix of the form

                  `[  op( Q1 )  op( Q2 ) ]`
             `Q =  [                     ],`
                  `[ -op( Q2 )  op( Q1 ) ]`

where Q1 and Q2 are N-by-N matrices. This residual can be used to
test wether Q is numerically an orthogonal symplectic matrix.



See the SLICOT documentation for details.
"""
function ma02jd! end

"""

Compute || Q^H Q - I ||_F for a complex matrix of the form

                  `[  op( Q1 )  op( Q2 ) ]`
             `Q =  [                     ],`
                  `[ -op( Q2 )  op( Q1 ) ]`

where Q1 and Q2 are N-by-N matrices. This residual can be used to
test wether Q is numerically a unitary symplectic matrix.



See the SLICOT documentation for details.
"""
function ma02jz! end

"""

Compute the value of the one norm, or the Frobenius norm, or
the infinity norm, or the element of largest absolute value
of a real skew-symmetric matrix.

Note that for this kind of matrices the infinity norm is equal
to the one norm.



See the SLICOT documentation for details.
"""
function ma02md! end

"""

Compute the value of the one norm, or the Frobenius norm, or
the infinity norm, or the element of largest absolute value
of a complex skew-Hermitian matrix.

Note that for this kind of matrices the infinity norm is equal
to the one norm.



See the SLICOT documentation for details.
"""
function ma02mz! end

"""

Permute two specified rows and corresponding columns of a
(skew-)symmetric/Hermitian complex matrix.



See the SLICOT documentation for details.
"""
function ma02nz! end

"""

Compute the number of zero rows (and zero columns) of a real
(skew-)Hamiltonian matrix,

      `(  A    D   )`
  `H = (           ).`
      `(  E  +/-A' )`



See the SLICOT documentation for details.
"""
function ma02od! end

"""

Compute the number of zero rows (and zero columns) of a complex
(skew-)Hamiltonian matrix,

      `(  A    D   )`
  `H = (           ).`
      `(  E  +/-A' )`



See the SLICOT documentation for details.
"""
function ma02oz! end

"""

Compute the number of zero rows and zero columns of a real
matrix.



See the SLICOT documentation for details.
"""
function ma02pd! end

"""

Compute the number of zero rows and zero columns of a complex
matrix.



See the SLICOT documentation for details.
"""
function ma02pz! end

"""

Perform one of the skew-symmetric rank 2k operations

    `C := alpha*A*B' - alpha*B*A' + beta*C,`

or

    `C := alpha*A'*B - alpha*B'*A + beta*C,`

where alpha and beta are scalars, C is a real N-by-N skew-
symmetric matrix and A, B are N-by-K matrices in the first case
and K-by-N matrices in the second case.

This is a modified version of the vanilla implemented BLAS
routine DSYR2K written by Jack Dongarra, Iain Duff,
Jeremy Du Croz and Sven Hammarling.



See the SLICOT documentation for details.
"""
function mb01kd! end

"""

Compute the matrix formula
   `_`
   `R = alpha*R + beta*op( A )*X*op( A )',`
                                            `_`
where alpha and beta are scalars, R, X, and R are skew-symmetric
matrices, A is a general matrix, and op( A ) is one of

   `op( A ) = A   or   op( A ) = A'.`

The result is overwritten on R.



See the SLICOT documentation for details.
"""
function mb01ld! end

"""

Perform the matrix-vector operation

   `y := alpha*A*x + beta*y,`

where alpha and beta are scalars, x and y are vectors of length
n and A is an n-by-n skew-symmetric matrix.

This is a modified version of the vanilla implemented BLAS
routine DSYMV written by Jack Dongarra, Jeremy Du Croz,
Sven Hammarling, and Richard Hanson.



See the SLICOT documentation for details.
"""
function mb01md! end

"""

Perform the skew-symmetric rank 2 operation

     `A := alpha*x*y' - alpha*y*x' + A,`

where alpha is a scalar, x and y are vectors of length n and A is
an n-by-n skew-symmetric matrix.

This is a modified version of the vanilla implemented BLAS
routine DSYR2 written by Jack Dongarra, Jeremy Du Croz,
Sven Hammarling, and Richard Hanson.



See the SLICOT documentation for details.
"""
function mb01nd! end

"""

Perform one of the special symmetric rank 2k operations

   `R := alpha*R + beta*H*X + beta*X*H',`
or
   `R := alpha*R + beta*H'*X + beta*X*H,`

where alpha and beta are scalars, R and X are N-by-N symmetric
matrices, and H is an N-by-N upper Hessenberg matrix.



See the SLICOT documentation for details.
"""
function mb01oc! end

"""

Compute the matrix formula

  `R := alpha*R + beta*( op( H )*X*op( E )' + op( E )*X*op( H )' ),`

where alpha and beta are scalars, R and X are symmetric matrices,
H is an upper Hessenberg matrix, E is an upper triangular matrix,
and op( M ) is one of

  `op( M ) = M   or   op( M ) = M'.`

The result is overwritten on R.



See the SLICOT documentation for details.
"""
function mb01od! end

"""

Compute one of the symmetric rank 2k operations

   `R := alpha*R + beta*H*E' + beta*E*H',`

or

   `R := alpha*R + beta*H'*E + beta*E'*H,`

where alpha and beta are scalars, R, E, and H are N-by-N matrices,
with H upper Hessenberg and E upper triangular.



See the SLICOT documentation for details.
"""
function mb01oe! end

"""

Compute one of the symmetric rank 2k operations

   `R := alpha*R + beta*H*A' + beta*A*H',`

or

   `R := alpha*R + beta*H'*A + beta*A'*H,`

where alpha and beta are scalars, R, A, and H are N-by-N matrices,
with A and H upper Hessenberg.



See the SLICOT documentation for details.
"""
function mb01oh! end

"""

Compute either P or P', with P defined by the matrix formula

   `P = op( H )*X*op( E )',`

where H is an upper Hessenberg matrix, X is a symmetric matrix,
E is an upper triangular matrix, and op( M ) is one of

   `op( M ) = M   or   op( M ) = M'.`



See the SLICOT documentation for details.
"""
function mb01oo! end

"""

Compute P = H*X or P = X*H, where H is an upper Hessenberg
matrix and X is a symmetric matrix.



See the SLICOT documentation for details.
"""
function mb01os! end

"""

Compute one of the symmetric rank 2k operations

   `R := alpha*R + beta*E*T' + beta*T*E',`

or

   `R := alpha*R + beta*E'*T + beta*T'*E,`

where alpha and beta are scalars, R, T, and E are N-by-N matrices,
with T and E upper triangular.



See the SLICOT documentation for details.
"""
function mb01ot! end

"""

Scale a matrix or undo scaling.  Scaling is performed, if
necessary, so that the matrix norm will be in a safe range of
representable numbers.



See the SLICOT documentation for details.
"""
function mb01pd! end

"""

Multiply the M by N real matrix A by the real scalar CTO/CFROM.
This is done without over/underflow as long as the final result
CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that
A may be full, (block) upper triangular, (block) lower triangular,
(block) upper Hessenberg, or banded.



See the SLICOT documentation for details.
"""
function mb01qd! end

"""

Compute either the upper or lower triangular part of one of the
matrix formulas
   `_`
   `R = alpha*R + beta*op( A )*B,                               (1)`
   `_`
   `R = alpha*R + beta*B*op( A ),                               (2)`
                                        `_`
where alpha and beta are scalars, R and R are m-by-m matrices,
op( A ) and B are m-by-n and n-by-m matrices for (1), or n-by-m
and m-by-n matrices for (2), respectively, and op( A ) is one of

   `op( A ) = A   or   op( A ) = A',  the transpose of A.`

The result is overwritten on R.



See the SLICOT documentation for details.
"""
function mb01rb! end

"""

Compute the matrix formula
   `_`
   `R = alpha*R + beta*op( A )*X*op( A )',`
                                            `_`
where alpha and beta are scalars, R, X, and R are symmetric
matrices, A is a general matrix, and op( A ) is one of

   `op( A ) = A   or   op( A ) = A'.`

The result is overwritten on R.



See the SLICOT documentation for details.
"""
function mb01rd! end

"""

Compute the matrix formula

   `R := alpha*R + beta*op( H )*X*op( H )',`

where alpha and beta are scalars, R and X are symmetric matrices,
H is an upper Hessenberg matrix, and op( H ) is one of

   `op( H ) = H   or   op( H ) = H'.`

The result is overwritten on R.



See the SLICOT documentation for details.
"""
function mb01rh! end

"""

Compute the matrix formula

   `R := alpha*R + beta*op( E )*X*op( E )',`

where alpha and beta are scalars, R and X are symmetric matrices,
E is an upper triangular matrix, and op( E ) is one of

   `op( E ) = E   or   op( E ) = E'.`

The result is overwritten on R.



See the SLICOT documentation for details.
"""
function mb01rt! end

"""

Compute the matrix formula
   `_`
   `R = alpha*R + beta*op( A )*X*op( A )',`
                                            `_`
where alpha and beta are scalars, R, X, and R are symmetric
matrices, A is a general matrix, and op( A ) is one of

   `op( A ) = A   or   op( A ) = A'.`

The result is overwritten on R.



See the SLICOT documentation for details.
"""
function mb01ru! end

"""

Compute the transformation of the symmetric matrix A by the
matrix Z in the form

   `A := op(Z)*A*op(Z)',`

where op(Z) is either Z or its transpose, Z'.



See the SLICOT documentation for details.
"""
function mb01rw! end

"""

Compute either the upper or lower triangular part of one of the
matrix formulas
   `_`
   `R = alpha*R + beta*op( A )*B,                               (1)`
   `_`
   `R = alpha*R + beta*B*op( A ),                               (2)`
                                        `_`
where alpha and beta are scalars, R and R are m-by-m matrices,
op( A ) and B are m-by-n and n-by-m matrices for (1), or n-by-m
and m-by-n matrices for (2), respectively, and op( A ) is one of

   `op( A ) = A   or   op( A ) = A',  the transpose of A.`

The result is overwritten on R.



See the SLICOT documentation for details.
"""
function mb01rx! end

"""

Compute either the upper or lower triangular part of one of the
matrix formulas
   `_`
   `R = alpha*R + beta*op( H )*B,                               (1)`
   `_`
   `R = alpha*R + beta*B*op( H ),                               (2)`
                                               `_`
where alpha and beta are scalars, H, B, R, and R are m-by-m
matrices, H is an upper Hessenberg matrix, and op( H ) is one of

   `op( H ) = H   or   op( H ) = H',  the transpose of H.`

The result is overwritten on R.



See the SLICOT documentation for details.
"""
function mb01ry! end

"""

Scale a general M-by-N matrix A using the row and column
scaling factors in the vectors R and C.



See the SLICOT documentation for details.
"""
function mb01sd! end

"""

Scale a symmetric N-by-N matrix A using the row and column
scaling factors stored in the vector D.



See the SLICOT documentation for details.
"""
function mb01ss! end

"""

Compute the matrix product A * B, where A and B are upper
quasi-triangular matrices (that is, block upper triangular with
1-by-1 or 2-by-2 diagonal blocks) with the same structure.
The result is returned in the array B.



See the SLICOT documentation for details.
"""
function mb01td! end

"""

Compute one of the matrix products

   `B = alpha*op( H ) * A, or B = alpha*A * op( H ),`

where alpha is a scalar, A and B are m-by-n matrices, H is an
upper Hessenberg matrix, and op( H ) is one of

   `op( H ) = H   or   op( H ) = H',  the transpose of H.`



See the SLICOT documentation for details.
"""
function mb01ud! end

"""

Compute one of the matrix products

   `A : = alpha*op( H ) * A, or A : = alpha*A * op( H ),`

where alpha is a scalar, A is an m-by-n matrix, H is an upper
Hessenberg matrix, and op( H ) is one of

   `op( H ) = H   or   op( H ) = H',  the transpose of H.`



See the SLICOT documentation for details.
"""
function mb01uw! end

"""

Compute one of the matrix products

  `A : = alpha*op( T ) * A, or A : = alpha*A * op( T ),`

where alpha is a scalar, A is an m-by-n matrix, T is a quasi-
triangular matrix, and op( T ) is one of

   `op( T ) = T   or   op( T ) = T',  the transpose of T.`



See the SLICOT documentation for details.
"""
function mb01ux! end

"""

Perform the following matrix operation

   `C = alpha*kron( op(A), op(B) ) + beta*C,`

where alpha and beta are real scalars, op(M) is either matrix M or
its transpose, M', and kron( X, Y ) denotes the Kronecker product
of the matrices X and Y.



See the SLICOT documentation for details.
"""
function mb01vd! end

"""

Compute the matrix formula
_
R = alpha*( op( A )'*op( T )'*op( T ) + op( T )'*op( T )*op( A ) )
    `+ beta*R,                                                  (1)`

if DICO = 'C', or
_
R = alpha*( op( A )'*op( T )'*op( T )*op( A ) -  op( T )'*op( T ))
    `+ beta*R,                                                  (2)`
                                                        `_`
if DICO = 'D', where alpha and beta are scalars, R, and R are
symmetric matrices, T is a triangular matrix, A is a general or
Hessenberg matrix, and op( M ) is one of

   `op( M ) = M   or   op( M ) = M'.`

The result is overwritten on R.



See the SLICOT documentation for details.
"""
function mb01wd! end

"""

Compute the matrix product U' * U or L * L', where U and L are
upper and lower triangular matrices, respectively, stored in the
corresponding upper or lower triangular part of the array A.

If UPLO = 'U' then the upper triangle of the result is stored,
overwriting the matrix U in A.
If UPLO = 'L' then the lower triangle of the result is stored,
overwriting the matrix L in A.



See the SLICOT documentation for details.
"""
function mb01xd! end

"""

Compute the matrix product U' * U or L * L', where U and L are
upper and lower triangular matrices, respectively, stored in the
corresponding upper or lower triangular part of the array A.

If UPLO = 'U' then the upper triangle of the result is stored,
overwriting the matrix U in A.
If UPLO = 'L' then the lower triangle of the result is stored,
overwriting the matrix L in A.



See the SLICOT documentation for details.
"""
function mb01xy! end

"""

Perform the symmetric rank k operations

   `C := alpha*op( A )*op( A )' + beta*C,`

where alpha and beta are scalars, C is an n-by-n symmetric matrix,
op( A ) is an n-by-k matrix, and op( A ) is one of

   `op( A ) = A   or   op( A ) = A'.`

The matrix A has l nonzero codiagonals, either upper or lower.



See the SLICOT documentation for details.
"""
function mb01yd! end

"""

Compute the matrix product

   `H := alpha*op( T )*H,   or   H := alpha*H*op( T ),`

where alpha is a scalar, H is an m-by-n upper or lower
Hessenberg-like matrix (with l nonzero subdiagonals or
superdiagonals, respectively), T is a unit, or non-unit,
upper or lower triangular matrix, and op( T ) is one of

   `op( T ) = T   or   op( T ) = T'.`



See the SLICOT documentation for details.
"""
function mb01zd! end

"""

Compute the Cholesky factor and the generator and/or the
Cholesky factor of the inverse of a symmetric positive definite
(s.p.d.) block Toeplitz matrix T, defined by either its first
block row, or its first block column, depending on the routine
parameter TYPET. Transformation information is stored.



See the SLICOT documentation for details.
"""
function mb02cd! end

"""

Bring the first blocks of a generator to proper form.
The positive part of the generator is contained in the arrays A1
and A2. The negative part of the generator is contained in B.
Transformation information will be stored and can be applied
via SLICOT Library routine MB02CV.



See the SLICOT documentation for details.
"""
function mb02cu! end

"""

Apply the transformations created by the SLICOT Library routine
MB02CU on other columns / rows of the generator, contained in the
arrays F1, F2 and G.



See the SLICOT documentation for details.
"""
function mb02cv! end

"""

Bring the first blocks of a generator in proper form.
The columns / rows of the positive and negative generators
are contained in the arrays A and B, respectively.
Transformation information will be stored and can be applied
via SLICOT Library routine MB02CY.



See the SLICOT documentation for details.
"""
function mb02cx! end

"""

Apply the transformations created by the SLICOT Library
routine MB02CX on other columns / rows of the generator,
contained in the arrays A and B of positive and negative
generators, respectively.



See the SLICOT documentation for details.
"""
function mb02cy! end

"""

Update the Cholesky factor and the generator and/or the
Cholesky factor of the inverse of a symmetric positive definite
(s.p.d.) block Toeplitz matrix T, given the information from
a previous factorization and additional blocks in TA of its first
block row, or its first block column, depending on the routine
parameter TYPET. Transformation information is stored.



See the SLICOT documentation for details.
"""
function mb02dd! end

"""

Solve a system of linear equations  T*X = B  or  X*T = B  with
a symmetric positive definite (s.p.d.) block Toeplitz matrix T.
T is defined either by its first block row or its first block
column, depending on the parameter TYPET.



See the SLICOT documentation for details.
"""
function mb02ed! end

"""

Compute the incomplete Cholesky (ICC) factor of a symmetric
positive definite (s.p.d.) block Toeplitz matrix T, defined by
either its first block row, or its first block column, depending
on the routine parameter TYPET.

By subsequent calls of this routine, further rows / columns of
the Cholesky factor can be added.
Furthermore, the generator of the Schur complement of the leading
(P+S)*K-by-(P+S)*K block in T is available, which can be used,
e.g., for measuring the quality of the ICC factorization.



See the SLICOT documentation for details.
"""
function mb02fd! end

"""

Compute the Cholesky factor of a banded symmetric positive
definite (s.p.d.) block Toeplitz matrix, defined by either its
first block row, or its first block column, depending on the
routine parameter TYPET.

By subsequent calls of this routine the Cholesky factor can be
computed block column by block column.



See the SLICOT documentation for details.
"""
function mb02gd! end

"""

Compute, for a banded K*M-by-L*N block Toeplitz matrix T with
block size (K,L), specified by the nonzero blocks of its first
block column TC and row TR, a LOWER triangular matrix R (in band
storage scheme) such that
                     `T          T`
                    `T  T  =  R R .                             (1)`

It is assumed that the first MIN(M*K, N*L) columns of T are
linearly independent.

By subsequent calls of this routine, the matrix R can be computed
block column by block column.



See the SLICOT documentation for details.
"""
function mb02hd! end

"""

Solve the overdetermined or underdetermined real linear systems
involving an M*K-by-N*L block Toeplitz matrix T that is specified
by its first block column and row. It is assumed that T has full
rank.
The following options are provided:

1. If JOB = 'O' or JOB = 'A' :  find the least squares solution of
   `an overdetermined system, i.e., solve the least squares problem`

             `minimize || B - T*X ||.                           (1)`

2. If JOB = 'U' or JOB = 'A' :  find the minimum norm solution of
   `the undetermined system`
              `T`
             `T * X = C.                                        (2)`



See the SLICOT documentation for details.
"""
function mb02id! end

"""

Compute a lower triangular matrix R and a matrix Q with
Q^T Q = I such that
                               `T`
                      `T  =  Q R ,`

where T is a K*M-by-L*N block Toeplitz matrix with blocks of size
(K,L). The first column of T will be denoted by TC and the first
row by TR. It is assumed that the first MIN(M*K, N*L) columns of T
have full rank.

By subsequent calls of this routine the factors Q and R can be
computed block column by block column.



See the SLICOT documentation for details.
"""
function mb02jd! end

"""

Compute a low rank QR factorization with column pivoting of a
K*M-by-L*N block Toeplitz matrix T with blocks of size (K,L);
specifically,
                                `T`
                      `T P =  Q R ,`

where R is lower trapezoidal, P is a block permutation matrix
and Q^T Q = I. The number of columns in R is equivalent to the
numerical rank of T with respect to the given tolerance TOL1.
Note that the pivoting scheme is local, i.e., only columns
belonging to the same block in T are permuted.



See the SLICOT documentation for details.
"""
function mb02jx! end

"""

Compute the matrix product

          `C = alpha*op( T )*B + beta*C,`

where alpha and beta are scalars and T is a block Toeplitz matrix
specified by its first block column TC and first block row TR;
B and C are general matrices of appropriate dimensions.



See the SLICOT documentation for details.
"""
function mb02kd! end

"""

Solve the Total Least Squares (TLS) problem using a Singular
Value Decomposition (SVD) approach.
The TLS problem assumes an overdetermined set of linear equations
AX = B, where both the data matrix A as well as the observation
matrix B are inaccurate. The routine also solves determined and
underdetermined sets of equations by computing the minimum norm
solution.
It is assumed that all preprocessing measures (scaling, coordinate
transformations, whitening, ... ) of the data have been performed
in advance.



See the SLICOT documentation for details.
"""
function mb02md! end

"""

Solve the Total Least Squares (TLS) problem using a Partial
Singular Value Decomposition (PSVD) approach.
The TLS problem assumes an overdetermined set of linear equations
AX = B, where both the data matrix A as well as the observation
matrix B are inaccurate. The routine also solves determined and
underdetermined sets of equations by computing the minimum norm
solution.
It is assumed that all preprocessing measures (scaling, coordinate
transformations, whitening, ... ) of the data have been performed
in advance.



See the SLICOT documentation for details.
"""
function mb02nd! end

"""

Separate a zero singular value of a bidiagonal submatrix of
order k, k <= p, of the bidiagonal matrix

          `|Q(1) E(1)  0    ...   0   |`
          `| 0   Q(2) E(2)        .   |`
      `J = | .                    .   |`
          `| .                  E(p-1)|`
          `| 0   ...  ...   ...  Q(p) |`

with p = MIN(M,N), by annihilating one or two superdiagonal
elements E(i-1) (if i > 1) and/or E(i) (if i < k).



See the SLICOT documentation for details.
"""
function mb02ny! end

"""

Solve (if well-conditioned) one of the matrix equations

   `op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,`

where alpha is a scalar, X and B are m-by-n matrices, A is a unit,
or non-unit, upper or lower triangular matrix and op( A ) is one
of

   `op( A ) = A   or   op( A ) = A'.`

An estimate of the reciprocal of the condition number of the
triangular matrix A, in either the 1-norm or the infinity-norm, is
also computed as

   `RCOND = 1 / ( norm(A) * norm(inv(A)) ).`

and the specified matrix equation is solved only if RCOND is
larger than a given tolerance TOL.  In that case, the matrix X is
overwritten on B.



See the SLICOT documentation for details.
"""
function mb02od! end

"""

Solve (if well-conditioned) the matrix equations

   `op( A )*X = B,`

where X and B are N-by-NRHS matrices, A is an N-by-N matrix and
op( A ) is one of

   `op( A ) = A   or   op( A ) = A'.`

Error bounds on the solution and a condition estimate are also
provided.



See the SLICOT documentation for details.
"""
function mb02pd! end

"""

Compute a solution, optionally corresponding to specified free
elements, to a real linear least squares problem:

    `minimize || A * X - B ||`

using a complete orthogonal factorization of the M-by-N matrix A,
which may be rank-deficient.

Several right hand side vectors b and solution vectors x can be
handled in a single call; they are stored as the columns of the
M-by-NRHS right hand side matrix B and the N-by-NRHS solution
matrix X.



See the SLICOT documentation for details.
"""
function mb02qd! end

"""

Determine the minimum-norm solution to a real linear least
squares problem:

    `minimize || A * X - B ||,`

using the rank-revealing QR factorization of a real general
M-by-N matrix  A,  computed by SLICOT Library routine  MB03OD.



See the SLICOT documentation for details.
"""
function mb02qy! end

"""

Solve a system of linear equations
   `H * X = B  or  H' * X = B`
with an upper Hessenberg N-by-N matrix H using the LU
factorization computed by MB02SD.



See the SLICOT documentation for details.
"""
function mb02rd! end

"""

Solve a system of linear equations
   `H * X = B,  H' * X = B  or  H**H * X = B`
with a complex upper Hessenberg N-by-N matrix H using the LU
factorization computed by MB02SZ.



See the SLICOT documentation for details.
"""
function mb02rz! end

"""

Compute an LU factorization of an n-by-n upper Hessenberg
matrix H using partial pivoting with row interchanges.



See the SLICOT documentation for details.
"""
function mb02sd! end

"""

Compute an LU factorization of a complex n-by-n upper
Hessenberg matrix H using partial pivoting with row interchanges.



See the SLICOT documentation for details.
"""
function mb02sz! end

"""

Estimate the reciprocal of the condition number of an upper
Hessenberg matrix H, in either the 1-norm or the infinity-norm,
using the LU factorization computed by MB02SD.



See the SLICOT documentation for details.
"""
function mb02td! end

"""

Estimate the reciprocal of the condition number of a complex
upper Hessenberg matrix H, in either the 1-norm or the
infinity-norm, using the LU factorization computed by MB02SZ.



See the SLICOT documentation for details.
"""
function mb02tz! end

"""

Compute the minimum norm least squares solution of one of the
following linear systems

   `op(R)*X = alpha*B,                                          (1)`
   `X*op(R) = alpha*B,                                          (2)`

where alpha is a real scalar, op(R) is either R or its transpose,
R', R is an L-by-L real upper triangular matrix, B is an M-by-N
real matrix, and L = M for (1), or L = N for (2). Singular value
decomposition, R = Q*S*P', is used, assuming that R is rank
deficient.



See the SLICOT documentation for details.
"""
function mb02ud! end

"""

Solve for x in A * x = scale * RHS, using the LU factorization
of the N-by-N matrix A computed by SLICOT Library routine MB02UV.
The factorization has the form A = P * L * U * Q, where P and Q
are permutation matrices, L is unit lower triangular and U is
upper triangular.



See the SLICOT documentation for details.
"""
function mb02uu! end

"""

Compute an LU factorization, using complete pivoting, of the
N-by-N matrix A. The factorization has the form A = P * L * U * Q,
where P and Q are permutation matrices, L is lower triangular with
unit diagonal elements and U is upper triangular.



See the SLICOT documentation for details.
"""
function mb02uv! end

"""

Solve a system of the form  A X = s B  or  A' X = s B  with
possible scaling ("s") and perturbation of A.  (A' means
A-transpose.)  A is an N-by-N real matrix, and X and B are
N-by-M matrices.  N may be 1 or 2.  The scalar "s" is a scaling
factor (.LE. 1), computed by this subroutine, which is so chosen
that X can be computed without overflow.  X is further scaled if
necessary to assure that norm(A)*norm(X) is less than overflow.



See the SLICOT documentation for details.
"""
function mb02uw! end

"""

Compute the solution to a real system of linear equations
   `X * op(A) = B,`
where op(A) is either A or its transpose, A is an N-by-N matrix,
and X and B are M-by-N matrices.
The LU decomposition with partial pivoting and row interchanges,
A = P * L * U, is used, where P is a permutation matrix, L is unit
lower triangular, and U is upper triangular.



See the SLICOT documentation for details.
"""
function mb02vd! end

"""

Determine a vector x which solves the system of linear
equations

      `A*x = b ,     D*x = 0 ,`

in the least squares sense, where A is an m-by-n matrix,
D is an n-by-n diagonal matrix, and b is an m-vector.
It is assumed that a QR factorization, with column pivoting, of A
is available, that is, A*P = Q*R, where P is a permutation matrix,
Q has orthogonal columns, and R is an upper triangular matrix
with diagonal elements of nonincreasing magnitude.
The routine needs the full upper triangle of R, the permutation
matrix P, and the first n components of Q'*b (' denotes the
transpose). The system A*x = b, D*x = 0, is then equivalent to

      `R*z = Q'*b ,  P'*D*P*z = 0 ,                             (1)`

where x = P*z. If this system does not have full rank, then a
least squares solution is obtained. On output, MB02YD also
provides an upper triangular matrix S such that

      `P'*(A'*A + D*D)*P = S'*S .`

The system (1) is equivalent to S*z = c , where c contains the
first n components of the vector obtained by applying to
[ (Q'*b)'  0 ]' the transformations which triangularized
[ R'  P'*D*P ]', getting S.



See the SLICOT documentation for details.
"""
function mb02yd! end

"""

Compute two Givens rotations (C1,S1) and (C2,S2) such that the
orthogonal matrix

          `[ Q  0 ]        [  C1  S1  0 ]   [ 1  0   0  ]`
      `Z = [      ],  Q := [ -S1  C1  0 ] * [ 0  C2  S2 ],`
          `[ 0  I ]        [  0   0   1 ]   [ 0 -S2  C2 ]`

makes the first column of the real Wilkinson double shift
polynomial of the product of matrices in periodic upper Hessenberg
form, stored in the array A, parallel to the first unit vector.
Only the rotation defined by C1 and S1 is needed for the real
Wilkinson single shift polynomial (see the SLICOT Library routines
MB03BE or MB03BF). The shifts are defined based on the eigenvalues
(computed externally by the SLICOT Library routine MB03BB) of the
trailing 2-by-2 submatrix of the matrix product. See the
definitions of the arguments W1 and W2.



See the SLICOT documentation for details.
"""
function mb03ab! end

"""

Compute two Givens rotations (C1,S1) and (C2,S2) such that the
orthogonal matrix

          `[ Q  0 ]        [  C1  S1  0 ]   [ 1  0   0  ]`
      `Z = [      ],  Q := [ -S1  C1  0 ] * [ 0  C2  S2 ],`
          `[ 0  I ]        [  0   0   1 ]   [ 0 -S2  C2 ]`

makes the first column of the real Wilkinson double shift
polynomial of the product of matrices in periodic upper Hessenberg
form, stored in the array A, parallel to the first unit vector.
Only the rotation defined by C1 and S1 is used for the real
Wilkinson single shift polynomial (see SLICOT Library routine
MB03BE).



See the SLICOT documentation for details.
"""
function mb03ad! end

"""

Compute two Givens rotations (C1,S1) and (C2,S2) such that the
orthogonal matrix

          `[ Q  0 ]        [  C1  S1  0 ]   [ 1  0   0  ]`
      `Z = [      ],  Q := [ -S1  C1  0 ] * [ 0  C2  S2 ],`
          `[ 0  I ]        [  0   0   1 ]   [ 0 -S2  C2 ]`

makes the first column of the real Wilkinson double shift
polynomial of the product of matrices in periodic upper Hessenberg
form, stored in the array A, parallel to the first unit vector.
Only the rotation defined by C1 and S1 is used for the real
Wilkinson single shift polynomial (see SLICOT Library routines
MB03BE or MB03BF). All factors whose exponents differ from that of 
the Hessenberg factor are assumed nonsingular. The trailing 2-by-2
submatrix and the five nonzero elements in the first two columns
of the matrix product are evaluated when a double shift is used.



See the SLICOT documentation for details.
"""
function mb03ae! end

"""

Compute two Givens rotations (C1,S1) and (C2,S2) such that the
orthogonal matrix

          `[ Q  0 ]        [  C1  S1  0 ]   [ 1  0   0  ]`
      `Z = [      ],  Q := [ -S1  C1  0 ] * [ 0  C2  S2 ],`
          `[ 0  I ]        [  0   0   1 ]   [ 0 -S2  C2 ]`

makes the first column of the real Wilkinson double shift
polynomial of the product of matrices in periodic upper Hessenberg
form, stored in the array A, parallel to the first unit vector.
Only the rotation defined by C1 and S1 is used for the real
Wilkinson single shift polynomial (see SLICOT Library routines
MB03BE or MB03BF).



See the SLICOT documentation for details.
"""
function mb03af! end

"""

Compute two Givens rotations (C1,S1) and (C2,S2) such that the
orthogonal matrix

          `[ Q  0 ]        [  C1  S1  0 ]   [ 1  0   0  ]`
      `Z = [      ],  Q := [ -S1  C1  0 ] * [ 0  C2  S2 ],`
          `[ 0  I ]        [  0   0   1 ]   [ 0 -S2  C2 ]`

makes the first column of the real Wilkinson double shift
polynomial of the product of matrices in periodic upper Hessenberg
form, stored in the array A, parallel to the first unit vector.
Only the rotation defined by C1 and S1 is used for the real
Wilkinson single shift polynomial (see SLICOT Library routines
MB03BE or MB03BF). All factors whose exponents differ from that of 
the Hessenberg factor are assumed nonsingular. The matrix product
is evaluated.



See the SLICOT documentation for details.
"""
function mb03ag! end

"""

Compute two Givens rotations (C1,S1) and (C2,S2) such that the
orthogonal matrix

          `[ Q  0 ]        [  C1  S1  0 ]   [ 1  0   0  ]`
      `Z = [      ],  Q := [ -S1  C1  0 ] * [ 0  C2  S2 ],`
          `[ 0  I ]        [  0   0   1 ]   [ 0 -S2  C2 ]`

makes the first column of the real Wilkinson double shift
polynomial of the product of matrices in periodic upper Hessenberg
form, stored in the array A, parallel to the first unit vector.
Only the rotation defined by C1 and S1 is used for the real
Wilkinson single shift polynomial (see SLICOT Library routines
MB03BE or MB03BF). All factors whose exponents differ from that of 
the Hessenberg factor are assumed nonsingular. The trailing 2-by-2
submatrix and the five nonzero elements in the first two columns
of the matrix product are evaluated when a double shift is used.



See the SLICOT documentation for details.
"""
function mb03ah! end

"""

Compute two Givens rotations (C1,S1) and (C2,S2)
such that the orthogonal matrix

          `[ Q  0 ]        [  C1  S1  0 ]   [ 1  0   0  ]`
      `Z = [      ],  Q := [ -S1  C1  0 ] * [ 0  C2  S2 ],`
          `[ 0  I ]        [  0   0   1 ]   [ 0 -S2  C2 ]`

makes the first column of the real Wilkinson double shift
polynomial of the product of matrices in periodic upper Hessenberg
form, stored in the array A, parallel to the first unit vector.
Only the rotation defined by C1 and S1 is used for the real
Wilkinson single shift polynomial (see SLICOT Library routines
MB03BE or MB03BF). All factors whose exponents differ from that of 
the Hessenberg factor are assumed nonsingular. The matrix product
is evaluated.



See the SLICOT documentation for details.
"""
function mb03ai! end

"""

Compute the suitable maps for Hessenberg index H and
signature array S. Auxiliary routine for the periodic QZ
algorithms.



See the SLICOT documentation for details.
"""
function mb03ba! end

"""

Compute the eigenvalues of a general 2-by-2 matrix product via
a complex single shifted periodic QZ algorithm.



See the SLICOT documentation for details.
"""
function mb03bb! end

"""

Compute the product singular value decomposition of the K-1
triangular factors corresponding to a 2-by-2 product of K
factors in upper Hessenberg-triangular form.
For a general product of 2-by-2 triangular matrices

                   `S(2)        S(3)            S(K)`
       `A = A(:,:,2)    A(:,:,3)    ... A(:,:,K),`

Givens rotations are computed so that
                                                     `S(i)`
  `[  CV(i-1) SV(i-1) ] [ A(1,1,i)(in)  A(1,2,i)(in) ]`
  `[ -SV(i-1) CV(i-1) ] [     0         A(2,2,i)(in) ]`
                                 `S(i)`
  `[ A(1,1,i)(out) A(1,2,i)(out) ]    [  CV(i) SV(i) ]`
= [     0         A(2,2,i)(out) ]    [ -SV(i) CV(i) ]

stays upper triangular and

  `[  CV(1) SV(1) ]       [ CV(K) -SV(K) ]`
  `[ -SV(1) CV(1) ] * A * [ SV(K)  CV(K) ]`

is diagonal.



See the SLICOT documentation for details.
"""
function mb03bc! end

"""

Find the eigenvalues of the generalized matrix product

             `S(1)           S(2)                 S(K)`
     `A(:,:,1)     * A(:,:,2)     * ... * A(:,:,K)`

where A(:,:,H) is upper Hessenberg and A(:,:,i), i <> H, is upper
triangular, using a double-shift version of the periodic
QZ method. In addition, A may be reduced to periodic Schur form:
A(:,:,H) is upper quasi-triangular and all the other factors
A(:,:,I) are upper triangular. Optionally, the 2-by-2 triangular
matrices corresponding to 2-by-2 diagonal blocks in A(:,:,H)
are so reduced that their product is a 2-by-2 diagonal matrix.

If COMPQ = 'U' or COMPQ = 'I', then the orthogonal factors are
computed and stored in the array Q so that for S(I) = 1,

                    `T`
        `Q(:,:,I)(in)   A(:,:,I)(in)   Q(:,:,MOD(I,K)+1)(in)`
                                                            `T  (1)`
    `=   Q(:,:,I)(out)  A(:,:,I)(out)  Q(:,:,MOD(I,K)+1)(out),`

and for S(I) = -1,

                             `T`
        `Q(:,:,MOD(I,K)+1)(in)   A(:,:,I)(in)   Q(:,:,I)(in)`
                                                            `T  (2)`
    `=   Q(:,:,MOD(I,K)+1)(out)  A(:,:,I)(out)  Q(:,:,I)(out).`

A partial generation of the orthogonal factors can be realized
via the array QIND.



See the SLICOT documentation for details.
"""
function mb03bd! end

"""

Apply at most 20 iterations of a real single shifted
periodic QZ algorithm to the 2-by-2 product of matrices stored
in the array A.



See the SLICOT documentation for details.
"""
function mb03be! end

"""

Apply at most 20 iterations of a real single shifted
periodic QZ algorithm to the 2-by-2 product of matrices stored
in the array A. The Hessenberg matrix is the last one of the
formal matrix product.



See the SLICOT documentation for details.
"""
function mb03bf! end

"""

Compute the eigenvalues of the 2-by-2 trailing submatrix of the
matrix product

             `S(1)           S(2)                 S(K)`
     `A(:,:,1)     * A(:,:,2)     * ... * A(:,:,K)`

where A(:,:,AMAP(K)) is upper Hessenberg and A(:,:,AMAP(i)),
1 <= i < K, is upper triangular. All factors to be inverted
(depending on S and SINV) are assumed nonsingular. Moreover,
AMAP(K) is either 1 or K.



See the SLICOT documentation for details.
"""
function mb03bg! end

"""

Find the eigenvalues of the complex generalized matrix product

             `S(1)           S(2)                 S(K)`
     `A(:,:,1)     * A(:,:,2)     * ... * A(:,:,K)    ,  S(1) = 1,`

where A(:,:,1) is upper Hessenberg and A(:,:,i) is upper
triangular, i = 2, ..., K, using a single-shift version of the
periodic QZ method. In addition, A may be reduced to periodic
Schur form by unitary transformations: all factors A(:,:,i) become
upper triangular.

If COMPQ = 'V' or COMPQ = 'I', then the unitary factors are
computed and stored in the array Q so that for S(I) = 1,

                    `H`
        `Q(:,:,I)(in)   A(:,:,I)(in)   Q(:,:,MOD(I,K)+1)(in)`
                     `H                                        (1)`
    `=   Q(:,:,I)(out)  A(:,:,I)(out)  Q(:,:,MOD(I,K)+1)(out),`

and for S(I) = -1,

                             `H`
        `Q(:,:,MOD(I,K)+1)(in)   A(:,:,I)(in)   Q(:,:,I)(in)`
                              `H                               (2)`
    `=   Q(:,:,MOD(I,K)+1)(out)  A(:,:,I)(out)  Q(:,:,I)(out).`



See the SLICOT documentation for details.
"""
function mb03bz! end

"""

Compute orthogonal matrices Q1, Q2, Q3 for a real 2-by-2,
3-by-3, or 4-by-4 regular block upper triangular pencil

               `( A11 A12 ) ( B11 B12 )     ( D11 D12 )`
  `aAB - bD = a (         ) (         ) - b (         ),        (1)`
               `(  0  A22 ) (  0  B22 )     (  0  D22 )`

such that the pencil a(Q3' A Q2 )(Q2' B Q1 ) - b(Q3' D Q1) is
still in block upper triangular form, but the eigenvalues in
Spec(A11 B11, D11), Spec(A22 B22, D22) are exchanged, where
Spec(X,Y) denotes the spectrum of the matrix pencil (X,Y), and M'
denotes the transpose of the matrix M.

Optionally, to upper triangularize the real regular pencil in
block lower triangular form

             `( A11  0  ) ( B11  0  )     ( D11  0  )`
aAB - bD = a (         ) (         ) - b (         ),          (2)
             `( A21 A22 ) ( B21 B22 )     ( D21 D22 )`

while keeping the eigenvalues in the same diagonal position.



See the SLICOT documentation for details.
"""
function mb03cd! end

"""

Compute unitary matrices Q1, Q2, and Q3 for a complex 2-by-2
regular pencil aAB - bD, with A, B, D upper triangular, such that
Q3' A Q2, Q2' B Q1, Q3' D Q1 are still upper triangular, but the
eigenvalues are in reversed order. The matrices Q1, Q2, and Q3 are
represented by

     `(  CO1  SI1  )       (  CO2  SI2  )       (  CO3  SI3  )`
Q1 = (            ), Q2 = (            ), Q3 = (            ).
     `( -SI1' CO1  )       ( -SI2' CO2  )       ( -SI3' CO3  )`

The notation M' denotes the conjugate transpose of the matrix M.



See the SLICOT documentation for details.
"""
function mb03cz! end

"""

Compute orthogonal matrices Q1 and Q2 for a real 2-by-2,
3-by-3, or 4-by-4 regular block upper triangular pencil

               `( A11 A12 )     ( B11 B12 )`
  `aA - bB =  a (         ) - b (         ),                    (1)`
               `(  0  A22 )     (  0  B22 )`

such that the pencil a(Q2' A Q1) - b(Q2' B Q1) is still in block
upper triangular form, but the eigenvalues in Spec(A11, B11),
Spec(A22, B22) are exchanged, where Spec(X,Y) denotes the spectrum
of the matrix pencil (X,Y) and the notation M' denotes the
transpose of the matrix M.

Optionally, to upper triangularize the real regular pencil in
block lower triangular form

              `( A11  0  )     ( B11  0  )`
  `aA - bB = a (         ) - b (         ),                     (2)`
              `( A21 A22 )     ( B21 B22 )`

while keeping the eigenvalues in the same diagonal position.



See the SLICOT documentation for details.
"""
function mb03dd! end

"""

Compute unitary matrices Q1 and Q2 for a complex 2-by-2 regular
pencil aA - bB with A, B upper triangular, such that
Q2' (aA - bB) Q1 is still upper triangular but the eigenvalues are
in reversed order. The matrices Q1 and Q2 are represented by

     `(  CO1  SI1  )       (  CO2  SI2  )`
Q1 = (            ), Q2 = (            ).
     `( -SI1' CO1  )       ( -SI2' CO2  )`

The notation M' denotes the conjugate transpose of the matrix M.



See the SLICOT documentation for details.
"""
function mb03dz! end

"""

Compute orthogonal matrices Q1, Q2, Q3 for a real 2-by-2 or
4-by-4 regular pencil

               `( A11  0  ) ( B11  0  )     (  0  D12 )`
  `aAB - bD = a (         ) (         ) - b (         ),        (1)`
               `(  0  A22 ) (  0  B22 )     ( D21  0  )`

such that Q3' A Q2 and Q2' B Q1 are upper triangular, Q3' D Q1 is
upper quasi-triangular, and the eigenvalues with negative real
parts (if there are any) are allocated on the top. The notation M'
denotes the transpose of the matrix M. The submatrices A11, A22,
B11, B22 and D12 are upper triangular. If D21 is 2-by-2, then all
other blocks are nonsingular and the product
   `-1        -1    -1        -1`
A22   D21 B11   A11   D12 B22   has a pair of complex conjugate
eigenvalues.



See the SLICOT documentation for details.
"""
function mb03ed! end

"""

Compute orthogonal matrices Q1 and Q2 for a real 2-by-2 or
4-by-4 regular pencil

              `( A11  0  )     (  0  B12 )`
  `aA - bB = a (         ) - b (         ),                     (1)`
              `(  0  A22 )     ( B21  0  )`

such that Q2' A Q1 is upper triangular, Q2' B Q1 is upper quasi-
triangular, and the eigenvalues with negative real parts (if there
are any) are allocated on the top. The notation M' denotes the
transpose of the matrix M. The submatrices A11, A22, and B12 are
upper triangular. If B21 is 2-by-2, then all the other blocks are
                               `-1        -1`
nonsingular and the product A11   B12 A22   B21 has a pair of
complex conjugate eigenvalues.



See the SLICOT documentation for details.
"""
function mb03fd! end

"""

Compute the eigenvalues of a complex N-by-N skew-Hamiltonian/
Hamiltonian pencil aS - bH, with

                        `(  B  F  )      (  Z11  Z12  )`
  `S = J Z' J' Z and H = (        ), Z = (            ),`
                        `(  G -B' )      (  Z21  Z22  )`
                                                              `(1)`
      `(  0  I  )`
  `J = (        ).`
      `( -I  0  )`

The structured Schur form of the embedded real skew-Hamiltonian/

skew-Hamiltonian pencil, a`B_S` - b`B_T`, with `B_S` = J `B_Z`' J' `B_Z`,

        `(  Re(Z11)  -Im(Z11)  |  Re(Z12)  -Im(Z12)  )`
        `(                     |                     )`
        `(  Im(Z11)   Re(Z11)  |  Im(Z12)   Re(Z12)  )`
        `(                     |                     )`
  `B_Z = (---------------------+---------------------) ,`
        `(                     |                     )`
        `(  Re(Z21)  -Im(Z21)  |  Re(Z22)  -Im(Z22)  )`
        `(                     |                     )`
        `(  Im(Z21)   Re(Z21)  |  Im(Z22)   Re(Z22)  )`
                                                               `(2)`
        `( -Im(B)  -Re(B)  | -Im(F)  -Re(F)  )`
        `(                 |                 )`
        `(  Re(B)  -Im(B)  |  Re(F)  -Im(F)  )`
        `(                 |                 )`
  `B_T = (-----------------+-----------------) ,  T = i*H,`
        `(                 |                 )`
        `( -Im(G)  -Re(G)  | -Im(B')  Re(B') )`
        `(                 |                 )`
        `(  Re(G)  -Im(G)  | -Re(B') -Im(B') )`

is determined and used to compute the eigenvalues. Optionally, if
COMPQ = 'C', an orthonormal basis of the right deflating subspace,
Def_-(S, H), of the pencil aS - bH in (1), corresponding to the
eigenvalues with strictly negative real part, is computed. Namely,
after transforming a`B_S` - b`B_H`, in the factored form, by unitary
matrices, we have `B_S`out = J `B_Z`out' J' `B_Z`out,

           `( BA  BD  )              ( BB  BF  )`
  `B_Zout = (         ) and B_Hout = (         ),               (3)`
           `(  0  BC  )              (  0 -BB' )`

and the eigenvalues with strictly negative real part of the
complex pencil a`B_S`out - b`B_H`out are moved to the top. The 
notation M' denotes the conjugate transpose of the matrix M.
Optionally, if COMPU = 'C', an orthonormal basis of the companion
subspace, range(P_U) [1], which corresponds to the eigenvalues
with negative real part, is computed. The embedding doubles the
multiplicities of the eigenvalues of the pencil aS - bH.



See the SLICOT documentation for details.
"""
function mb03fz! end

"""

Compute an orthogonal matrix Q and an orthogonal symplectic
matrix U for a real regular 2-by-2 or 4-by-4 skew-Hamiltonian/
Hamiltonian pencil a J B' J' B - b D with

      `( B11  B12 )      (  D11  D12  )      (  0  I  )`
  `B = (          ), D = (            ), J = (        ),`
      `(  0   B22 )      (   0  -D11' )      ( -I  0  )`

such that J Q' J' D Q and U' B Q keep block triangular form, but
the eigenvalues are reordered. The notation M' denotes the
transpose of the matrix M.



See the SLICOT documentation for details.
"""
function mb03gd! end

"""

Compute a unitary matrix Q and a unitary symplectic matrix U
for a complex regular 2-by-2 skew-Hamiltonian/Hamiltonian pencil
aS - bH with S = J Z' J' Z, where

       `(  Z11  Z12  )         (  H11  H12  )`
   `Z = (            ) and H = (            ),`
       `(   0   Z22  )         (   0  -H11' )`

such that U' Z Q, (J Q J' )' H Q are both upper triangular, but the  
eigenvalues of (J Q J')' ( aS - bH ) Q are in reversed order.
The matrices Q and U are represented by

       `(  CO1  SI1  )         (  CO2  SI2  )`
   `Q = (            ) and U = (            ), respectively.`
       `( -SI1' CO1  )         ( -SI2' CO2  )`

The notation M' denotes the conjugate transpose of the matrix M.



See the SLICOT documentation for details.
"""
function mb03gz! end

"""

Determine an orthogonal matrix Q, for a real regular 2-by-2 or
4-by-4 skew-Hamiltonian/Hamiltonian pencil

                `( A11 A12  )     ( B11  B12  )`
    `aA - bB = a (          ) - b (           )`
                `(  0  A11' )     (  0  -B11' )`

in structured Schur form, such that  J Q' J' (aA - bB) Q  is still
in structured Schur form but the eigenvalues are exchanged. The
notation M' denotes the transpose of the matrix M.



See the SLICOT documentation for details.
"""
function mb03hd! end

"""

Compute a unitary matrix Q for a complex regular 2-by-2 
skew-Hamiltonian/Hamiltonian pencil aS - bH with

    `(  S11  S12  )        (  H11  H12  )`
S = (            ),   H = (            ),
    `(   0   S11' )        (   0  -H11' )`

such that J Q' J' (aS - bH) Q is upper triangular but the
eigenvalues are in reversed order. The matrix Q is represented by

    `(  CO  SI  )`
Q = (          ).
    `( -SI' CO  )`

The notation M' denotes the conjugate transpose of the matrix M.



See the SLICOT documentation for details.
"""
function mb03hz! end

"""

Move the eigenvalues with strictly negative real parts of an
N-by-N real skew-Hamiltonian/Hamiltonian pencil aS - bH in
structured Schur form, with

                     `(  0  I  )      (  A  D  )      (  B  F  )`
  `S = J Z' J' Z, J = (        ), Z = (        ), H = (        ),`
                     `( -I  0  )      (  0  C  )      (  0 -B' )`

to the leading principal subpencil, while keeping the triangular
form. Above, A is upper triangular, B is upper quasi-triangular,
and C is lower triangular.
The matrices Z and H are transformed by an orthogonal symplectic
matrix U and an orthogonal matrix Q such that

                  `(  Aout  Dout  )`
  `Zout = U' Z Q = (              ), and`
                  `(    0   Cout  )`
                                                               `(1)`
                       `(  Bout  Fout  )`
  `Hout = J Q' J' H Q = (              ),`
                       `(    0  -Bout' )`

where Aout, Bout and Cout remain in triangular form. The notation
M' denotes the transpose of the matrix M.
Optionally, if COMPQ = 'I' or COMPQ = 'U', the orthogonal matrix Q
that fulfills (1) is computed.
Optionally, if COMPU = 'I' or COMPU = 'U', the orthogonal
symplectic matrix

      `(  U1  U2  )`
  `U = (          )`
      `( -U2  U1  )`

that fulfills (1) is computed.



See the SLICOT documentation for details.
"""
function mb03id! end

"""

Move the eigenvalues with strictly negative real parts of an
N-by-N complex skew-Hamiltonian/Hamiltonian pencil aS - bH in
structured Schur form, with

                           `(  0  I  )`
  `S = J Z' J' Z, where J = (        ),`
                           `( -I  0  )`

to the leading principal subpencil, while keeping the triangular
form. On entry, we have

      `(  A  D  )      (  B  F  )`
  `Z = (        ), H = (        ),`
      `(  0  C  )      (  0 -B' )`

where A and B are upper triangular and C is lower triangular.
Z and H are transformed by a unitary symplectic matrix U and a
unitary matrix Q such that

                  `(  Aout  Dout  )`
  `Zout = U' Z Q = (              ), and`
                  `(    0   Cout  )`
                                                               `(1)`
                       `(  Bout  Fout  )`
  `Hout = J Q' J' H Q = (              ), `
                       `(    0  -Bout' )`

where Aout, Bout and Cout remain in triangular form. The notation
M' denotes the conjugate transpose of the matrix M.
Optionally, if COMPQ = 'I' or COMPQ = 'U', the unitary matrix Q
that fulfills (1) is computed.
Optionally, if COMPU = 'I' or COMPU = 'U', the unitary symplectic
matrix 

      `(  U1  U2  )`
  `U = (          )`
      `( -U2  U1  )   `

that fulfills (1) is computed.



See the SLICOT documentation for details.
"""
function mb03iz! end

"""

Move the eigenvalues with strictly negative real parts of an
N-by-N real skew-Hamiltonian/Hamiltonian pencil aS - bH in
structured Schur form,

      `(  A  D  )      (  B  F  )`
  `S = (        ), H = (        ),`
      `(  0  A' )      (  0 -B' )`

with A upper triangular and B upper quasi-triangular, to the
leading principal subpencil, while keeping the triangular form.
The notation M' denotes the transpose of the matrix M.
The matrices S and H are transformed by an orthogonal matrix Q
such that

                       `(  Aout  Dout  )  `
  `Sout = J Q' J' S Q = (              ),`
                       `(    0   Aout' )  `
                                                               `(1)`
                       `(  Bout  Fout  )           (  0  I  )`
  `Hout = J Q' J' H Q = (              ), with J = (        ),`
                       `(  0    -Bout' )           ( -I  0  )`

where Aout is upper triangular and Bout is upper quasi-triangular.
Optionally, if COMPQ = 'I' or COMPQ = 'U', the orthogonal matrix Q
that fulfills (1), is computed.



See the SLICOT documentation for details.
"""
function mb03jd! end

"""

Move the eigenvalues with strictly negative real parts of an
N-by-N real skew-Hamiltonian/Hamiltonian pencil aS - bH in
structured Schur form,

      `(  A  D  )      (  B  F  )`
  `S = (        ), H = (        ),`
      `(  0  A' )      (  0 -B' )`

with A upper triangular and B upper quasi-triangular, to the
leading principal subpencil, while keeping the triangular form.
The notation M' denotes the transpose of the matrix M.
The matrices S and H are transformed by an orthogonal matrix Q
such that

                       `(  Aout  Dout  )  `
  `Sout = J Q' J' S Q = (              ),`
                       `(    0   Aout' )  `
                                                               `(1)`
                       `(  Bout  Fout  )           (  0  I  )`
  `Hout = J Q' J' H Q = (              ), with J = (        ),`
                       `(  0    -Bout' )           ( -I  0  )`

where Aout is upper triangular and Bout is upper quasi-triangular.
Optionally, if COMPQ = 'I' or COMPQ = 'U', the orthogonal matrix Q
that fulfills (1), is computed.



See the SLICOT documentation for details.
"""
function mb03jp! end

"""

Move the eigenvalues with strictly negative real parts of an
N-by-N complex skew-Hamiltonian/Hamiltonian pencil aS - bH in
structured Schur form to the leading principal subpencil, while
keeping the triangular form. On entry, we have

      `(  A  D  )      (  B  F  )`
  `S = (        ), H = (        ),`
      `(  0  A' )      (  0 -B' )`

where A and B are upper triangular.
S and H are transformed by a unitary matrix Q such that

                       `(  Aout  Dout  )`
  `Sout = J Q' J' S Q = (              ), and`
                       `(    0   Aout' )`
                                                               `(1)`
                       `(  Bout  Fout  )           (  0  I  )`
  `Hout = J Q' J' H Q = (              ), with J = (        ),`
                       `(    0  -Bout' )           ( -I  0  )`

where Aout and Bout remain in upper triangular form. The notation
M' denotes the conjugate transpose of the matrix M.
Optionally, if COMPQ = 'I' or COMPQ = 'U', the unitary matrix Q
that fulfills (1) is computed.



See the SLICOT documentation for details.
"""
function mb03jz! end

"""

Reorder the diagonal blocks of the formal matrix product

   `T22_K^S(K) * T22_K-1^S(K-1) * ... * T22_1^S(1),             (1)`

of length K, in the generalized periodic Schur form

         `[  T11_k  T12_k  T13_k  ]`
   `T_k = [    0    T22_k  T23_k  ],    k = 1, ..., K,          (2)`
         `[    0      0    T33_k  ]`

where

- the submatrices T11_k are NI(k+1)-by-NI(k), if S(k) = 1, or
  `NI(k)-by-NI(k+1), if S(k) = -1, and contain dimension-induced`
  `infinite eigenvalues,`
- the submatrices T22_k are NC-by-NC and contain core eigenvalues,
  `which are generically neither zero nor infinite,`
- the submatrices T33_k contain dimension-induced zero
  `eigenvalues,`

such that the block with starting row index IFST in (1) is moved
to row index ILST. The indices refer to the T22_k submatrices.

Optionally, the transformation matrices `Q_1`,...,`Q_K` from the
reduction into generalized periodic Schur form are updated with
respect to the performed reordering.



See the SLICOT documentation for details.
"""
function mb03ka! end

"""

Reorder the diagonal blocks of the formal matrix product

   `T22_K^S(K) * T22_K-1^S(K-1) * ... * T22_1^S(1)              (1)`

of length K in the generalized periodic Schur form

         `[  T11_k  T12_k  T13_k  ]`
   `T_k = [    0    T22_k  T23_k  ],    k = 1, ..., K,          (2)`
         `[    0      0    T33_k  ]`

where

- the submatrices T11_k are NI(k+1)-by-NI(k), if S(k) = 1, or
  `NI(k)-by-NI(k+1), if S(k) = -1, and contain dimension-induced`
  `infinite eigenvalues,`
- the submatrices T22_k are NC-by-NC and contain core eigenvalues,
  `which are generically neither zero nor infinite,`
- the submatrices T33_k contain dimension-induced zero
  `eigenvalues,`

such that pairs of adjacent diagonal blocks of sizes 1 and/or 2 in
the product (1) are swapped.

Optionally, the transformation matrices `Q_1`,...,`Q_K` from the
reduction into generalized periodic Schur form are updated with
respect to the performed reordering.



See the SLICOT documentation for details.
"""
function mb03kb! end

"""

Reduce a 2-by-2 general, formal matrix product A of length K,

   `A_K^s(K) * A_K-1^s(K-1) * ... * A_1^s(1),`

to the periodic Hessenberg-triangular form using a K-periodic
sequence of elementary reflectors (Householder matrices). The
matrices A_k, k = 1, ..., K, are stored in the N-by-N-by-K array A
starting in the R-th row and column, and N can be 3 or 4.

Each elementary reflector H_k is represented as

   `H_k = I - tau_k * v_k * v_k',                               (1)`

where I is the 2-by-2 identity, tau_k is a real scalar, and v_k is
a vector of length 2, k = 1,...,K, and it is constructed such that
the following holds for k = 1,...,K:

       `H_{k+1} * A_k * H_k = T_k, if s(k) = 1,`
                                                               `(2)`
       `H_k * A_k * H_{k+1} = T_k, if s(k) = -1,`

with H_{K+1} = `H_1` and all `T_k` upper triangular except for
T_{khess} which is full. Clearly,

   `T_K^s(K) *...* T_1^s(1) = H_1 * A_K^s(K) *...* A_1^s(1) * H_1.`

The reflectors are suitably applied to the whole, extended N-by-N
matrices `Ae_k`, not only to the submatrices `A_k`, k = 1, ..., K.



See the SLICOT documentation for details.
"""
function mb03kc! end

"""

Reorder the diagonal blocks of the formal matrix product

   `T22_K^S(K) * T22_K-1^S(K-1) * ... * T22_1^S(1),             (1)`

of length K, in the generalized periodic Schur form,

         `[  T11_k  T12_k  T13_k  ]`
   `T_k = [    0    T22_k  T23_k  ],    k = 1, ..., K,          (2)`
         `[    0      0    T33_k  ]`

where

- the submatrices T11_k are NI(k+1)-by-NI(k), if S(k) = 1, or
  `NI(k)-by-NI(k+1), if S(k) = -1, and contain dimension-induced`
  `infinite eigenvalues,`
- the submatrices T22_k are NC-by-NC and contain core eigenvalues,
  `which are generically neither zero nor infinite,`
- the submatrices T33_k contain dimension-induced zero
  `eigenvalues,`

such that the M selected eigenvalues pointed to by the logical
vector SELECT end up in the leading part of the matrix sequence
T22_k.

Given that N(k) = N(k+1) for all k where S(k) = -1, the T11_k are
void and the first M columns of the updated orthogonal
transformation matrix sequence `Q_1`, ..., `Q_K` span a periodic
deflating subspace corresponding to the same eigenvalues.



See the SLICOT documentation for details.
"""
function mb03kd! end

"""

Solve small periodic Sylvester-like equations (PSLE)

 `op(A(i))*X( i ) + isgn*X(i+1)*op(B(i)) = -scale*C(i), S(i) =  1,`
 `op(A(i))*X(i+1) + isgn*X( i )*op(B(i)) = -scale*C(i), S(i) = -1.`

i = 1, ..., K, where op(A) means A or A**T, for the K-periodic
matrix sequence X(i) = X(i+K), where A, B and C are K-periodic
matrix sequences and A and B are in periodic real Schur form. The
matrices A(i) are M-by-M and B(i) are N-by-N, with 1 <= M, N <= 2.



See the SLICOT documentation for details.
"""
function mb03ke! end

"""

Compute the relevant eigenvalues of a real N-by-N skew-
Hamiltonian/Hamiltonian pencil aS - bH, with

      `(  A  D  )         (  B  F  )`
  `S = (        ) and H = (        ),                           (1)`
      `(  E  A' )         (  G -B' )`

where the notation M' denotes the transpose of the matrix M.
Optionally, if COMPQ = 'C', an orthogonal basis of the right
deflating subspace of aS - bH corresponding to the eigenvalues
with strictly negative real part is computed.



See the SLICOT documentation for details.
"""
function mb03ld! end

"""

Compute the relevant eigenvalues of a real N-by-N skew-
Hamiltonian/Hamiltonian pencil aS - bH, with

                              `(  B  F  )      (  0  I  )`
  `S = T Z = J Z' J' Z and H = (        ), J = (        ),      (1)`
                              `(  G -B' )      ( -I  0  )`

where the notation M' denotes the transpose of the matrix M.
Optionally, if COMPQ = 'C', an orthogonal basis of the right
deflating subspace of aS - bH corresponding to the eigenvalues
with strictly negative real part is computed. Optionally, if
COMPU = 'C', an orthonormal basis of the companion subspace,
range(P_U) [1], which corresponds to the eigenvalues with strictly
negative real part, is computed.



See the SLICOT documentation for details.
"""
function mb03lf! end

"""

Compute the relevant eigenvalues of a real N-by-N skew-
Hamiltonian/Hamiltonian pencil aS - bH, with

      `(  A  D  )         (  B  F  )`
  `S = (        ) and H = (        ),                           (1)`
      `(  E  A' )         (  G -B' )`

where the notation M' denotes the transpose of the matrix M.
Optionally, if COMPQ = 'C', an orthogonal basis of the right
deflating subspace of aS - bH corresponding to the eigenvalues
with strictly negative real part is computed.



See the SLICOT documentation for details.
"""
function mb03lp! end

"""

Compute the eigenvalues of a complex N-by-N skew-Hamiltonian/
Hamiltonian pencil aS - bH, with

      `(  A  D  )         (  B  F  )`
  `S = (        ) and H = (        ).                           (1)`
      `(  E  A' )         (  G -B' )`

The structured Schur form of the embedded real skew-Hamiltonian/
skew-Hamiltonian pencil a`B_S` - b`B_T`, defined as

        `(  Re(A)  -Im(A)  |  Re(D)  -Im(D)  )`
        `(                 |                 )`
        `(  Im(A)   Re(A)  |  Im(D)   Re(D)  )`
        `(                 |                 )`
  `B_S = (-----------------+-----------------) , and`
        `(                 |                 )`
        `(  Re(E)  -Im(E)  |  Re(A')  Im(A') )`
        `(                 |                 )`
        `(  Im(E)   Re(E)  | -Im(A')  Re(A') )`
                                                               `(2)`
        `( -Im(B)  -Re(B)  | -Im(F)  -Re(F)  )`
        `(                 |                 )`
        `(  Re(B)  -Im(B)  |  Re(F)  -Im(F)  )`
        `(                 |                 )`
  `B_T = (-----------------+-----------------) ,  T = i*H,`
        `(                 |                 )`
        `( -Im(G)  -Re(G)  | -Im(B')  Re(B') )`
        `(                 |                 )`
        `(  Re(G)  -Im(G)  | -Re(B') -Im(B') )`

is determined and used to compute the eigenvalues. The notation M'
denotes the conjugate transpose of the matrix M. Optionally,
if COMPQ = 'C', an orthonormal basis of the right deflating
subspace of the pencil aS - bH, corresponding to the eigenvalues
with strictly negative real part, is computed. Namely, after
transforming a`B_S` - b`B_H` by unitary matrices, we have

           `( BA  BD  )              ( BB  BF  )`
  `B_Sout = (         ) and B_Hout = (         ),               (3)`
           `(  0  BA' )              (  0 -BB' )`

and the eigenvalues with strictly negative real part of the
complex pencil a`B_S`out - b`B_H`out are moved to the top. The
embedding doubles the multiplicities of the eigenvalues of the
pencil aS - bH.



See the SLICOT documentation for details.
"""
function mb03lz! end

"""

Compute an upper bound THETA using a bisection method such that
the bidiagonal matrix

         `|q(1) e(1)  0    ...   0   |`
         `| 0   q(2) e(2)        .   |`
     `J = | .                    .   |`
         `| .                  e(N-1)|`
         `| 0   ...        ...  q(N) |`

has precisely L singular values less than or equal to THETA plus
a given tolerance TOL.

This routine is mainly intended to be called only by other SLICOT
routines.



See the SLICOT documentation for details.
"""
function mb03md! end

"""

Compute the absolute minimal value of NX elements in an array.
The function returns the value zero if NX < 1.



See the SLICOT documentation for details.
"""
function mb03my! end

"""

Find the number of singular values of the bidiagonal matrix

         `|q(1) e(1)  .    ...    0   |`
         `| 0   q(2) e(2)         .   |`
     `J = | .                     .   |`
         `| .                   e(N-1)|`
         `| 0   ...     ...   0  q(N) |`

which are less than or equal to a given bound THETA.

This routine is intended to be called only by other SLICOT
routines.



See the SLICOT documentation for details.
"""
function mb03nd! end

"""

Compute the smallest singular value of A - jwI.



See the SLICOT documentation for details.
"""
function mb03ny! end

"""

Compute (optionally) a rank-revealing QR factorization of a
real general M-by-N matrix  A,  which may be rank-deficient,
and estimate its effective rank using incremental condition
estimation.

The routine uses a QR factorization with column pivoting:
   `A * P = Q * R,  where  R = [ R11 R12 ],`
                              `[  0  R22 ]`
with R11 defined as the largest leading submatrix whose estimated
condition number is less than 1/RCOND.  The order of R11, RANK,
is the effective rank of A.

MB03OD  does not perform any scaling of the matrix A.



See the SLICOT documentation for details.
"""
function mb03od! end

"""

Compute a rank-revealing QR factorization of a real general
M-by-N matrix  A,  which may be rank-deficient, and estimate its
effective rank using incremental condition estimation.

The routine uses a truncated QR factorization with column pivoting
                              `[ R11 R12 ]`
   `A * P = Q * R,  where  R = [         ],`
                              `[  0  R22 ]`
with R11 defined as the largest leading upper triangular submatrix
whose estimated condition number is less than 1/RCOND.  The order
of R11, RANK, is the effective rank of A.  Condition estimation is
performed during the QR factorization process.  Matrix R22 is full
(but of small norm), or empty.

MB03OY  does not perform any scaling of the matrix A.



See the SLICOT documentation for details.
"""
function mb03oy! end

"""

Compute (optionally) a rank-revealing RQ factorization of a
real general M-by-N matrix  A,  which may be rank-deficient,
and estimate its effective rank using incremental condition
estimation.

The routine uses an RQ factorization with row pivoting:
   `P * A = R * Q,  where  R = [ R11 R12 ],`
                              `[  0  R22 ]`
with R22 defined as the largest trailing submatrix whose estimated
condition number is less than 1/RCOND.  The order of R22, RANK,
is the effective rank of A.

MB03PD  does not perform any scaling of the matrix A.



See the SLICOT documentation for details.
"""
function mb03pd! end

"""

Compute a rank-revealing RQ factorization of a real general
M-by-N matrix  A,  which may be rank-deficient, and estimate its
effective rank using incremental condition estimation.

The routine uses a truncated RQ factorization with row pivoting:
                              `[ R11 R12 ]`
   `P * A = R * Q,  where  R = [         ],`
                              `[  0  R22 ]`
with R22 defined as the largest trailing upper triangular
submatrix whose estimated condition number is less than 1/RCOND.
The order of R22, RANK, is the effective rank of A.  Condition
estimation is performed during the RQ factorization process.
Matrix R11 is full (but of small norm), or empty.

MB03PY  does not perform any scaling of the matrix A.



See the SLICOT documentation for details.
"""
function mb03py! end

"""

Reorder the diagonal blocks of a principal submatrix of an
upper quasi-triangular matrix A together with their eigenvalues by
constructing an orthogonal similarity transformation UT.
After reordering, the leading block of the selected submatrix of A
has eigenvalues in a suitably defined domain of interest, usually
related to stability/instability in a continuous- or discrete-time
sense.



See the SLICOT documentation for details.
"""
function mb03qd! end

"""

Reorder the diagonal blocks of a principal subpencil of an
upper quasi-triangular matrix pencil A-lambda*E together with
their generalized eigenvalues, by constructing orthogonal
similarity transformations UT and VT.
After reordering, the leading block of the selected subpencil of
A-lambda*E has generalized eigenvalues in a suitably defined
domain of interest, usually related to stability/instability in a
continuous- or discrete-time sense.



See the SLICOT documentation for details.
"""
function mb03qg! end

"""

Compute the eigenvalues of an upper quasi-triangular matrix
pencil.



See the SLICOT documentation for details.
"""
function mb03qv! end

"""

Compute the eigenvalues of a selected 2-by-2 diagonal block
pair of an upper quasi-triangular pencil, to reduce the selected
block pair to the standard form and to split it in the case of
real eigenvalues, by constructing orthogonal matrices UT and VT.
The transformations UT and VT are applied to the pair (A,E) by
computing (UT'*A*VT, UT'*E*VT ), to the matrices U and V,
by computing U*UT and V*VT.



See the SLICOT documentation for details.
"""
function mb03qw! end

"""

Compute the eigenvalues of an upper quasi-triangular matrix.



See the SLICOT documentation for details.
"""
function mb03qx! end

"""

Compute the eigenvalues of a selected 2-by-2 diagonal block
of an upper quasi-triangular matrix, to reduce the selected block
to the standard form and to split the block in the case of real
eigenvalues by constructing an orthogonal transformation UT.
This transformation is applied to A (by similarity) and to
another matrix U from the right.



See the SLICOT documentation for details.
"""
function mb03qy! end

"""

Reduce a matrix A in real Schur form to a block-diagonal form
using well-conditioned non-orthogonal similarity transformations.
The condition numbers of the transformations used for reduction
are roughly bounded by PMAX*PMAX, where PMAX is a given value.
The transformations are optionally postmultiplied in a given
matrix X. The real Schur form is optionally ordered, so that
clustered eigenvalues are grouped in the same block.



See the SLICOT documentation for details.
"""
function mb03rd! end

"""

Reorder the diagonal blocks of the principal submatrix between
the indices KL and KU (KU >= KL) of a real Schur form matrix A
together with their eigenvalues, using orthogonal similarity
transformations, such that the block specified by KU is moved in
the position KL. The transformations are optionally postmultiplied
in a given matrix X.



See the SLICOT documentation for details.
"""
function mb03rx! end

"""

Solve the Sylvester equation -AX + XB = C, where A and B are
M-by-M and N-by-N matrices, respectively, in real Schur form.

This routine is intended to be called only by SLICOT Library
routine MB03RD. For efficiency purposes, the computations are
aborted when the infinity norm of an elementary submatrix of X is
greater than a given value PMAX.



See the SLICOT documentation for details.
"""
function mb03ry! end

"""

Compute the eigenvalues of an N-by-N square-reduced Hamiltonian
matrix

           `( A'   G'  )`
    `H'  =  (        T ).                                       (1)`
           `( Q'  -A'  )`

Here, A' is an N-by-N matrix, and G' and Q' are symmetric N-by-N
matrices.  It is assumed without a check that H' is square-
reduced, i.e., that

      `2    ( A''   G'' )`
    `H'  =  (         T )    with A'' upper Hessenberg.         (2)`
           `( 0    A''  )`

                       `T                2`
(Equivalently, Q'A'- A' Q' = 0, A'' = A' + G'Q', and for i > j+1,
 `A''(i,j) = 0.)  Ordinarily, H' is the output from SLICOT Library`
routine MB04ZD. The eigenvalues of H' are computed as the square
roots of the eigenvalues of A''.



See the SLICOT documentation for details.
"""
function mb03sd! end

"""

Reorder a matrix X in skew-Hamiltonian Schur form:

              `[  A   G  ]          T`
        `X  =  [       T ],   G = -G,`
              `[  0   A  ]`

or in Hamiltonian Schur form:

              `[  A   G  ]          T`
        `X  =  [       T ],   G =  G,`
              `[  0  -A  ]`

where A is in upper quasi-triangular form, so that a selected
cluster of eigenvalues appears in the leading diagonal blocks
of the matrix A (in X) and the leading columns of [ U1; -U2 ] form
an orthonormal basis for the corresponding right invariant
subspace.

If X is skew-Hamiltonian, then each eigenvalue appears twice; one
copy corresponds to the j-th diagonal element and the other to the
(n+j)-th diagonal element of X. The logical array LOWER controls
which copy is to be reordered to the leading part of A.

If X is Hamiltonian then the eigenvalues appear in pairs
(lambda,-lambda); lambda corresponds to the j-th diagonal
element and -lambda to the (n+j)-th diagonal element of X.
The logical array LOWER controls whether lambda or -lambda is to
be reordered to the leading part of A.

The matrix A must be in Schur canonical form (as returned by the
LAPACK routine DHSEQR), that is, block upper triangular with
1-by-1 and 2-by-2 diagonal blocks; each 2-by-2 diagonal block has
its diagonal elements equal and its off-diagonal elements of
opposite sign.



See the SLICOT documentation for details.
"""
function mb03td! end

"""

Swap diagonal blocks A11 and A22 of order 1 or 2 in the upper
quasi-triangular matrix A contained in a skew-Hamiltonian matrix

              `[  A   G  ]          T`
        `X  =  [       T ],   G = -G,`
              `[  0   A  ]`

or in a Hamiltonian matrix

              `[  A   G  ]          T`
        `X  =  [       T ],   G =  G.`
              `[  0  -A  ]`

This routine is a modified version of the LAPACK subroutine
DLAEX2.

The matrix A must be in Schur canonical form (as returned by the
LAPACK routine DHSEQR), that is, block upper triangular with
1-by-1 and 2-by-2 diagonal blocks; each 2-by-2 diagonal block has
its diagonal elements equal and its off-diagonal elements of
opposite sign.



See the SLICOT documentation for details.
"""
function mb03ts! end

"""

Compute all, or part, of the singular value decomposition of a
real upper triangular matrix.

The N-by-N upper triangular matrix A is factored as  A = Q*S*P',
where Q and P are N-by-N orthogonal matrices and S is an
N-by-N diagonal matrix with non-negative diagonal elements,
SV(1), SV(2), ..., SV(N), ordered such that

   `SV(1) >= SV(2) >= ... >= SV(N) >= 0.`

The columns of Q are the left singular vectors of A, the diagonal
elements of S are the singular values of A and the columns of P
are the right singular vectors of A.

Either or both of Q and P' may be requested.
When P' is computed, it is returned in A.



See the SLICOT documentation for details.
"""
function mb03ud! end

"""

Reduce a product of p real general matrices A = `A_1`*`A_2`*...*`A_p`
to upper Hessenberg form, H = `H_1`*`H_2`*...*`H_p`, where `H_1` is
upper Hessenberg, and `H_2`, ..., `H_p` are upper triangular, by using
orthogonal similarity transformations on A,

        `Q_1' * A_1 * Q_2 = H_1,`
        `Q_2' * A_2 * Q_3 = H_2,`
               `...`
        `Q_p' * A_p * Q_1 = H_p.`



See the SLICOT documentation for details.
"""
function mb03vd! end

"""

Generate the real orthogonal matrices `Q_1`, `Q_2`, ..., `Q_p`,
which are defined as the product of ihi-ilo elementary reflectors
of order n, as returned by SLICOT Library routine MB03VD:

   `Q_j = H_j(ilo) H_j(ilo+1) . . . H_j(ihi-1).`



See the SLICOT documentation for details.
"""
function mb03vy! end

"""

Swap adjacent diagonal blocks A11*B11 and A22*B22 of size
1-by-1 or 2-by-2 in an upper (quasi) triangular matrix product
A*B by an orthogonal equivalence transformation.

(A, B) must be in periodic real Schur canonical form (as returned
by SLICOT Library routine MB03XP), i.e., A is block upper
triangular with 1-by-1 and 2-by-2 diagonal blocks, and B is upper
triangular.

Optionally, the matrices Q and Z of generalized Schur vectors are
updated.

    `Q(in) * A(in) * Z(in)' = Q(out) * A(out) * Z(out)',`
    `Z(in) * B(in) * Q(in)' = Z(out) * B(out) * Q(out)'.`

This routine is largely based on the LAPACK routine DTGEX2
developed by Bo Kagstrom and Peter Poromaa.



See the SLICOT documentation for details.
"""
function mb03wa! end

"""

Compute the Schur decomposition and the eigenvalues of a
product of matrices, H = `H_1`*`H_2`*...*`H_p`, with `H_1` an upper
Hessenberg matrix and `H_2`, ..., `H_p` upper triangular matrices,
without evaluating the product. Specifically, the matrices Z_i
are computed, such that

        `Z_1' * H_1 * Z_2 = T_1,`
        `Z_2' * H_2 * Z_3 = T_2,`
               `...`
        `Z_p' * H_p * Z_1 = T_p,`

where `T_1` is in real Schur form, and `T_2`, ..., `T_p` are upper
triangular.

The routine works primarily with the Hessenberg and triangular
submatrices in rows and columns ILO to IHI, but optionally applies
the transformations to all the rows and columns of the matrices
H_i, i = 1,...,p. The transformations can be optionally
accumulated.



See the SLICOT documentation for details.
"""
function mb03wd! end

"""

Compute the eigenvalues of a product of matrices,
T = `T_1`*`T_2`*...*`T_p`, where `T_1` is an upper quasi-triangular
matrix and `T_2`, ..., `T_p` are upper triangular matrices.



See the SLICOT documentation for details.
"""
function mb03wx! end

"""

Compute the eigenvalues of a Hamiltonian matrix,

              `[  A   G  ]         T        T`
        `H  =  [       T ],   G = G,   Q = Q,                  (1)`
              `[  Q  -A  ]`

where A, G and Q are real n-by-n matrices.

Due to the structure of H all eigenvalues appear in pairs
(lambda,-lambda). This routine computes the eigenvalues of H
using an algorithm based on the symplectic URV and the periodic
Schur decompositions as described in [1],

      `T       [  T   G  ]`
     `U H V =  [       T ],                                    (2)`
              `[  0   S  ]`

where U and V are 2n-by-2n orthogonal symplectic matrices,
S is in real Schur form and T is upper triangular.

The algorithm is backward stable and preserves the eigenvalue
pairings in finite precision arithmetic.

Optionally, a symplectic balancing transformation to improve the
conditioning of eigenvalues is computed (see MB04DD). In this
case, the matrix H in decomposition (2) must be replaced by the
balanced matrix.

The SLICOT Library routine MB03ZD can be used to compute invariant
subspaces of H from the output of this routine.



See the SLICOT documentation for details.
"""
function mb03xd! end

"""

Compute the periodic Schur decomposition and the eigenvalues of
a product of matrices, H = A*B, with A upper Hessenberg and B
upper triangular without evaluating any part of the product.
Specifically, the matrices Q and Z are computed, so that

     `Q' * A * Z = S,    Z' * B * Q = T`

where S is in real Schur form, and T is upper triangular.



See the SLICOT documentation for details.
"""
function mb03xp! end

"""

Compute the eigenvalues and real skew-Hamiltonian Schur form of
a skew-Hamiltonian matrix,

              `[  A   G  ]`
        `W  =  [       T ],`
              `[  Q   A  ]`

where A is an N-by-N matrix and G, Q are N-by-N skew-symmetric
matrices. Specifically, an orthogonal symplectic matrix U is
computed so that

          `T       [  Aout  Gout  ]`
         `U W U =  [            T ] ,`
                  `[    0   Aout  ]`

where Aout is in Schur canonical form (as returned by the LAPACK
routine DHSEQR). That is, Aout is block upper triangular with
1-by-1 and 2-by-2 diagonal blocks; each 2-by-2 diagonal block has
its diagonal elements equal and its off-diagonal elements of
opposite sign.

Optionally, the matrix U is returned in terms of its first N/2
rows

                 `[  U1   U2 ]`
             `U = [          ].`
                 `[ -U2   U1 ]`



See the SLICOT documentation for details.
"""
function mb03xs! end

"""

Reduce 2*nb columns and rows of a real (k+2n)-by-(k+2n)
matrix H:

        `[ op(A)   G   ]`
    `H = [             ],`
        `[  Q    op(B) ]`

so that elements in the first nb columns below the k-th
subdiagonal of the (k+n)-by-n matrix op(A), in the first nb
columns and rows of the n-by-n matrix Q and in the first nb rows
above the diagonal of the n-by-(k+n) matrix op(B) are zero.
The reduction is performed by orthogonal symplectic
transformations UU'*H*VV and matrices U, V, YA, YB, YG, YQ, XA,
XB, XG, and XQ are returned so that

               `[ op(Aout)+U*YA'+XA*V'     G+U*YG'+XG*V'    ]`
    `UU' H VV = [                                           ].`
               `[   Qout+U*YQ'+XQ*V'   op(Bout)+U*YB'+XB*V' ]`

This is an auxiliary routine called by MB04TB.



See the SLICOT documentation for details.
"""
function mb03xu! end

"""

Compute the eigenvalues of a Hamiltonian matrix,

              `[  A   G  ]         H        H`
        `H  =  [       H ],   G = G ,  Q = Q ,                  (1)`
              `[  Q  -A  ]`

where A, G and Q are complex n-by-n matrices.

Due to the structure of H, if lambda is an eigenvalue, then
-conjugate(lambda) is also an eigenvalue. This does not mean that
purely imaginary eigenvalues are necessarily multiple. The routine
computes the eigenvalues of H using an embedding to a real skew-
Hamiltonian matrix He,

               `[  Ae   Ge  ]            T            T`
        `He  =  [         T ],   Ge = -Ge ,   Qe = -Qe ,        (2)`
               `[  Qe   Ae  ]`

where Ae, Ge, and Qe are real 2*n-by-2*n matrices, defined by

               `[   Im(A)   Re(A)  ]`
        `Ae  =  [                  ],`
               `[  -Re(A)   Im(A)  ]`

               `[  triu(Im(G))     Re(G)     ]`
   `triu(Ge) =  [                            ],`
               `[       0       triu(Im(G))  ]`

               `[  tril(Im(Q))       0       ]`
   `tril(Qe) =  [                            ], `
               `[     -Re(Q)    tril(Im(Q))  ]`

and triu and tril denote the upper and lower triangle,
respectively. Then, an orthogonal symplectic matrix Ue is used to
reduce He to the structured real Schur form

      `T          [  Se   De ]            T`
     `Ue He Ue =  [        T ],   De = -De ,                    (3)`
                 `[  0    Se ]`

where Ue is a 4n-by-4n real symplectic matrix, and Se is upper
quasi-triangular (real Schur form).

Optionally, if JOB = 'S', or JOB = 'G', the matrix i*He is further
transformed to the structured complex Schur form

      `H            [  Sc  Gc ]           H`
     `U (i*He) U =  [       H ],   Gc = Gc ,                    (4)`
                   `[  0  -Sc ]`

where U is a 4n-by-4n unitary symplectic matrix, and Sc is upper
triangular (Schur form).

The algorithm is backward stable and preserves the spectrum
structure in finite precision arithmetic.

Optionally, a symplectic balancing transformation to improve the
conditioning of eigenvalues is computed (see the SLICOT Library
routine MB04DZ). In this case, the matrix He in decompositions (3)
and (4) must be replaced by the balanced matrix.



See the SLICOT documentation for details.
"""
function mb03xz! end

"""

Annihilate one or two entries on the subdiagonal of the
Hessenberg matrix A for dealing with zero elements on the diagonal
of the triangular matrix B.

MB03YA is an auxiliary routine called by SLICOT Library routines
MB03XP and MB03YD.



See the SLICOT documentation for details.
"""
function mb03ya! end

"""

Deal with small subtasks of the product eigenvalue problem.

MB03YD is an auxiliary routine called by SLICOT Library routine
MB03XP.



See the SLICOT documentation for details.
"""
function mb03yd! end

"""

Compute the periodic Schur factorization of a real 2-by-2
matrix pair (A,B) where B is upper triangular. This routine
computes orthogonal (rotation) matrices given by CSL, SNL and CSR,
SNR such that

1) if the pair (A,B) has two real eigenvalues, then

   `[ a11 a12 ] := [  CSL  SNL ] [ a11 a12 ] [  CSR -SNR ]`
   `[  0  a22 ]    [ -SNL  CSL ] [ a21 a22 ] [  SNR  CSR ]`

   `[ b11 b12 ] := [  CSR  SNR ] [ b11 b12 ] [  CSL -SNL ]`
   `[  0  b22 ]    [ -SNR  CSR ] [  0  b22 ] [  SNL  CSL ],`

2) if the pair (A,B) has a pair of complex conjugate eigenvalues,
   `then`

   `[ a11 a12 ] := [  CSL  SNL ] [ a11 a12 ] [  CSR -SNR ]`
   `[ a21 a22 ]    [ -SNL  CSL ] [ a21 a22 ] [  SNR  CSR ]`

   `[ b11  0  ] := [  CSR  SNR ] [ b11 b12 ] [  CSL -SNL ]`
   `[  0  b22 ]    [ -SNR  CSR ] [  0  b22 ] [  SNL  CSL ].`

This is a modified version of the LAPACK routine DLAGV2 for
computing the real, generalized Schur decomposition of a
two-by-two matrix pencil.



See the SLICOT documentation for details.
"""
function mb03yt! end

"""

1. To compute, for a given matrix pair (A,B) in periodic Schur
   `form, orthogonal matrices Ur and Vr so that`

       `T           [ A11  A12 ]     T           [ B11  B12 ]`
     `Vr * A * Ur = [          ],  Ur * B * Vr = [          ], (1)`
                   `[  0   A22 ]                 [  0   B22 ]`

   `is in periodic Schur form, and the eigenvalues of A11*B11`
   `form a selected cluster of eigenvalues.`

2. To compute an orthogonal matrix W so that

              `T  [  0  -A11 ]       [  R11   R12 ]`
             `W * [          ] * W = [            ],           (2)`
                 `[ B11   0  ]       [   0    R22 ]`

   `where the eigenvalues of R11 and -R22 coincide and have`
   `positive real part.`

Optionally, the matrix C is overwritten by Ur'*C*Vr.

All eigenvalues of A11*B11 must either be complex or real and
negative.



See the SLICOT documentation for details.
"""
function mb03za! end

"""

Compute the stable and unstable invariant subspaces for a
Hamiltonian matrix with no eigenvalues on the imaginary axis,
using the output of the SLICOT Library routine MB03XD.



See the SLICOT documentation for details.
"""
function mb03zd! end

"""

Compute the eigenvalues of a real N-by-N skew-Hamiltonian/
Hamiltonian pencil aS - bH with

                                 `(  0  I  )`
  `S = T Z = J Z' J' Z, where J = (        ),                   (1)`
                                 `( -I  0  )`

via generalized symplectic URV decomposition. That is, orthogonal
matrices Q1 and Q2 and orthogonal symplectic matrices U1 and U2
are computed such that

                              `( T11  T12 )`
  `Q1' T U1 = Q1' J Z' J' U1 = (          ) = Tout,`
                              `(  0   T22 )`

             `( Z11  Z12 )`
  `U2' Z Q2 = (          ) = Zout,                              (2)`
             `(  0   Z22 )`

             `( H11  H12 )`
  `Q1' H Q2 = (          ) = Hout,`
             `(  0   H22 )`

where T11, T22', Z11, Z22', H11 are upper triangular and H22' is
upper quasi-triangular. The notation M' denotes the transpose of
the matrix M.
Optionally, if COMPQ1 = 'I' or COMPQ1 = 'U', the orthogonal
transformation matrix Q1 will be computed.
Optionally, if COMPQ2 = 'I' or COMPQ2 = 'U', the orthogonal
transformation matrix Q2 will be computed.
Optionally, if COMPU1 = 'I' or COMPU1 = 'U', the orthogonal
symplectic transformation matrix

       `(  U11  U12  )`
  `U1 = (            )`
       `( -U12  U11  )`

will be computed.
Optionally, if COMPU2 = 'I' or COMPU2 = 'U', the orthogonal
symplectic transformation matrix

       `(  U21  U22  )`
  `U2 = (            )`
       `( -U22  U21  )`

will be computed.



See the SLICOT documentation for details.
"""
function mb04ad! end

"""

Compute the eigenvalues of a complex N-by-N skew-Hamiltonian/
Hamiltonian pencil aS - bH, with

         `H  T           (  B  F  )       (  Z11  Z12  )`
  `S = J Z  J  Z and H = (      H ), Z =: (            ).       (1)`
                        `(  G -B  )       (  Z21  Z22  )`

The structured Schur form of the embedded real skew-Hamiltonian/
                                                      `H  T`
skew-Hamiltonian pencil, a`B_S` - b`B_T`, with `B_S` = J `B_Z`  J  `B_Z`,

        `(  Re(Z11)  -Im(Z11)  |  Re(Z12)  -Im(Z12)  )`
        `(                     |                     )`
        `(  Im(Z11)   Re(Z11)  |  Im(Z12)   Re(Z12)  )`
        `(                     |                     )`
  `B_Z = (---------------------+---------------------) ,`
        `(                     |                     )`
        `(  Re(Z21)  -Im(Z21)  |  Re(Z22)  -Im(Z22)  )`
        `(                     |                     )`
        `(  Im(Z21)   Re(Z21)  |  Im(Z22)   Re(Z22)  )`
                                                               `(2)`
        `( -Im(B)  -Re(B)  | -Im(F)  -Re(F)  )`
        `(                 |                 )`
        `(  Re(B)  -Im(B)  |  Re(F)  -Im(F)  )`
        `(                 |                 )`
  `B_T = (-----------------+-----------------) ,  T = i*H,`
        `(                 |      T       T  )`
        `( -Im(G)  -Re(G)  | -Im(B )  Re(B ) )`
        `(                 |      T       T  )`
        `(  Re(G)  -Im(G)  | -Re(B ) -Im(B ) )`

is determined and used to compute the eigenvalues. Optionally,
if JOB = 'T', the pencil a`B_S` - b`B_H` is transformed by a unitary
matrix Q and a unitary symplectic matrix U to the structured Schur
                                              `H  T`
form a`B_S`out - b`B_H`out, with `B_S`out = J `B_Z`out  J  `B_Z`out,

           `( BA  BD  )              ( BB  BF  )`
  `B_Zout = (         ) and B_Hout = (       H ),               (3)`
           `(  0  BC  )              (  0 -BB  )`

where BA and BB are upper triangular, BC is lower triangular,
and BF is Hermitian. `B_H` above is defined as `B_H` = -i*`B_T`.
The embedding doubles the multiplicities of the eigenvalues of
the pencil aS - bH.
Optionally, if COMPQ = 'C', the unitary matrix Q is computed.
Optionally, if COMPU = 'C', the unitary symplectic matrix U is
computed.



See the SLICOT documentation for details.
"""
function mb04az! end

"""

Compute the eigenvalues of a real N-by-N skew-Hamiltonian/
Hamiltonian pencil aS - bH with

      `(  A  D  )         (  C  V  )`
  `S = (        ) and H = (        ).                           (1)`
      `(  E  A' )         (  W -C' )`

Optionally, if JOB = 'T', decompositions of S and H will be
computed via orthogonal transformations Q1 and Q2 as follows:

                  `(  Aout  Dout  )`
  `Q1' S J Q1 J' = (              ),`
                  `(   0    Aout' )`

                  `(  Bout  Fout  )`
  `J' Q2' J S Q2 = (              ) =: T,                       (2)`
                  `(   0    Bout' )`

             `(  C1out  Vout  )            (  0  I  )`
  `Q1' H Q2 = (               ), where J = (        )`
             `(  0     C2out' )            ( -I  0  )`

and Aout, Bout, C1out are upper triangular, C2out is upper quasi-
triangular and Dout and Fout are skew-symmetric. The notation M'
denotes the transpose of the matrix M.
Optionally, if COMPQ1 = 'I' or COMPQ1 = 'U', then the orthogonal
transformation matrix Q1 will be computed.
Optionally, if COMPQ2 = 'I' or COMPQ2 = 'U', then the orthogonal
transformation matrix Q2 will be computed.



See the SLICOT documentation for details.
"""
function mb04bd! end

"""

Compute the eigenvalues of a real N-by-N skew-Hamiltonian/
Hamiltonian pencil aS - bH with

      `(  A  D  )         (  C  V  )`
  `S = (        ) and H = (        ).                           (1)`
      `(  E  A' )         (  W -C' )`

Optionally, if JOB = 'T', decompositions of S and H will be
computed via orthogonal transformations Q1 and Q2 as follows:

                  `(  Aout  Dout  )`
  `Q1' S J Q1 J' = (              ),`
                  `(   0    Aout' )`

                  `(  Bout  Fout  )`
  `J' Q2' J S Q2 = (              ) =: T,                       (2)`
                  `(   0    Bout' )`

             `(  C1out  Vout  )            (  0  I  )`
  `Q1' H Q2 = (               ), where J = (        )`
             `(  0     C2out' )            ( -I  0  )`

and Aout, Bout, C1out are upper triangular, C2out is upper quasi-
triangular and Dout and Fout are skew-symmetric. The notation M'
denotes the transpose of the matrix M.
Optionally, if COMPQ1 = 'I' or COMPQ1 = 'U', then the orthogonal
transformation matrix Q1 will be computed.
Optionally, if COMPQ2 = 'I' or COMPQ2 = 'U', then the orthogonal
transformation matrix Q2 will be computed.



See the SLICOT documentation for details.
"""
function mb04bp! end

"""

Compute the eigenvalues of a complex N-by-N skew-Hamiltonian/
Hamiltonian pencil aS - bH, with

      `(  A  D  )         (  B  F  )`
  `S = (      H ) and H = (      H ).                           (1)`
      `(  E  A  )         (  G -B  )`

This routine computes the eigenvalues using an embedding to a real
skew-Hamiltonian/skew-Hamiltonian pencil a`B_S` - b`B_T`, defined as

        `(  Re(A)  -Im(A)  |  Re(D)  -Im(D)  )`
        `(                 |                 )`
        `(  Im(A)   Re(A)  |  Im(D)   Re(D)  )`
        `(                 |                 )`
  `B_S = (-----------------+-----------------) , and`
        `(                 |      T       T  )`
        `(  Re(E)  -Im(E)  |  Re(A )  Im(A ) )`
        `(                 |      T       T  )`
        `(  Im(E)   Re(E)  | -Im(A )  Re(A ) )`
                                                               `(2)`
        `( -Im(B)  -Re(B)  | -Im(F)  -Re(F)  )`
        `(                 |                 )`
        `(  Re(B)  -Im(B)  |  Re(F)  -Im(F)  )`
        `(                 |                 )`
  `B_T = (-----------------+-----------------) ,  T = i*H.`
        `(                 |      T       T  )`
        `( -Im(G)  -Re(G)  | -Im(B )  Re(B ) )`
        `(                 |      T       T  )`
        `(  Re(G)  -Im(G)  | -Re(B ) -Im(B ) )`

Optionally, if JOB = 'T', the pencil a`B_S` - b`B_H` (`B_H` = -i*`B_T`) is
transformed by a unitary matrix Q to the structured Schur form

           `( BA  BD  )              ( BB  BF  )`
  `B_Sout = (       H ) and B_Hout = (       H ),               (3)`
           `(  0  BA  )              (  0 -BB  )`

where BA and BB are upper triangular, BD is skew-Hermitian, and
BF is Hermitian. The embedding doubles the multiplicities of the
eigenvalues of the pencil aS - bH. Optionally, if COMPQ = 'C', the
unitary matrix Q is computed.



See the SLICOT documentation for details.
"""
function mb04bz! end

"""

Compute the transformed matrices A, B and D, using orthogonal
matrices Q1, Q2 and Q3 for a real N-by-N regular pencil

                `( A11   0  ) ( B11   0  )     (  0   D12 )`
  `aA*B - bD = a (          ) (          ) - b (          ),    (1)`
                `(  0   A22 ) (  0   B22 )     ( D21   0  )`

where A11, A22, B11, B22 and D12 are upper triangular, D21 is
upper quasi-triangular and the generalized matrix product 
   `-1        -1    -1        -1`
A11   D12 B22   A22   D21 B11   is upper quasi-triangular, such
that Q3' A Q2, Q2' B Q1 are upper triangular, Q3' D Q1 is upper
quasi-triangular and the transformed pencil
a(Q3' A B Q1) - b(Q3' D Q1) is in generalized Schur form. The
notation M' denotes the transpose of the matrix M.



See the SLICOT documentation for details.
"""
function mb04cd! end

"""

Apply from the left the inverse of a balancing transformation,
computed by the SLICOT Library routine MB04DP, to the matrix

          `[   V1   ]`
          `[        ],`
          `[ sgn*V2 ]`

where sgn is either +1 or -1.



See the SLICOT documentation for details.
"""
function mb04db! end

"""

Balance a real Hamiltonian matrix,

              `[  A   G  ]`
         `H =  [       T ] ,`
              `[  Q  -A  ]`

where A is an N-by-N matrix and G, Q are N-by-N symmetric
matrices. This involves, first, permuting H by a symplectic
similarity transformation to isolate eigenvalues in the first
1:ILO-1 elements on the diagonal of A; and second, applying a
diagonal similarity transformation to rows and columns
ILO:N, N+ILO:2*N to make the rows and columns as close in 1-norm
as possible. Both steps are optional.



See the SLICOT documentation for details.
"""
function mb04dd! end

"""

Apply the inverse of a balancing transformation, computed by
the SLICOT Library routines MB04DD or MB04DS, to a 2*N-by-M matrix

          `[   V1   ]`
          `[        ],`
          `[ sgn*V2 ]`

where sgn is either +1 or -1.



See the SLICOT documentation for details.
"""
function mb04di! end

"""

Balance a pair of N-by-N real matrices (A,B). This involves,
first, permuting A and B by equivalence transformations to isolate
eigenvalues in the first 1 to ILO-1 and last IHI+1 to N elements
on the diagonal of A and B; and second, applying a diagonal
equivalence transformation to rows and columns ILO to IHI to make
the rows and columns as close in 1-norm as possible. Both steps
are optional. Balancing may reduce the 1-norms of the matrices,
and improve the accuracy of the computed eigenvalues and/or
eigenvectors in the generalized eigenvalue problem
A*x = lambda*B*x.

This routine may optionally improve the conditioning of the
scaling transformation compared to the LAPACK routine DGGBAL.



See the SLICOT documentation for details.
"""
function mb04dl! end

"""

Balance the 2*N-by-2*N skew-Hamiltonian/Hamiltonian pencil
aS - bH, with

      `(  A  D  )         (  C  V  )`
  `S = (        ) and H = (        ),  A, C N-by-N,             (1)`
      `(  E  A' )         (  W -C' )`

where D and E are skew-symmetric, and V and W are symmetric
matrices. This involves, first, permuting aS - bH by a symplectic
equivalence transformation to isolate eigenvalues in the first
1:ILO-1 elements on the diagonal of A and C; and second, applying
a diagonal equivalence transformation to make the pairs of rows
and columns ILO:N and N+ILO:2*N as close in 1-norm as possible.
Both steps are optional. Balancing may reduce the 1-norms of the
matrices S and H.



See the SLICOT documentation for details.
"""
function mb04dp! end

"""

Balance a real skew-Hamiltonian matrix

              `[  A   G  ]`
         `S =  [       T ] ,`
              `[  Q   A  ]`

where A is an N-by-N matrix and G, Q are N-by-N skew-symmetric
matrices. This involves, first, permuting S by a symplectic
similarity transformation to isolate eigenvalues in the first
1:ILO-1 elements on the diagonal of A; and second, applying a
diagonal similarity transformation to rows and columns
ILO:N, N+ILO:2*N to make the rows and columns as close in 1-norm
as possible. Both steps are optional.



See the SLICOT documentation for details.
"""
function mb04ds! end

"""

Perform a symplectic scaling on the Hamiltonian matrix

         `( A    G  )`
     `H = (       T ),                                          (1)`
         `( Q   -A  )`

i.e., perform either the symplectic scaling transformation

                              `-1`
            `( A'   G'  )   ( D   0 ) ( A   G  ) ( D  0   )`
     `H' <-- (        T ) = (       ) (      T ) (     -1 ),    (2)`
            `( Q'  -A'  )   ( 0   D ) ( Q  -A  ) ( 0  D   )`

where D is a diagonal scaling matrix, or the symplectic norm
scaling transformation

             `( A''   G''  )    1  (   A   G/tau )`
     `H'' <-- (          T ) = --- (           T ),             (3)`
             `( Q''  -A''  )   tau ( tau Q   -A  )`

where tau is a real scalar.  Note that if tau is not equal to 1,
then (3) is NOT a similarity transformation.  The eigenvalues
of H are then tau times the eigenvalues of H''.

For symplectic scaling (2), D is chosen to give the rows and
columns of A' approximately equal 1-norms and to give Q' and G'
approximately equal norms.  (See METHOD below for details.) For
norm scaling, tau = MAX(1, ||A||, ||G||, ||Q||) where ||.||
denotes the 1-norm (column sum norm).



See the SLICOT documentation for details.
"""
function mb04dy! end

"""

Balance a complex Hamiltonian matrix,

              `[  A   G  ]`
         `H =  [       H ] ,`
              `[  Q  -A  ]`

where A is an N-by-N matrix and G, Q are N-by-N Hermitian
matrices. This involves, first, permuting H by a symplectic
similarity transformation to isolate eigenvalues in the first
1:ILO-1 elements on the diagonal of A; and second, applying a
diagonal similarity transformation to rows and columns
ILO:N, N+ILO:2*N to make the rows and columns as close in 1-norm
as possible. Both steps are optional. Assuming ILO = 1, let D be a
diagonal matrix of order N with the scaling factors on the
diagonal. The scaled Hamiltonian is defined by

               `[  D**-1*A*D   D**-1*G*D**-1  ]`
         `Hs =  [                   H         ] .`
               `[    D*Q*D      -D*A *D**-1   ]`




See the SLICOT documentation for details.
"""
function mb04dz! end

"""

Compute the eigenvalues of a real N-by-N skew-Hamiltonian/
skew-Hamiltonian pencil aS - bT with

                        `(  B  F  )            (  0  I  )`
  `S = J Z' J' Z and T = (        ), where J = (        ).      (1)`
                        `(  G  B' )            ( -I  0  )`

Optionally, if JOB = 'T', the pencil aS - bT will be transformed
to the structured Schur form: an orthogonal transformation matrix
Q and an orthogonal symplectic transformation matrix U are
computed, such that

           `(  Z11  Z12  )`
  `U' Z Q = (            ) = Zout, and`
           `(   0   Z22  )`
                                                               `(2)`
                `(  Bout  Fout  )`
  `J Q' J' T Q = (              ),`
                `(   0    Bout' )`

where Z11 and Z22' are upper triangular and Bout is upper quasi-
triangular. The notation M' denotes the transpose of the matrix M.
Optionally, if COMPQ = 'I', the orthogonal transformation matrix Q
will be computed.
Optionally, if COMPU = 'I' or COMPU = 'U', the orthogonal
symplectic transformation matrix

      `(  U1  U2  )`
  `U = (          )`
      `( -U2  U1  )`

will be computed.



See the SLICOT documentation for details.
"""
function mb04ed! end

"""

Compute the eigenvalues of a real N-by-N skew-Hamiltonian/
skew-Hamiltonian pencil aS - bT with

      `(  A  D  )         (  B  F  )`
  `S = (        ) and T = (        ).                           (1)`
      `(  E  A' )         (  G  B' )`

Optionally, if JOB = 'T', the pencil aS - bT will be transformed
to the structured Schur form: an orthogonal transformation matrix
Q is computed such that

                `(  Aout  Dout  )`
  `J Q' J' S Q = (              ), and`
                `(   0    Aout' )`
                                                               `(2)`
                `(  Bout  Fout  )            (  0  I  )`
  `J Q' J' T Q = (              ), where J = (        ),`
                `(   0    Bout' )            ( -I  0  )`

Aout is upper triangular, and Bout is upper quasi-triangular. The
notation M' denotes the transpose of the matrix M.
Optionally, if COMPQ = 'I' or COMPQ = 'U', the orthogonal
transformation matrix Q will be computed.



See the SLICOT documentation for details.
"""
function mb04fd! end

"""

Compute the eigenvalues of a real N-by-N skew-Hamiltonian/
skew-Hamiltonian pencil aS - bT with

      `(  A  D  )         (  B  F  )`
  `S = (        ) and T = (        ).                           (1)`
      `(  E  A' )         (  G  B' )`

Optionally, if JOB = 'T', the pencil aS - bT will be transformed
to the structured Schur form: an orthogonal transformation matrix
Q is computed such that

                `(  Aout  Dout  )`
  `J Q' J' S Q = (              ), and`
                `(   0    Aout' )`
                                                               `(2)`
                `(  Bout  Fout  )            (  0  I  )`
  `J Q' J' T Q = (              ), where J = (        ),`
                `(   0    Bout' )            ( -I  0  )`

Aout is upper triangular, and Bout is upper quasi-triangular. The
notation M' denotes the transpose of the matrix M.
Optionally, if COMPQ = 'I' or COMPQ = 'U', the orthogonal
transformation matrix Q will be computed.



See the SLICOT documentation for details.
"""
function mb04fp! end

"""

Compute an RQ factorization with row pivoting of a
real m-by-n matrix A: P*A = R*Q.



See the SLICOT documentation for details.
"""
function mb04gd! end

"""

Compute the transformed matrices A and B, using orthogonal
matrices Q1 and Q2 for a real N-by-N regular pencil

              `( A11   0  )     (  0   B12 )`
  `aA - bB = a (          ) - b (          ),                   (1)`
              `(  0   A22 )     ( B21   0  )`

where A11, A22 and B12 are upper triangular, B21 is upper
quasi-triangular and the generalized matrix product
   `-1        -1`
A11   B12 A22   B21 is in periodic Schur form, such that the
matrix Q2' A Q1 is upper triangular, Q2' B Q1 is upper
quasi-triangular and the transformed pencil
a(Q2' A Q1) - b(Q2' B Q1) is in generalized Schur form. The
notation M' denotes the transpose of the matrix M.



See the SLICOT documentation for details.
"""
function mb04hd! end

"""

Compute a QR factorization of an n-by-m matrix A (A = Q * R),
having a p-by-min(p,m) zero triangle in the lower left-hand side
corner, as shown below, for n = 8, m = 7, and p = 2:

       `[ x x x x x x x ]`
       `[ x x x x x x x ]`
       `[ x x x x x x x ]`
       `[ x x x x x x x ]`
   `A = [ x x x x x x x ],`
       `[ x x x x x x x ]`
       `[ 0 x x x x x x ]`
       `[ 0 0 x x x x x ]`

and optionally apply the transformations to an n-by-l matrix B
(from the left). The problem structure is exploited. This
computation is useful, for instance, in combined measurement and
time update of one iteration of the time-invariant Kalman filter
(square root information filter).



See the SLICOT documentation for details.
"""
function mb04id! end

"""

Overwrite the real n-by-m matrix  C  with  Q' * C,  Q * C,
C * Q',  or  C * Q,  according to the following table

                `SIDE = 'L'     SIDE = 'R'`
TRANS = 'N':      Q * C          C * Q
TRANS = 'T':      Q'* C          C * Q'

where  Q  is a real orthogonal matrix defined as the product of
k elementary reflectors

   `Q = H(1) H(2) . . . H(k)`

as returned by SLICOT Library routine MB04ID.  Q  is of order n
if  SIDE = 'L'  and of order m if  SIDE = 'R'.



See the SLICOT documentation for details.
"""
function mb04iy! end

"""

Compute a QR factorization of an n-by-m matrix A (A = Q * R),
having a p-by-min(p,m) zero triangle in the lower left-hand side
corner, as shown below, for n = 8, m = 7, and p = 2:

       `[ x x x x x x x ]`
       `[ x x x x x x x ]`
       `[ x x x x x x x ]`
       `[ x x x x x x x ]`
   `A = [ x x x x x x x ],`
       `[ x x x x x x x ]`
       `[ 0 x x x x x x ]`
       `[ 0 0 x x x x x ]`

and optionally apply the transformations to an n-by-l matrix B
(from the left). The problem structure is exploited. This
computation is useful, for instance, in combined measurement and
time update of one iteration of the time-invariant Kalman filter
(square root information filter).



See the SLICOT documentation for details.
"""
function mb04iz! end

"""

Compute an LQ factorization of an n-by-m matrix A (A = L * Q),
having a min(n,p)-by-p zero triangle in the upper right-hand side
corner, as shown below, for n = 8, m = 7, and p = 2:

       `[ x x x x x 0 0 ]`
       `[ x x x x x x 0 ]`
       `[ x x x x x x x ]`
       `[ x x x x x x x ]`
   `A = [ x x x x x x x ],`
       `[ x x x x x x x ]`
       `[ x x x x x x x ]`
       `[ x x x x x x x ]`

and optionally apply the transformations to an l-by-m matrix B
(from the right). The problem structure is exploited. This
computation is useful, for instance, in combined measurement and
time update of one iteration of the time-invariant Kalman filter
(square root covariance filter).



See the SLICOT documentation for details.
"""
function mb04jd! end

"""

Calculate a QR factorization of the first block column and
apply the orthogonal transformations (from the left) also to the
second block column of a structured matrix, as follows
                     `_`
       `[ R   0 ]   [ R   C ]`
  `Q' * [       ] = [       ]`
       `[ A   B ]   [ 0   D ]`
            `_`
where R and R are upper triangular. The matrix A can be full or
upper trapezoidal/triangular. The problem structure is exploited.
This computation is useful, for instance, in combined measurement
and time update of one iteration of the Kalman filter (square
root information filter).



See the SLICOT documentation for details.
"""
function mb04kd! end

"""

Calculate an LQ factorization of the first block row and apply
the orthogonal transformations (from the right) also to the second
block row of a structured matrix, as follows
                   `_`
   `[ L   A ]     [ L   0 ]`
   `[       ]*Q = [       ]`
   `[ 0   B ]     [ C   D ]`
            `_`
where L and L are lower triangular. The matrix A can be full or
lower trapezoidal/triangular. The problem structure is exploited.
This computation is useful, for instance, in combined measurement
and time update of one iteration of the Kalman filter (square
root covariance filter).



See the SLICOT documentation for details.
"""
function mb04ld! end

"""

Reduce the 1-norm of a general real matrix A by balancing.
This involves diagonal similarity transformations applied
iteratively to A to make the rows and columns as close in norm as
possible.

This routine can be used instead LAPACK Library routine DGEBAL,
when no reduction of the 1-norm of the matrix is possible with
DGEBAL, as for upper triangular matrices. LAPACK Library routine
DGEBAK, with parameters ILO = 1, IHI = N, and JOB = 'S', should
be used to apply the backward transformation.



See the SLICOT documentation for details.
"""
function mb04md! end

"""

Calculate an RQ factorization of the first block row and
apply the orthogonal transformations (from the right) also to the
second block row of a structured matrix, as follows
                         `_`
  `[ A   R ]        [ 0   R ]`
  `[       ] * Q' = [ _   _ ]`
  `[ C   B ]        [ C   B ]`
            `_`
where R and R are upper triangular. The matrix A can be full or
upper trapezoidal/triangular. The problem structure is exploited.



See the SLICOT documentation for details.
"""
function mb04nd! end

"""

Apply a real elementary reflector H to a real m-by-(n+1)
matrix C = [ A  B ], from the right, where A has one column. H is
represented in the form
                                   `( 1 )`
      `H = I - tau * u *u',    u  = (   ),`
                                   `( v )`
where tau is a real scalar and v is a real n-vector.

If tau = 0, then H is taken to be the unit matrix.

In-line code is used if H has order < 11.



See the SLICOT documentation for details.
"""
function mb04ny! end

"""

Calculate a QR factorization of the first block column and
apply the orthogonal transformations (from the left) also to the
second block column of a structured matrix, as follows
                     `_   _`
       `[ R   B ]   [ R   B ]`
  `Q' * [       ] = [     _ ]`
       `[ A   C ]   [ 0   C ]`
            `_`
where R and R are upper triangular. The matrix A can be full or
upper trapezoidal/triangular. The problem structure is exploited.



See the SLICOT documentation for details.
"""
function mb04od! end

"""

Perform the QR factorization

   `( U  ) = Q*( R ),  where  U = ( U1  U2 ),  R = ( R1  R2 ),`
   `( x' )     ( 0 )              ( 0   T  )       ( 0   R3 )`

where U and R are (m+n)-by-(m+n) upper triangular matrices, x is
an m+n element vector, U1 is m-by-m, T is n-by-n, stored
separately, and Q is an (m+n+1)-by-(m+n+1) orthogonal matrix.

The matrix ( U1 U2 ) must be supplied in the m-by-(m+n) upper
trapezoidal part of the array A and this is overwritten by the
corresponding part ( R1 R2 ) of R. The remaining upper triangular
part of R, R3, is overwritten on the array T.

The transformations performed are also applied to the (m+n+1)-by-p
matrix ( B' C' d )' (' denotes transposition), where B, C, and d'
are m-by-p, n-by-p, and 1-by-p matrices, respectively.



See the SLICOT documentation for details.
"""
function mb04ow! end

"""

Perform the QR factorization

   `(U ) = Q*(R),`
   `(x')     (0)`

where U and R are n-by-n upper triangular matrices, x is an
n element vector and Q is an (n+1)-by-(n+1) orthogonal matrix.

U must be supplied in the n-by-n upper triangular part of the
array A and this is overwritten by R.



See the SLICOT documentation for details.
"""
function mb04ox! end

"""

Apply a real elementary reflector H to a real (m+1)-by-n
matrix C = [ A ], from the left, where A has one row. H is
           `[ B ]`
represented in the form
                                   `( 1 )`
      `H = I - tau * u *u',    u  = (   ),`
                                   `( v )`
where tau is a real scalar and v is a real m-vector.

If tau = 0, then H is taken to be the unit matrix.

In-line code is used if H has order < 11.



See the SLICOT documentation for details.
"""
function mb04oy! end

"""

Reduce a Hamiltonian like matrix

              `[  A   G  ]           T          T`
         `H =  [       T ] ,    G = G ,    Q = Q,`
              `[  Q  -A  ]`

or a skew-Hamiltonian like matrix

              `[  A   G  ]            T          T`
         `W =  [       T ] ,    G = -G ,   Q = -Q,`
              `[  Q   A  ]`

so that elements below the (k+1)-th subdiagonal in the first nb
columns of the (k+n)-by-n matrix A, and offdiagonal elements
in the first nb columns and rows of the n-by-n matrix Q are zero.

The reduction is performed by an orthogonal symplectic
transformation UU'*H*UU and matrices U, XA, XG, XQ, and YA are
returned so that

               `[ Aout + U*XA'+ YA*U'   Gout + U*XG'+ XG*U' ]`
    `UU'*H*UU = [                                           ].`
               `[ Qout + U*XQ'+ XQ*U'  -Aout'- XA*U'- U*YA' ]`

Similarly,

               `[ Aout + U*XA'+ YA*U'   Gout + U*XG'- XG*U' ]`
    `UU'*W*UU = [                                           ].`
               `[ Qout + U*XQ'- XQ*U'   Aout'+ XA*U'+ U*YA' ]`

This is an auxiliary routine called by MB04PB.



See the SLICOT documentation for details.
"""
function mb04pa! end

"""

Reduce a Hamiltonian matrix,

              `[  A   G  ]`
         `H =  [       T ] ,`
              `[  Q  -A  ]`

where A is an N-by-N matrix and G,Q are N-by-N symmetric matrices,
to Paige/Van Loan (PVL) form. That is, an orthogonal symplectic U
is computed so that

          `T       [  Aout   Gout  ]`
         `U H U =  [             T ] ,`
                  `[  Qout  -Aout  ]`

where Aout is upper Hessenberg and Qout is diagonal.
Blocked version.



See the SLICOT documentation for details.
"""
function mb04pb! end

"""

Reduce a Hamiltonian matrix,

              `[  A   G  ]`
         `H =  [       T ] ,`
              `[  Q  -A  ]`

where A is an N-by-N matrix and G,Q are N-by-N symmetric matrices,
to Paige/Van Loan (PVL) form. That is, an orthogonal symplectic U
is computed so that

          `T       [  Aout   Gout  ]`
         `U H U =  [             T ] ,`
                  `[  Qout  -Aout  ]`

where Aout is upper Hessenberg and Qout is diagonal.
Unblocked version.



See the SLICOT documentation for details.
"""
function mb04pu! end

"""

Apply a real elementary reflector H to a real m-by-n matrix
C, from either the left or the right. H is represented in the form
                                   `( 1 )`
      `H = I - tau * u *u',    u  = (   ),`
                                   `( v )`
where tau is a real scalar and v is a real vector.

If tau = 0, then H is taken to be the unit matrix.

In-line code is used if H has order < 11.



See the SLICOT documentation for details.
"""
function mb04py! end

"""

Overwrite general real m-by-n matrices C and D, or their
transposes, with

          `[ op(C) ]`
    `Q  *  [       ]   if TRANQ = 'N', or`
          `[ op(D) ]`

     `T    [ op(C) ]`
    `Q  *  [       ]   if TRANQ = 'T',`
          `[ op(D) ]`

where Q is defined as the product of symplectic reflectors and
Givens rotations,

    `Q = diag( H(1),H(1) ) G(1) diag( F(1),F(1) )`
        `diag( H(2),H(2) ) G(2) diag( F(2),F(2) )`
                          `....`
        `diag( H(k),H(k) ) G(k) diag( F(k),F(k) ).`

Blocked version.



See the SLICOT documentation for details.
"""
function mb04qb! end

"""

Apply the orthogonal symplectic block reflector

         `[  I+V*T*V'  V*R*S*V'  ]`
    `Q =  [                      ]`
         `[ -V*R*S*V'  I+V*T*V'  ]`

or its transpose to a real 2m-by-n matrix [ op(A); op(B) ] from
the left.
The k-by-k upper triangular blocks of the matrices

                            `[ S1 ]       [ T11 T12 T13 ]`
    `R  = [ R1 R2 R3 ],  S = [ S2 ],  T = [ T21 T22 T23 ],`
                            `[ S3 ]       [ T31 T32 T33 ]`

with R2 unit and S1, R3, T21, T31, T32 strictly upper triangular,
are stored rowwise in the arrays RS and T, respectively.



See the SLICOT documentation for details.
"""
function mb04qc! end

"""

Form the triangular block factors R, S and T of a symplectic
block reflector SH, which is defined as a product of 2k
concatenated Householder reflectors and k Givens rotations,

    `SH = diag( H(1),H(1) ) G(1) diag( F(1),F(1) )`
         `diag( H(2),H(2) ) G(2) diag( F(2),F(2) )`
                           `....`
         `diag( H(k),H(k) ) G(k) diag( F(k),F(k) ).`

The upper triangular blocks of the matrices

                            `[ S1 ]       [ T11 T12 T13 ]`
    `R  = [ R1 R2 R3 ],  S = [ S2 ],  T = [ T21 T22 T23 ],`
                            `[ S3 ]       [ T31 T32 T33 ]`

with R2 unit and S1, R3, T21, T31, T32 strictly upper triangular,
are stored rowwise in the arrays RS and T, respectively.



See the SLICOT documentation for details.
"""
function mb04qf! end

"""

Overwrites general real m-by-n/n-by-m matrices C and D with

          `[ op(C) ]`
    `U  *  [       ]   if TRANU = 'N', or`
          `[ op(D) ]`

     `T    [ op(C) ]`
    `U  *  [       ]   if TRANU = 'T',`
          `[ op(D) ]`

where U is defined as the product of symplectic reflectors and
Givens rotations,

    `U = diag( H(1),H(1) ) G(1) diag( F(1),F(1) )`
        `diag( H(2),H(2) ) G(2) diag( F(2),F(2) )`
                          `....`
        `diag( H(k),H(k) ) G(k) diag( F(k),F(k) ),`

with k = m-1, as returned by the SLICOT Library routines MB04PU
or MB04RU.



See the SLICOT documentation for details.
"""
function mb04qs! end

"""

Overwrite general real m-by-n matrices C and D, or their
transposes, with

          `[ op(C) ]`
    `Q  *  [       ]   if TRANQ = 'N', or`
          `[ op(D) ]`

     `T    [ op(C) ]`
    `Q  *  [       ]   if TRANQ = 'T',`
          `[ op(D) ]`

where Q is defined as the product of symplectic reflectors and
Givens rotations,

    `Q = diag( H(1),H(1) ) G(1) diag( F(1),F(1) )`
        `diag( H(2),H(2) ) G(2) diag( F(2),F(2) )`
                          `....`
        `diag( H(k),H(k) ) G(k) diag( F(k),F(k) ).`

Unblocked version.



See the SLICOT documentation for details.
"""
function mb04qu! end

"""

Reduce a skew-Hamiltonian matrix,

              `[  A   G  ]`
        `W  =  [       T ] ,`
              `[  Q   A  ]`

where A is an N-by-N matrix and G, Q are N-by-N skew-symmetric
matrices, to Paige/Van Loan (PVL) form. That is, an orthogonal
symplectic matrix U is computed so that

          `T       [  Aout  Gout  ]`
         `U W U =  [            T ] ,`
                  `[    0   Aout  ]`

where Aout is in upper Hessenberg form.
Blocked version.



See the SLICOT documentation for details.
"""
function mb04rb! end

"""

Reduce a skew-Hamiltonian matrix,

              `[  A   G  ]`
        `W  =  [       T ] ,`
              `[  Q   A  ]`

where A is an N-by-N matrix and G, Q are N-by-N skew-symmetric
matrices, to Paige/Van Loan (PVL) form. That is, an orthogonal
symplectic matrix U is computed so that

          `T       [  Aout  Gout  ]`
         `U W U =  [            T ] ,`
                  `[    0   Aout  ]`

where Aout is in upper Hessenberg form.
Unblocked version.



See the SLICOT documentation for details.
"""
function mb04ru! end

"""

Compute a symplectic QR decomposition of a real 2M-by-N matrix
[A; B],

          `[ A ]             [ R11  R12 ]`
          `[   ] = Q * R = Q [          ],`
          `[ B ]             [ R21  R22 ]`

where Q is a symplectic orthogonal matrix, R11 is upper triangular
and R21 is strictly upper triangular.
If [A; B] is symplectic then, theoretically, R21 = 0 and
R22 = inv(R11)^T. Unblocked version.



See the SLICOT documentation for details.
"""
function mb04su! end

"""

Compute a symplectic URV (SURV) decomposition of a real
2N-by-2N matrix H,

       `[ op(A)   G   ]                 [ op(R11)   R12   ]`
   `H = [             ] = U R V'  = U * [                 ] * V' ,`
       `[  Q    op(B) ]                 [   0     op(R22) ]`

where A, B, G, Q, R12 are real N-by-N matrices, op(R11) is a real
N-by-N upper triangular matrix, op(R22) is a real N-by-N lower
Hessenberg matrix and U, V are 2N-by-2N orthogonal symplectic
matrices. Blocked version.



See the SLICOT documentation for details.
"""
function mb04tb! end

"""

Compute a symplectic URV (SURV) decomposition of a real
2N-by-2N matrix H:

        `[ op(A)   G   ]        T       [ op(R11)   R12   ]    T`
    `H = [             ] = U R V  = U * [                 ] * V ,`
        `[  Q    op(B) ]                [   0     op(R22) ]`

where A, B, G, Q, R12 are real N-by-N matrices, op(R11) is a real
N-by-N upper triangular matrix, op(R22) is a real N-by-N lower
Hessenberg matrix and U, V are 2N-by-2N orthogonal symplectic
matrices. Unblocked version.



See the SLICOT documentation for details.
"""
function mb04ts! end

"""

Let A and E be M-by-N matrices with E in column echelon form.
Let AA and EE be the following submatrices of A and E:
  `AA := A(IFIRA : M ; IFICA : N)`
  `EE := E(IFIRA : M ; IFICA : N).`
Let Aj and Ej be the following submatrices of AA and EE:
  `Aj := A(IFIRA : M ; IFICA : IFICA + NCA - 1) and`
  `Ej := E(IFIRA : M ; IFICA + NCA : N).`

Transform (AA,EE) such that Aj is row compressed while keeping
matrix Ej in column echelon form (which may be different from the
form on entry).
In fact the routine performs the j-th step of Algorithm 3.2.1 in
[1]. Furthermore, it determines the rank RANK of the submatrix Ej,
which is equal to the number of corner points in submatrix Ej.



See the SLICOT documentation for details.
"""
function mb04tt! end

"""

Perform the Givens transformation, defined by C (cos) and S
(sin), and interchange the vectors involved, i.e.

   `|X(i)|    | 0   1 |   | C   S |   |X(i)|`
   `|    | := |       | x |       | x |    |, i = 1,...N.`
   `|Y(i)|    | 1   0 |   |-S   C |   |Y(i)|`

REMARK. This routine is a modification of DROT from BLAS.
        `This routine is called only by the SLICOT routines MB04TX`
        `and MB04VX.`

NUMERICAL ASPECTS

The algorithm is backward stable.



See the SLICOT documentation for details.
"""
function mb04tu! end

"""

Reduce a submatrix A(k) of A to upper triangular form by column
Givens rotations only.
Here A(k) = A(IFIRA:ma,IFICA:na) where ma = IFIRA - 1 + NRA,
na = IFICA - 1 + NCA.
Matrix A(k) is assumed to have full row rank on entry. Hence, no
pivoting is done during the reduction process. See Algorithm 2.3.1
and Remark 2.3.4 in [1].
The constructed column transformations are also applied to matrix
E(k) = E(1:IFIRA-1,IFICA:na).
Note that in E columns are transformed with the same column
indices as in A, but with row indices different from those in A.



See the SLICOT documentation for details.
"""
function mb04tv! end

"""

Reduce a submatrix E(k) of E to upper triangular form by row
Givens rotations only.
Here E(k) = E(IFIRE:me,IFICE:ne), where me = IFIRE - 1 + NRE,
                                        `ne = IFICE - 1 + NCE.`
Matrix E(k) is assumed to have full column rank on entry. Hence,
no pivoting is done during the reduction process. See Algorithm
2.3.1 and Remark 2.3.4 in [1].
The constructed row transformations are also applied to matrix
A(k) = A(IFIRE:me,IFICA:N).
Note that in A(k) rows are transformed with the same row indices
as in E but with column indices different from those in E.



See the SLICOT documentation for details.
"""
function mb04tw! end

"""

Separate the pencils s*E(eps)-A(eps) and s*E(inf)-A(inf) in
s*E(eps,inf)-A(eps,inf) using Algorithm 3.3.3 in [1].

On entry, it is assumed that the M-by-N matrices A and E have
been obtained after applying the Algorithms 3.2.1 and 3.3.1 to
the pencil s*E - A as described in [1], i.e.

                   `| s*E(eps,inf)-A(eps,inf) |      X      |`
   `Q'(s*E - A)Z  = |-------------------------|-------------|`
                   `|             0           | s*E(r)-A(r) |`

Here the pencil s*E(eps,inf)-A(eps,inf) is in staircase form.
This pencil contains all Kronecker column indices and infinite
elementary divisors of the pencil s*E - A.
The pencil s*E(r)-A(r) contains all Kronecker row indices and
finite elementary divisors of s*E - A.
Furthermore, the submatrices having full row and column rank in
the pencil s*E(eps,inf)-A(eps,inf) are assumed to be
triangularized.

On exit, the result then is

                   `Q'(s*E - A)Z =`

     `| s*E(eps)-A(eps) |        X        |      X      |`
     `|-----------------|-----------------|-------------|`
     `|        0        | s*E(inf)-A(inf) |      X      |`
     `|===================================|=============|`
     `|                                   |             |`
     `|                 0                 | s*E(r)-A(r) |`

Note that the pencil s*E(r)-A(r) is not reduced further.



See the SLICOT documentation for details.
"""
function mb04tx! end

"""

Perform the triangularization of the submatrices having full
row and column rank in the pencil s*E(eps,inf)-A(eps,inf) below

               `| s*E(eps,inf)-A(eps,inf) |     X       |`
     `s*E - A = |-------------------------|-------------| ,`
               `|            0            | s*E(r)-A(r) |`

using Algorithm 3.3.1 in [1].
On entry, it is assumed that the M-by-N matrices A and E have
been transformed to generalized Schur form by unitary
transformations (see Algorithm 3.2.1 in [1]), and that the pencil
s*E(eps,inf)-A(eps,inf) is in staircase form.
This pencil contains all Kronecker column indices and infinite
elementary divisors of the pencil s*E - A.
The pencil s*E(r)-A(r) contains all Kronecker row indices and
finite elementary divisors of s*E - A.



See the SLICOT documentation for details.
"""
function mb04ty! end

"""

Compute orthogonal transformations Q and Z such that the
transformed pencil Q'(sE-A)Z has the E matrix in column echelon
form, where E and A are M-by-N matrices.



See the SLICOT documentation for details.
"""
function mb04ud! end

"""

Compute orthogonal transformations Q and Z such that the
transformed pencil Q'(sE-A)Z is in upper block triangular form,
where E is an M-by-N matrix in column echelon form (see SLICOT
Library routine MB04UD) and A is an M-by-N matrix.

If MODE = 'B', then the matrices A and E are transformed into the
following generalized Schur form by unitary transformations Q1
and Z1 :

                 `| sE(eps,inf)-A(eps,inf) |      X     |`
   `Q1'(sE-A)Z1 = |------------------------|------------|.   (1)`
                 `|            O           | sE(r)-A(r) |`

The pencil sE(eps,inf)-A(eps,inf) is in staircase form, and it
contains all Kronecker column indices and infinite elementary
divisors of the pencil sE-A. The pencil sE(r)-A(r) contains all
Kronecker row indices and elementary divisors of sE-A.
Note: X is a pencil.

If MODE = 'T', then the submatrices having full row and column
rank in the pencil sE(eps,inf)-A(eps,inf) in (1) are
triangularized by applying unitary transformations Q2 and Z2 to
Q1'*(sE-A)*Z1.

If MODE = 'S', then the pencil sE(eps,inf)-A(eps,inf) in (1) is
separated into sE(eps)-A(eps) and sE(inf)-A(inf) by applying
unitary transformations Q3 and Z3 to Q2'*Q1'*(sE-A)*Z1*Z2.

This gives

           `| sE(eps)-A(eps) |        X       |      X     |`
           `|----------------|----------------|------------|`
           `|        O       | sE(inf)-A(inf) |      X     |`
Q'(sE-A)Z =|=================================|============| (2)
           `|                                 |            |`
           `|                O                | sE(r)-A(r) |`

where Q = Q1*Q2*Q3 and Z = Z1*Z2*Z3.
Note: the pencil sE(r)-A(r) is not reduced further.



See the SLICOT documentation for details.
"""
function mb04vd! end

"""

Separate the pencils s*E(eps)-A(eps) and s*E(inf)-A(inf) in
s*E(eps,inf)-A(eps,inf) using Algorithm 3.3.3 in [1].

On entry, it is assumed that the M-by-N matrices A and E have
been obtained after applying the Algorithms 3.2.1 and 3.3.1 to
the pencil s*E - A as described in [1], i.e.

                   `| s*E(eps,inf)-A(eps,inf) |      X      |`
   `Q'(s*E - A)Z  = |-------------------------|-------------|`
                   `|             0           | s*E(r)-A(r) |`

Here the pencil s*E(eps,inf)-A(eps,inf) is in staircase form.
This pencil contains all Kronecker column indices and infinite
elementary divisors of the pencil s*E - A.
The pencil s*E(r)-A(r) contains all Kronecker row indices and
finite elementary divisors of s*E - A.
Furthermore, the submatrices having full row and column rank in
the pencil s*E(eps,inf)-A(eps,inf) are assumed to be
triangularized.

On exit, the result then is

                   `Q'(s*E - A)Z =`

     `| s*E(eps)-A(eps) |        X        |      X      |`
     `|-----------------|-----------------|-------------|`
     `|        0        | s*E(inf)-A(inf) |      X      |`
     `|===================================|=============|`
     `|                                   |             |`
     `|                 0                 | s*E(r)-A(r) |`

Note that the pencil s*E(r)-A(r) is not reduced further.



See the SLICOT documentation for details.
"""
function mb04vx! end

"""

Generate a matrix Q with orthogonal columns (spanning an
isotropic subspace), which is defined as the first n columns
of a product of symplectic reflectors and Givens rotations,

    `Q = diag( H(1),H(1) ) G(1) diag( F(1),F(1) )`
        `diag( H(2),H(2) ) G(2) diag( F(2),F(2) )`
                          `....`
        `diag( H(k),H(k) ) G(k) diag( F(k),F(k) ).`

The matrix Q is returned in terms of its first 2*M rows

                 `[  op( Q1 )   op( Q2 ) ]`
             `Q = [                      ].`
                 `[ -op( Q2 )   op( Q1 ) ]`

Blocked version of the SLICOT Library routine MB04WU.



See the SLICOT documentation for details.
"""
function mb04wd! end

"""

Generate an orthogonal symplectic matrix U, which is defined as
a product of symplectic reflectors and Givens rotations

U = diag( H(1),H(1) )      G(1)  diag( F(1),F(1) )
    `diag( H(2),H(2) )      G(2)  diag( F(2),F(2) )`
                           `....`
    `diag( H(n-1),H(n-1) ) G(n-1) diag( F(n-1),F(n-1) ).`

as returned by MB04PU. The matrix U is returned in terms of its
first N rows

                 `[  U1   U2 ]`
             `U = [          ].`
                 `[ -U2   U1 ]`



See the SLICOT documentation for details.
"""
function mb04wp! end

"""

Generate orthogonal symplectic matrices U or V, defined as
products of symplectic reflectors and Givens rotations

U = diag( HU(1),HU(1) )  GU(1)  diag( FU(1),FU(1) )
    `diag( HU(2),HU(2) )  GU(2)  diag( FU(2),FU(2) )`
                         `....`
    `diag( HU(n),HU(n) )  GU(n)  diag( FU(n),FU(n) ),`

V = diag( HV(1),HV(1) )       GV(1)   diag( FV(1),FV(1) )
    `diag( HV(2),HV(2) )       GV(2)   diag( FV(2),FV(2) )`
                              `....`
    `diag( HV(n-1),HV(n-1) )  GV(n-1)  diag( FV(n-1),FV(n-1) ),`

as returned by the SLICOT Library routines MB04TS or MB04TB. The
matrices U and V are returned in terms of their first N/2 rows:

            `[  U1   U2 ]           [  V1   V2 ]`
        `U = [          ],      V = [          ].`
            `[ -U2   U1 ]           [ -V2   V1 ]`



See the SLICOT documentation for details.
"""
function mb04wr! end

"""

Generate a matrix Q with orthogonal columns (spanning an
isotropic subspace), which is defined as the first n columns
of a product of symplectic reflectors and Givens rotations,

    `Q = diag( H(1),H(1) ) G(1) diag( F(1),F(1) )`
        `diag( H(2),H(2) ) G(2) diag( F(2),F(2) )`
                          `....`
        `diag( H(k),H(k) ) G(k) diag( F(k),F(k) ).`

The matrix Q is returned in terms of its first 2*M rows

                 `[  op( Q1 )   op( Q2 ) ]`
             `Q = [                      ].`
                 `[ -op( Q2 )   op( Q1 ) ]`



See the SLICOT documentation for details.
"""
function mb04wu! end

"""

Compute a basis for the left and/or right singular subspace of
an M-by-N matrix A corresponding to its smallest singular values.



See the SLICOT documentation for details.
"""
function mb04xd! end

"""

Apply the Householder transformations Pj stored in factored
form into the columns of the array X, to the desired columns of
the matrix U by premultiplication, and/or the Householder
transformations Qj stored in factored form into the rows of the
array X, to the desired columns of the matrix V by
premultiplication. The Householder transformations Pj and Qj
are stored as produced by LAPACK Library routine DGEBRD.



See the SLICOT documentation for details.
"""
function mb04xy! end

"""

Partially diagonalize the bidiagonal matrix

          `|q(1) e(1)  0    ...       0      |`
          `| 0   q(2) e(2)            .      |`
      `J = | .                        .      |                  (1)`
          `| .                  e(MIN(M,N)-1)|`
          `| 0   ...        ...  q(MIN(M,N)) |`

using QR or QL iterations in such a way that J is split into
unreduced bidiagonal submatrices whose singular values are either
all larger than a given bound or are all smaller than (or equal
to) this bound. The left- and right-hand Givens rotations
performed on J (corresponding to each QR or QL iteration step) may
be optionally accumulated in the arrays U and V.



See the SLICOT documentation for details.
"""
function mb04yd! end

"""

Perform either one QR or QL iteration step onto the unreduced
bidiagonal submatrix Jk:

         `|D(l) E(l)    0  ...    0   |`
         `| 0   D(l+1) E(l+1)     .   |`
    `Jk = | .                     .   |`
         `| .                     .   |`
         `| .                   E(k-1)|`
         `| 0   ...        ...   D(k) |`

with k <= p and l >= 1, p = MIN(M,N), of the bidiagonal matrix J:

         `|D(1) E(1)  0    ...   0   |`
         `| 0   D(2) E(2)        .   |`
     `J = | .                    .   |.`
         `| .                    .   |`
         `| .                  E(p-1)|`
         `| 0   ...        ...  D(p) |`

Hereby, Jk is transformed to  S' Jk T with S and T products of
Givens rotations. These Givens rotations S (respectively, T) are
postmultiplied into U (respectively, V), if UPDATU (respectively,
UPDATV) is .TRUE..



See the SLICOT documentation for details.
"""
function mb04yw! end

"""

Transform a Hamiltonian matrix

          `( A   G  )`
      `H = (      T )                                           (1)`
          `( Q  -A  )`

into a square-reduced Hamiltonian matrix

           `( A'  G'  )`
      `H' = (       T )                                         (2)`
           `( Q' -A'  )`
                                                            `T`
by an orthogonal symplectic similarity transformation H' = U H U,
where
          `(  U1   U2 )`
      `U = (          ).                                        (3)`
          `( -U2   U1 )`
                                                         `T`
The square-reduced Hamiltonian matrix satisfies Q'A' - A' Q' = 0,
and

      `2       T     2     ( A''   G''  )`
    `H'  :=  (U  H U)   =  (          T ).`
                          `( 0     A''  )`

In addition, A'' is upper Hessenberg and G'' is skew symmetric.
The square roots of the eigenvalues of A'' = A'*A' + G'*Q' are the
eigenvalues of H.



See the SLICOT documentation for details.
"""
function mb04zd! end

"""

Compute exp(A*delta) where A is a real N-by-N non-defective
matrix with real or complex eigenvalues and delta is a scalar
value. The routine also returns the eigenvalues and eigenvectors
of A as well as (if all eigenvalues are real) the matrix product
exp(Lambda*delta) times the inverse of the eigenvector matrix
of A, where Lambda is the diagonal matrix of eigenvalues.
Optionally, the routine computes a balancing transformation to
improve the conditioning of the eigenvalues and eigenvectors.



See the SLICOT documentation for details.
"""
function mb05md! end

"""

Compute, for an N-by-N real nonsymmetric matrix A, the
orthogonal matrix Q reducing it to real Schur form T, the
eigenvalues, and the right eigenvectors of T.

The right eigenvector r(j) of T satisfies
                 `T * r(j) = lambda(j) * r(j)`
where lambda(j) is its eigenvalue.

The matrix of right eigenvectors R is upper triangular, by
construction.



See the SLICOT documentation for details.
"""
function mb05my! end

"""

Compute

(a)    F(delta) =  exp(A*delta) and

(b)    H(delta) =  Int[F(s) ds] from s = 0 to s = delta,

where A is a real N-by-N matrix and delta is a scalar value.



See the SLICOT documentation for details.
"""
function mb05nd! end

"""

Compute exp(A*delta) where A is a real N-by-N matrix and delta
is a scalar value. The routine also returns the minimal number of
accurate digits in the 1-norm of exp(A*delta) and the number of
accurate digits in the 1-norm of exp(A*delta) at 95% confidence
level.



See the SLICOT documentation for details.
"""
function mb05od! end

"""

Restore a matrix after it has been transformed by applying
balancing transformations (permutations and scalings), as
determined by LAPACK Library routine DGEBAL.



See the SLICOT documentation for details.
"""
function mb05oy! end

"""

Move the eigenvalues with strictly negative real parts of an
N-by-N complex skew-Hamiltonian/Hamiltonian pencil aS - bH in
structured Schur form to the leading principal subpencil, while
keeping the triangular form. On entry, we have

      `(  A  D  )      (  B  F  )`
  `S = (        ), H = (        ),`
      `(  0  A' )      (  0 -B' )`

where A and B are upper triangular.
S and H are transformed by a unitary matrix Q such that

                       `(  Aout  Dout  )`
  `Sout = J Q' J' S Q = (              ), and`
                       `(    0   Aout' )`
                                                               `(1)`
                       `(  Bout  Fout  )           (  0  I  )`
  `Hout = J Q' J' H Q = (              ), with J = (        ),`
                       `(    0  -Bout' )           ( -I  0  )`

where Aout and Bout remain in upper triangular form. The notation
M' denotes the conjugate transpose of the matrix M.
Optionally, if COMPQ = 'I' or COMPQ = 'U', the unitary matrix Q
that fulfills (1) is computed.



See the SLICOT documentation for details.
"""
function mb3jzp! end

"""

Compute the eigenvalues of a complex N-by-N skew-Hamiltonian/
Hamiltonian pencil aS - bH, with

      `(  A  D  )         (  B  F  )`
  `S = (        ) and H = (        ).                           (1)`
      `(  E  A' )         (  G -B' )`

The structured Schur form of the embedded real skew-Hamiltonian/
skew-Hamiltonian pencil a`B_S` - b`B_T`, defined as

        `(  Re(A)  -Im(A)  |  Re(D)  -Im(D)  )`
        `(                 |                 )`
        `(  Im(A)   Re(A)  |  Im(D)   Re(D)  )`
        `(                 |                 )`
  `B_S = (-----------------+-----------------) , and`
        `(                 |                 )`
        `(  Re(E)  -Im(E)  |  Re(A')  Im(A') )`
        `(                 |                 )`
        `(  Im(E)   Re(E)  | -Im(A')  Re(A') )`
                                                               `(2)`
        `( -Im(B)  -Re(B)  | -Im(F)  -Re(F)  )`
        `(                 |                 )`
        `(  Re(B)  -Im(B)  |  Re(F)  -Im(F)  )`
        `(                 |                 )`
  `B_T = (-----------------+-----------------) ,  T = i*H,`
        `(                 |                 )`
        `( -Im(G)  -Re(G)  | -Im(B')  Re(B') )`
        `(                 |                 )`
        `(  Re(G)  -Im(G)  | -Re(B') -Im(B') )`

is determined and used to compute the eigenvalues. The notation M'
denotes the conjugate transpose of the matrix M. Optionally,
if COMPQ = 'C', an orthonormal basis of the right deflating
subspace of the pencil aS - bH, corresponding to the eigenvalues
with strictly negative real part, is computed. Namely, after
transforming a`B_S` - b`B_H` by unitary matrices, we have

           `( BA  BD  )              ( BB  BF  )`
  `B_Sout = (         ) and B_Hout = (         ),               (3)`
           `(  0  BA' )              (  0 -BB' )`

and the eigenvalues with strictly negative real part of the
complex pencil a`B_S`out - b`B_H`out are moved to the top. The
embedding doubles the multiplicities of the eigenvalues of the
pencil aS - bH.



See the SLICOT documentation for details.
"""
function mb3lzp! end

"""

Compute a rank-revealing QR factorization of a complex general
M-by-N matrix  A,  which may be rank-deficient, and estimate its
effective rank using incremental condition estimation.

The routine uses a truncated QR factorization with column pivoting
                              `[ R11 R12 ]`
   `A * P = Q * R,  where  R = [         ],`
                              `[  0  R22 ]`
with R11 defined as the largest leading upper triangular submatrix
whose estimated condition number is less than 1/RCOND.  The order
of R11, RANK, is the effective rank of A.  Condition estimation is
performed during the QR factorization process.  Matrix R22 is full
(but of small norm), or empty.

MB3OYZ  does not perform any scaling of the matrix A.



See the SLICOT documentation for details.
"""
function mb3oyz! end

"""

Compute a rank-revealing RQ factorization of a complex general
M-by-N matrix  A,  which may be rank-deficient, and estimate its
effective rank using incremental condition estimation.

The routine uses a truncated RQ factorization with row pivoting:
                              `[ R11 R12 ]`
   `P * A = R * Q,  where  R = [         ],`
                              `[  0  R22 ]`
with R22 defined as the largest trailing upper triangular
submatrix whose estimated condition number is less than 1/RCOND.
The order of R22, RANK, is the effective rank of A.  Condition
estimation is performed during the RQ factorization process.
Matrix R11 is full (but of small norm), or empty.

MB3PYZ  does not perform any scaling of the matrix A.



See the SLICOT documentation for details.
"""
function mb3pyz! end

"""

Apply from the left the inverse of a balancing transformation,
computed by the SLICOT Library routine MB4DPZ, to the complex
matrix

          `[   V1   ]`
          `[        ],`
          `[ sgn*V2 ]`

where sgn is either +1 or -1.



See the SLICOT documentation for details.
"""
function mb4dbz! end

"""

Balance a pair of N-by-N complex matrices (A,B). This involves,
first, permuting A and B by equivalence transformations to isolate
eigenvalues in the first 1 to ILO-1 and last IHI+1 to N elements
on the diagonal of A and B; and second, applying a diagonal
equivalence transformation to rows and columns ILO to IHI to make
the rows and columns as close in 1-norm as possible. Both steps
are optional. Balancing may reduce the 1-norms of the matrices,
and improve the accuracy of the computed eigenvalues and/or
eigenvectors in the generalized eigenvalue problem
A*x = lambda*B*x.

This routine may optionally improve the conditioning of the
scaling transformation compared to the LAPACK routine ZGGBAL.



See the SLICOT documentation for details.
"""
function mb4dlz! end

"""

Balance the 2*N-by-2*N complex skew-Hamiltonian/Hamiltonian
pencil aS - bH, with

      `(  A  D  )         (  C  V  )`
  `S = (        ) and H = (        ),  A, C N-by-N,             (1)`
      `(  E  A' )         (  W -C' )`

where D and E are skew-Hermitian, V and W are Hermitian matrices,
and ' denotes conjugate transpose. This involves, first, permuting
aS - bH by a symplectic equivalence transformation to isolate
eigenvalues in the first 1:ILO-1 elements on the diagonal of A
and C; and second, applying a diagonal equivalence transformation
to make the pairs of rows and columns ILO:N and N+ILO:2*N as close
in 1-norm as possible. Both steps are optional. Balancing may
reduce the 1-norms of the matrices S and H.



See the SLICOT documentation for details.
"""
function mb4dpz! end

"""

Calculate, for a given real polynomial P(x) and a real scalar
alpha, the leading K coefficients of the shifted polynomial
                                                          `K-1`
   `P(x) = q(1) + q(2) * (x-alpha) + ... + q(K) * (x-alpha)   + ...`

using Horner's algorithm.



See the SLICOT documentation for details.
"""
function mc01md! end

"""

Compute the value of the real polynomial P(x) at a given
complex point x = x0 using Horner's algorithm.



See the SLICOT documentation for details.
"""
function mc01nd! end

"""

Compute the coefficients of a complex polynomial P(x) from its
zeros.



See the SLICOT documentation for details.
"""
function mc01od! end

"""

Compute the coefficients of a real polynomial P(x) from its
zeros.



See the SLICOT documentation for details.
"""
function mc01pd! end

"""

Compute the coefficients of a real polynomial P(x) from its
zeros. The coefficients are stored in decreasing order of the
powers of x.



See the SLICOT documentation for details.
"""
function mc01py! end

"""

Compute, for two given real polynomials A(x) and B(x), the
quotient polynomial Q(x) and the remainder polynomial R(x) of
A(x) divided by B(x).

The polynomials Q(x) and R(x) satisfy the relationship

   `A(x) = B(x) * Q(x) + R(x),`

where the degree of R(x) is less than the degree of B(x).



See the SLICOT documentation for details.
"""
function mc01qd! end

"""

Compute the coefficients of the polynomial

   `P(x) = P1(x) * P2(x) + alpha * P3(x),`

where P1(x), P2(x) and P3(x) are given real polynomials and alpha
is a real scalar.

Each of the polynomials P1(x), P2(x) and P3(x) may be the zero
polynomial.



See the SLICOT documentation for details.
"""
function mc01rd! end

"""

Scale the coefficients of the real polynomial P(x) such that
the coefficients of the scaled polynomial Q(x) = sP(tx) have
minimal variation, where s and t are real scalars.



See the SLICOT documentation for details.
"""
function mc01sd! end

"""

Find the mantissa M and the exponent E of a real number A such
that
   `A = M * B**E`
   `1 <= ABS( M ) < B`
if A is non-zero. If A is zero, then M and E are set to 0.



See the SLICOT documentation for details.
"""
function mc01sw! end

"""

Compute the variation V of the exponents of a series of
non-zero floating-point numbers: a(j) = MANT(j) * beta**(E(j)),
where beta is the base of the machine representation of
floating-point numbers, i.e.,
V = max(E(j)) - min(E(j)), j = LB,...,UB and MANT(j) non-zero.



See the SLICOT documentation for details.
"""
function mc01sx! end

"""

Find a real number A from its mantissa M and its exponent E,
i.e.,
   `A = M * B**E.`
M and E need not be the standard floating-point values.
If ABS(A) < B**(EMIN-1), i.e. the smallest positive model number,
then the routine returns A = 0.
If M = 0, then the routine returns A = 0 regardless of the value
of E.



See the SLICOT documentation for details.
"""
function mc01sy! end

"""

Determine whether or not a given polynomial P(x) with real
coefficients is stable, either in the continuous-time or discrete-
time case.

A polynomial is said to be stable in the continuous-time case
if all its zeros lie in the left half-plane, and stable in the
discrete-time case if all its zeros lie inside the unit circle.



See the SLICOT documentation for details.
"""
function mc01td! end

"""

Compute the roots of a quadratic equation with real
coefficients.



See the SLICOT documentation for details.
"""
function mc01vd! end

"""

Compute, for a given real polynomial P(x) and a quadratic
polynomial B(x), the quotient polynomial Q(x) and the linear
remainder polynomial R(x) such that

   `P(x) = B(x) * Q(x) + R(x),`

                            `2`
where B(x) = u1 + u2 * x + x , R(x) = q(1) + q(2) * (u2 + x)
and u1, u2, q(1) and q(2) are real scalars.



See the SLICOT documentation for details.
"""
function mc01wd! end

"""

Compute the roots of the polynomial

    `P(t) = ALPHA + BETA*t + GAMMA*t^2 + DELTA*t^3 .`



See the SLICOT documentation for details.
"""
function mc01xd! end

"""

Compute the coefficients of the real polynomial matrix

   `P(x) = P1(x) * P2(x) + alpha * P3(x),`

where P1(x), P2(x) and P3(x) are given real polynomial matrices
and alpha is a real scalar.

Each of the polynomial matrices P1(x), P2(x) and P3(x) may be the
zero matrix.



See the SLICOT documentation for details.
"""
function mc03md! end

"""

Compute the coefficients of a minimal polynomial basis
                                            `DK`
    `K(s) = K(0) + K(1) * s + ... + K(DK) * s`

for the right nullspace of the MP-by-NP polynomial matrix of
degree DP, given by
                                            `DP`
    `P(s) = P(0) + P(1) * s + ... + P(DP) * s  ,`

which corresponds to solving the polynomial matrix equation
P(s) * K(s) = 0.



See the SLICOT documentation for details.
"""
function mc03nd! end

"""

Given an MP-by-NP polynomial matrix of degree dp
                               `dp-1            dp`
P(s) = P(0) + ... + P(dp-1) * s     + P(dp) * s            (1)

the routine composes the related pencil s*E-A where

    `| I              |           | O          -P(dp) |`
    `|   .            |           | I .           .   |`
A = |     .          |  and  E = |   . .         .   |.    (2)
    `|       .        |           |     . O       .   |`
    `|         I      |           |       I  O -P(2)  |`
    `|           P(0) |           |          I -P(1)  |`

==================================================================
REMARK: This routine is intended to be called only from the SLICOT
        `routine MC03ND.`
==================================================================



See the SLICOT documentation for details.
"""
function mc03nx! end

"""

Determine a minimal basis of the right nullspace of the
subpencil s*E(eps)-A(eps) using the method given in [1] (see
Eqs.(4.6.8), (4.6.9)).
This pencil only contains Kronecker column indices, and it must be
in staircase form as supplied by SLICOT Library Routine MB04VD.
The basis vectors are represented by matrix V(s) having the form

           `| V11(s) V12(s) V13(s)   . .   V1n(s) |`
           `|        V22(s) V23(s)         V2n(s) |`
           `|               V33(s)           .    |`
    `V(s) = |                  .             .    |`
           `|                      .         .    |`
           `|                          .     .    |`
           `|                              Vnn(s) |`

where n is the number of full row rank blocks in matrix A(eps) and

                                          `k               j-i`
    `Vij(s) = Vij,0 + Vij,1*s +...+ Vij,k*s +...+ Vij,j-i*s   . (1)`

In other words, Vij,k is the coefficient corresponding to degree k
in the matrix polynomial Vij(s).
Vij,k has dimensions mu(i)-by-(mu(j)-nu(j)).
The coefficients Vij,k are stored in the matrix VEPS as follows
(for the case n = 3):

    `sizes      m1-n1    m2-n2   m2-n2    m3-n3   m3-n3   m3-n3`

        `m1 { | V11,0 || V12,0 | V12,1 || V13,0 | V13,1 | V13,2 ||`
             `|       ||       |       ||       |       |       ||`
 `VEPS = m2 { |       || V22,0 |       || V23,0 | V23,1 |       ||`
             `|       ||       |       ||       |       |       ||`
        `m3 { |       ||       |       || V33,0 |       |       ||`

where mi = mu(i), ni = nu(i).
Matrix VEPS has dimensions nrv-by-ncv where
  `nrv = Sum(i=1,...,n) mu(i)`
  `ncv = Sum(i=1,...,n) i*(mu(i)-nu(i))`

==================================================================
REMARK: This routine is intended to be called only from the SLICOT
        `routine MC03ND.`
==================================================================



See the SLICOT documentation for details.
"""
function mc03ny! end

"""

Compute the QR factorization with column pivoting of an
m-by-n Jacobian matrix J (m >= n), that is, J*P = Q*R, where Q is
a matrix with orthogonal columns, P a permutation matrix, and
R an upper trapezoidal matrix with diagonal elements of
nonincreasing magnitude, and to apply the transformation Q' on
the error vector e (in-situ). The 1-norm of the scaled gradient
is also returned.

This routine is an interface to SLICOT Library routine MD03BX,
for solving standard nonlinear least squares problems using SLICOT
routine MD03BD.



See the SLICOT documentation for details.
"""
function md03ba! end

"""

Determine a value for the parameter PAR such that if x solves
the system

      `A*x = b ,     sqrt(PAR)*D*x = 0 ,`

in the least squares sense, where A is an m-by-n matrix, D is an
n-by-n nonsingular diagonal matrix, and b is an m-vector, and if
DELTA is a positive number, DXNORM is the Euclidean norm of D*x,
then either PAR is zero and

      `( DXNORM - DELTA ) .LE. 0.1*DELTA ,`

or PAR is positive and

      `ABS( DXNORM - DELTA ) .LE. 0.1*DELTA .`

It is assumed that a QR factorization, with column pivoting, of A
is available, that is, A*P = Q*R, where P is a permutation matrix,
Q has orthogonal columns, and R is an upper triangular matrix
with diagonal elements of nonincreasing magnitude.
The routine needs the full upper triangle of R, the permutation
matrix P, and the first n components of Q'*b (' denotes the
transpose). On output, MD03BB also provides an upper triangular
matrix S such that

      `P'*(A'*A + PAR*D*D)*P = S'*S .`

Matrix S is used in the solution process.

This routine is an interface to SLICOT Library routine MD03BY,
for solving standard nonlinear least squares problems using SLICOT
routine MD03BD.



See the SLICOT documentation for details.
"""
function md03bb! end

"""

This is the FCN routine for solving a standard nonlinear least
squares problem using SLICOT Library routine MD03BD. See the
parameter FCN in the routine MD03BD for the description of
parameters.

The example programmed in this routine is adapted from that
accompanying the MINPACK routine LMDER.

******************************************************************



See the SLICOT documentation for details.
"""
function md03bf! end

"""

Compute the QR factorization with column pivoting of an
m-by-n matrix J (m >= n), that is, J*P = Q*R, where Q is a matrix
with orthogonal columns, P a permutation matrix, and R an upper
trapezoidal matrix with diagonal elements of nonincreasing
magnitude, and to apply the transformation Q' on the error
vector e (in-situ). The 1-norm of the scaled gradient is also
returned. The matrix J could be the Jacobian of a nonlinear least
squares problem.



See the SLICOT documentation for details.
"""
function md03bx! end

"""

Determine a value for the parameter PAR such that if x solves
the system

      `A*x = b ,     sqrt(PAR)*D*x = 0 ,`

in the least squares sense, where A is an m-by-n matrix, D is an
n-by-n nonsingular diagonal matrix, and b is an m-vector, and if
DELTA is a positive number, DXNORM is the Euclidean norm of D*x,
then either PAR is zero and

      `( DXNORM - DELTA ) .LE. 0.1*DELTA ,`

or PAR is positive and

      `ABS( DXNORM - DELTA ) .LE. 0.1*DELTA .`

It is assumed that a QR factorization, with column pivoting, of A
is available, that is, A*P = Q*R, where P is a permutation matrix,
Q has orthogonal columns, and R is an upper triangular matrix
with diagonal elements of nonincreasing magnitude.
The routine needs the full upper triangle of R, the permutation
matrix P, and the first n components of Q'*b (' denotes the
transpose). On output, MD03BY also provides an upper triangular
matrix S such that

      `P'*(A'*A + PAR*D*D)*P = S'*S .`

Matrix S is used in the solution process.



See the SLICOT documentation for details.
"""
function md03by! end

