/* -*- c++ -*- */
#ifndef LAPACKPROTOMOL_H
#define LAPACKPROTOMOL_H

//Lapack/Blas Fortran routines for Normal Mode routines
//Includes:
//dgemv   : y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y, A mxn matrix.
//dlamch  : determines double precision machine parameters. 
//dsyevr  : computes selected eigenvalues and eigenvectors of a real symmetric matrix A.
//dgemm   : C := alpha*op( A )*op( B ) + beta*C.
//ddot    : dodt product of vectors
//dnrm2   : euclidean norm of vector.

/******************************************************************************************************************************************************************
*      SUBROUTINE DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
*     .. Scalar Arguments .. DOUBLE PRECISION ALPHA,BETA, INTEGER INCX,INCY,LDA,M,N, CHARACTER TRANS
*     .. Array Arguments .. DOUBLE PRECISION A(LDA,*),X(*),Y(*)
*  Purpose
*  =======
*  DGEMV  performs one of the matrix-vector operations
*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y, where alpha and beta are scalars, x and y are vectors and A is an m by n matrix.
*  Arguments
*  ==========
*  TRANS  - CHARACTER*1. On entry, TRANS specifies the operation to be performed as follows:
*              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
*              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
*              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y. Unchanged on exit.
*  M      - INTEGER. On entry, M specifies the number of rows of the matrix A. M must be at least zero. Unchanged on exit.
*  N      - INTEGER. On entry, N specifies the number of columns of the matrix A. N must be at least zero. Unchanged on exit.
*  ALPHA  - DOUBLE PRECISION. On entry, ALPHA specifies the scalar alpha. Unchanged on exit.
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ). Before entry, the leading m by n part of the array A must contain the matrix of coefficients. Unchanged on exit.
*  LDA    - INTEGER. On entry, LDA specifies the first dimension of A as declared in the calling (sub) program. LDA must be at least max( 1, m ). Unchanged on exit.
*  X      - DOUBLE PRECISION array of DIMENSION at least ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n' and at least ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
*           Before entry, the incremented array X must contain the vector x. Unchanged on exit.
*  INCX   - INTEGER. On entry, INCX specifies the increment for the elements of X. INCX must not be zero. Unchanged on exit.
*  BETA   - DOUBLE PRECISION. On entry, BETA specifies the scalar beta. When BETA is supplied as zero then Y need not be set on input. Unchanged on exit.
*  Y      - DOUBLE PRECISION array of DIMENSION at least ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n' and at least ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
*           Before entry with BETA non-zero, the incremented array Y must contain the vector y. On exit, Y is overwritten by the updated vector y.
*  INCY   - INTEGER. On entry, INCY specifies the increment for the elements of Y. INCY must not be zero. Unchanged on exit.
*/
extern "C" void dgemv_(char *transA, int *m, int *n, double *alpha, double *A, int *lda, double *x, int *incx, 
                       double *beta, double *Y, int *incY);

/******************************************************************************************************************************************************************
*       SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
*
*     .. Scalar Arguments ..CHARACTER   JOBZ, UPLO,INTEGER  INFO, LDA, LWORK, N
*     .. Array Arguments  ..DOUBLE PRECISION   A( LDA, * ), W( * ), WORK( * )
*
*  Purpose
*  =======
*
*  DSYEV computes all eigenvalues and, optionally, eigenvectors of a
*  real symmetric matrix A.
*
*  Arguments
*  =========
*
*  JOBZ    (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only;
*          = 'V':  Compute eigenvalues and eigenvectors.
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
*          On entry, the symmetric matrix A.  If UPLO = 'U', the
*          leading N-by-N upper triangular part of A contains the
*          upper triangular part of the matrix A.  If UPLO = 'L',
*          the leading N-by-N lower triangular part of A contains
*          the lower triangular part of the matrix A.
*          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
*          orthonormal eigenvectors of the matrix A.
*          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
*          or the upper triangle (if UPLO='U') of A, including the
*          diagonal, is destroyed.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  W       (output) DOUBLE PRECISION array, dimension (N)
*          If INFO = 0, the eigenvalues in ascending order.
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The length of the array WORK.  LWORK >= max(1,3*N-1).
*          For optimal efficiency, LWORK >= (NB+2)*N,
*          where NB is the blocksize for DSYTRD returned by ILAENV.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, the algorithm failed to converge; i
*                off-diagonal elements of an intermediate tridiagonal
*                form did not converge to zero.
*/
extern "C" void dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda,  
                       double *w, double *work, int *lwork, int *info);


/******************************************************************************************************************************************************************
*      SUBROUTINE DSYEVR( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, IWORK, LIWORK, INFO )
*  -- LAPACK driver routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*     .. Scalar Arguments .. CHARACTER JOBZ, RANGE, UPLO, INTEGER IL, INFO, IU, LDA, LDZ, LIWORK, LWORK, M, N, DOUBLE PRECISION   ABSTOL, VL, VU
*     .. Array Arguments .. INTEGER ISUPPZ( * ), IWORK( * ), DOUBLE PRECISION   A( LDA, * ), W( * ), WORK( * ), Z( LDZ, * )
*  Purpose
*  =======
*  DSYEVR computes selected eigenvalues and, optionally, eigenvectors of a real symmetric matrix A.  Eigenvalues and eigenvectors can be selected by specifying either a range of values or a range of
*  indices for the desired eigenvalues. DSYEVR first reduces the matrix A to tridiagonal form T with a call to DSYTRD.  Then, whenever possible, DSYEVR calls DSTEMR to compute
*  the eigenspectrum using Relatively Robust Representations.  DSTEMR computes eigenvalues by the dqds algorithm, while orthogonal eigenvectors are computed from various "good" L D L^T representations
*  (also known as Relatively Robust Representations). Gram-Schmidt orthogonalization is avoided as far as possible. More specifically, the various steps of the algorithm are as follows.
*  For each unreduced block (submatrix) of T,
*     (a) Compute T - sigma I  = L D L^T, so that L and D define all the wanted eigenvalues to high relative accuracy. This means that small relative changes in the entries of D and L
*         cause only small relative changes in the eigenvalues and eigenvectors. The standard (unfactored) representation of the tridiagonal matrix T does not have this property in general.
*     (b) Compute the eigenvalues to suitable accuracy. If the eigenvectors are desired, the algorithm attains full accuracy of the computed eigenvalues only right before
*         the corresponding vectors have to be computed, see steps c) and d).
*     (c) For each cluster of close eigenvalues, select a new shift close to the cluster, find a new factorization, and refine the shifted eigenvalues to suitable accuracy.
*     (d) For each eigenvalue with a large enough relative separation compute the corresponding eigenvector by forming a rank revealing twisted factorization. Go back to (c) for any clusters that remain.
*  The desired accuracy of the output can be specified by the input parameter ABSTOL.
*  For more details, see DSTEMR's documentation.
*  Note 1 : DSYEVR calls DSTEMR when the full spectrum is requested on machines which conform to the ieee-754 floating point standard.
*  DSYEVR calls DSTEBZ and SSTEIN on non-ieee machines and when partial spectrum requests are made.
*  Normal execution of DSTEMR may create NaNs and infinities and hence may abort due to a floating point exception in environments
*  which do not handle NaNs and infinities in the ieee standard default manner.
*  Arguments
*  =========
*  JOBZ    (input) CHARACTER*1 = 'N':  Compute eigenvalues only; = 'V':  Compute eigenvalues and eigenvectors.
*  RANGE   (input) CHARACTER*1 = 'A': all eigenvalues will be found. = 'V': all eigenvalues in the half-open interval (VL,VU] will be found. = 'I': the IL-th through IU-th eigenvalues will be found.
********** For RANGE = 'V' or 'I' and IU - IL < N - 1, DSTEBZ and DSTEIN are called
*  UPLO    (input) CHARACTER*1 = 'U':  Upper triangle of A is stored; = 'L':  Lower triangle of A is stored.
*  N       (input) INTEGER The order of the matrix A.  N >= 0.
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N) On entry, the symmetric matrix A.  If UPLO = 'U', the leading N-by-N upper triangular part of A contains the
*          upper triangular part of the matrix A.  If UPLO = 'L', the leading N-by-N lower triangular part of A contains the lower triangular part of the matrix A.
*          On exit, the lower triangle (if UPLO='L') or the upper triangle (if UPLO='U') of A, including the diagonal, is destroyed.
*  LDA     (input) INTEGER The leading dimension of the array A.  LDA >= max(1,N).
*  VL      (input) DOUBLE PRECISION
*  VU      (input) DOUBLE PRECISION: If RANGE='V', the lower and upper bounds of the interval to be searched for eigenvalues. VL < VU. Not referenced if RANGE = 'A' or 'I'.
*  IL      (input) INTEGER
*  IU      (input) INTEGER: If RANGE='I', the indices (in ascending order) of the smallest and largest eigenvalues to be returned. 1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
*          Not referenced if RANGE = 'A' or 'V'.
*  ABSTOL  (input) DOUBLE PRECISION The absolute error tolerance for the eigenvalues. An approximate eigenvalue is accepted as converged when it is determined to lie in an interval [a,b]
*          of width less than or equal to ABSTOL + EPS *   max( |a|,|b| ) , where EPS is the machine precision.  If ABSTOL is less than or equal to zero, then  EPS*|T|  will be used in its place,
*          where |T| is the 1-norm of the tridiagonal matrix obtained by reducing A to tridiagonal form.
*          If high relative accuracy is important, set ABSTOL to DLAMCH( 'Safe minimum' ).  Doing so will guarantee that eigenvalues are computed to high relative accuracy when
*          possible in future releases.  The current code does not make any guarantees about high relative accuracy, but future releases will.
*  M       (output) INTEGER The total number of eigenvalues found.  0 <= M <= N. If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.
*  W       (output) DOUBLE PRECISION array, dimension (N) The first M elements contain the selected eigenvalues in ascending order.
*  Z       (output) DOUBLE PRECISION array, dimension (LDZ, max(1,M)) If JOBZ = 'V', then if INFO = 0, the first M columns of Z contain the orthonormal eigenvectors of the matrix A
*          corresponding to the selected eigenvalues, with the i-th column of Z holding the eigenvector associated with W(i). If JOBZ = 'N', then Z is not referenced.
*          Note: the user must ensure that at least max(1,M) columns are supplied in the array Z; if RANGE = 'V', the exact value of M is not known in advance and an upper bound must be used.
*          Supplying N columns is always safe.
*  LDZ     (input) INTEGER The leading dimension of the array Z.  LDZ >= 1, and if JOBZ = 'V', LDZ >= max(1,N).
*  ISUPPZ  (output) INTEGER array, dimension ( 2*max(1,M) ) The support of the eigenvectors in Z, i.e., the indices indicating the nonzero elements in Z. The i-th eigenvector
*          is nonzero only in elements ISUPPZ( 2*i-1 ) through ISUPPZ( 2*i ).
********** Implemented only for RANGE = 'A' or 'I' and IU - IL = N - 1
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK)) On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*  LWORK   (input) INTEGER The dimension of the array WORK.  LWORK >= max(1,26*N). For optimal efficiency, LWORK >= (NB+6)*N, where NB is the max of the blocksize for DSYTRD and DORMTR
*          returned by ILAENV. If LWORK = -1, then a workspace query is assumed; the routine only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error message related to LWORK is issued by XERBLA.
*  IWORK   (workspace/output) INTEGER array, dimension (MAX(1,LIWORK)) On exit, if INFO = 0, IWORK(1) returns the optimal LWORK.
*  LIWORK  (input) INTEGER The dimension of the array IWORK.  LIWORK >= max(1,10*N). If LIWORK = -1, then a workspace query is assumed; the
*          routine only calculates the optimal size of the IWORK array, returns this value as the first entry of the IWORK array, and no error message related to LIWORK is issued by XERBLA.
*  INFO    (output) INTEGER = 0:  successful exit < 0:  if INFO = -i, the i-th argument had an illegal value > 0:  Internal error
*/
extern "C" void dsyevr_(char *jobz, char *range, char *uplo, int *n, double *a, int *lda, 
                    double *vl, double *vu, int *il, int *iu, double *abstol, int *m, double *w, double *z, 
                    int *ldz, int *isuppz, double *work, int *lwork, int *iwork, int *liwork, int *info);

/******************************************************************************************************************************************************************
*      DOUBLE PRECISION FUNCTION DLAMCH( CMACH )
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*     .. Scalar Arguments .. CHARACTER          CMACH
*  Purpose
*  =======
*  DLAMCH determines double precision machine parameters.
*  Arguments
*  =========
*  CMACH   (input) CHARACTER*1 Specifies the value to be returned by DLAMCH:
*          = 'E' or 'e',   DLAMCH := eps;  = 'S' or 's ,   DLAMCH := sfmin;
*          = 'B' or 'b',   DLAMCH := base; = 'P' or 'p',   DLAMCH := eps*base;
*          = 'N' or 'n',   DLAMCH := t;    = 'R' or 'r',   DLAMCH := rnd;
*          = 'M' or 'm',   DLAMCH := emin; = 'U' or 'u',   DLAMCH := rmin;
*          = 'L' or 'l',   DLAMCH := emax; = 'O' or 'o',   DLAMCH := rmax;
*          where
*          eps   = relative machine precision; sfmin = safe minimum, such that 1/sfmin does not overflow; base  = base of the machine;
*          prec  = eps*base; t = number of (base) digits in the mantissa; rnd = 1.0 when rounding occurs in addition, 0.0 otherwise;
*          emin  = minimum exponent before (gradual) underflow; rmin  = underflow threshold - base**(emin-1);
*          emax  = largest exponent before overflow; rmax  = overflow threshold  - (base**emax)*(1-eps)
*/
extern "C" double  dlamch_(char *cmach);

/***************************************************************************************
*     SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
*     .. Scalar Arguments .. DOUBLE PRECISION ALPHA,BETA, INTEGER K,LDA,LDB,LDC,M,N, CHARACTER TRANSA,TRANSB
*     .. Array Arguments .. DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
*  Purpose
*  =======
*  DGEMM  performs one of the matrix-matrix operations C := alpha*op( A )*op( B ) + beta*C, where  op( X ) is one of op( X ) = X   or   op( X ) = X',
*  alpha and beta are scalars, and A, B and C are matrices, with op( A ) an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
*  Arguments
*  ==========
*  TRANSA - CHARACTER*1. On entry, TRANSA specifies the form of op( A ) to be used in the matrix multiplication as follows:
*              TRANSA = 'N' or 'n',  op( A ) = A.
*              TRANSA = 'T' or 't',  op( A ) = A'.
*              TRANSA = 'C' or 'c',  op( A ) = A'. Unchanged on exit.
*  TRANSB - CHARACTER*1. On entry, TRANSB specifies the form of op( B ) to be used in the matrix multiplication as follows:
*              TRANSB = 'N' or 'n',  op( B ) = B.
*              TRANSB = 'T' or 't',  op( B ) = B'.
*              TRANSB = 'C' or 'c',  op( B ) = B'. Unchanged on exit.
*  M      - INTEGER. On entry,  M  specifies  the number  of rows  of the  matrix op( A )  and of the  matrix  C.  M  must  be at least  zero. Unchanged on exit.
*  N      - INTEGER. On entry,  N  specifies the number  of columns of the matrix op( B ) and the number of columns of the matrix C. N must be at least zero. Unchanged on exit.
*  K      - INTEGER. On entry,  K  specifies  the number of columns of the matrix op( A ) and the number of rows of the matrix op( B ). K must be at least  zero. Unchanged on exit.
*  ALPHA  - DOUBLE PRECISION. On entry, ALPHA specifies the scalar alpha. Unchanged on exit.
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is k  when  TRANSA = 'N' or 'n',  and is  m  otherwise. Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
*           part of the array  A  must contain the matrix  A,  otherwise the leading  k by m  part of the array  A  must contain  the matrix A. Unchanged on exit.
*  LDA    - INTEGER. On entry, LDA specifies the first dimension of A as declared in the calling (sub) program. When  TRANSA = 'N' or 'n' then
*           LDA must be at least  max( 1, m ), otherwise  LDA must be at least  max( 1, k ). Unchanged on exit.
*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is n  when  TRANSB = 'N' or 'n',  and is  k  otherwise. Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
*           part of the array  B  must contain the matrix  B,  otherwise the leading  n by k  part of the array  B  must contain  the matrix B. Unchanged on exit.
*  LDB    - INTEGER. On entry, LDB specifies the first dimension of B as declared in the calling (sub) program. When  TRANSB = 'N' or 'n' then
*           LDB must be at least  max( 1, k ), otherwise  LDB must be at least  max( 1, n ). Unchanged on exit.
*  BETA   - DOUBLE PRECISION. On entry,  BETA  specifies the scalar  beta.  When  BETA  is supplied as zero then C need not be set on input. Unchanged on exit.
*  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ). Before entry, the leading  m by n  part of the array  C must contain the matrix  C,  except when  beta  is zero, in which
*           case C need not be set on entry. On exit, the array  C  is overwritten by the  m by n  matrix ( alpha*op( A )*op( B ) + beta*C ).
*  LDC    - INTEGER. On entry, LDC specifies the first dimension of C as declared in  the  calling  (sub)  program.   LDC  must  be  at  least max( 1, m ). Unchanged on exit.
*/
extern "C" void dgemm_ (char *transA, char *transB, int *m, int *n, int *k, double *alpha, 
                            double *A, int *lda, double *B, int *ldb, double *beta, double *C, int *l);

/***************************************************************************************
      DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
*     .. Scalar Arguments  INTEGER INCX,INCY,N
*/
extern "C" double ddot_(int *n, double *x, int *incx, double *y, int *incy);

/***************************************************************************************
      DOUBLE PRECISION FUNCTION DNRM2(N,X,INCX)
*     .. Scalar Arguments .. INTEGER INCX,N
*/
extern "C" double dnrm2_(int *n, double *x, int *incx);

/***************************************************************************************
      SUBROUTINE DPOTRI( UPLO, N, A, LDA, INFO )
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
*  Purpose
*  =======
*  DPOTRI computes the inverse of a real symmetric positive definite
*  matrix A using the Cholesky factorization A = U**T*U or A = L*L**T
*  computed by DPOTRF.
*  Arguments
*  =========
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the triangular factor U or L from the Cholesky
*          factorization A = U**T*U or A = L*L**T, as computed by
*          DPOTRF.
*          On exit, the upper or lower triangle of the (symmetric)
*          inverse of A, overwriting the input factor U or L.
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, the (i,i) element of the factor U or L is
*                zero, and the inverse could not be computed.
*/
extern "C" void dpotri_(char *transA, int *n, double *A, int *lda, int *info);

/***************************************************************************************
     SUBROUTINE DPOTRF( UPLO, N, A, LDA, INFO )
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
*     ..
*  Purpose
*  =======
*  DPOTRF computes the Cholesky factorization of a real symmetric
*  positive definite matrix A.
*  The factorization has the form
*     A = U**T * U,  if UPLO = 'U', or
*     A = L  * L**T,  if UPLO = 'L',
*  where U is an upper triangular matrix and L is lower triangular.
*  This is the block version of the algorithm, calling Level 3 BLAS.
*  Arguments
*  =========
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
*          N-by-N upper triangular part of A contains the upper
*          triangular part of the matrix A, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading N-by-N lower triangular part of A contains the lower
*          triangular part of the matrix A, and the strictly upper
*          triangular part of A is not referenced.
*
*          On exit, if INFO = 0, the factor U or L from the Cholesky
*          factorization A = U**T*U or A = L*L**T.
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, the leading minor of order i is not
*                positive definite, and the factorization could not be
*                completed.
*/
extern "C" void dpotrf_(char *transA, int *n, double *A, int *lda, int *info);

/***************************************************************************************
      SUBROUTINE DPOSV( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
*  -- LAPACK driver routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*  Purpose
*  =======
*  DPOSV computes the solution to a real system of linear equations
*     A * X = B,
*  where A is an N-by-N symmetric positive definite matrix and X and B
*  are N-by-NRHS matrices.
*  The Cholesky decomposition is used to factor A as
*     A = U**T* U,  if UPLO = 'U', or
*     A = L * L**T,  if UPLO = 'L',
*  where U is an upper triangular matrix and L is a lower triangular
*  matrix.  The factored form of A is then used to solve the system of
*  equations A * X = B.
*  Arguments
*  =========
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*  N       (input) INTEGER
*          The number of linear equations, i.e., the order of the
*          matrix A.  N >= 0.
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
*          N-by-N upper triangular part of A contains the upper
*          triangular part of the matrix A, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading N-by-N lower triangular part of A contains the lower
*          triangular part of the matrix A, and the strictly upper
*          triangular part of A is not referenced.
*
*          On exit, if INFO = 0, the factor U or L from the Cholesky
*          factorization A = U**T*U or A = L*L**T.
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
*          On entry, the N-by-NRHS right hand side matrix B.
*          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, the leading minor of order i of A is not
*                positive definite, so the factorization could not be
*                completed, and the solution has not been computed.
*/
extern "C" void dposv_(char *transA, int *n, int *nrhs, double *a, int *lda, double *b, int *ldb,int *info);

/***************************************************************************************
      SUBROUTINE DTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA
      INTEGER LDA,LDB,M,N
      CHARACTER DIAG,SIDE,TRANSA,UPLO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),B(LDB,*)
*     ..
*  Purpose
*  =======
*  DTRMM  performs one of the matrix-matrix operations
*     B := alpha*op( A )*B,   or   B := alpha*B*op( A ),
*  where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
*  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
*     op( A ) = A   or   op( A ) = A'.
*  Arguments
*  ==========
*  SIDE   - CHARACTER*1.
*           On entry,  SIDE specifies whether  op( A ) multiplies B from
*           the left or right as follows:
*              SIDE = 'L' or 'l'   B := alpha*op( A )*B.
*              SIDE = 'R' or 'r'   B := alpha*B*op( A ).
*           Unchanged on exit.
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix A is an upper or
*           lower triangular matrix as follows:
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*           Unchanged on exit.
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*              TRANSA = 'N' or 'n'   op( A ) = A.
*              TRANSA = 'T' or 't'   op( A ) = A'.
*              TRANSA = 'C' or 'c'   op( A ) = A'.
*           Unchanged on exit.
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit triangular
*           as follows:
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*           Unchanged on exit.
*  M      - INTEGER.
*           On entry, M specifies the number of rows of B. M must be at
*           least zero.
*           Unchanged on exit.
*  N      - INTEGER.
*           On entry, N specifies the number of columns of B.  N must be
*           at least zero.
*           Unchanged on exit.
*  ALPHA  - DOUBLE PRECISION.
*           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
*           zero then  A is not referenced and  B need not be set before
*           entry.
*           Unchanged on exit.
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m
*           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
*           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
*           upper triangular part of the array  A must contain the upper
*           triangular matrix  and the strictly lower triangular part of
*           A is not referenced.
*           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
*           lower triangular part of the array  A must contain the lower
*           triangular matrix  and the strictly upper triangular part of
*           A is not referenced.
*           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
*           A  are not referenced either,  but are assumed to be  unity.
*           Unchanged on exit.
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
*           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
*           then LDA must be at least max( 1, n ).
*           Unchanged on exit.
*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
*           Before entry,  the leading  m by n part of the array  B must
*           contain the matrix  B,  and  on exit  is overwritten  by the
*           transformed matrix.
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in  the  calling  (sub)  program.   LDB  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*/
extern "C" void dtrmm_(char *sideA, char *ulA, char *transA, char *diagA, int *m, int *n, double *alpha, double *A, int *lda, double *B, int *ldb);

/***************************************************************************************
*	  SUBROUTINE DTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
*     .. Scalar Arguments DOUBLE PRECISION ALPHA    INTEGER LDA,LDB,M,N    CHARACTER DIAG,SIDE,TRANSA,UPLO
*     .. Array Arguments  DOUBLE PRECISION A(LDA,*),B(LDB,*)
*  Purpose
*  =======
*  DTRSM  solves one of the matrix equations
*     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
*  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
*  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
*     op( A ) = A   or   op( A ) = A'.
*  The matrix X is overwritten on B.
*  Arguments
*  ==========
*  SIDE   - CHARACTER*1.
*           On entry, SIDE specifies whether op( A ) appears on the left
*           or right of X as follows:
*              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
*              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
*           Unchanged on exit.
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix A is an upper or
*           lower triangular matrix as follows:
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*           Unchanged on exit.
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*              TRANSA = 'N' or 'n'   op( A ) = A.
*              TRANSA = 'T' or 't'   op( A ) = A'.
*              TRANSA = 'C' or 'c'   op( A ) = A'.
*           Unchanged on exit.
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit triangular
*           as follows:
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*           Unchanged on exit.
*  M      - INTEGER.
*           On entry, M specifies the number of rows of B. M must be at
*           least zero.
*           Unchanged on exit.
*  N      - INTEGER.
*           On entry, N specifies the number of columns of B.  N must be
*           at least zero.
*           Unchanged on exit.
*  ALPHA  - DOUBLE PRECISION.
*           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
*           zero then  A is not referenced and  B need not be set before
*           entry.
*           Unchanged on exit.
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m
*           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
*           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
*           upper triangular part of the array  A must contain the upper
*           triangular matrix  and the strictly lower triangular part of
*           A is not referenced.
*           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
*           lower triangular part of the array  A must contain the lower
*           triangular matrix  and the strictly upper triangular part of
*           A is not referenced.
*           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
*           A  are not referenced either,  but are assumed to be  unity.
*           Unchanged on exit.
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
*           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
*           then LDA must be at least max( 1, n ).
*           Unchanged on exit.
*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
*           Before entry,  the leading  m by n part of the array  B must
*           contain  the  right-hand  side  matrix  B,  and  on exit  is
*           overwritten by the solution matrix  X.
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in  the  calling  (sub)  program.   LDB  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*/
extern "C" void dtrsm_(char *sideA, char *ulA, char *transA, char *diagA, int *m, int *n, double *alpha, double *A, int *lda, double *B, int *ldb);


#endif /* LAPACKPROTOMOL_H */
