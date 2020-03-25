!>    @brief performs one of the matrix-vector operations
!>    @param[in] trans - CHARACTER*1.
!>    @details
!>    On entry, TRANS specifies the operation to be performed as\n
!>    follows:\n
!>    = 'N' or 'n'   y := alpha*A*x + beta*y.\n
!>    = 'T' or 't'   y := alpha*A'*x + beta*y.\n
!>    = 'C' or 'c'   y := alpha*A'*x + beta*y.\n
!>    @param[in] m - INTEGER.\n
!>    On entry, M specifies the number of rows of the matrix A.\n
!>    Must be at least zero.\n
!>    @param[in] n - INTEGER.\n
!>    On entry, N specifies the number of columns of the matrix A.\n
!>    Must be at least zero.\n
!>    @param[in] alpha - DOUBLE PRECISION.\n
!>    On entry, ALPHA specifies the scalar alpha.\n
!>    @param[in] a - DOUBLE PRECISION array of DIMENSION ( LDA, n ).\n
!>    Before entry, the leading m by n part of the array A must\n
!>    contain the matrix of coefficients.\n
!>    @param[in] lda - INTEGER.\n
!>    On entry, LDA specifies the first dimension of A as declared\n
!>    in the calling (sub) program. LDA must be at least\n
!>    max( 1, m ).\n
!>    @param[in] x - DOUBLE PRECISION array of DIMENSION at least\n
!>    ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'\n
!>    and at least\n
!>    ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.\n
!>    Before entry, the incremented array X must contain the\n
!>    vector x.\n
!>    @param[in] incx - INTEGER.\n
!>    On entry, INCX specifies the increment for the elements of X.\n
!>    Must not be zero.\n
!>    @param[in] beta - DOUBLE PRECISION.\n
!>    On entry, BETA specifies the scalar beta. When BETA is\n
!>    supplied as zero then Y need not be set on input.\n
!>    @param[in] y - DOUBLE PRECISION array of DIMENSION at least\n
!>    ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'\n
!>    and at least\n
!>    ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.\n
!>    Before entry with BETA non-zero, the incremented array Y\n
!>    must contain the vector y. On exit, Y is overwritten by the\n
!>    updated vector y.\n
!>    @param[in] incy - INTEGER.\n
!>    On entry, INCY specifies the increment for the elements of Y.\n
!>    Must not be zero.\n
!>    @param[out] y - DOUBLE PRECISION array of DIMENSION at least\n
!>    ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'\n
!>    and at least\n
!>    ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.\n
!>    Before entry with BETA non-zero, the incremented array Y\n
!>    must contain the vector y. On exit, Y is overwritten by the\n
!>    updated vector y.\n
!>    @details\n
!> g_DGEMV  performs one of the matrix-vector operations\n
!>       g_y := alpha*g_A*x + alpha*A*g_x + beta*g_y,\n
!>  or   g_y := alpha*g_A'*x + alpha*A'*g_x + beta*g_y,\n
!> where alpha and beta are scalars, x, g_x and y, g_y are vectors \n
!> and A, g_A is an m by n matrix.\n
!>    This procedure implements a special differentiated version of the\n
!>    BLAS2 routine dgemv(). Based on original BLAS2 and changed\n
!>    to a differentiated handcoded version.\n
!>    Technique:\n
!> Level 2 Blas routine.\n
!> -- Written on 22-October-1986.\n
!>    Jack Dongarra, Argonne National Lab.\n
!>    Jeremy Du Croz, Nag Central Office.\n
!>    Sven Hammarling, Nag Central Office.\n
!>    Richard Hanson, Sandia National Labs.\n
!>    Original BLAS description:\n
!>    .. Scalar Arguments ..\n
      SUBROUTINE g_dgemv(trans,m,n,alpha,a,g_a,lda,x,g_x,incx,beta,y, &
          g_y,incy)
        DOUBLE PRECISION alpha, beta
        INTEGER incx, incy, lda, m, n
        CHARACTER*1 trans
!     .. Array Arguments ..
        DOUBLE PRECISION a(lda,*), x(*), y(*)
!     ..



!     .. Executable Statements ..

!     Test the input parameters.

        DOUBLE PRECISION g_y(*), g_x(*), g_a(lda,*)


!     g_y := alpha*g_A(')*x + alpha*A(')*g_x + beta*g_y :

!     g_y_t := alpha*A(')*g_x + beta*g_y
        CALL dgemv(trans,m,n,alpha,a,lda,g_x,incx,beta,g_y,incy)

!     g_y := alpha*g_A(')*x + g_y_t
        CALL dgemv(trans,m,n,alpha,g_a,lda,x,incx,1.0D0,g_y,incy)


!     original BLAS2 call for evaluating the original function values
        CALL dgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)

        RETURN

!     End of g_DGEMV .

      END
