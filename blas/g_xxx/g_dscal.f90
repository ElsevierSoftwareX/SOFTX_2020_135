!>    @brief This procedure implements the differentiated version of the
!>    @param[in] n vector length
!>    @param[in] dx vector x
!>    @param[in] g_dx derivative object associated to dx
!>    @param[in] da scalar a = alpha
!>    @param[in] incx stride in vector x
!>    @param[out] dx vector x
!>    @param[out] g_dx derivative object associated to dx
!>    @details
!>    This procedure implements the differentiated version of the\n
!>    BLAS1 routine dscal(). Based on original BLAS1 and changed\n
!>    to a differentiated handcoded version.\n
!>    Technique:\n
!>    Changes: add support for (incx .le. 0)\n
!>    Original BLAS description:\n
!>    scales a vector by a constant.\n
!>    uses unrolled loops for increment equal to one.\n
!>    jack dongarra, linpack, 3/11/78.\n
!>    modified 3/93 to return if incx .le. 0.\n
!>    modified 12/3/93, array(1) declarations changed to array(*)\n
      SUBROUTINE g_dscal(n,da,dx,g_dx,incx)
        DOUBLE PRECISION da, dx(*), g_dx(*)
        INTEGER incx, n

!     CheCk out input:

!hAW  CorreCt inCx .le. 0
!hAW  if( n.le.0 .or. inCx.le.0 )return
        IF (n<=0) RETURN

!     Start of Computation :

!hAW  seems a little bit slower than the Code for not equal to 1
        CALL dscal(n,da,g_dx,incx)
!     Call orig. funCtion after
        CALL dscal(n,da,dx,incx)
        RETURN
      END
