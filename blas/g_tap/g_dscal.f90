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
      SUBROUTINE g_dscal(n,da,g_da,dx,g_dx,incx)
        DOUBLE PRECISION da, g_da, dx(*), g_dx(*)
        INTEGER incx, n
        INTRINSIC abs

!     check out input:
!hAW  correct incx .le. 0
!hAW  if( n.le.0 .or. incx.le.0 )return
        IF (n<=0) RETURN

!     Start of computation :
        IF (abs(g_da)<=1.0d-300) THEN
!hAW      seems a little bit slower than the code for not equal to 1
          CALL dscal(n,da,g_dx,incx)
        ELSE
          CALL dscal(n,da,g_dx,incx)
          CALL daxpy(n,g_da,dx,incx,g_dx,incx)
        END IF
!       call orig. function after
        CALL dscal(n,da,dx,incx)

        RETURN
      END
