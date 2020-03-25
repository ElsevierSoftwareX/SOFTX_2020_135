!>    @brief This procedure implements the differentiated version of the
!>    @param[in] n vector length
!>    @param[in] dx vector x
!>    @param[in] g_dx derivative object associated to dx
!>    @param[in] incx stride in vector x
!>    @param[in] incy stride in vector y
!>    @param[out] dy vector y
!>    @param[out] g_dy derivative object associated to dy
!>    @details
!>    This procedure implements the differentiated version of the\n
!>    BLAS1 routine dcopy(). Based on original BLAS1 and changed\n
!>    to a differentiated handcoded version.\n
!>    Technique:\n
!>    Changes:\n
!>    Original BLAS description:\n
!>    copies a vector, x, to a vector, y.\n
!>    uses unrolled loops for increments equal to one.\n
!>    jack dongarra, linpack, 3/11/78.\n
!>    modified 12/3/93, array(1) declarations changed to array(*)\n
      SUBROUTINE g_dcopy(n,dx,g_dx,incx,dy,g_dy,incy)
        DOUBLE PRECISION dx(*), dy(*)
        INTEGER incx, incy, n
        DOUBLE PRECISION g_dy(*), g_dx(*)

!     Check out input:
        IF (n<=0) RETURN

!     Start of computation :
!       call orig. dcopy to copy dx to dy
        CALL dcopy(n,dx,incx,dy,incy)
        CALL dcopy(n,g_dx,incx,g_dy,incy)

        RETURN
      END
