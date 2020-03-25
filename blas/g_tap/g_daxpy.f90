!>    @brief This procedure implements the differentiated version of the
!>    @param[in] n vector length
!>    @param[in] dx vector x
!>    @param[in] g_dx derivative object associated to dx
!>    @param[in] dy vector y
!>    @param[in] g_dy derivative object associated to dy
!>    @param[in] da scalar a = alpha
!>    @param[in] incx stride in vector x
!>    @param[in] incy stride in vector y
!>    @param[out] dy vector y
!>    @param[out] g_dy derivative object associated to dy
!>    @details
!>    This procedure implements the differentiated version of the\n
!>    BLAS1 routine daxpy(). Based on original BLAS1 and changed\n
!>    to a differentiated handcoded version.\n
!>    Technique:\n
!>    Changes: use of blocking technics (bl_size=blocksize)\n
!>    Original BLAS description:\n
!>    constant times a vector plus a vector.\n
!>    uses unrolled loops for increments equal to one.\n
!>    jack dongarra, linpack, 3/11/78.\n
!>    modified 12/3/93, array(1) declarations changed to array(*)\n
      SUBROUTINE g_daxpy(n,da,g_da,dx,g_dx,incx,dy,g_dy,incy)
        DOUBLE PRECISION da, g_da, dx(*), dy(*), g_dx(*), g_dy(*)
        INTEGER incx, incy, n

!     check out input:
        IF (n<=0) RETURN

!     Start of computation :
        IF (abs(g_da)<=1.0d-300) THEN
!         for more spec. performance code
          CALL daxpy(n,da,g_dx,incx,g_dy,incy)
        ELSE
          CALL daxpy(n,da,g_dx,incx,g_dy,incy)
          CALL daxpy(n,g_da,dx,incx,g_dy,incy)
        END IF
!       call orig. function for none-differentiated computation
        CALL daxpy(n,da,dx,incx,dy,incy)

        RETURN
      END
