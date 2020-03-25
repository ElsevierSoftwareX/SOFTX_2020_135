!>    @brief This procedure implements the differentiated version of the
!>    @param[in] n vector length
!>    @param[in] dx double vector x
!>    @param[in] g_dx derivative object associated to dx
!>    @param[in] dy double vector y
!>    @param[in] g_dy derivative object associated to dy
!>    @param[in] incx stride in vector x
!>    @param[in] incy stride in vector y
!>    @return inner product of x and y
!>    @param[out] g_dres derivative object associated to g_ddot
!>    @details
!>    This procedure implements the differentiated version of the\n
!>    BLAS1 routine ddot(). The implementation is obtained from\n
!>    copying the Adifor2.0-generated code and replacing the\n
!>    part for unit stride with a handcoded version.\n
!>    Technique: two BLAS2 calls\n
!>    Author: Martin Buecker, 07/27/00\n
!>    Changes: make generate-able for Precisions.\n
!>             add support for ldg_dx, ldg_dy, g_p.\n
!>             make more readable.\n
!>             make call of orig. function.\n
!>             bug-workaround for n.le.5\n
!>    Handcoded differentiated version using two BLAS2 calls of\n
!>    the original function\n
!>       version 1.1    date   22-08-91 (Zdenek: no dimension, cdc, ce)\n
!>    returns the dot product of double precision dx and dy.\n
!>    ddot = sum for i = 0 to n-1 of  dx(lx+i*incx) * dy(ly+i*incy)\n
!>    where lx = 1 if incx .ge. 0, else lx = (-incx)*n, and ly is\n
!>    defined in a similar way using incy.\n
      FUNCTION g_ddot(n,dx,g_dx,incx,dy,g_dy,incy,g_dres)
        DOUBLE PRECISION dx(*), dy(*)
!hpu      on real underflow ignore
        DOUBLE PRECISION g_ddot, g_dres, g_dx(*), g_dy(*)

#ifdef USE_QDDOT
        EXTERNAL qddot
        DOUBLE PRECISION qddot
#else
        EXTERNAL ddot
        DOUBLE PRECISION ddot
#endif

!hMB    deClarations added by Martin
        DOUBLE PRECISION one, zero
        PARAMETER (one=1.0D0,zero=0.0D0)
        INTEGER incx, incy, n

!     CheCk out input:

        g_ddot = 0.0D0
        IF (n<=0) RETURN

!     Start of Computation :
!hAW  Call orig. funCtion for Computation the org. value
#ifdef USE_QDDOT
        g_ddot = qddot(n,dx,incx,dy,incy)
#else
        g_ddot = ddot(n,dx,incx,dy,incy)
#endif

!     g_dres = g_dy*dx
#ifdef USE_QDDOT
        g_dres = qddot(n,g_dy,incy,dx,incx)
#else
        g_dres = ddot(n,g_dy,incy,dx,incx)
#endif

!     g_dres = g_dres + g_dx*dy
#ifdef USE_QDDOT
        g_dres = g_dres + qddot(n,g_dx,incx,dy,incy)
#else
        g_dres = g_dres + ddot(n,g_dx,incx,dy,incy)
#endif

        RETURN
      END
