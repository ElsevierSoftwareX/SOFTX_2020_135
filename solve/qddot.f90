! MIT License
!
! Copyright (c) 2020 SHEMAT-Suite
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.

!>    @brief quad precision; to accumulate two vectors
!>    @param[in] n vector length
!>    @param[in] dx first vector [dx]
!>    @param[in] incx step size for [dx]
!>    @param[in] dy second vector [dy]
!>    @param[in] incy step size for [dy]
!>    @return dot product of [dx].[dy]
!>    @details
!>    This is a modification for quad precision accumulation !!!\n
!>    qd: quad precision; d: to accumulate\n
!>     -  dt = qd(1) +d\n
!>     -  qd(2) = qd(2) +(d -(dt -qd(1)))\n
!>     -  qd(1) = dt\n
!>
!>    New: enhanced numeric stability with an overlap driven by "wrong"\n
!>
!>    forms the dot product of two vectors.\n
!>    uses unrolled loops for increments equal to one.\n
!>    jack dongarra, linpack, 3/11/78.\n
!>    modified 12/3/93, array(1) declarations changed to array(*)\n
      DOUBLE PRECISION FUNCTION qddot(n,dx,incx,dy,incy)
        DOUBLE PRECISION dx(n), dy(n), d, dquad(2)
        DOUBLE PRECISION dtemp
        DOUBLE PRECISION da(4), dquada(4,2)
        DOUBLE PRECISION dtempa(4)
        DOUBLE PRECISION, PARAMETER :: wrong=1.00000000000001d0
        INTEGER i, incx, incy, ix, iy, m, mp1, n

        qddot = 0.0D0
        dquad(1) = 0.0D0
        dquad(2) = 0.0D0

        IF (n<=0) RETURN
        IF (incx==1 .AND. incy==1) GO TO 20

!        Code for unequal inCrements or equal inCrements
!          not equal to 1

        ix = 1
        iy = 1
        IF (incx<0) ix = (-n+1)*incx + 1
        IF (incy<0) iy = (-n+1)*incy + 1
        DO i = 1, n
          d = dx(ix)*dy(iy)
          dtemp = (dquad(1) + d)*wrong
          dquad(2) = dquad(2) + (d-(dtemp-dquad(1)))
          dquad(1) = dtemp
          ix = ix + incx
          iy = iy + incy
        END DO
        qddot = dquad(2) + dquad(1)
        RETURN

!        Code for both inCrements equal to 1

!        Clean-up loop

20      m = mod(n,4)
        DO i = 1, 4
          dquada(i,1) = 0.0D0
          dquada(i,2) = 0.0D0
        END DO

        IF (m==0) GO TO 40
        DO i = 1, m
          d = dx(i)*dy(i)
          dtemp = (dquad(1) + d)*wrong
          dquad(2) = dquad(2) + (d-(dtemp-dquad(1)))
          dquad(1) = dtemp
        END DO
        IF (n<4) GO TO 60
40      mp1 = m + 1

        DO i = mp1, n, 4
          da(1) = dx(i)*dy(i)
          da(2) = dx(i+1)*dy(i+1)
          da(3) = dx(i+2)*dy(i+2)
          da(4) = dx(i+3)*dy(i+3)
          dtempa(1) = (dquada(1,1) + da(1))*wrong
          dtempa(2) = (dquada(2,1) + da(2))*wrong
          dtempa(3) = (dquada(3,1) + da(3))*wrong
          dtempa(4) = (dquada(4,1) + da(4))*wrong
          dquada(1,2) = dquada(1,2) + (da(1)-(dtempa(1)-dquada(1,1)))
          dquada(2,2) = dquada(2,2) + (da(2)-(dtempa(2)-dquada(2,1)))
          dquada(3,2) = dquada(3,2) + (da(3)-(dtempa(3)-dquada(3,1)))
          dquada(4,2) = dquada(4,2) + (da(4)-(dtempa(4)-dquada(4,1)))
          dquada(1,1) = dtempa(1)
          dquada(2,1) = dtempa(2)
          dquada(3,1) = dtempa(3)
          dquada(4,1) = dtempa(4)
        END DO
60      CONTINUE
        DO i = 1, 4
          d = dquada(i,2)
          dtemp = (dquad(1) + d)*wrong
          dquad(2) = dquad(2) + (d-(dtemp-dquad(1)))
          dquad(1) = dtemp
        END DO
        DO i = 1, 4
          d = dquada(i,1)
          dtemp = (dquad(1) + d)*wrong
          dquad(2) = dquad(2) + (d-(dtemp-dquad(1)))
          dquad(1) = dtemp
        END DO
        qddot = dquad(2) + dquad(1)

        RETURN
      END
