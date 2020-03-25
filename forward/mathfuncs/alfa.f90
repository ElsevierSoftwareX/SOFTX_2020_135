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

!>    @brief truncated series expansion of cosh(x) for x < 0.1
!>    @param[in] x input value
!>    @return truncated function value
      double precision FUNCTION alfa(x)
        IMPLICIT NONE
        DOUBLE PRECISION x
!      double precision x2, x3, x5, x7
        alfa = 0.D0
        IF (abs(x)>1.D-30) alfa = (x/tanh(x)-1.D0)/x

!      if (abs(x).lt.1.d-1) then
!         x2 = x*x
!         x3 = x2*x
!         x5 = x3*x2
!         x7 = x5*x2
!         alfa = x/3.d0 - x3/4.5d1 + 2.d0*x5/9.45d2 -
!     & x7/4.725d3
!         alfa = x*3.333333333333333d-01
!     &        - x3*2.222222222222222d-02
!     &        + x5*2.116402116402117d-03
!     &        - x7*2.116402116402116e-04
!      else
!         alfa=(x/tanh(x) - 1.d0)/x
!      endif
        RETURN
      END
