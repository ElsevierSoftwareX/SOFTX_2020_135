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

!>    @brief rhof(i,j,k,ismpl) calculates the viscosity in (in Pa s) of  pure water,
!>    @param[in] i cell index, direction I0
!>    @param[in] j cell index, direction J0
!>    @param[in] k cell index, direction K0
!>    @param[in] ismpl local sample index
!>    @return visf  [Pa s]
!>    @details
!>    rhof(i,j,k,ismpl) calculates the viscosity in (in Pa s) of  pure water,\n
!>    given temperature (t, in C), and pressure (p,in Pa) at node(i,j,k)\n
!>    derived from the formulation given in:\n
!>          zylkovskij et al: models and methods summary for\n
!>          the fehmn application,\n
!>           ecd 22, la-ur-94-3787, los alamos nl, 1994.\n
!>    Speedy, R.J. (1987) Thermodynamic properties of supercooled water\n
!>          at 1 atm. Journal of Physical Chemistry, 91: 3354–3358.
!>    range of validity:\n
!>      pressures   0.01 - 110 mpa,\n
!>      temperature   15 - 350 °c and -46°c - 0°c\n
!>    input:\n
!>      pressure                            plocal [Pa]\n
!>      temperature                         tlocal in [C]\n
      DOUBLE PRECISION FUNCTION visf(i,j,k,ismpl)
        use arrays
        use mod_flow
        IMPLICIT NONE

        INTEGER i, j, k, ismpl
        DOUBLE PRECISION cf(20), bf(6)
        DOUBLE PRECISION ta, tb, tlocal, plocal, t, t2, t3, tred, p, &
          p2, p3, p4, tp, t2p, tp2

        DATA cf/0.17409149D-02, 0.18894882D-04, -0.66439332D-07, &
          -0.23122388D-09, -0.31534914D-05, 0.11120716D-07, &
          -0.48576020D-10, 0.28006861D-07, 0.23225035D-09, &
          0.47180171D-10, 0.10000000D+01, 0.10523153D-01, &
          -0.22658391D-05, -0.31796607D-06, 0.29869141D-01, &
          0.21844248D-03, -0.87658855D-06, 0.41690362D-03, &
          -0.25147022D-05, 0.22144660D-05/
!     new: after Speedy (1987) for T < 0 to -46 C
        DATA bf/26.312D0, -144.565D0, 1239.075D0, -8352.579D0, &
          31430.760, -48576.798D0/
!     end new

        plocal = pres(i,j,k,ismpl)*pa_conv1
        tlocal = temp(i,j,k,ismpl)
        IF (tlocal<-45D0) tlocal = -45.D0

        IF (tlocal<0.D0) THEN
! tloCal = 0.d0
!     new: after Speedy (1987) for T < 0 to -46 C
          tred = (tlocal+273.15D0-227.15D0)/227.15D0
          visf = bf(1)/sqrt(tred) + bf(2) + bf(3)*tred + &
            bf(4)*tred*tred + bf(5)*tred*tred*tred + &
            bf(6)*tred*tred*tred*tred
          visf = visf*0.001
        ELSE
!     end new
          IF (tlocal>300.D0) tlocal = 300.D0

          p = plocal
          t = tlocal
          p2 = p*p
          p3 = p2*p
          p4 = p3*p
          t2 = t*t
          t3 = t2*t
          tp = p*t
          t2p = t2*p
          tp2 = t*p2

          ta = cf(1) + cf(2)*p + cf(3)*p2 + cf(4)*p3 + cf(5)*t + &
            cf(6)*t2 + cf(7)*t3 + cf(8)*tp + cf(10)*t2p + cf(9)*tp2
          tb = cf(11) + cf(12)*p + cf(13)*p2 + cf(14)*p3 + cf(15)*t + &
            cf(16)*t2 + cf(17)*t3 + cf(18)*tp + cf(20)*t2p + &
            cf(19)*tp2
          visf = ta/tb
        END IF

        RETURN
      END
