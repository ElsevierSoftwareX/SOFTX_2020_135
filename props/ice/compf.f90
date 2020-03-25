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

!>    @brief compf calculates compressibility of pure water
!>    @param[in] i cell index, direction I0
!>    @param[in] j cell index, direction J0
!>    @param[in] k cell index, direction K0
!>    @param[in] ismpl local sample index
!>    @return  compressibility                     compf  [1./Pa]
!>    @details
!>    compf calculates compressibility of pure water \n
!>    given temperature (t, in C), and pressure (p,in Pa)\n
!>    at node(plocal,tlocal).\n
!>    method: compf = 1/rhof d/dP rhof, rhof= fluid density.\n
!>    derived from the formulation given in:\n
!>          zylkovskij et al: models and methods summary for\n
!>          the fehmn application,\n
!>           ecd 22, la-ur-94-3787, los alamos nl, 1994.\n
!>    range of validity:\n
!>                   pressures      0.01 - 110 MPa,\n
!>                   temperature   0.001 - 350 °C and -46°C - 0°C\n
!>    input:\n
!>      pressure                               plocal [Pa]\n
!>      temperature                         tlocal in [C]\n
      DOUBLE PRECISION FUNCTION compf(i,j,k,ismpl)
        use arrays
        use mod_flow
        IMPLICIT NONE

        INTEGER i, j, k, ismpl
        DOUBLE PRECISION cf(20), bf(6)
        DOUBLE PRECISION ta, tb, da, db, b2, rhof_loc, drhodp, t, t2, t3, &
          tlocal, tred, p, p2, p3, p4, plocal, tp, t2p, tp2

        DATA cf/0.10000000D+01, 0.17472599D-01, -0.20443098D-04, &
          -0.17442012D-06, 0.49564109D-02, -0.40757664D-04, &
          0.50676664D-07, 0.50330978D-04, 0.33914814D-06, &
          -0.18383009D-06, 0.10009476D-02, 0.16812589D-04, &
          -0.24582622D-07, -0.17014984D-09, 0.48841156D-05, &
          -0.32967985D-07, 0.28619380D-10, 0.53249055D-07, &
          0.30456698D-09, -0.12221899D-09/

!     new: after Speedy (1987) for T < 0 to -46 C
        DATA bf/20.D0, 4.12D0, -1.13D0, 77.817D0, -78.143D0, 54.29D0/
!     end new

        plocal = pres(i,j,k,ismpl)*pa_conv1
        tlocal = temp(i,j,k,ismpl)
        IF (tlocal<-45D0) tlocal = -45.D0

        IF (tlocal<0.D0) THEN
!     new: after Speedy (1987) for T < 0 to -46 C
          tred = (tlocal+273.15D0-227.15D0)/227.15D0
          compf = bf(1)/sqrt(tred) + bf(2) + bf(3)*tred + &
            bf(4)*tred*tred + bf(5)*tred*tred*tred + &
            bf(6)*tred*tred*tred*tred
          compf = compf*1.D-11
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

! liquid density
          ta = cf(1) + cf(2)*p + cf(3)*p2 + cf(4)*p3 + cf(5)*t + &
            cf(6)*t2 + cf(7)*t3 + cf(8)*tp + cf(10)*t2p + cf(9)*tp2
          tb = cf(11) + cf(12)*p + cf(13)*p2 + cf(14)*p3 + cf(15)*t + &
            cf(16)*t2 + cf(17)*t3 + cf(18)*tp + cf(20)*t2p + &
            cf(19)*tp2
          rhof_loc = ta/tb

! derivative C   C2+2*C3*p+3*C4*p^2+C8*t+C10*t^2+2*C9*t*p
          da = cf(2) + 2.D0*cf(3)*p + 3.D0*cf(4)*p2 + cf(8)*t + &
            2.D0*cf(9)*tp + cf(10)*t2
! derivative C12+2*C13*p+3*C14*p^2+C18*t+C20*t^2+2*C19*t*p
          db = cf(12) + 2.D0*cf(13)*p + 3.D0*cf(14)*p2 + cf(18)*t + &
            2.0*cf(19)*tp + cf(20)*t2

          b2 = tb*tb
          drhodp = (da*tb-ta*db)/b2
!      Compf=rhof_loc/drhodp
!      Compf=rhof_loc
          compf = 1.E-6*drhodp/rhof_loc
        END IF

        RETURN
      END
