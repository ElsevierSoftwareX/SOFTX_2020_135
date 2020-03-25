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

!>    @brief rhof(i,j,k,ismpl) calculates the density in (in kg/m^3) of pure water,
!>    @param[in] i cell index, direction I0
!>    @param[in] j cell index, direction J0
!>    @param[in] k cell index, direction K0
!>    @param[in] ismpl local sample index
!>    @return rhof  [kg/m^3]
!>    @details
!>    rhof(i,j,k,ismpl) calculates the density in (in kg/m^3) of pure water,\n
!>    given temperature (t, in c), and pressure (p,in Mpa) at node(i,j,k)\n
!>    derived from the formulation given in:\n
!>          zylkovskij et al: models and methods summary for\n
!>          the fehmn application,\n
!>           ecd 22, la-ur-94-3787, los alamos nl, 1994.\n
!>          Speedy, R.J. (1987) Thermodynamic properties of supercooled water\n
!>          at 1 atm. Journal of Physical Chemistry, 91: 3354–3358.
!>    range of validity:\n
!>      pressures   0.01 - 110 MPa,\n
!>      temperature   15 - 350 °C\n
!>    input:\n
!>      pressure                            plocal [Pa]\n
!>      temperature                         tlocal in [C]\n
      DOUBLE PRECISION FUNCTION rhof(i,j,k,ismpl)
        use arrays
        use mod_flow
        IMPLICIT NONE

        INTEGER i, j, k, ismpl
        DOUBLE PRECISION cf(20), bf(6)
        DOUBLE PRECISION ta, tb, tlocal, plocal, t, t2, t3, tred, p, &
          p2, p3, p4, tp, t2p, tp2, contfac

        DATA cf/0.10000000D+01, 0.17472599D-01, -0.20443098D-04, &
          -0.17442012D-06, 0.49564109D-02, -0.40757664D-04, &
          0.50676664D-07, 0.50330978D-04, 0.33914814D-06, &
          -0.18383009D-06, 0.10009476D-02, 0.16812589D-04, &
          -0.24582622D-07, -0.17014984D-09, 0.48841156D-05, &
          -0.32967985D-07, 0.28619380D-10, 0.53249055D-07, &
          0.30456698D-09, -0.12221899D-09/
!     new: after Speedy (1987) for T < 0 to -46 C (see also Grant, S. A., Physical \n
!          and Chemical Factors Affecting Contaminant Hydrology in Cold Environments, \n
!          Technical Report ERDC/CRREL TR-00-21, US Army Corps of Engineers, 2000) \n
        DATA bf/901.5328593D0, -0.0011761652D0, 0.0038442382D0, &
          -0.0157270761D0, 0.0744064614D0, -0.1406432653D0/
!     end new

        contfac = 0.999195706402050D0
        plocal = pres(i,j,k,ismpl)*pa_conv1
        tlocal = temp(i,j,k,ismpl)
        IF (tlocal<-45.D0) tlocal = -45.D0

        IF (tlocal<0.D0) THEN
!     new: after Speedy (1987) for T < 0 to -46 C
          tred = (tlocal+273.15D0-227.15D0)/227.15D0
          rhof = contfac*bf(1)*exp(-227.15D0*(bf(3)*tred+0.5D0*bf(4)* &
            tred*tred+0.333D0*bf(5)*tred*tred*tred+0.25D0*bf( &
            6)*tred*tred*tred*tred+2.D0*bf(2)*sqrt(tred)))
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
          rhof = ta/tb
        END IF

        RETURN
      END
