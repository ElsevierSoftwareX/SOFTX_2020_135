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

!>    @brief cpf(i,j,k,ismpl) calculates the isobaric heat capacity in (in J/kg/K)
!>    @param[in] i cell index, direction I0
!>    @param[in] j cell index, direction J0
!>    @param[in] k cell index, direction K0
!>    @param[in] ismpl local sample index
!>    @return cpf [J/kg/K]
!>    @details
!>    cpf(i,j,k,ismpl) calculates the isobaric heat capacity in (in J/kg/K)\n
!>    of pure water, given temperature (t, in C), and pressure (p,in Pa)\n
!>    at node(i,j,k).\n
!>    method: c_p = d/dT E, E= fluid enthaply.\n
!>    derived from the formulation given in:\n
!>          zylkovskij et al: models and methods summary for\n
!>          the fehmn application,\n
!>          ecd 22, la-ur-94-3787, los alamos nl, 1994.\n
!>          Speedy, R.J. (1987) Thermodynamic properties of supercooled water \n
!>          at 1 atm. Journal of Physical Chemistry, 91: 3354–3358. \n
!>    range of validity:\n
!>                   pressures      0.01 - 110 MPa,\n
!>                   temperature   0.001 - 350 °c\n
!>    input:\n
!>      pressure                            plocal [Pa]\n
!>      temperature                         tlocal in [C]\n
      DOUBLE PRECISION FUNCTION cpf(i,j,k,ismpl)
        use arrays
        use mod_flow
        IMPLICIT NONE

        INTEGER i, j, k, ismpl
        DOUBLE PRECISION plocal, tlocal, enth, denthdt, tred, p, p1, &
          p2, p3, p4, t, t1, t2, t3, tp, t2p, tp2, ta, tb, da, db, b2
        DOUBLE PRECISION cf(20), bf(6)
        DATA cf/0.25623465D-3, 0.10184405D-2, 0.22554970D-4, &
          0.34836663D-7, 0.41769866D-2, -0.21244879D-4, 0.25493516D-7, &
          0.89557885D-4, 0.10855046D-6, -0.21720560D-6, 0.10000000D+1, &
          0.23513278D-1, 0.48716386D-4, -0.19935046D-8, &
          -0.50770309D-2, 0.57780287D-5, 0.90972916D-9, &
          -0.58981537D-4, -0.12990752D-7, 0.45872518D-8/

!     new: after Speedy (1987) for T < 0 to -46 C
        DATA bf/14.2D0, 25.952D0, 128.281D0, -221.405D0, 196.894D0, &
          -64.812D0/
!     end new


        plocal = pres(i,j,k,ismpl)*pa_conv1
        tlocal = temp(i,j,k,ismpl)
        IF (tlocal<-45D0) tlocal = -45.D0

        IF (tlocal<0.D0) THEN
!     new: after Speedy (1987) for T < 0 to -46 C
          tred = (tlocal+273.15D0-227.15D0)/227.15D0
          cpf = bf(1)/sqrt(tred) + bf(2) + bf(3)*tred + &
            bf(4)*tred*tred + bf(5)*tred*tred*tred + &
            bf(6)*tred*tred*tred*tred
          cpf = cpf*0.99048992406520D0*1000.D0/18.D0
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
! enthalpy
          ta = cf(1) + cf(2)*p + cf(3)*p2 + cf(4)*p3 + cf(5)*t + &
            cf(6)*t2 + cf(7)*t3 + cf(8)*tp + cf(10)*t2p + cf(9)*tp2
          tb = cf(11) + cf(12)*p + cf(13)*p2 + cf(14)*p3 + cf(15)*t + &
            cf(16)*t2 + cf(17)*t3 + cf(18)*tp + cf(20)*t2p + &
            cf(19)*tp2
          enth = ta/tb

! derivative
          da = cf(5) + 2.D0*cf(6)*t + 3.D0*cf(7)*t2 + cf(8)*p + &
            2.D0*cf(10)*tp + cf(9)*p2
          db = cf(15) + 2.D0*cf(16)*t + 3D0*cf(17)*t2 + cf(18)*p + &
            2.D0*cf(20)*tp + cf(19)*p2

          b2 = tb*tb
          denthdt = da/tb - ta*db/b2
          cpf = denthdt*1.D6
        END IF

        RETURN
      END
