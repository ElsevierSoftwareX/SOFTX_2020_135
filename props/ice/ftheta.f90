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

!>    @brief calculates the fluid/ice partition function
!>    @param[in] tlocal = temperature
!>    @param[in] tliq = phase boundary temperature
!>    @param[in] wmushy =  scaling factor for "mushy" region = Tsolidus-tliqidus
!>    @param[in] ismpl local sample index
!>    @param[out] theta = value of partition function  (0<Theta<1)
!>    @param[out] dtheta = value of derivative of partition function with respect to temperature
!#>   returns the fluid/ice partition function theta and its derivative dtheta
!>    @details
!>calculates the fluid/ice partition function used for the \n
!>apparent heat capacity approach to phase change.\n
      SUBROUTINE ftheta(tlocal,tliq,wmushy,theta,dtheta,ismpl)
        use arrays
                use ice
        IMPLICIT NONE
        INTEGER ismpl
        DOUBLE PRECISION tlocal, theta, dtheta, val, tliq, wmushy

        theta = 1.D0
        dtheta = 0.D0

        IF (tlocal<tliq) THEN
          val = (tlocal-tliq)/wmushy
!         Theta=exp((-(T-Tf)/w)^2)
          theta = exp(-val*val)
!         dTheta=-2*(T-Tf)/w^2*Theta;
          dtheta = -2.D0*val*theta/wmushy
!       IF (theta.lt.liqmin) THEN
!         theta  =  liqmin
!         dtheta =  0.
!       END IF

        END IF
        RETURN
      END


!>    @brief calculates the fluid/ice partition function used for the
!>    @param[in] tlocal = temperature
!>    @param[in] tliq = phase boundary temperature
!>    @param[in] wmushy = scaling factor for "mushy" region = tsolidus-tliqidus
!>    @param[in] ismpl local sample index
!>    @param[out] theta = value of partition function (0<Theta<1)
!>    @param[out] dtheta = value of derivative of partition function with respect to temperature
!>    @details
!>calculates the fluid/ice partition function used for the \n
!>apparent heat capacity approach to phase change (galushkin)\n
      SUBROUTINE ftheta_gal(tlocal,tliq,wmushy,theta,dtheta,ismpl)
        use arrays
                use ice
        IMPLICIT NONE
        INTEGER ismpl
        DOUBLE PRECISION solid, tlocal, tref, tlimit, km, theta, &
          dtheta, a0, a1, a2, a3, a4, a5, a6, a7, a8, tliq, wmushy
        DOUBLE PRECISION t, t1, t2, t3, t4, t5, t6, t7, t8
        PARAMETER (a0=1.0D0,a1=0.60152823763179D0, &
          a2=0.23218232347212D0,a3=0.04669792788297D0, &
          a4=0.00535597924776D0,a5=0.00036415588418D0, &
          a6=0.00001450956751D0,a7=0.00000031279149D0, &
          a8=0.00000000281502D0)

        theta = 1.D0
        dtheta = 0.0

        IF (tlocal<tliq) THEN
          t1 = tlocal
          t2 = t1*t1
          t3 = t2*t1
          t4 = t3*t1
          t5 = t4*t1
          t6 = t5*t1
          t7 = t6*t1
          t8 = t7*t1
          theta = a0 + a1*t1 + a2*t2 + a3*t3 + a4*t4 + a5*t5 + a6*t6 + &
            a7*t7 + a8*t8
          dtheta = a1 + 2.D0*a2*t1 + 3.D0*a3*t2 + 4.D0*a4*t3 + &
            5.D0*a5*t4 + 6.D0*a6*t5 + 7.D0*a7*t6 + 8.D0*a8*t7
        END IF
        RETURN
      END

!>    @brief calculates the fluid/ice partition function used for the 
!>    @param[in] tlocal = temperature
!>    @param[in] tliq = phase boundary temperature
!>    @param[in] wpara = scaling factor for "mushy" region = tsolidus-tliqidus
!>    @param[in] ismpl local sample index
!>    @param[out] theta = value of partition function (0<Theta<1)
!>    @param[out] dtheta = value of derivative of partition function with respect to temperature
!>    @details
!>calculates the fluid/ice partition function used for the \n
!>apparent heat capacity approach to phase change (nikolsky2009)\n
      SUBROUTINE ftheta_nik(tlocal,tliq,wpara,theta,dtheta,ismpl)
        use arrays
                use ice
        IMPLICIT NONE
        INTEGER ismpl
        DOUBLE PRECISION solid,tlocal,tref,tlimit,km,theta,dtheta
        DOUBLE PRECISION tliq,wpara

        theta    =   1.d0
        dtheta   =   0.0d0

        IF (tlocal<tliq .and. wpara>=0.5d0 .and. wpara<=0.8d0) THEN
          theta = (abs(tliq)**wpara)*(abs(tlocal)**(-wpara))

!AW
          dtheta= (abs(tliq)**wpara)* &
            (-wpara)*abs(tlocal)**(-wpara-1.d0)
!AW oben stehende Zwischenloesung muss mit unten stehendem Original korrigiert werden <- SIGN ist aber falsch!!!
          WRITE(*,*) 'error: "ftheta_nik" needs to be modified!!!'
          STOP
!AW-original          dtheta= (abs(tliq)**wpara)* &
!AW-original            (-wpara)*abs(tlocal)**(-wpara-1.d0)*sign(tlocal)

        END IF
        RETURN
      END
