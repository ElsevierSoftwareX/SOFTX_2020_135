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

!> @brief rhof(i,j,k,ismpl) calculates the density in (in kg/m^3) of pure water,
!> @param[in] i cell index, direction I0
!> @param[in] j cell index, direction J0
!> @param[in] k cell index, direction K0
!> @param[in] ismpl local sample index
!> @return rho  [kg/m^3]
!> @details
!> rhof(i,j,k,ismpl) calculates the density in (in kg/m^3) of pure
!> water, given temperature (t, in degC), and pressure (p,in Pa) at
!> node(i,j,k)\n \n
!>
!>  Main source Zyvoloski1997: \n
!>
!> Zyvoloski, G.A., Robinson, B.A., Dash, Z.V., & Trease, L.L. Summary
!> of the models and methods for the FEHM application - a
!> finite-element heat- and mass-transfer code. United
!> States. doi:10.2172/565545. \n \n
!>
!> See Section 8.4.3. of Zyvoloski1997 for an explanation of the
!> "Rational function approximation" used in this subroutine. \n \n
!> The approximation uses the table of coefficients in Appendix 10 of
!> Zyvoloski1997.\n
!>
!> Alternative source (same text, more modern, without doi): \n
!> https://fehm.lanl.gov/orgs/ees/fehm/pdfs/fehm_mms.pdf \n \n
!>
!> The table of coefficients from Zyvoloski1997 describes the physical
!> values found in Haar1984: \n
!>
!> Lester Haar, John Gallagher, George Kell, NBS/NRC Steam Tables:
!> Thermodynamic and Transport Properties and Computer Programs for
!> Vapor and Liquid States of Water in SI Units, Hemisphere Publishing
!> Corporation, Washington, 1984. \n \n
!>
!>    range of validity:\n
!>    - pressures   0.001 - 110 MPa,\n
!>    - temperature   15 - 360 degC\n
      double precision function rhof(i,j,k,ismpl)
        use arrays
        use mod_flow

        implicit none

        ! Location indices
        integer, intent (in) :: i
        integer, intent (in) :: j
        integer, intent (in) :: k

        ! Sample index
        integer :: ismpl

        ! Temperature (degC)
        double precision :: tlocal

        ! Pressure (MPa)
        double precision :: plocal

        ! Monomials of temperature and pressure
        double precision :: t, t2, t3
        double precision :: p, p2, p3, p4
        double precision :: tp, t2p, tp2

        ! Coefficients of numerator of rational function approximation
        double precision, parameter :: Y0 = 0.10000000D+01
        double precision, parameter :: Y1 = 0.17472599D-01
        double precision, parameter :: Y2 = -0.20443098D-04
        double precision, parameter :: Y3 = -0.17442012D-06
        double precision, parameter :: Y4 = 0.49564109D-02
        double precision, parameter :: Y5 = -0.40757664D-04
        double precision, parameter :: Y6 = 0.50676664D-07
        double precision, parameter :: Y7 = 0.50330978D-04
        double precision, parameter :: Y8 = 0.33914814D-06
        double precision, parameter :: Y9 = -0.18383009D-06

        ! Coefficients of denominator of rational function approximation
        double precision, parameter :: Z0 = 0.10009476D-02
        double precision, parameter :: Z1 = 0.16812589D-04
        double precision, parameter :: Z2 = -0.24582622D-07
        double precision, parameter :: Z3 = -0.17014984D-09
        double precision, parameter :: Z4 = 0.48841156D-05
        double precision, parameter :: Z5 = -0.32967985D-07
        double precision, parameter :: Z6 = 0.28619380D-10
        double precision, parameter :: Z7 = 0.53249055D-07
        double precision, parameter :: Z8 = 0.30456698D-09
        double precision, parameter :: Z9 = -0.12221899D-09

        ! Numerator and denominator of rational function approximation
        double precision :: ta, tb


        ! Local Pressure in MPa
        plocal = pres(i,j,k,ismpl)*pa_conv1

        ! Local Temperature in degC
        tlocal = temp(i,j,k,ismpl)

        ! Temperature out of bounds
        if (tlocal > 360.0d0) then
          write (*,*) "[E1]: Error: Temperature (",&
              tlocal,") out of bounds (> 360 degC) at ", i,j,k
          stop
        end if
        if (tlocal < 0.0d0) then
          ! Relax table boundary of 15degC to error boundary 0degC
          write (*,*) "[E2]: Error: Temperature (",&
              tlocal,") out of bounds (< 0 degC) at ", i,j,k
          stop
        end if

        ! Pressure out of bounds
        if (plocal > 110.0d0) then
          write (*,*) "[E3]: Error: Pressure (",&
              plocal,") out of bounds (> 110 MPa) at ", i,j,k
          stop
        end if
        if (plocal < 0.001d0) then
          write (*,*) "[E4]: Error: Pressure (",&
              plocal,") out of bounds (< 0.001 MPa) at ", i,j,k
          stop
        end if

        ! Compute monomials in pressure and temperature
        p = plocal
        t = tlocal
        p2 = p*p
        p3 = p2*p
        p4 = p3*p
        t2 = t*t
        t3 = t2*t
        tp = p*t
        tp2 = t*p2
        t2p = t2*p

        ! Numerator of rational function approximation
        ta = Y0 + Y1*p + Y2*p2 + Y3*p3 + Y4*t + &
            Y5*t2 + Y6*t3 + Y7*tp + Y8*tp2 + Y9*t2p

        ! Denominator of rational function approximation
        tb = Z0 + Z1*p + Z2*p2 + Z3*p3 + Z4*t + &
            Z5*t2 + Z6*t3 + Z7*tp + Z8*tp2 + Z9*t2p

        ! Water density
        rhof = ta/tb

        return

      end function rhof
