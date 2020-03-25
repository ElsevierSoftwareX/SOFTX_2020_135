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

!> @brief cpf(i,j,k,ismpl) calculates the isobaric heat capacity in (in J/kg/K)
!> @param[in] i cell index, direction I0
!> @param[in] j cell index, direction J0
!> @param[in] k cell index, direction K0
!> @param[in] ismpl local sample index
!> @return cpf  [J/kg/K]
!> @details
!> cpf(i,j,k,ismpl) calculates the isobaric heat capacity in (in
!> J/kg/K)\n of pure water, given temperature (t, in C), and
!> pressure (p,in Pa)\n at node(i,j,k).\n \n
!>
!> method: c_p = d/dT E, E= fluid enthalpy.\n \n
!>
!>  Main source Zyvoloski1997: \n
!>
!> Zyvoloski, G.A., Robinson, B.A., Dash, Z.V., & Trease, L.L. Summary
!> of the models and methods for the FEHM application - a
!> finite-element heat- and mass-transfer code. United
!> States. doi:10.2172/565545. \n \n
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
!> range of validity:\n
!> pressures      0.001 - 110 MPa,\n
!> temperature   15 - 350 degC\n
      double precision function cpf(i,j,k,ismpl)
        use arrays, only: pres, temp
        use mod_flow, only: pa_conv1

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

        ! Enthalpy (J/kg)
        double precision :: enth

        ! Derivative of enthalpy wrt T (J/kg/K)
        double precision :: denthdt
        
        ! Monomials of temperature and pressure
        double precision :: t, t2, t3
        double precision :: p, p2, p3, p4
        double precision :: tp, t2p, tp2

        ! Coefficients of numerator of rational function approximation
        double precision, parameter :: Y0 = 0.25623465D-3
        double precision, parameter :: Y1 = 0.10184405D-2
        double precision, parameter :: Y2 = 0.22554970D-4
        double precision, parameter :: Y3 = 0.34836663D-7
        double precision, parameter :: Y4 = 0.41769866D-2 
        double precision, parameter :: Y5 = -0.21244879D-4
        double precision, parameter :: Y6 = 0.25493516D-7
        double precision, parameter :: Y7 = 0.89557885D-4
        double precision, parameter :: Y8 = 0.10855046D-6
        double precision, parameter :: Y9 = -0.21720560D-6

        ! Coefficients of denominator of rational function approximation
        double precision, parameter :: Z0 = 0.10000000D+1
        double precision, parameter :: Z1 = 0.23513278D-1
        double precision, parameter :: Z2 = 0.48716386D-4
        double precision, parameter :: Z3 = -0.19935046D-8
        double precision, parameter :: Z4 = -0.50770309D-2
        double precision, parameter :: Z5 = 0.57780287D-5
        double precision, parameter :: Z6 = 0.90972916D-9
        double precision, parameter :: Z7 = -0.58981537D-4
        double precision, parameter :: Z8 = -0.12990752D-7
        double precision, parameter :: Z9 = 0.45872518D-8

        ! Numerator and denominator of rational function approximation
        double precision :: ta, tb

        ! Derivative of numerator wrt T
        double precision :: da

        ! Derivative of denominator wrt T
        double precision :: db

        ! Denominator squared
        double precision :: b2


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

        ! Enthalpy
        enth = ta/tb

        ! Derivative of numerator
        da = Y4 + 2.0d0*Y5*t + 3.0d0*Y6*t2 + Y7*p + &
            Y8*p2 + 2.0d0*Y9*tp

        ! Derivative of denominator
        db = Z4 + 2.0d0*Z5*t + 3.0d0*Z6*t2 + Z7*p + &
            Z8*p2 + 2.0d0*Z9*tp

        ! Denominator squared
        b2 = tb*tb

        ! Derivative, quotient rule
        denthdt = da/tb - ta*db/b2

        ! Isobaric heat capacity (J/kg/K)
        cpf = denthdt*1.0d6

        return

      end function cpf
