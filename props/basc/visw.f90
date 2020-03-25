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

!> @brief rhof(i,j,k,ismpl) calculates the viscosity in (in Pa s) of pure water
!> @param[in] p cell pressure [MPa]
!> @param[in] t temperature [degC]
!> @return visw  [Pa s]
!> @details
!> rhof(i,j,k,ismpl) calculates the viscosity in (in Pa s) of  pure water,\n
!> given temperature (t, in C), and pressure (p,in Pa) at node(i,j,k)\n\n
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
      double precision function visw(p,t)

        implicit none

        ! Temperature (degC)
        double precision, intent (in) :: t

        ! Pressure (MPa)
        double precision, intent (in) :: p

        ! Coefficients of numerator of rational function approximation
        double precision, parameter :: Y0 = 0.17409149D-02
        double precision, parameter :: Y1 = 0.18894882D-04
        double precision, parameter :: Y2 = -0.66439332D-07
        double precision, parameter :: Y3 = -0.23122388D-09
        double precision, parameter :: Y4 = -0.31534914D-05
        double precision, parameter :: Y5 = 0.11120716D-07
        double precision, parameter :: Y6 = -0.48576020D-10
        double precision, parameter :: Y7 = 0.28006861D-07
        double precision, parameter :: Y8 = 0.23225035D-09
        double precision, parameter :: Y9 = 0.47180171D-10

        ! Coefficients of denominator of rational function approximation
        double precision, parameter :: Z0 =  0.10000000D+01
        double precision, parameter :: Z1 = 0.10523153D-01
        double precision, parameter :: Z2 = -0.22658391D-05
        double precision, parameter :: Z3 = -0.31796607D-06
        double precision, parameter :: Z4 = 0.29869141D-01
        double precision, parameter :: Z5 = 0.21844248D-03
        double precision, parameter :: Z6 = -0.87658855D-06
        double precision, parameter :: Z7 = 0.41690362D-03
        double precision, parameter :: Z8 = -0.25147022D-05
        double precision, parameter :: Z9 = 0.22144660D-05
        
        ! Monomials of temperature and pressure
        double precision :: t2, t3
        double precision :: p2, p3, p4
        double precision :: tp, t2p, tp2

        ! Numerator and denominator of rational function approximation
        double precision :: ta, tb


        ! Compute monomials in pressure and temperature
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

        ! Viscosity
        visw = ta/tb

        return

      end function visw
