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

!> @brief compw calculates compressibility of pure water
!> @param[in] p_h pressure [MPa]
!> @param[in] t_h temperature [degC]
!> @return  compressibility                     compw  [1./Pa]
!> @details
!> compw calculates compressibility of pure water [1/Pa]
!> given temperature (t, in C), and pressure (p,in MPa)
!> at pressure/temperature (p_h,t).\n \n
!>
!> Method: \n
!>
!> compw = 1/rhow d/dP rhow, \n
!>
!> where rhow= water density.\n \n
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
!>    - temperature   15 - 360 degC\n \n
!>
!> input:\n
!>   pressure                               p [MPa]\n
!>   temperature                         t in [degC]\n
      double precision function compw(p_h,t_h)
        use mod_flow

        implicit none

        ! Input Pressure (MPa)
        DOUBLE PRECISION  p_h

        ! Input temperature (degc)
        DOUBLE PRECISION  t_h

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

        ! Derivative of numerator wrt P
        double precision :: da

        ! Derivative of denominator wrt P
        double precision :: db

        ! Denominator squared
        double precision :: b2

        ! Water density (local)
        double precision ::  rhow_loc

        ! Derivative of water density wrt P
        double precision :: drhodp

        ! Compressibiliy in Mpa
        double precision :: compw_mpa


        ! pressure [MPa]
        p = p_h

        ! temperature [degC]
        t = t_h

        ! Compute monomials in pressure and temperature
        p2 = p*p
        p3 = p2*p
        p4 = p3*p
        t2 = t*t
        t3 = t2*t
        tp = p*t
        t2p = t2*p
        tp2 = t*p2

        ! Numerator of rational function approximation
        ta = Y0 + Y1*p + Y2*p2 + Y3*p3 + Y4*t + &
            Y5*t2 + Y6*t3 + Y7*tp + Y8*tp2 + Y9*t2p

        ! Denominator of rational function approximation
        tb = Z0 + Z1*p + Z2*p2 + Z3*p3 + Z4*t + &
            Z5*t2 + Z6*t3 + Z7*tp + Z8*tp2 + Z9*t2p

        ! Water density
        rhow_loc = ta/tb

        ! Derivative of numerator
        da = Y1 + 2.D0*Y2*p + 3.D0*Y3*p2 + Y7*t + &
            2.D0*Y8*tp + Y9*t2

        ! Derivative of denominator
        db = Z1 + 2.D0*Z2*p + 3.D0*Z3*p2 + Z7*t + &
            2.0*Z8*tp + Z9*t2

        ! Denominator squared
        b2 = tb*tb

        ! Derivative, quotient rule
        drhodp = (da*tb-ta*db)/b2

        ! Compressibility: (1/rhow_loc) * drhodp [1/MPa]
        compw_mpa = drhodp/rhow_loc

        ! Compressibility [1/Pa]
        compw = compw_mpa / pa_conv


        return

      end function compw
