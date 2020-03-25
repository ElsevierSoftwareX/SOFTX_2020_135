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

!> @brief rhow calculates the density of pure water [kg/m3]
!> @param[in] p pressure [MPa]
!> @param[in] t temperature [degC]
!> @return rhow density of pure water [kg/m3]
!> @details
!> rhow(i,j,k,ismpl) calculates the density in (in kg/m^3) of pure
!> water, given temperature (t, in degC), and pressure (p,in Pa)\n \n
!>
!> Source: \n
!>
!> Magri, F. (2005), Mechanismus und Fluiddynamik der
!> Salzwasserzirkulation im Norddeutschen Becken: Ergebnisse
!> thermohaliner numerischer Simulationen, Doktorarbeit,
!> Geoforschungszentrum Potsdam \n \n
!>
!> Specifically: Equation (2.12) and Table 2-1
      double precision function rhow(p,t)

        implicit none

        ! Local temperature (degC)
        double precision, intent (in) :: t

        ! Local pressure (MPa)
        double precision, intent (in) :: p

        ! Pressure in kPa
        double precision :: pk

        ! Monomials in pressure and temperature
        double precision :: pk2, t2, t4, t6

        ! Coefficients from Magri2005 Equation (2.12)
        double precision a, b, c, d, e, f, g

        ! Coefficients from Magri2005 Table 2-1
        double precision, parameter :: a0 = 9.99792877961606D02
        double precision, parameter :: a1 = 5.07605113140940D-04
        double precision, parameter :: a2 = -5.28425478164183D-10
        double precision, parameter :: b0 = 5.13864847162196D-02
        double precision, parameter :: b1 = -3.61991396354483D-06
        double precision, parameter :: b2 = 7.97204102509724D-12
        double precision, parameter :: c0 = -7.53557031774437D-03
        double precision, parameter :: c1 = 6.32712093275576D-08
        double precision, parameter :: c2 = -1.66203631393248D-13
        double precision, parameter :: d0 = 4.60380647957350D-05
        double precision, parameter :: d1 = -5.61299059722121D-10
        double precision, parameter :: d2 = 1.80924436489400D-15
        double precision, parameter :: e0 = -2.26651454175013D-07
        double precision, parameter :: e1 = 3.36874416675978D-12
        double precision, parameter :: e2 = -1.30352149261326D-17
        double precision, parameter :: f0 = 6.14889851856743D-10
        double precision, parameter :: f1 = -1.06165223196756D-14
        double precision, parameter :: f2 = 4.75014903737416D-20
        double precision, parameter :: g0 = -7.39221950969522D-13
        double precision, parameter :: g1 = 1.42790422913922D-17
        double precision, parameter :: g2 = -7.13130230531541D-23


        ! Unit [kPa] is needed after Magri(2005); the forwarded p is in [MPa]
        pk = p*1000

        ! Monomials in pressure and temperature
        pk2 = pk*pk
        t2 = t*t
        t4 = t2*t2
        t6 = t4*t2

        ! Coefficients from Magri2005, Equation (2.12)
        a = a0 + a1*pk + a2*pk2
        b = b0 + b1*pk + b2*pk2
        c = c0 + c1*pk + c2*pk2
        d = d0 + d1*pk + d2*pk2
        e = e0 + e1*pk + e2*pk2
        f = f0 + f1*pk + f2*pk2
        g = g0 + g1*pk + g2*pk2

        ! Water density [kg/m3]
        rhow = a + b*t + c*t2 + d*t2*t + e*t4 + f*t4*t + g*t6

        return

      end function rhow
