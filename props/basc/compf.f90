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

!> @brief compf(i,j,k,ismpl) calculates the compressibility of the fluid
!> @param[in] i cell index, direction I0
!> @param[in] j cell index, direction J0
!> @param[in] k cell index, direction K0
!> @param[in] ismpl local sample index
!> @return compf  [1 / Pa]
!> @details
!> compf(i,j,k,ismpl) calculates the compressibility in (in 1/Pa) of
!> salin water, given temperature (t, in c) pressure (p,in pa), and
!> salinity (s, in mol/L) at node(i,j,k)\n\n
!>
!> Sources:\n
!> Driesner & Heinrich, The system H2O-NaCl. Part I: Correlation
!>    formulae for phase relations in temperature-pressure-composition
!>    space from 0 to 1000 C, 0 to 5000 bar, and 0 to 1 XNaCl
!>    Geochimica et Cosmochimica Acta 71 (2007) 4880-4901\n\n
!>
!> Driesner, The system H2O-NaCl. Part II: Correlations for molar
!>    volume, enthalpy, and isobaric heat capacity from 0 to 1000 C, 1
!>    to 5000 bar, and 0 to 1 XNaCl Geochimica et Cosmochimica Acta 71
!>    (2007) 4902-4919\n\n
      double precision function compf(i,j,k,ismpl)
        use arrays
        use mod_flow

        implicit none

        ! Location indices
        integer, intent (in) :: i
        integer, intent (in) :: j
        integer, intent (in) :: k

        ! Sample index
        integer :: ismpl

        ! Molar mass of NaCl [g/mol]
        double precision, parameter :: mmnacl = 58.44277d0
        
        ! Mass fraction of NaCl in solution [-]
        double precision :: fracnacl
        
        ! Molar mass of water, H2O [g/mol]
        double precision, parameter :: mmwater = 18.01528d0
        
        ! Temperature, [degC]
        double precision :: t

        ! Pressure [MPa]
        double precision :: p

        ! Monomials of pressure [bar]
        double precision :: pb, pb2, pb3

        ! Pure water compressibility [1/Pa]
        double precision cw
        double precision, external :: compw

        ! Pure water density [kg/m3]
        double precision rw
        double precision, external :: rhow

        ! Local salinity [mol/kg]
        double precision :: s

        ! Mole fraction of NaCl [-]
        double precision :: xnacl
        
        ! Temperature at which pure water has the same molar volume
        ! as the solution [K]
        double precision :: tv

        ! Driesner2007: Coefficients

        ! Equation (9)
        double precision :: n1, n10, n11, n12

        ! Equation (10)
        double precision :: n2, n20, n21, n22, n23

        ! Equation (11), (12)
        double precision :: n1x, n2x

        ! Table 4
        double precision :: n30, n300, n301, n302
        double precision :: n31, n310, n311, n312

        ! Deviation function (13)
        double precision :: dt


        ! Local temperature [degC]
        t = temp(i,j,k,ismpl)

        ! Local pressure [MPa]
        p = pres(i,j,k,ismpl)*pa_conv1

        ! Local salinity [mol/L]
        s = tsal(i,j,k,ismpl)

        ! pure water compressibility [1 / Pa]
        cw = compw(p,t)

        if (s<=0.0d0) then

          compf = cw

        else
          ! pure water density [g/L = kg/m3]
          rw = rhow(p,t)

          ! Pressure [bar] and pressure monomials
          pb = p*10.0d0
          pb2 = pb*pb
          pb3 = pb*pb2

          ! Mass fraction of NaCl, mol/L > (g/L) / (g/L) mass fraction
          fracnacl = s*mmnacl/(rw+s*mmnacl)

          ! Mole fraction of NaCl, (mol/L) / (mol/L)
          xnacl = (fracnacl/mmnacl) / (fracnacl/mmnacl+(1-fracnacl)/mmwater)

          ! Driesner2007 parameters, Table 4
          n11 = -54.2958d0 - 45.7623d0*exp(-9.44785d-4*pb)
          n21 = -2.6142d0 - 0.000239092d0*pb
          n22 = 0.0356828d0 + 4.37235d-6*pb + 2.0566d-9*pb2
          n300 = 7.60664d6 / (pb+472.051d0)**2
          n301 = -50.0d0 - 86.1446d0*exp(-6.21128d-4*pb)
          n302 = 294.318d0 * exp(-5.66735d-3*pb)
          n310 = -0.0732761d0 * exp(-2.3772d-3*pb) - 5.2948d-5*pb
          n311 = -47.2747d0 + 24.3653d0*exp(-1.25533D-3*pb)
          n312 = -0.278529 - 0.00081381*pb

          ! Driesner2007: n1 and n2 parameters for liquid NaCl
          n1x = 330.47d0 + 0.942876d0*sqrt(pb) + 0.0817193*pb - &
            2.47556D-8*pb2 + 3.45052D-10*pb3
          n2x = -0.0370751 + 0.00237723*sqrt(pb) + 5.42049D-5*pb + &
            5.84709D-9*pb2 - 5.99373D-13*pb3

          ! Driesner2007: Set xnacl=1 in (9), => n1 = n1x
          n10 = n1x

          ! Driesner2007: Set xnacl=0 in (9) => n1 = 0.0d0
          n12 = -n10 - n11

          ! Driesner2007: Set xnacl=0 in (10) => n2 = 1.0d0
          n20 = 1.0d0 - n21*sqrt(n22)

          ! Driesner2007: Set xnacl=1 in (10) => n2 = n2x
          n23 = n2x - n20 - n21*sqrt(1.0d0+n22)

          ! Driesner2007: Equation (9)
          n1 = n10 + n11*(1.0d0-xnacl) + n12*(1.0d0-xnacl)**2

          ! Driesner2007: Equation (10)
          n2 = n20 + n21*sqrt(xnacl+n22) + n23*xnacl

          ! Driesner2007: deviation function D(T) Equations (14)-(16)
          n30 = n300 * (exp(n301*xnacl) - 1.0d0) + n302 * xnacl
          n31 = n310 * exp(n311*xnacl) + n312 * xnacl
          dt = n30 * exp(n31 * t)

          ! Temperature at which pure water has the same molar volume
          ! as the solution
          tv = n1 + n2*t+dt

          ! Pure water compressibility at tv
          compf = compw(p,tv)
          
        end if
        
        return
        
      end function compf
