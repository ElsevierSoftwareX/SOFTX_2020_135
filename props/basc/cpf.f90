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

!> @brief cpf(i,j,k,ismpl) calculates the isobaric heat capacity of the fluid [J/kg/K]
!> @param[in] i cell index, direction I0
!> @param[in] j cell index, direction J0
!> @param[in] k cell index, direction K0
!> @param[in] ismpl local sample index
!> @return cpf  [J/kg/K]
!> @details
!> cpf(i,j,k,ismpl) calculates the isobaric heat capacity in (in
!> J/kg/K) of brine, given temperature (t, in degC) pressure (p,in MPa),
!> and salinity (s, in mol/L) at node(i,j,k)\\n \n
!>
!> Source:\n
!> Driesner & Heinrich, The system H2O-NaCl. Part I: Correlation
!>    formulae for phase relations in temperature-pressure-composition
!>    space from 0 to 1000 C, 0 to 5000 bar, and 0 to 1 XNaCl
!>    Geochimica et Cosmochimica Acta 71 (2007) 4880-4901\n\n
!>
!> Driesner, The system H2O-NaCl. Part II: Correlations for molar
!>    volume, enthalpy, and isobaric heat capacity from 0 to 1000 C, 1
!>    to 5000 bar, and 0 to 1 XNaCl Geochimica et Cosmochimica Acta 71
!>    (2007) 4902-4919\n \n
      double precision function cpf(i,j,k,ismpl)
        use arrays, only: temp, pres, tsal
        use mod_flow, only: pa_conv1

        implicit none

        ! Location indices
        integer, intent (in) :: i
        integer, intent (in) :: j
        integer, intent (in) :: k

        ! Sample index
        integer :: ismpl

        ! Local Temperature, [degC]
        double precision :: t

        ! Local Pressure [MPa]
        double precision :: p

        ! Local salinity [mol/kg]
        double precision :: s

        ! Monomials of pressure [bar]
        double precision :: pb, pb2

        ! Pure water compressibility [1/Pa]
        double precision cw
        double precision, external :: compw

        ! Pure water density [kg/m3]
        double precision rw
        double precision, external :: rhow

        ! Pure water thermal conductivity [J/kg/K]
        double precision, external :: cpw

        ! Molar mass of NaCl [g/mol]
        double precision, parameter :: mmnacl = 58.44277d0
        
        ! Mass fraction of NaCl in solution [-]
        double precision :: fracnacl
        
        ! Molar mass of water, H2O [g/mol]
        double precision, parameter :: mmwater = 18.01528d0
        
        ! Mole fraction of NaCl [-]
        double precision :: xnacl
        
        ! Temperature at which pure water has the same specific
        ! enthalpy as the solution [K]
        double precision :: tv

        ! Driesner2007: Coefficients

        ! Equation (23)
        double precision :: q1, q10, q11, q12

        ! Equation (24)
        double precision :: q2, q20, q21, q22, q23

        ! Equation (25), (26)
        double precision :: q1x, q2x


        ! Local temperature [degC]
        t = temp(i,j,k,ismpl)

        ! Local pressure [MPa]
        p = pres(i,j,k,ismpl)*pa_conv1

        ! Local salinity [mol/L]
        s = tsal(i,j,k,ismpl)

        if (s<=0.0d0) then

          ! Pure water heat capacity
          cpf = cpw(p,t)

        else

          ! pure water density [g/L = kg/m3]
          rw = rhow(p,t)

          ! Pressure monomials [bar]
          pb = p*10.0d0
          pb2 = pb*pb

          ! Mass fraction of NaCl, mol/L > (g/L) / (g/L) mass fraction
          fracnacl = s*mmnacl/(rw+s*mmnacl)

          ! Mole fraction of NaCl, (mol/L) / (mol/L)
          xnacl = (fracnacl/mmnacl) / (fracnacl/mmnacl+(1-fracnacl)/mmwater)

          ! Driesner2007 parameters, Table 5
          q11 = -32.1724d0 + 0.0621255d0*pb
          q21 = -1.69513d0 - 4.52781d-4*pb - 6.04279d-8*pb2
          q22 = 0.0612567d0 + 1.88082d-5*pb

          ! Driesner2007: q1 and q2 parameters for liquid NaCl
          q1x = 47.9048d0 - 9.36994d-3*pb + 6.51059d-6*pb2
          q2x = 0.241022d0 + 3.45087d-5*pb - 4.28356d-9*pb2

          ! Driesner2007: Set xnacl=1 in (23) => q1 = q1x
          q10 = q1x

          ! Driesner2007: Set xnacl=0 in (23) => q1 = 0.0d0
          q12 = -q10 - q11

          ! Driesner2007: Set xnacl=0 in (24) => q2 = 1.0d0
          q20 = 1.0d0 - q21*sqrt(q22)

          ! Driesner2007: Set xnacl=1 in (10) => q2 = q2x
          q23 = q2x - q20 - q21*sqrt(1.0d0+q22)

          ! Driesner2007: Equation (23)
          q1 = q10 + q11*(1.0d0-xnacl) + q12*(1.0d0-xnacl)**2

          ! Driesner2007: Equation (24)
          q2 = q20 + q21*sqrt(xnacl+q22) + q23*xnacl

          ! Temperature at which pure water has the same specific
          ! enthalpy as the solution
          tv = q1 + q2*t

          ! Isobaric heat capacity of solution
          cpf = q2*cpw(p,tv)

        end if

        return

      end function cpf
