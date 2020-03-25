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

!> @brief rhof(i,j,k,ismpl) calculates the density of the fluid (in kg/m^3),
!> @param[in] i cell index, direction I0
!> @param[in] j cell index, direction J0
!> @param[in] k cell index, direction K0
!> @param[in] ismpl local sample index
!> @return rho  [kg/m^3]
!> @details
!> rhow(i,j,k,ismpl) calculates the density in (in kg/m^3) of brine,
!> given temperature (t, in c) pressure (p,in pa), and salinity (s, in
!> mol/L) at node(i,j,k)\n
!>
!> Source:\n
!> Batzle, M., & Wang, Z., Seismic properties of pore fluids,
!> GEOPHYSICS, 57(11), 1396–1408 (1992).
!> http://dx.doi.org/10.1190/1.1443207 \n\n
!>
!>   Pressures 5-100 MPa, Temperature 20-350°C, Salinity <=320 g/L\n \n
!>
!>  CODE VERIFICATION:\n
!>   INPUT:  TEMP = 298.15K  P =0.1013 MPa  S = 0.25 g/g      OUTPUT: RHO = 1187.35 kg/m3\n
!>   INPUT:  TEMP = 393.15K  P =   30 MPa   S = 0.10 g/g      OUTPUT: RHO = 1027.06 kg/m3\n\n
!>
!> ARGUMENTS NAME    TYPE    UNITS           DESCRIPTION\n
!>     INPUT:        Temp  Real    T                 C       Temperature \n
!>                         Real    P         Pa      Pressure\n
!>                         Real    S         g/g     Salinity in mass fraction\n
!>    OUTPUT:        LABEL            RHO            Real     kg/m3  Density of brine \n
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
        double precision :: t

        ! Pressure (MPa)
        double precision :: p

        ! Salinity (mol/L)
        double precision :: s

        ! Salinity fraction (g/L)
        double precision :: sr

        ! Molar mass of NaCl [g/mol]
        ! double precision, parameter :: mmnacl = 58.44277d0
        double precision, parameter :: mmnacl = 58.44d0

        ! Pure water density (kg/m3)
        double precision :: rw
        double precision, external :: rhow

        ! Pure water density (g/cm3)
        double precision :: rw_gcm3

        ! Fluid density (g/cm3)
        double precision :: rhof_gcm3


        ! Local Temperature (degC)
        t = temp(i,j,k,ismpl)

        ! Local Pressure [MPa]
        p = pres(i,j,k,ismpl)*pa_conv1

        ! Local salinity [mol/L]
        s = tsal(i,j,k,ismpl)

        ! Pure water density [kg/m3]
        rw = rhow(p,t)


        if (s<=0.0d0) then

          rhof = rw

        else

          ! mol/L (Molarity) > g/g (Mass fraCtion)
          sr = s*mmnacl/(rw+s*mmnacl)

          ! Pure water density [g/cm3]
          rw_gcm3 = rw / 1.0d3

          ! Batzle, Equation (27b), densities in g/cm3
          rhof_gcm3 = rw_gcm3 + sr*(0.668d0 + 0.44d0*sr &
              + 1.0D-6*(3.0d2*p - 2.4d3*p*sr &
              + t * (80.0d0 + 3.0d0*t - 3.3d3*sr - 13.0d0*p + 47.0d0*p*sr)))

          ! Fluid density [kg/m3]
          rhof = rhof_gcm3 * 1.0d3
        end if

        return

      end function rhof
