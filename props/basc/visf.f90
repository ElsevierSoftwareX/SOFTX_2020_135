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

!> @brief calculates viscosity of aqueous NaCl solutions
!> @param[in] i cell index, direction I0
!> @param[in] j cell index, direction J0
!> @param[in] k cell index, direction K0
!> @param[in] ismpl local sample index
!> @return  visc
!> @details
!> Source:\n
!> Batzle, M., & Wang, Z., Seismic properties of pore fluids,
!> GEOPHYSICS, 57(11), 1396–1408 (1992).
!> http://dx.doi.org/10.1190/1.1443207 \n\n
!>\n
!>   Pressures 5-100 MPa, Temperature 20-350°C, Salinity <=320 g/L\n
!>\n
      double precision function visf(i,j,k,ismpl)
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

        ! Mass fraction [-]
        double precision :: sr

        ! Pure water density (kg/m3)
        double precision :: rw
        double precision, external :: rhow

        ! Pure water viscosity [Pa s]
        double precision, external :: visw

        ! Molar mass of NaCl [g/mol]
        ! double precision, parameter :: mmnacl = 58.44277d0
        double precision, parameter :: mmnacl = 58.44d0

        ! Viscosity [cP]
        double precision :: visf_cp


        ! Pressure [MPa]
        p = pres(i,j,k,ismpl)*pa_conv1

        ! Temperature [degC]
        t = temp(i,j,k,ismpl)

        ! Salinity [mol/L]
        s = tsal(i,j,k,ismpl)

        
        if (s <= 0.D0) then

          visf = visw(p,t)
          
        else

          ! Pure water density [kg/m3]
          rw = rhow(p,t)
          
          ! mol/L (Molarity) > [g/L] = [kg/m3; ]g/g (mass fraction)      
          sr = s*mmnacl/(rw+s*mmnacl)

          ! Viscosity formula after Batzle & Wang [cP] 
          visf_cp = 0.1d0 + 0.333d0 * sr  &
              + (1.65d0 + 91.9d0 * sr**3) * &
              exp(-(0.42d0 * (sr**(0.8d0) - 0.17d0)**2 + 0.045d0) * t**(0.8d0))

          ! Conversion of viscosity from [cP] to [Pa s]
          visf = visf_cp * 1.0d-3
          
        end if

        return

      end function visf
