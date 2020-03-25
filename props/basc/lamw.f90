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

!> @brief calculate the thermal conductivity kf in W/(m*K) of water
!> @param[in] p cell pressure [MPa]
!> @param[in] t temperature [degC]
!> @return thermal conductivity [W/(m*K)]
!> @details
!> Calculate the thermal conductivity kf in W/(m*K) of freshwater,
!> given temperature in degC. Thermal conductivity of freshwater, kfw
!> is calculated using the Phillips (1981) formulation (page 8). \n\n
!>
!> Source:\n\n
!>
!> Phillips, S., Igbene, A., Fair, J., Ozbek, H., & Tavana, M.,
!> Technical databook for geothermal energy utilization (1981).
!> http://dx.doi.org/10.2172/6301274 \n\n
!>
!> Range of validity:  20 to 330 degC\n\n
!>
!>      temperature tlocal in [C]\n
      double precision function lamw(p,t)
        use arrays, only: temp

        implicit none

        ! Pressure p [MPa]
        double precision, intent (in) :: p

        ! Temperature t [degC]
        double precision, intent (in) :: t

        ! Temperature (degC)
        double precision :: tlocal

        ! Monomials of temperatures quotient
        double precision  :: tr, tr2, tr3, tr4

        ! Coefficients of approximation
        double precision, parameter :: c0 = -0.92247d0
        double precision, parameter :: c1 = 2.8395d0
        double precision, parameter :: c2 = 1.8007d0
        double precision, parameter :: c3 = 0.52577d0
        double precision, parameter :: c4 = 0.07344d0


        ! Monomials of temperature quotient
        tr = (t+273.15d0) / 273.15d0
        tr2 = tr*tr
        tr3 = tr2*tr
        tr4 = tr3*tr

        ! Thermal conductivity [W/(m*K)]
        lamw = (c0 + c1*tr - c2*tr2 + c3*tr3 - c4*tr4)

        return

      end function lamw
