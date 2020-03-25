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

!> @brief calculate the thermal conductivity kf of fluid [W/(m*K)]
!> @param[in] i cell index, direction I0
!> @param[in] j cell index, direction J0
!> @param[in] k cell index, direction K0
!> @param[in] ismpl local sample index
!> @return  thermal conductivity                lamf[W/(m*K)]
!> @details
!> calculate the thermal conductivity kf in W/(m*K) of saline water,
!> given temperature in degC, and salinity in mass fraction (g/g)of
!> NaCl. Thermal conductivity of freshwater, kfw is calculated using
!> the Phillips (1981) formulation (page8). \n \n
!>
!>      C = S./(1 + S)*1.d2;C2=C.*C;    % C=salinity in mol/kg \n\n
!>    kf = kfw.*(1.d0 - (2.3434d-3 - 7.924d-6*T + 3.924d-8*T2).*C ... \n
!>                         + (1.06d-5 - 2.d-8*T - 1.2d-10*T2).*C2) \n\n
!> Source:\n\n
!>
!> Phillips, S., Igbene, A., Fair, J., Ozbek, H., & Tavana, M.,
!> Technical databook for geothermal energy utilization (1981).
!> http://dx.doi.org/10.2172/6301274 \n\n
!>
!> Range of validity:  20 to 330degC and up to 4 molal NaCl\n\n
!> input:\n
!> pressure                             p [MPa]\n
!> temperature                          t in [C]\n
!> salinity                              s in [mol/L]\n
      double precision function lamf(i,j,k,ismpl)
        use arrays
        use mod_flow
        IMPLICIT NONE

        ! Location indices
        integer, intent (in) :: i
        integer, intent (in) :: j
        integer, intent (in) :: k

        ! Sample index
        integer :: ismpl

        ! Local Temperature (degC)
        double precision :: t

        ! Local Pressure (MPa)
        double precision :: p

        ! Local salinity [mol/kg / mol/L]
        double precision :: s

        ! Monomials of temperature
        double precision :: t2, t3, t4

        ! Salinity from Phillips1981 [-]
        double precision :: sr

        ! Monomial of salinity from Phillips1981 [-]
        double precision :: sr2

        ! Factor for thermal conductivity, Phillips1981 (2)
        double precision :: lamfac

        ! Pure water thermal conductivity [W/(m*K)]
        double precision, external :: lamw


        ! Local pressure [MPa]
        p = pres(i,j,k,ismpl)*pa_conv1

        ! Local temperature [degC]
        t = temp(i,j,k,ismpl)

        ! Local salinity [mol/L / mol/kg]
        s = tsal(i,j,k,ismpl)

        if (s<=0.0d0) then

          ! Pure water conductivity
          lamf = lamw(p,t)

        else

          ! Salinity according to Phillips (1981) between (2) and (3)
          sr =5844.3d0*s/(1.0d3+58.443d0*s)   

          ! Monomials in salinity and temperature
          sr2 = sr*sr
          t2 = t*t
          t3 = t2*t
          t4 = t3*t

          ! Factor lamf/lamw from Phillips1981, eq (2)
          lamfac = 1.0d0 - (2.3434d-3 - 7.924d-6*t + 3.924d-8*t2) * sr + &
              (1.06d-5 - 2.0d-8*t + 1.2d-10*t2) * sr2

          ! Thermal conductivity of fluid
          lamf = lamfac*lamw(p,t)

        end if

        return

      end function lamf
