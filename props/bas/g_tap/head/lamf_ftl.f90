!        Generated by TAPENADE     (INRIA, Ecuador team)
!  tapenade 3.x
!
!  Differentiation of lamf in forward (tangent) mode:
!   variations   of useful results: lamf
!   with respect to varying inputs: *temp
!   Plus diff mem management of: temp:in
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
!> @param[in] i cell index, direction I0
!> @param[in] j cell index, direction J0
!> @param[in] k cell index, direction K0
!> @param[in] ismpl local sample index
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
DOUBLE PRECISION FUNCTION g_LAMF(i, j, k, ismpl, lamf)
  USE ARRAYS

  USE g_ARRAYS

  IMPLICIT NONE
! Location indices
  INTEGER, INTENT(IN) :: i
  INTEGER, INTENT(IN) :: j
  INTEGER, INTENT(IN) :: k
! Sample index
  INTEGER :: ismpl
! Temperature (degC)
  DOUBLE PRECISION :: tlocal
  DOUBLE PRECISION :: g_tlocal
! Monomials of temperatures quotient
  DOUBLE PRECISION :: tr, tr2, tr3, tr4
  DOUBLE PRECISION :: g_tr, g_tr2, g_tr3, g_tr4
! Coefficients of approximation
  DOUBLE PRECISION, PARAMETER :: c0=-0.92247d0
  DOUBLE PRECISION, PARAMETER :: c1=2.8395d0
  DOUBLE PRECISION, PARAMETER :: c2=1.8007d0
  DOUBLE PRECISION, PARAMETER :: c3=0.52577d0
  DOUBLE PRECISION, PARAMETER :: c4=0.07344d0
  DOUBLE PRECISION :: lamf
! Local Temperature in degC
  g_tlocal = g_temp(i, j, k, ismpl)
  tlocal = temp(i, j, k, ismpl)
! Monomials of temperature quotient
  g_tr = g_tlocal/273.15d0
  tr = (tlocal+273.15d0)/273.15d0
  g_tr2 = 2*tr*g_tr
  tr2 = tr*tr
  g_tr3 = tr*g_tr2 + tr2*g_tr
  tr3 = tr2*tr
  g_tr4 = tr*g_tr3 + tr3*g_tr
  tr4 = tr3*tr
! Thermal conductivity [W/(m*K)]
  g_lamf = c1*g_tr - c2*g_tr2 + c3*g_tr3 - c4*g_tr4
  lamf = c0 + c1*tr - c2*tr2 + c3*tr3 - c4*tr4
  RETURN
END FUNCTION g_LAMF

