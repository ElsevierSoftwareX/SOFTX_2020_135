!        Generated by TAPENADE     (INRIA, Ecuador team)
!  tapenade 3.x
!
!  Differentiation of rhocm in forward (tangent) mode:
!   variations   of useful results: rhocm
!   with respect to varying inputs: *temp *propunit
!   Plus diff mem management of: temp:in propunit:in
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
!> @brief calculates  heat capacity*density of rock.
!> @param[in] i cell index, direction I0
!> @param[in] j cell index, direction J0
!> @param[in] k cell index, direction K0
!> @param[in] ismpl local sample index
!> @return  rhoc [W/(m*K)]
!> @details
!> temperature tlocal in [C]\n\n
!>
!> Under input file "# rhocm", the temperature variation coefficients
!> cma1, cma2, cma3 can be set. \n
!> Default: cma1 = 1.0d0, cma2 = cma3 = 0.0d0
DOUBLE PRECISION FUNCTION g_RHOCM(i, j, k, ismpl, rhocm)
  USE ARRAYS

  USE g_ARRAYS

  USE MOD_TEMP
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
  DOUBLE PRECISION :: temp0
  DOUBLE PRECISION :: temp1
  DOUBLE PRECISION :: rhocm
! Local Temperature in degC
  g_tlocal = g_temp(i, j, k, ismpl)
  tlocal = temp(i, j, k, ismpl)
! Volumetric heat capacity from input file [J/(kg*m3)]
  temp0 = cma1 + cma2*tlocal + cma3*(tlocal*tlocal)
  temp1 = propunit(uindex(i, j, k), idx_rc, ismpl)
  g_rhocm = temp0*g_propunit(uindex(i, j, k), idx_rc, ismpl) + temp1&
&   *(cma2+cma3*2*tlocal)*g_tlocal
  rhocm = temp1*temp0
  RETURN
END FUNCTION g_RHOCM

