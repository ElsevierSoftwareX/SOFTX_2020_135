!        Generated by TAPENADE     (INRIA, Ecuador team)
!  tapenade 3.x
!
!  Differentiation of rhoceff in forward (tangent) mode:
!   variations   of useful results: rhoceff
!   with respect to varying inputs: *temp *propunit *tsal *pres
!   Plus diff mem management of: temp:in propunit:in tsal:in pres:in
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
!> @brief calculates volumetric heat capacity of the cell
!> @param[in] i cell index, direction I0
!> @param[in] j cell index, direction J0
!> @param[in] k cell index, direction K0
!> @param[in] ismpl local sample index
!> @return volumetric heat capacity
!> @details
!> calculates volumetric heat capacity of the system
!> matrix-porosity [J/(K*m3)].\n
DOUBLE PRECISION FUNCTION g_RHOCEFF(i, j, k, ismpl, rhoceff)
  USE ARRAYS

  USE g_ARRAYS

  IMPLICIT NONE
! Location indices
  INTEGER, INTENT(IN) :: i
  INTEGER, INTENT(IN) :: j
  INTEGER, INTENT(IN) :: k
! Sample index
  INTEGER :: ismpl
! Local temperature [degC]
  DOUBLE PRECISION :: tlocal
! Local porosity [-]
  DOUBLE PRECISION :: porlocal
  DOUBLE PRECISION :: g_porlocal
  DOUBLE PRECISION, EXTERNAL :: POR
  DOUBLE PRECISION, EXTERNAL :: g_POR
! Matrix fraction in cell
  DOUBLE PRECISION :: fm
  DOUBLE PRECISION :: g_fm
! Fluid fraction in cell
  DOUBLE PRECISION :: ff
  DOUBLE PRECISION :: g_ff
! Heat capacity of the matrix
  DOUBLE PRECISION, EXTERNAL :: RHOCM
! Heat capacity of the fluid
  DOUBLE PRECISION, EXTERNAL :: RHOCF
	double precision, external :: g_rhocm
	double precision, external :: g_rhocf
  DOUBLE PRECISION :: result1
  DOUBLE PRECISION :: g_result1
  DOUBLE PRECISION :: result2
  DOUBLE PRECISION :: g_result2
  DOUBLE PRECISION :: rhoceff
! Local Temperature in degC
  tlocal = temp(i, j, k, ismpl)
! Local porosity
  g_porlocal = g_POR(i, j, k, ismpl, porlocal)
! Matrix fraction
  g_fm = -g_porlocal
  fm = 1.d0 - porlocal
! Fluid fraction
  g_ff = g_porlocal
  ff = porlocal
! Heat capacity in cell, arithmetic mean
  g_result1 = g_RHOCF(i, j, k, ismpl, result1)
  g_result2 = g_RHOCM(i, j, k, ismpl, result2)
  g_rhoceff = result1*g_ff + ff*g_result1 + result2*g_fm + fm*&
&   g_result2
  rhoceff = ff*result1 + fm*result2
  RETURN
END FUNCTION g_RHOCEFF

