!        Generated by TAPENADE     (INRIA, Ecuador team)
!  tapenade 3.x
!
!  Differentiation of rhoceff in reverse (adjoint) mode (with options noISIZE i8):
!   gradient     of useful results: *temp *propunit *pres rhoceff
!   with respect to varying inputs: *temp *propunit *pres
!   Plus diff mem management of: temp:in propunit:in pres:in
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
SUBROUTINE RHOCEFF_AD(i, j, k, ismpl, rhoceff_adv)
  use arrays

  USE ARRAYS_AD

  IMPLICIT NONE
  double precision :: rhoceff_adv
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
  DOUBLE PRECISION :: porlocal_ad
  DOUBLE PRECISION, EXTERNAL :: POR
! Matrix fraction in cell
  DOUBLE PRECISION :: fm
  DOUBLE PRECISION :: fm_ad
! Fluid fraction in cell
  DOUBLE PRECISION :: ff
  DOUBLE PRECISION :: ff_ad
! Heat capacity of the matrix
  DOUBLE PRECISION, EXTERNAL :: RHOCM
! Heat capacity of the fluid
  DOUBLE PRECISION, EXTERNAL :: RHOCF
  DOUBLE PRECISION :: result1
  DOUBLE PRECISION :: result1_ad
  DOUBLE PRECISION :: result2
  DOUBLE PRECISION :: result2_ad
  DOUBLE PRECISION :: rhoceff
! Local Temperature in degC
! Local porosity
  porlocal = POR(i, j, k, ismpl)
! Matrix fraction
  fm = 1.d0 - porlocal
! Fluid fraction
  ff = porlocal
! Heat capacity in cell, arithmetic mean
  result1 = RHOCF(i, j, k, ismpl)
  result2 = RHOCM(i, j, k, ismpl)
  ff_ad = result1*rhoceff_adv
  result1_ad = ff*rhoceff_adv
  fm_ad = result2*rhoceff_adv
  result2_ad = fm*rhoceff_adv
  CALL RHOCM_AD0(i, j, k, ismpl, result2_ad)
  CALL RHOCF_AD(i, j, k, ismpl, result1_ad)
  porlocal_ad = ff_ad - fm_ad
  CALL POR_AD0(i, j, k, ismpl, porlocal_ad)
END SUBROUTINE RHOCEFF_AD

