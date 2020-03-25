!        Generated by TAPENADE     (INRIA, Ecuador team)
!  tapenade 3.x
!
!  Differentiation of buoy in reverse (adjoint) mode (with options noISIZE i8):
!   gradient     of useful results: *propunit buoy
!   with respect to varying inputs: *propunit
!   Plus diff mem management of: propunit:in
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
!>    @brief calculate buoyancy for head equation
!>    @param[in] i cell index, direction I0
!>    @param[in] j cell index, direction J0
!>    @param[in] k cell index, direction K0
!>    @param[in] ismpl local sample index
!>    @return buoyancy
!>    @details
!>calculate buoyancy for head equation\n
!>sign convention: negative for positive buoyancy\n
SUBROUTINE BUOY_AD(i, j, k, ismpl, buoy_adv)
  use arrays

  USE ARRAYS_AD

  USE MOD_GENRL
  USE MOD_FLOW
  IMPLICIT NONE
  double precision :: buoy_adv
  INTEGER :: ismpl
  INTEGER :: i, j, k
  DOUBLE PRECISION :: rhor, rhav, hh, h0, h1, prod, summ
  DOUBLE PRECISION :: hh_ad, h0_ad, h1_ad, prod_ad, summ_ad
  EXTERNAL RHOF, KZ, VISF
  EXTERNAL KZ_AD
  DOUBLE PRECISION :: RHOF, KZ, VISF
  INTEGER :: arg1
  DOUBLE PRECISION :: result1
  DOUBLE PRECISION :: result1_ad
  DOUBLE PRECISION :: result2
  DOUBLE PRECISION :: result3
  INTEGER :: arg2
  INTEGER :: arg3
  DOUBLE PRECISION :: temporary_ad
  INTEGER :: branch
  DOUBLE PRECISION :: buoy
  result1 = RHOF(i, j, arg1, ismpl)
  result2 = RHOF(i, j, k, ismpl)
  rhav = 0.5d0*(result1+result2)
  rhor = (rhav-rref)/rref
  result1 = KZ(i, j, k, ismpl)
  result2 = RHOF(i, j, k, ismpl)
  result3 = VISF(i, j, k, ismpl)
  h0 = result1*result2*grav/result3
  arg1 = k + 1
  result1 = KZ(i, j, arg1, ismpl)
  CALL PUSHREAL8(result2)
  result2 = RHOF(i, j, arg2, ismpl)
  CALL PUSHREAL8(result3)
  result3 = VISF(i, j, arg3, ismpl)
  h1 = result1*result2*grav/result3
  summ = h0 + h1
  prod = h0*h1
  IF (summ .GT. 0.d0) THEN
    CALL PUSHCONTROL1B(0)
  ELSE
    CALL PUSHCONTROL1B(1)
  END IF
  hh_ad = rhor*buoy_adv
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 0) THEN
    temporary_ad = 2.0d0*hh_ad/summ
    prod_ad = temporary_ad
    summ_ad = -(prod*temporary_ad/summ)
  ELSE
    prod_ad = 0.D0
    summ_ad = 0.D0
  END IF
  h0_ad = h1*prod_ad + summ_ad
  h1_ad = h0*prod_ad + summ_ad
  result1_ad = result2*grav*h1_ad/result3
  CALL POPREAL8(result3)
  CALL POPREAL8(result2)
  CALL KZ_AD(i, j, arg1, ismpl, result1_ad)
  result1_ad = result2*grav*h0_ad/result3
  CALL KZ_AD(i, j, k, ismpl, result1_ad)
END SUBROUTINE BUOY_AD

