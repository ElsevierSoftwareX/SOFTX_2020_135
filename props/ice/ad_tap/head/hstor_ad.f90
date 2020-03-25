!        Generated by TAPENADE     (INRIA, Ecuador team)
!  tapenade 3.x
!
!  Differentiation of hstor in reverse (adjoint) mode (with options noISIZE i8):
!   gradient     of useful results: *temp *propunit *pres hstor
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
!>    @brief calculates the bulk storativity
!>    @param[in] i cell index, direction I0
!>    @param[in] j cell index, direction J0
!>    @param[in] k cell index, direction K0
!>    @param[in] ismpl local sample index
!>    @return bulk storativity
!>    @details
!>    storb(i,j,k,ismpl) calculates the bulk storativity \n
!>    at node(i,j,k).\n
SUBROUTINE HSTOR_AD0(i, j, k, ismpl, hstor_ad)
  use arrays

  USE ARRAYS_AD

  USE MOD_FLOW
  IMPLICIT NONE
  INTEGER :: i, j, k, ismpl
  EXTERNAL RHOF, COMPM, COMPF, POR
  EXTERNAL RHOF_AD, COMPM_AD, COMPF_AD0, POR_AD0
  DOUBLE PRECISION :: RHOF, COMPM, COMPF, &
& POR
  DOUBLE PRECISION :: result1
  DOUBLE PRECISION :: result1_ad
  DOUBLE PRECISION :: result2
  DOUBLE PRECISION :: result2_ad
  DOUBLE PRECISION :: result3
  DOUBLE PRECISION :: result3_ad
  DOUBLE PRECISION :: result4
  DOUBLE PRECISION :: result4_ad
  DOUBLE PRECISION :: temporary_ad
  DOUBLE PRECISION :: hstor_ad
  DOUBLE PRECISION :: hstor
  result1 = RHOF(i, j, k, ismpl)
  result2 = COMPM(i, j, k, ismpl)
  result3 = POR(i, j, k, ismpl)
  result4 = COMPF(i, j, k, ismpl)
  result1_ad = (result2+result3*result4)*grav*hstor_ad
  temporary_ad = result1*grav*hstor_ad
  result2_ad = temporary_ad
  result3_ad = result4*temporary_ad
  result4_ad = result3*temporary_ad
  CALL COMPF_AD0(i, j, k, ismpl, result4_ad)
  CALL POR_AD0(i, j, k, ismpl, result3_ad)
  CALL COMPM_AD(i, j, k, ismpl, result2_ad)
  CALL RHOF_AD(i, j, k, ismpl, result1_ad)
END SUBROUTINE HSTOR_AD0

