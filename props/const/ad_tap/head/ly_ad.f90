!        Generated by TAPENADE     (INRIA, Ecuador team)
!  tapenade 3.x
!
!  Differentiation of ly in reverse (adjoint) mode (with options noISIZE i8):
!   gradient     of useful results: *propunit ly
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
!> @brief calculates effective thermal conductivity of the cell
!> @param[in] i cell index, direction I0
!> @param[in] j cell index, direction J0
!> @param[in] k cell index, direction K0
!> @param[in] ismpl local sample index
!> @return  thermal conductivity ly[W/(m*K)]
!> @details
!> calculates effective thermal conductivity of the two phase system
!> matrix-porosity .\n\n
!>
!> input:\n
!> porosity                            porlocal [-]\n
!> thermal conductivity of fluid           lamf [W/(m*K)]\n
!> thermal conductivity of matrix       lammref [W/(m*K)]\n
SUBROUTINE LY_AD(i, j, k, ismpl, ly_adv)
  use arrays

  USE ARRAYS_AD

  IMPLICIT NONE
  double precision :: ly_adv
! Location indices
  INTEGER, INTENT(IN) :: i
  INTEGER, INTENT(IN) :: j
  INTEGER, INTENT(IN) :: k
! Sample index
  INTEGER :: ismpl
! Local fluid thermal conductivity [W/(m*K)]
  DOUBLE PRECISION :: lamfluid
  DOUBLE PRECISION, EXTERNAL :: LAMF
! Matrix thermal conductivity [W/(m*K)]
  DOUBLE PRECISION :: lammref
  DOUBLE PRECISION :: lammref_ad
! Local porosity [-]
  DOUBLE PRECISION, EXTERNAL :: POR
  DOUBLE PRECISION :: result1
  DOUBLE PRECISION :: result1_ad
  DOUBLE PRECISION :: result2
  DOUBLE PRECISION :: result2_ad
  DOUBLE PRECISION :: temp0
  DOUBLE PRECISION :: temp1
  DOUBLE PRECISION :: ly
! Local fluid thermal conductivity [W/(m*K)]
  lamfluid = LAMF(i, j, k, ismpl)
! Reference matrix thermal conductivity [W/(m*K)]
  lammref = propunit(uindex(i, j, k), idx_lz, ismpl)*propunit(uindex(i, &
&   j, k), idx_an_ly, ismpl)
  result1 = POR(i, j, k, ismpl)
  result2 = POR(i, j, k, ismpl)
  temp0 = lamfluid**result2
  temp1 = lammref**(-result1+1.0d0)
  IF (lammref .LE. 0.0 .AND. (1.0d0 - result1 .EQ. 0.0 .OR. 1.0d0 - &
&     result1 .NE. INT(1.0d0 - result1))) THEN
    lammref_ad = 0.D0
  ELSE
    lammref_ad = (1.0d0-result1)*lammref**(-result1)*temp0*ly_adv
  END IF
  IF (lammref .LE. 0.0) THEN
    result1_ad = 0.D0
  ELSE
    result1_ad = -(temp1*LOG(lammref)*temp0*ly_adv)
  END IF
  IF (lamfluid .LE. 0.0) THEN
    result2_ad = 0.D0
  ELSE
    result2_ad = temp0*LOG(lamfluid)*temp1*ly_adv
  END IF
  CALL POR_AD0(i, j, k, ismpl, result2_ad)
  CALL POR_AD0(i, j, k, ismpl, result1_ad)
  propunit_ad(uindex(i, j, k), idx_lz, ismpl) = propunit_ad(uindex(i, j&
&   , k), idx_lz, ismpl) + propunit(uindex(i, j, k), idx_an_ly, ismpl)*&
&   lammref_ad
  propunit_ad(uindex(i, j, k), idx_an_ly, ismpl) = propunit_ad(uindex(i&
&   , j, k), idx_an_ly, ismpl) + propunit(uindex(i, j, k), idx_lz, ismpl&
&   )*lammref_ad
END SUBROUTINE LY_AD

