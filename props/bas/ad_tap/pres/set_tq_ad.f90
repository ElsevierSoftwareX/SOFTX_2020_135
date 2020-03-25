!        Generated by TAPENADE     (INRIA, Ecuador team)
!  tapenade 3.x
!
!  Differentiation of set_tq in reverse (adjoint) mode (with options noISIZE i8):
!   gradient     of useful results: *w *propunit
!   with respect to varying inputs: *w *propunit
!   Plus diff mem management of: w:in propunit:in simtime:in
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
!>    @brief modify coefficents for the temperature equation
!>    @param[in] ismpl local sample index
!>    @details
!> modify coefficents for the temperature equation according to the prescribed sources and sinks\n
!> rhs stored in w.\n
SUBROUTINE SET_TQ_AD(ismpl)
  use arrays

  USE ARRAYS_AD

  USE MOD_GENRL
  use mod_time

  USE MOD_TIME_AD

  IMPLICIT NONE
  INTEGER :: ismpl
  INTEGER :: i, j, k
  EXTERNAL DELTAT, QT
  EXTERNAL QT_AD
  DOUBLE PRECISION :: DELTAT, deltt, QT
  DOUBLE PRECISION :: result1
  DOUBLE PRECISION :: result1_ad
! rhs: sources
  IF (transient .AND. tr_switch(ismpl)) THEN
    DO k=k0,1,-1
      DO j=j0,1,-1
        DO i=i0,1,-1
          result1_ad = -w_ad(i, j, k, ismpl)
          CALL QT_AD(i, j, k, ismpl, result1_ad)
        END DO
      END DO
    END DO
  ELSE
    DO k=k0,1,-1
      DO j=j0,1,-1
        DO i=i0,1,-1
          result1_ad = -w_ad(i, j, k, ismpl)
          CALL QT_AD(i, j, k, ismpl, result1_ad)
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE SET_TQ_AD

