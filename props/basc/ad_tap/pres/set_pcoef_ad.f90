!        Generated by TAPENADE     (INRIA, Ecuador team)
!  tapenade 3.x
!
!  Differentiation of set_pcoef in reverse (adjoint) mode (with options noISIZE i8):
!   gradient     of useful results: *d *e *f *g *temp *propunit
!                *tsal *pres *a *b *c
!   with respect to varying inputs: *d *e *f *g *temp *propunit
!                *tsal *pres *a *b *c
!   Plus diff mem management of: d:in e:in f:in g:in temp:in propunit:in
!                tsal:in pres:in a:in b:in c:in
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
!>    @brief calculate coefficents for the head equation
!>    @param[in] ismpl local sample index
!>    @details
!> calculate coefficents for the head equation\n
!> coefficients are stored as vectors in the diagonals a-g (d center) and rhs in w.\n
SUBROUTINE SET_PCOEF_AD(ismpl)
  use arrays

  USE ARRAYS_AD

  USE MOD_GENRL
  USE MOD_GENRLC
  USE MOD_FLOW
  use mod_time

  USE MOD_TIME_AD

  USE MOD_LINFOS
  IMPLICIT NONE
  INTEGER :: ismpl
  INTEGER :: i, j, k
  EXTERNAL FI, FJ, FK
  EXTERNAL FI_AD, FJ_AD, FK_AD
  DOUBLE PRECISION :: FI, FJ, FK
  DOUBLE PRECISION :: result1
  DOUBLE PRECISION :: result1_ad
  INTEGER :: arg1
  INTEGER :: branch
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  DO k=1,k0
    DO j=1,j0
      DO i=1,i0
        IF (i0 .GT. 1) THEN
          IF (i .LT. i0) THEN
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
          IF (i .GT. 1) THEN
            CALL PUSHCONTROL2B(0)
          ELSE
            CALL PUSHCONTROL2B(1)
          END IF
        ELSE
          CALL PUSHCONTROL2B(2)
        END IF
        IF (j0 .GT. 1) THEN
          IF (j .LT. j0) THEN
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
          IF (j .GT. 1) THEN
            CALL PUSHCONTROL2B(0)
          ELSE
            CALL PUSHCONTROL2B(1)
          END IF
        ELSE
          CALL PUSHCONTROL2B(2)
        END IF
        IF (k0 .GT. 1) THEN
          IF (k .LT. k0) THEN
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
          IF (k .GT. 1) THEN
            CALL PUSHCONTROL2B(0)
          ELSE
            CALL PUSHCONTROL2B(1)
          END IF
        ELSE
          CALL PUSHCONTROL2B(2)
        END IF
      END DO
    END DO
  END DO
  DO k=k0,1,-1
    DO j=j0,1,-1
      DO 120 i=i0,1,-1
        e_ad(i, j, k, ismpl) = e_ad(i, j, k, ismpl) - d_ad(i, j, k, &
&         ismpl)
        c_ad(i, j, k, ismpl) = c_ad(i, j, k, ismpl) - d_ad(i, j, k, &
&         ismpl)
        f_ad(i, j, k, ismpl) = f_ad(i, j, k, ismpl) - d_ad(i, j, k, &
&         ismpl)
        b_ad(i, j, k, ismpl) = b_ad(i, j, k, ismpl) - d_ad(i, j, k, &
&         ismpl)
        g_ad(i, j, k, ismpl) = g_ad(i, j, k, ismpl) - d_ad(i, j, k, &
&         ismpl)
        a_ad(i, j, k, ismpl) = a_ad(i, j, k, ismpl) - d_ad(i, j, k, &
&         ismpl)
        d_ad(i, j, k, ismpl) = 0.D0
        CALL POPCONTROL2B(branch)
        IF (branch .EQ. 0) THEN
          result1_ad = a_ad(i, j, k, ismpl)/delz(k)
          a_ad(i, j, k, ismpl) = 0.D0
          arg1 = k - 1
          CALL FK_AD(i, j, arg1, ismpl, result1_ad)
        ELSE IF (branch .NE. 1) THEN
          GOTO 100
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          result1_ad = g_ad(i, j, k, ismpl)/delz(k)
          g_ad(i, j, k, ismpl) = 0.D0
          CALL FK_AD(i, j, k, ismpl, result1_ad)
        END IF
 100    CALL POPCONTROL2B(branch)
        IF (branch .EQ. 0) THEN
          result1_ad = b_ad(i, j, k, ismpl)/dely(j)
          b_ad(i, j, k, ismpl) = 0.D0
          arg1 = j - 1
          CALL FJ_AD(i, arg1, k, ismpl, result1_ad)
        ELSE IF (branch .NE. 1) THEN
          GOTO 110
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          result1_ad = f_ad(i, j, k, ismpl)/dely(j)
          f_ad(i, j, k, ismpl) = 0.D0
          CALL FJ_AD(i, j, k, ismpl, result1_ad)
        END IF
 110    CALL POPCONTROL2B(branch)
        IF (branch .EQ. 0) THEN
          result1_ad = c_ad(i, j, k, ismpl)/delx(i)
          c_ad(i, j, k, ismpl) = 0.D0
          arg1 = i - 1
          CALL FI_AD(arg1, j, k, ismpl, result1_ad)
        ELSE IF (branch .NE. 1) THEN
          GOTO 120
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          result1_ad = e_ad(i, j, k, ismpl)/delx(i)
          e_ad(i, j, k, ismpl) = 0.D0
          CALL FI_AD(i, j, k, ismpl, result1_ad)
        END IF
 120  CONTINUE
    END DO
  END DO
END SUBROUTINE SET_PCOEF_AD

