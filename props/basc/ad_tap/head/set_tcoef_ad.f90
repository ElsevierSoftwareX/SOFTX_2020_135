!        Generated by TAPENADE     (INRIA, Ecuador team)
!  tapenade 3.x
!
!  Differentiation of set_tcoef in reverse (adjoint) mode (with options noISIZE i8):
!   gradient     of useful results: *d *e *f *g *temp *head *propunit
!                *tsal *pres *a *b *c
!   with respect to varying inputs: *d *e *f *g *temp *head *propunit
!                *tsal *pres *a *b *c
!   Plus diff mem management of: d:in e:in f:in g:in temp:in head:in
!                propunit:in tsal:in pres:in a:in b:in c:in
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
!>    @brief calculate coefficents for the heat equation
!>    @param[in] ismpl local sample index
!>    @details
!> calculate coefficents for the heat equation\n
!> coefficients are stored as vectors in the diagonals a-g (d center) and rhs in w.\n
SUBROUTINE SET_TCOEF_AD(ismpl)
  use arrays

  USE ARRAYS_AD

  USE MOD_GENRL
  USE MOD_GENRLC
  USE MOD_TEMP
  use mod_time

  USE MOD_TIME_AD

  USE MOD_LINFOS
  IMPLICIT NONE
  INTEGER :: ismpl
  INTEGER :: i, j, k
  EXTERNAL LI, LJ, LK, RHOCF, VX, &
&     VY, VZ, ALFA, AMEAN
  EXTERNAL LI_AD, LJ_AD, LK_AD, RHOCF_AD, VX_AD, VY_AD, VZ_AD, ALFA_AD, &
&     AMEAN_AD
  DOUBLE PRECISION :: LI, LJ, LK, RHOCF, &
& VX, VY, VZ, ALFA, AMEAN
  DOUBLE PRECISION :: rijk, ra, rc, va, la, alf, p2
  DOUBLE PRECISION :: rijk_ad, ra_ad, rc_ad, va_ad, la_ad, alf_ad, p2_ad
  INTEGER :: arg1
  DOUBLE PRECISION :: result1
  DOUBLE PRECISION :: result1_ad
  DOUBLE PRECISION :: temporary_ad
  INTEGER :: branch
! initialize coefficients for sparse solvers
! inner points of grid - - - - - - - - - - - - - - - - - - - - - - - - -
  DO k=1,k0
    DO j=1,j0
      DO i=1,i0
        rijk = RHOCF(i, j, k, ismpl)
        IF (i0 .GT. 1) THEN
          IF (i .LT. i0) THEN
            CALL PUSHREAL8(la)
            la = LI(i, j, k, ismpl)
            arg1 = i + 1
            result1 = RHOCF(arg1, j, k, ismpl)
            CALL PUSHREAL8(ra)
            ra = AMEAN(result1, rijk)
            CALL PUSHREAL8(va)
            va = VX(i, j, k, ismpl)
            rc = 0.5*ra*va
            IF (la .GT. 0.d0) THEN
              p2 = rc/la
              CALL PUSHREAL8(alf)
              alf = ALFA(p2)
              CALL PUSHCONTROL1B(0)
            ELSE
              CALL PUSHREAL8(alf)
              alf = 0.d0
              IF (va .LT. 0.d0) alf = -1.d0
              IF (va .GT. 0.d0) THEN
                CALL PUSHCONTROL1B(1)
                alf = 1.d0
              ELSE
                CALL PUSHCONTROL1B(1)
              END IF
            END IF
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
          IF (i .GT. 1) THEN
            arg1 = i - 1
            CALL PUSHREAL8(la)
            la = LI(arg1, j, k, ismpl)
            arg1 = i - 1
            result1 = RHOCF(arg1, j, k, ismpl)
            CALL PUSHREAL8(ra)
            ra = AMEAN(result1, rijk)
            arg1 = i - 1
            CALL PUSHREAL8(va)
            va = VX(arg1, j, k, ismpl)
            rc = 0.5*ra*va
            CALL PUSHREAL8(alf)
            alf = 0.d0
            IF (va .EQ. 0.d0) THEN
              CALL PUSHCONTROL2B(0)
              alf = 0.d0
            ELSE IF (la .GT. 0.d0) THEN
              p2 = rc/la
              alf = ALFA(p2)
              CALL PUSHCONTROL2B(1)
            ELSE
              IF (va .LT. 0.d0) alf = -1.d0
              IF (va .GT. 0.d0) THEN
                CALL PUSHCONTROL2B(2)
                alf = 1.d0
              ELSE
                CALL PUSHCONTROL2B(2)
              END IF
            END IF
            CALL PUSHCONTROL2B(0)
          ELSE
            CALL PUSHCONTROL2B(1)
          END IF
        ELSE
          CALL PUSHCONTROL2B(2)
        END IF
        IF (j0 .GT. 1) THEN
          IF (j .LT. j0) THEN
            CALL PUSHREAL8(la)
            la = LJ(i, j, k, ismpl)
            arg1 = j + 1
            result1 = RHOCF(i, arg1, k, ismpl)
            CALL PUSHREAL8(ra)
            ra = AMEAN(result1, rijk)
            CALL PUSHREAL8(va)
            va = VY(i, j, k, ismpl)
            rc = 0.5d0*ra*va
            CALL PUSHREAL8(alf)
            alf = 0.d0
            IF (va .EQ. 0.d0) THEN
              CALL PUSHCONTROL2B(0)
              alf = 0.d0
            ELSE IF (la .GT. 0.d0) THEN
              p2 = rc/la
              alf = ALFA(p2)
              CALL PUSHCONTROL2B(1)
            ELSE
              IF (va .LT. 0.d0) alf = -1.d0
              IF (va .GT. 0.d0) THEN
                CALL PUSHCONTROL2B(2)
                alf = 1.d0
              ELSE
                CALL PUSHCONTROL2B(2)
              END IF
            END IF
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
          IF (j .GT. 1) THEN
            arg1 = j - 1
            CALL PUSHREAL8(la)
            la = LJ(i, arg1, k, ismpl)
            arg1 = j - 1
            result1 = RHOCF(i, arg1, k, ismpl)
            CALL PUSHREAL8(ra)
            ra = AMEAN(result1, rijk)
            arg1 = j - 1
            CALL PUSHREAL8(va)
            va = VY(i, arg1, k, ismpl)
            rc = 0.5*ra*va
            CALL PUSHREAL8(alf)
            alf = 0.d0
            IF (va .EQ. 0.d0) THEN
              CALL PUSHCONTROL2B(0)
              alf = 0.d0
            ELSE IF (la .GT. 0.d0) THEN
              p2 = rc/la
              alf = ALFA(p2)
              CALL PUSHCONTROL2B(1)
            ELSE
              IF (va .LT. 0.d0) alf = -1.d0
              IF (va .GT. 0.d0) THEN
                CALL PUSHCONTROL2B(2)
                alf = 1.d0
              ELSE
                CALL PUSHCONTROL2B(2)
              END IF
            END IF
            CALL PUSHCONTROL2B(0)
          ELSE
            CALL PUSHCONTROL2B(1)
          END IF
        ELSE
          CALL PUSHCONTROL2B(2)
        END IF
        IF (k0 .GT. 1) THEN
          IF (k .LT. k0) THEN
            CALL PUSHREAL8(la)
            la = LK(i, j, k, ismpl)
            arg1 = k + 1
            result1 = RHOCF(i, j, arg1, ismpl)
            CALL PUSHREAL8(ra)
            ra = AMEAN(result1, rijk)
            CALL PUSHREAL8(va)
            va = VZ(i, j, k, ismpl)
            rc = 0.5*ra*va
            CALL PUSHREAL8(alf)
            alf = 0.d0
            IF (va .EQ. 0.d0) THEN
              CALL PUSHCONTROL2B(0)
              alf = 0.d0
            ELSE IF (la .GT. 0.d0) THEN
              p2 = rc/la
              alf = ALFA(p2)
              CALL PUSHCONTROL2B(1)
            ELSE
              IF (va .LT. 0.d0) alf = -1.d0
              IF (va .GT. 0.d0) THEN
                CALL PUSHCONTROL2B(2)
                alf = 1.d0
              ELSE
                CALL PUSHCONTROL2B(2)
              END IF
            END IF
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHCONTROL1B(1)
          END IF
          IF (k .GT. 1) THEN
            arg1 = k - 1
            CALL PUSHREAL8(la)
            la = LK(i, j, arg1, ismpl)
            arg1 = k - 1
            result1 = RHOCF(i, j, arg1, ismpl)
            CALL PUSHREAL8(ra)
            ra = AMEAN(result1, rijk)
            arg1 = k - 1
            CALL PUSHREAL8(va)
            va = VZ(i, j, arg1, ismpl)
            rc = 0.5*ra*va
            CALL PUSHREAL8(alf)
            alf = 0.d0
            IF (va .EQ. 0.d0) THEN
              CALL PUSHCONTROL2B(0)
              alf = 0.d0
            ELSE IF (la .GT. 0.d0) THEN
              p2 = rc/la
              alf = ALFA(p2)
              CALL PUSHCONTROL2B(1)
            ELSE
              IF (va .LT. 0.d0) alf = -1.d0
              IF (va .GT. 0.d0) THEN
                CALL PUSHCONTROL2B(2)
                alf = 1.d0
              ELSE
                CALL PUSHCONTROL2B(2)
              END IF
            END IF
            CALL PUSHCONTROL2B(2)
          ELSE
            CALL PUSHCONTROL2B(1)
          END IF
        ELSE
          CALL PUSHCONTROL2B(0)
        END IF
      END DO
    END DO
  END DO
  DO k=k0,1,-1
    DO j=j0,1,-1
      DO i=i0,1,-1
        CALL POPCONTROL2B(branch)
        IF (branch .EQ. 0) THEN
          rijk_ad = 0.D0
        ELSE
          IF (branch .EQ. 1) THEN
            rijk_ad = 0.D0
          ELSE
            rc = 0.5*ra*va
            temporary_ad = -(d_ad(i, j, k, ismpl)/delz(k))
            la_ad = temporary_ad
            alf_ad = rc*temporary_ad
            rc_ad = -((1.d0-alf)*temporary_ad)
            temporary_ad = a_ad(i, j, k, ismpl)/delz(k)
            a_ad(i, j, k, ismpl) = 0.D0
            la_ad = la_ad + temporary_ad
            alf_ad = alf_ad + rc*temporary_ad
            rc_ad = rc_ad + (alf+1.d0)*temporary_ad
            CALL POPCONTROL2B(branch)
            IF (branch .NE. 0) THEN
              IF (branch .EQ. 1) THEN
                p2 = rc/la
                CALL ALFA_AD(p2, p2_ad, alf_ad)
                rc_ad = rc_ad + p2_ad/la
                la_ad = la_ad - rc*p2_ad/la**2
              END IF
            END IF
            CALL POPREAL8(alf)
            ra_ad = va*0.5*rc_ad
            va_ad = ra*0.5*rc_ad
            arg1 = k - 1
            CALL POPREAL8(va)
            CALL VZ_AD(i, j, arg1, ismpl, va_ad)
            CALL POPREAL8(ra)
            rijk_ad = 0.D0
            CALL AMEAN_AD(result1, result1_ad, rijk, rijk_ad, ra_ad)
            arg1 = k - 1
            CALL RHOCF_AD(i, j, arg1, ismpl, result1_ad)
            arg1 = k - 1
            CALL POPREAL8(la)
            CALL LK_AD(i, j, arg1, ismpl, la_ad)
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            rc = 0.5*ra*va
            temporary_ad = -(d_ad(i, j, k, ismpl)/delz(k))
            la_ad = temporary_ad
            alf_ad = rc*temporary_ad
            rc_ad = (alf+1.d0)*temporary_ad
            temporary_ad = g_ad(i, j, k, ismpl)/delz(k)
            g_ad(i, j, k, ismpl) = 0.D0
            la_ad = la_ad + temporary_ad
            alf_ad = alf_ad + rc*temporary_ad
            rc_ad = rc_ad - (1.d0-alf)*temporary_ad
            CALL POPCONTROL2B(branch)
            IF (branch .NE. 0) THEN
              IF (branch .EQ. 1) THEN
                p2 = rc/la
                CALL ALFA_AD(p2, p2_ad, alf_ad)
                rc_ad = rc_ad + p2_ad/la
                la_ad = la_ad - rc*p2_ad/la**2
              END IF
            END IF
            CALL POPREAL8(alf)
            ra_ad = va*0.5*rc_ad
            va_ad = ra*0.5*rc_ad
            CALL POPREAL8(va)
            CALL VZ_AD(i, j, k, ismpl, va_ad)
            CALL POPREAL8(ra)
            CALL AMEAN_AD(result1, result1_ad, rijk, rijk_ad, ra_ad)
            arg1 = k + 1
            CALL RHOCF_AD(i, j, arg1, ismpl, result1_ad)
            CALL POPREAL8(la)
            CALL LK_AD(i, j, k, ismpl, la_ad)
          END IF
        END IF
        CALL POPCONTROL2B(branch)
        IF (branch .EQ. 0) THEN
          rc = 0.5*ra*va
          temporary_ad = -(d_ad(i, j, k, ismpl)/dely(j))
          la_ad = temporary_ad
          alf_ad = rc*temporary_ad
          rc_ad = -((1.d0-alf)*temporary_ad)
          temporary_ad = b_ad(i, j, k, ismpl)/dely(j)
          b_ad(i, j, k, ismpl) = 0.D0
          la_ad = la_ad + temporary_ad
          alf_ad = alf_ad + rc*temporary_ad
          rc_ad = rc_ad + (alf+1.d0)*temporary_ad
          CALL POPCONTROL2B(branch)
          IF (branch .NE. 0) THEN
            IF (branch .EQ. 1) THEN
              p2 = rc/la
              CALL ALFA_AD(p2, p2_ad, alf_ad)
              rc_ad = rc_ad + p2_ad/la
              la_ad = la_ad - rc*p2_ad/la**2
            END IF
          END IF
          CALL POPREAL8(alf)
          ra_ad = va*0.5*rc_ad
          va_ad = ra*0.5*rc_ad
          arg1 = j - 1
          CALL POPREAL8(va)
          CALL VY_AD(i, arg1, k, ismpl, va_ad)
          CALL POPREAL8(ra)
          CALL AMEAN_AD(result1, result1_ad, rijk, rijk_ad, ra_ad)
          arg1 = j - 1
          CALL RHOCF_AD(i, arg1, k, ismpl, result1_ad)
          arg1 = j - 1
          CALL POPREAL8(la)
          CALL LJ_AD(i, arg1, k, ismpl, la_ad)
        ELSE IF (branch .NE. 1) THEN
          GOTO 100
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          rc = 0.5d0*ra*va
          temporary_ad = -(d_ad(i, j, k, ismpl)/dely(j))
          la_ad = temporary_ad
          alf_ad = rc*temporary_ad
          rc_ad = (alf+1.d0)*temporary_ad
          temporary_ad = f_ad(i, j, k, ismpl)/dely(j)
          f_ad(i, j, k, ismpl) = 0.D0
          la_ad = la_ad + temporary_ad
          alf_ad = alf_ad + rc*temporary_ad
          rc_ad = rc_ad - (1.d0-alf)*temporary_ad
          CALL POPCONTROL2B(branch)
          IF (branch .NE. 0) THEN
            IF (branch .EQ. 1) THEN
              p2 = rc/la
              CALL ALFA_AD(p2, p2_ad, alf_ad)
              rc_ad = rc_ad + p2_ad/la
              la_ad = la_ad - rc*p2_ad/la**2
            END IF
          END IF
          CALL POPREAL8(alf)
          ra_ad = va*0.5d0*rc_ad
          va_ad = ra*0.5d0*rc_ad
          CALL POPREAL8(va)
          CALL VY_AD(i, j, k, ismpl, va_ad)
          CALL POPREAL8(ra)
          CALL AMEAN_AD(result1, result1_ad, rijk, rijk_ad, ra_ad)
          arg1 = j + 1
          CALL RHOCF_AD(i, arg1, k, ismpl, result1_ad)
          CALL POPREAL8(la)
          CALL LJ_AD(i, j, k, ismpl, la_ad)
        END IF
 100    CALL POPCONTROL2B(branch)
        IF (branch .EQ. 0) THEN
          rc = 0.5*ra*va
          temporary_ad = -(d_ad(i, j, k, ismpl)/delx(i))
          la_ad = temporary_ad
          alf_ad = rc*temporary_ad
          rc_ad = -((1.d0-alf)*temporary_ad)
          temporary_ad = c_ad(i, j, k, ismpl)/delx(i)
          c_ad(i, j, k, ismpl) = 0.D0
          la_ad = la_ad + temporary_ad
          alf_ad = alf_ad + rc*temporary_ad
          rc_ad = rc_ad + (alf+1.d0)*temporary_ad
          CALL POPCONTROL2B(branch)
          IF (branch .NE. 0) THEN
            IF (branch .EQ. 1) THEN
              p2 = rc/la
              CALL ALFA_AD(p2, p2_ad, alf_ad)
              rc_ad = rc_ad + p2_ad/la
              la_ad = la_ad - rc*p2_ad/la**2
            END IF
          END IF
          CALL POPREAL8(alf)
          ra_ad = va*0.5*rc_ad
          va_ad = ra*0.5*rc_ad
          arg1 = i - 1
          CALL POPREAL8(va)
          CALL VX_AD(arg1, j, k, ismpl, va_ad)
          CALL POPREAL8(ra)
          CALL AMEAN_AD(result1, result1_ad, rijk, rijk_ad, ra_ad)
          arg1 = i - 1
          CALL RHOCF_AD(arg1, j, k, ismpl, result1_ad)
          arg1 = i - 1
          CALL POPREAL8(la)
          CALL LI_AD(arg1, j, k, ismpl, la_ad)
        ELSE IF (branch .NE. 1) THEN
          GOTO 110
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          rc = 0.5*ra*va
          temporary_ad = -(d_ad(i, j, k, ismpl)/delx(i))
          la_ad = temporary_ad
          alf_ad = rc*temporary_ad
          rc_ad = (alf+1.d0)*temporary_ad
          temporary_ad = e_ad(i, j, k, ismpl)/delx(i)
          e_ad(i, j, k, ismpl) = 0.D0
          la_ad = la_ad + temporary_ad
          alf_ad = alf_ad + rc*temporary_ad
          rc_ad = rc_ad - (1.d0-alf)*temporary_ad
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            p2 = rc/la
            CALL POPREAL8(alf)
            CALL ALFA_AD(p2, p2_ad, alf_ad)
            rc_ad = rc_ad + p2_ad/la
            la_ad = la_ad - rc*p2_ad/la**2
          ELSE
            CALL POPREAL8(alf)
          END IF
          ra_ad = va*0.5*rc_ad
          va_ad = ra*0.5*rc_ad
          CALL POPREAL8(va)
          CALL VX_AD(i, j, k, ismpl, va_ad)
          CALL POPREAL8(ra)
          CALL AMEAN_AD(result1, result1_ad, rijk, rijk_ad, ra_ad)
          arg1 = i + 1
          CALL RHOCF_AD(arg1, j, k, ismpl, result1_ad)
          CALL POPREAL8(la)
          CALL LI_AD(i, j, k, ismpl, la_ad)
        END IF
 110    CALL RHOCF_AD(i, j, k, ismpl, rijk_ad)
      END DO
    END DO
  END DO
END SUBROUTINE SET_TCOEF_AD
