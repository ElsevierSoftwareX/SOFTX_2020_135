!        Generated by TAPENADE     (INRIA, Ecuador team)
!  tapenade 3.x
!
!  Differentiation of set_tcoef in forward (tangent) mode:
!   variations   of useful results: *d *e *f *g *a *b *c
!   with respect to varying inputs: *d *e *f *g *temp *propunit
!                *pres *a *b *c
!   Plus diff mem management of: d:in e:in f:in g:in temp:in propunit:in
!                pres:in a:in b:in c:in
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
SUBROUTINE g_SET_TCOEF(ismpl)
  USE ARRAYS

  USE g_ARRAYS

  USE MOD_GENRL
  USE MOD_GENRLC
  USE MOD_TEMP
  USE MOD_TIME

  USE g_MOD_TIME

  USE MOD_LINFOS
  IMPLICIT NONE
  INTEGER :: ismpl
  INTEGER :: i, j, k
  EXTERNAL LI, LJ, LK, RHOCF, VX, &
&     VY, VZ, ALFA, AMEAN
  EXTERNAL g_LI, g_LJ, g_LK, g_RHOCF, g_VX, g_VY, g_VZ, &
&     g_ALFA, g_AMEAN
  DOUBLE PRECISION :: LI, LJ, LK, RHOCF, &
& VX, VY, VZ, ALFA, AMEAN
  DOUBLE PRECISION :: g_LI, g_LJ, g_LK, g_RHOCF, g_VX, g_VY&
& , g_VZ, g_ALFA, g_AMEAN
  DOUBLE PRECISION :: rijk, ra, rc, va, la, alf, p2
  DOUBLE PRECISION :: g_rijk, g_ra, g_rc, g_va, g_la, g_alf&
& , g_p2
  EXTERNAL DUMMY
  DOUBLE PRECISION :: result1
  DOUBLE PRECISION :: g_result1
!$OMP master
  IF (linfos(3) .GE. 2) WRITE(*, *) ' ... tcoef'
!$OMP end master
! initialize coefficients for sparse solvers
! inner points of grid - - - - - - - - - - - - - - - - - - - - - - - - -
!$OMP do schedule(static) collapse(3)
  DO k=1,k0
    DO j=1,j0
      DO i=1,i0
        g_rijk = g_RHOCF(i, j, k, ismpl, rijk)
        IF (i0 .GT. 1) THEN
          IF (i .LT. i0) THEN
            g_la = g_LI(i, j, k, ismpl, la)
            g_result1 = g_RHOCF(i + 1, j, k, ismpl, result1)
            g_ra = g_AMEAN(result1, g_result1, rijk, g_rijk, ra)
            g_va = g_VX(i, j, k, ismpl, va)
            g_rc = 0.5*(va*g_ra+ra*g_va)
            rc = 0.5*ra*va
            IF (la .GT. 0.d0) THEN
              g_p2 = (g_rc-rc*g_la/la)/la
              p2 = rc/la
              g_alf = g_ALFA(p2, g_p2, alf)
            ELSE
              alf = 0.d0
              IF (va .LT. 0.d0) alf = -1.d0
              IF (va .GT. 0.d0) THEN
                alf = 1.d0
                g_alf = 0.D0
              ELSE
                g_alf = 0.D0
              END IF
            END IF
            g_e(i, j, k, ismpl) = (g_la-(1.d0-alf)*g_rc+rc*g_alf&
&             )/delx(i)
            e(i, j, k, ismpl) = (la-(1.d0-alf)*rc)/delx(i)
            g_d(i, j, k, ismpl) = g_d(i, j, k, ismpl) - (g_la+rc*&
&             g_alf+(alf+1.d0)*g_rc)/delx(i)
            d(i, j, k, ismpl) = d(i, j, k, ismpl) - (la+(1.d0+alf)*rc)/&
&             delx(i)
          END IF
          IF (i .GT. 1) THEN
            g_la = g_LI(i - 1, j, k, ismpl, la)
            g_result1 = g_RHOCF(i - 1, j, k, ismpl, result1)
            g_ra = g_AMEAN(result1, g_result1, rijk, g_rijk, ra)
            g_va = g_VX(i - 1, j, k, ismpl, va)
            g_rc = 0.5*(va*g_ra+ra*g_va)
            rc = 0.5*ra*va
            alf = 0.d0
            IF (va .EQ. 0.d0) THEN
              alf = 0.d0
              g_alf = 0.D0
            ELSE IF (la .GT. 0.d0) THEN
              g_p2 = (g_rc-rc*g_la/la)/la
              p2 = rc/la
              g_alf = g_ALFA(p2, g_p2, alf)
            ELSE
              IF (va .LT. 0.d0) alf = -1.d0
              IF (va .GT. 0.d0) THEN
                alf = 1.d0
                g_alf = 0.D0
              ELSE
                g_alf = 0.D0
              END IF
            END IF
            g_c(i, j, k, ismpl) = (g_la+rc*g_alf+(alf+1.d0)*g_rc&
&             )/delx(i)
            c(i, j, k, ismpl) = (la+(1.d0+alf)*rc)/delx(i)
            g_d(i, j, k, ismpl) = g_d(i, j, k, ismpl) - (g_la-(&
&             1.d0-alf)*g_rc+rc*g_alf)/delx(i)
            d(i, j, k, ismpl) = d(i, j, k, ismpl) - (la-(1.d0-alf)*rc)/&
&             delx(i)
          END IF
        END IF
        IF (j0 .GT. 1) THEN
          IF (j .LT. j0) THEN
            g_la = g_LJ(i, j, k, ismpl, la)
            g_result1 = g_RHOCF(i, j + 1, k, ismpl, result1)
            g_ra = g_AMEAN(result1, g_result1, rijk, g_rijk, ra)
            g_va = g_VY(i, j, k, ismpl, va)
            g_rc = 0.5d0*(va*g_ra+ra*g_va)
            rc = 0.5d0*ra*va
            alf = 0.d0
            IF (va .EQ. 0.d0) THEN
              alf = 0.d0
              g_alf = 0.D0
            ELSE IF (la .GT. 0.d0) THEN
              g_p2 = (g_rc-rc*g_la/la)/la
              p2 = rc/la
              g_alf = g_ALFA(p2, g_p2, alf)
            ELSE
              IF (va .LT. 0.d0) alf = -1.d0
              IF (va .GT. 0.d0) THEN
                alf = 1.d0
                g_alf = 0.D0
              ELSE
                g_alf = 0.D0
              END IF
            END IF
            g_f(i, j, k, ismpl) = (g_la-(1.d0-alf)*g_rc+rc*g_alf&
&             )/dely(j)
            f(i, j, k, ismpl) = (la-(1.d0-alf)*rc)/dely(j)
            g_d(i, j, k, ismpl) = g_d(i, j, k, ismpl) - (g_la+rc*&
&             g_alf+(alf+1.d0)*g_rc)/dely(j)
            d(i, j, k, ismpl) = d(i, j, k, ismpl) - (la+(1.d0+alf)*rc)/&
&             dely(j)
          END IF
          IF (j .GT. 1) THEN
            g_la = g_LJ(i, j - 1, k, ismpl, la)
            g_result1 = g_RHOCF(i, j - 1, k, ismpl, result1)
            g_ra = g_AMEAN(result1, g_result1, rijk, g_rijk, ra)
            g_va = g_VY(i, j - 1, k, ismpl, va)
            g_rc = 0.5*(va*g_ra+ra*g_va)
            rc = 0.5*ra*va
            alf = 0.d0
            IF (va .EQ. 0.d0) THEN
              alf = 0.d0
              g_alf = 0.D0
            ELSE IF (la .GT. 0.d0) THEN
              g_p2 = (g_rc-rc*g_la/la)/la
              p2 = rc/la
              g_alf = g_ALFA(p2, g_p2, alf)
            ELSE
              IF (va .LT. 0.d0) alf = -1.d0
              IF (va .GT. 0.d0) THEN
                alf = 1.d0
                g_alf = 0.D0
              ELSE
                g_alf = 0.D0
              END IF
            END IF
            g_b(i, j, k, ismpl) = (g_la+rc*g_alf+(alf+1.d0)*g_rc&
&             )/dely(j)
            b(i, j, k, ismpl) = (la+(1.d0+alf)*rc)/dely(j)
            g_d(i, j, k, ismpl) = g_d(i, j, k, ismpl) - (g_la-(&
&             1.d0-alf)*g_rc+rc*g_alf)/dely(j)
            d(i, j, k, ismpl) = d(i, j, k, ismpl) - (la-(1.d0-alf)*rc)/&
&             dely(j)
          END IF
        END IF
        IF (k0 .GT. 1) THEN
          IF (k .LT. k0) THEN
            g_la = g_LK(i, j, k, ismpl, la)
            g_result1 = g_RHOCF(i, j, k + 1, ismpl, result1)
            g_ra = g_AMEAN(result1, g_result1, rijk, g_rijk, ra)
            g_va = g_VZ(i, j, k, ismpl, va)
            g_rc = 0.5*(va*g_ra+ra*g_va)
            rc = 0.5*ra*va
            alf = 0.d0
            IF (va .EQ. 0.d0) THEN
              alf = 0.d0
              g_alf = 0.D0
            ELSE IF (la .GT. 0.d0) THEN
              g_p2 = (g_rc-rc*g_la/la)/la
              p2 = rc/la
              g_alf = g_ALFA(p2, g_p2, alf)
            ELSE
              IF (va .LT. 0.d0) alf = -1.d0
              IF (va .GT. 0.d0) THEN
                alf = 1.d0
                g_alf = 0.D0
              ELSE
                g_alf = 0.D0
              END IF
            END IF
            g_g(i, j, k, ismpl) = (g_la-(1.d0-alf)*g_rc+rc*g_alf&
&             )/delz(k)
            g(i, j, k, ismpl) = (la-(1.d0-alf)*rc)/delz(k)
            g_d(i, j, k, ismpl) = g_d(i, j, k, ismpl) - (g_la+rc*&
&             g_alf+(alf+1.d0)*g_rc)/delz(k)
            d(i, j, k, ismpl) = d(i, j, k, ismpl) - (la+(1.d0+alf)*rc)/&
&             delz(k)
          END IF
          IF (k .GT. 1) THEN
            g_la = g_LK(i, j, k - 1, ismpl, la)
            g_result1 = g_RHOCF(i, j, k - 1, ismpl, result1)
            g_ra = g_AMEAN(result1, g_result1, rijk, g_rijk, ra)
            g_va = g_VZ(i, j, k - 1, ismpl, va)
            g_rc = 0.5*(va*g_ra+ra*g_va)
            rc = 0.5*ra*va
            alf = 0.d0
            IF (va .EQ. 0.d0) THEN
              alf = 0.d0
              g_alf = 0.D0
            ELSE IF (la .GT. 0.d0) THEN
              g_p2 = (g_rc-rc*g_la/la)/la
              p2 = rc/la
              g_alf = g_ALFA(p2, g_p2, alf)
            ELSE
              IF (va .LT. 0.d0) alf = -1.d0
              IF (va .GT. 0.d0) THEN
                alf = 1.d0
                g_alf = 0.D0
              ELSE
                g_alf = 0.D0
              END IF
            END IF
            g_a(i, j, k, ismpl) = (g_la+rc*g_alf+(alf+1.d0)*g_rc&
&             )/delz(k)
            a(i, j, k, ismpl) = (la+(1.d0+alf)*rc)/delz(k)
            g_d(i, j, k, ismpl) = g_d(i, j, k, ismpl) - (g_la-(&
&             1.d0-alf)*g_rc+rc*g_alf)/delz(k)
            d(i, j, k, ismpl) = d(i, j, k, ismpl) - (la-(1.d0-alf)*rc)/&
&             delz(k)
          END IF
        END IF
      END DO
    END DO
  END DO
!$OMP end do nowait
  RETURN
END SUBROUTINE g_SET_TCOEF

