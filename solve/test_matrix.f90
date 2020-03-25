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

!>    @brief test the numerical stability of the system matrix
!>    @param[in] ismpl local sample index
      SUBROUTINE test_matrix(ismpl)
        use arrays
        use mod_genrl
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k
        DOUBLE PRECISION k1, k2, k3, mindiff, max_val
!     to small koeff.
        PARAMETER (mindiff=1.0D-30)
        PARAMETER (max_val=1.0D+150)
        INTEGER fehler, nijk
        INTRINSIC dabs

        nijk = i0*j0*k0
!     test about math stability
        fehler = 0
        DO k = 1, k0
          DO j = 1, j0
            DO i = 1, i0
!      equal matrix lines, but different right sides ?
              IF (i<i0) THEN
                k1 = max_val
                k2 = max_val
                k3 = max_val
                IF (d(i,j,k,ismpl)/=0.0D0) k1 = c(i+1,j,k,ismpl)/ &
                  d(i,j,k,ismpl)
                IF (e(i,j,k,ismpl)/=0.0D0) k2 = d(i+1,j,k,ismpl)/ &
                  e(i,j,k,ismpl)
                IF (w(i,j,k,ismpl)/=0.0D0) k3 = w(i+1,j,k,ismpl)/ &
                  w(i,j,k,ismpl)
                IF ((a(i,j,k,ismpl)==0.0D0) .AND. (b(i,j,k, &
                    ismpl)==0.0D0) .AND. (c(i,j,k, &
                    ismpl)==0.0D0) .AND. (0.0D0==a(i+1,j,k, &
                    ismpl)) .AND. (0.0D0==b(i+1,j,k,ismpl)) .AND. & !aw  *       (d(i,j,k,ismpl).eq.c(i+1,j,k,ismpl)).and.
!aw  *       (e(i,j,k,ismpl).eq.d(i+1,j,k,ismpl)).and.
                    (f(i,j,k,ismpl)==0.0D0) .AND. (g(i,j,k, &
                    ismpl)==0.0D0) .AND. (0.0D0==e(i+1,j,k, &
                    ismpl)) .AND. (0.0D0==f(i+1,j,k, &
                    ismpl)) .AND. (0.0D0==g(i+1,j,k,ismpl)) .AND. & !aw  *       (W(i,j,k,ismpl).ne.W(i+1,j,k,ismpl))
                    ((dabs(k1-k2)<mindiff) .AND. (dabs(k1-k3)> &
                    mindiff))) THEN
                  fehler = fehler + 1
                  WRITE(*,*) ' ', i, j, k, ' i+1'
                  WRITE(*,*) 'D            E            | W'
                  WRITE(*,*) d(i,j,k,ismpl), e(i,j,k,ismpl), &
                    w(i,j,k,ismpl)
                  WRITE(*,*) 'C            D            | W'
                  WRITE(*,*) c(i+1,j,k,ismpl), d(i+1,j,k,ismpl), &
                    w(i+1,j,k,ismpl)
                END IF
              END IF
!
              IF (j<j0) THEN
                k1 = max_val
                k2 = max_val
                k3 = max_val
                IF (d(i,j,k,ismpl)/=0.0D0) k1 = b(i,j+1,k,ismpl)/ &
                  d(i,j,k,ismpl)
                IF (f(i,j,k,ismpl)/=0.0D0) k2 = d(i,j+1,k,ismpl)/ &
                  f(i,j,k,ismpl)
                IF (w(i,j,k,ismpl)/=0.0D0) k3 = w(i,j+1,k,ismpl)/ &
                  w(i,j,k,ismpl)
                IF ((a(i,j,k,ismpl)==0.0D0) .AND. (b(i,j,k, &
                    ismpl)==0.0D0) .AND. (c(i,j,k, &
                    ismpl)==0.0D0) .AND. (0.0D0==a(i,j+1,k, &
                    ismpl)) .AND. (0.0D0==c(i,j+1,k,ismpl)) .AND. & !aw  *       (d(i,j,k,ismpl).eq.b(i,j+1,k,ismpl)).and.
!aw  *       (f(i,j,k,ismpl).eq.d(i,j+1,k,ismpl)).and.
                    (e(i,j,k,ismpl)==0.0D0) .AND. (g(i,j,k, &
                    ismpl)==0.0D0) .AND. (0.0D0==e(i,j+1,k, &
                    ismpl)) .AND. (0.0D0==f(i,j+1,k, &
                    ismpl)) .AND. (0.0D0==g(i,j+1,k,ismpl)) .AND. & !aw  *       (W(i,j,k,ismpl).ne.W(i+1,j,k,ismpl))
                    ((dabs(k1-k2)<mindiff) .AND. (dabs(k1-k3)> &
                    mindiff))) THEN
                  fehler = fehler + 1
                  WRITE(*,*) ' ', i, j, k, ' j+1'
                  WRITE(*,*) 'D            F            | W'
                  WRITE(*,*) d(i,j,k,ismpl), f(i,j,k,ismpl), &
                    w(i,j,k,ismpl)
                  WRITE(*,*) 'B            D            | W'
                  WRITE(*,*) b(i,j+1,k,ismpl), d(i,j+1,k,ismpl), &
                    w(i,j+1,k,ismpl)
                END IF
              END IF
!
              IF (k<k0) THEN
                k1 = max_val
                k2 = max_val
                k3 = max_val
                IF (d(i,j,k,ismpl)/=0.0D0) k1 = a(i,j,k+1,ismpl)/ &
                  d(i,j,k,ismpl)
                IF (g(i,j,k,ismpl)/=0.0D0) k2 = d(i,j,k+1,ismpl)/ &
                  g(i,j,k,ismpl)
                IF (w(i,j,k,ismpl)/=0.0D0) k3 = w(i,j,k+1,ismpl)/ &
                  w(i,j,k,ismpl)
                IF ((a(i,j,k,ismpl)==0.0D0) .AND. (b(i,j,k, &
                    ismpl)==0.0D0) .AND. (c(i,j,k, &
                    ismpl)==0.0D0) .AND. (0.0D0==c(i,j,k+1, &
                    ismpl)) .AND. (0.0D0==b(i,j,k+1,ismpl)) .AND. & !aw  *       (d(i,j,k,ismpl).eq.a(i,j,k+1,ismpl)).and.
!aw  *       (g(i,j,k,ismpl).eq.d(i,j,k+1,ismpl)).and.
                    (f(i,j,k,ismpl)==0.0D0) .AND. (e(i,j,k, &
                    ismpl)==0.0D0) .AND. (0.0D0==e(i,j,k+1, &
                    ismpl)) .AND. (0.0D0==f(i,j,k+1, &
                    ismpl)) .AND. (0.0D0==g(i,j,k+1,ismpl)) .AND. & !aw  *       (W(i,j,k,ismpl).ne.W(i+1,j,k,ismpl))
                    ((dabs(k1-k2)<mindiff) .AND. (dabs(k1-k3)> &
                    mindiff))) THEN
                  fehler = fehler + 1
                  WRITE(*,*) ' ', i, j, k, ' K+1'
                  WRITE(*,*) 'D            G            | W'
                  WRITE(*,*) d(i,j,k,ismpl), g(i,j,k,ismpl), &
                    w(i,j,k,ismpl)
                  WRITE(*,*) 'A            D            | W'
                  WRITE(*,*) a(i,j,k+1,ismpl), d(i,j,k+1,ismpl), &
                    w(i,j,k+1,ismpl)
                END IF
              END IF
            END DO
          END DO
        END DO
        IF (fehler>0) THEN
          WRITE(*,*) &
            'error: inconsistent system, can not be solved ! (', &
            fehler, &
            ' linear lines with different right hand sides)'
          STOP
        END IF
!
        RETURN
      END
