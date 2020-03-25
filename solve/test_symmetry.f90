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

!>    @brief test the symmetry of the system matrix
!>    @param[in] ismpl local sample index
!>    @return true: when the matrix is symmetric
      LOGICAL FUNCTION test_symmetry(ismpl)
        use arrays
        use mod_genrl
        use mod_linfos
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k
!     to small koeff.
        INTEGER fehler
        LOGICAL isit

!     test about symmetry
        fehler = 0
        DO k = 1, k0
          DO j = 1, j0
            DO i = 1, i0
              IF (d(i,j,k,ismpl)/=0.0D0) THEN
                IF (i>1) THEN
                  IF ((c(i,j,k,ismpl)/=e(i-1,j,k,ismpl)) .AND. (d(i-1, &
                    j,k,ismpl)/=0.0D0)) fehler = fehler + 1
                END IF
                IF (j>1) THEN
                  IF ((b(i,j,k,ismpl)/=f(i,j-1,k,ismpl)) .AND. (d(i, &
                    j-1,k,ismpl)/=0.0D0)) fehler = fehler + 1
                END IF
                IF (k>1) THEN
                  IF ((a(i,j,k,ismpl)/=g(i,j,k-1,ismpl)) .AND. (d(i,j, &
                    k-1,ismpl)/=0.0D0)) fehler = fehler + 1
                END IF
              END IF
            END DO
          END DO
        END DO
!
        IF (fehler==0) THEN
          IF (linfos(4)>=2) WRITE(*,*) ' Symmetric matrix !'
          isit = .TRUE.
        ELSE
          IF (linfos(4)>=2) WRITE(*,*) ' Unsymmetric matrix with ', &
            fehler, ' unsym. points !'
          isit = .FALSE.
        END IF
!
        test_symmetry = isit
!
        RETURN
      END
