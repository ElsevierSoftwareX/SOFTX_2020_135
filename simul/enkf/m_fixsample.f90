! MIT License
!
! Copyright (c) 2019 Geir Evensen
!
! Permission is hereby granted, free of charge, to any person
! obtaining a copy of this software and associated documentation
! files (the "Software"), to deal in the Software without
! restriction, including without limitation the rights to use, copy,
! modify, merge, publish, distribute, sublicense, and/or sell copies
! of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be
! included in all copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
! EXPRESS OR IMPLIED, INCLUDING BUT NOT NO EVENT SHALL THE AUTHORS OR
! COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
! ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE
! OR OTHER DEALINGS IN THE SOFTWARE.

      MODULE m_fixsample

      CONTAINS
        SUBROUTINE fixsample(e,n,m)
          integer, INTENT (IN) :: m
          integer, INTENT (IN) :: n
          double precision, INTENT (INOUT) :: e(n,m)

          integer :: iens, i
          double precision, ALLOCATABLE :: average(:), variance(:)

          ALLOCATE(average(n),variance(n))

          average = 0.0D0
          DO iens = 1, m
            average(:) = average(:) + e(:,iens)
          END DO
          average = average/dble(float(m))

          DO iens = 1, m
            e(:,iens) = e(:,iens) - average(:)
          END DO

          variance = 0.0D0
          DO iens = 1, m
            variance(:) = variance(:) + e(:,iens)**2.
          END DO

          DO i = 1, n
            variance(i) = 1.0D0/dsqrt(variance(i)/dble(float(m)))
          END DO

          DO iens = 1, m
            DO i = 1, n
              e(i,iens) = variance(i)*e(i,iens)
            END DO
          END DO

          DEALLOCATE(average,variance)

        END SUBROUTINE
      END MODULE m_fixsample






