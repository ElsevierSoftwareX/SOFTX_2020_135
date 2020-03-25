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

      MODULE m_random

      CONTAINS
        SUBROUTINE random(work1,n)
!  Returns a vector of random values N(variance=1,mean=0)
          IMPLICIT NONE
          integer, INTENT (IN) :: n
          double precision, INTENT (OUT), dimension (n) :: work1
          double precision, ALLOCATABLE, dimension (:) :: work2
          double precision, PARAMETER :: pi = 3.141592653589D0

          ALLOCATE(work2(n))

          CALL random_number(work1)
          CALL random_number(work2)
          work1 = dsqrt(-2.0D0*dlog(work1+dble(tiny(1.0))))* &
            dcos(2.0D0*pi*work2)

          DEALLOCATE(work2)
        END SUBROUTINE random
      END MODULE m_random
