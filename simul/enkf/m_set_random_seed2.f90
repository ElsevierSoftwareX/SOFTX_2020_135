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

      MODULE m_set_random_seed2

      CONTAINS
        SUBROUTINE set_random_seed2
! Sets a random seed based on the system and wall clock time
          IMPLICIT NONE

          integer, DIMENSION (8) :: val
          integer ::  cnt
          integer ::  sze
          integer, ALLOCATABLE, DIMENSION (:) :: pt

          CALL date_and_time(values=val)
          CALL system_clock(count=cnt)
          CALL random_seed(size=sze)
          ALLOCATE(pt(max(sze,2)))
          pt(1) = val(8)*val(3)
          pt(2) = cnt
          CALL random_seed(put=pt)
          DEALLOCATE(pt)
        END SUBROUTINE set_random_seed2
      END MODULE m_set_random_seed2
