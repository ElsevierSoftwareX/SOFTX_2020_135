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

! ### This file contains an OpenMP workaround ###
! ### for a bug found in the intel compiler ! ###
!    16.05.2008: intel compiler build 10.1.015


!>    @brief extended OpenMP ORDERED alternative
!>    @details
!> < this is an alternative workaround for OpenMP >\n
      MODULE omp_ordered
        IMPLICIT NONE
!     locking vector for ordering
        INTEGER, ALLOCATABLE :: ordered_locks(:)
!     number of locks
        INTEGER nordered_locks
      END MODULE omp_ordered

!>    @brief create the locking vector (for ORDERED)
!>    @param[in] llen max loop index
!>    @details
!> alternative to "!$omp do ordered"\n
      SUBROUTINE omp_ordered_create(llen)
        USE omp_ordered
        IMPLICIT NONE
        INTEGER i
!     max loop index
        INTEGER llen

!$OMP master
!       initialisation
        ALLOCATE(ordered_locks(llen))
!       clear the locking vector
        DO i = 1, llen
          ordered_locks(i) = 0
        END DO
        nordered_locks = llen
!$OMP end master
!     flush the locking vector
!$OMP barrier
        RETURN
      END

!>    @brief free the locking vector (for ORDERED)
!>    @details
!> additionally to "!$omp end do" (ORDERED at begin)\n
      SUBROUTINE omp_ordered_delete()
        USE omp_ordered
        IMPLICIT NONE

!     first, wait for all threads
!$OMP barrier
!$OMP master
!       initialisation
        DEALLOCATE(ordered_locks)
        nordered_locks = 0
!$OMP end master
        RETURN
      END

!>    @brief begin ordered section (independent levels)
!>    @param[in] idx "ordered" loop index
!>    @details
!> alternative to "!$omp ordered"\n
!> Allows more than one occurrence during a parllel loop\n
!> -> more than one level (different occurrence - different level).\n 
!> Each level is independent from the one before.\n
      SUBROUTINE omp_ordered_begin(idx)
        USE omp_ordered
        IMPLICIT NONE
!     in the right order?
        LOGICAL ok
!     temp. loop index
        INTEGER i
!     "ordered" loop index
        INTEGER idx

!     polling over the state of the previous iterations
10      ok = .TRUE.
        DO i = 1, idx - 1
          IF (ordered_locks(i)<=ordered_locks(idx)) ok = .FALSE.
        END DO
!     can entering the ordered section ?
        IF (ok) RETURN
!     flush the state
        CALL omp_ordered_flush()
        GO TO 10
      END

!>    @brief begin ordered section (chunk size dependent levels)
!>    @param[in] idx "ordered" loop index
!>    @details
!> alternative to "!$omp ordered"\n
!> Allows more than one occurrence during a parllel loop\n
!> -> more than one level (different occurrence - different level).\n 
!> Each level depends on the end of the one before,\n
!> but inside a chunk size (index interval).\n
!> All runs in the chunks before needs to have an equal level.\n
      SUBROUTINE omp_ordered_dep_begin(idx, ichunk)
        USE omp_ordered
        IMPLICIT NONE
!     in the right order?
        LOGICAL ok
!     temp. loop index
        INTEGER i
!     "ordered" loop index
        INTEGER idx
!     chunk size for dependency, start/end index of the chunk
        INTEGER ichunk, schunk, echunk
        INTRINSIC min

        schunk = ((idx -1) /ichunk) *ichunk +1
        echunk = min(schunk +ichunk -1, nordered_locks)
!     polling over the state of the previous iterations
10      ok = .TRUE.
!       test all chunks before
        DO i = 1, schunk -1
          IF (ordered_locks(i)/=ordered_locks(1)) ok = .FALSE.
        END DO
!       test runs before (same chunk)
        DO i = schunk, idx -1
          IF (ordered_locks(i)<=ordered_locks(idx)) ok = .FALSE.
        END DO
!       test runs after (same chunk)
        DO i = idx +1, echunk
          IF (ordered_locks(i)<ordered_locks(idx)) ok = .FALSE.
        END DO
!     can entering the ordered section ?
        IF (ok) RETURN
!     flush the state
        CALL omp_ordered_flush()
        GO TO 10
      END

!>    @brief finish ordered section
!>    @param[in] idx "ordered" loop index
!>    @details
!> alternative to "!$omp end ordered"\n
      SUBROUTINE omp_ordered_end(idx)
        USE omp_ordered
        IMPLICIT NONE
!     "ordered" loop index
        INTEGER idx

!     increase my state
        ordered_locks(idx) = ordered_locks(idx) + 1
!     flush the state
        CALL omp_ordered_flush()
        RETURN
      END

!>    @brief flush the locking vector (workaround)
      SUBROUTINE omp_ordered_flush()
        USE omp_ordered
        IMPLICIT NONE

!     flush the state
!$OMP flush(ordered_locks)
        RETURN
      END
