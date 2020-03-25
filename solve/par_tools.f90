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

!>    @brief deallocate special solver staff, [proza]
      SUBROUTINE par_end2()
        use arrays
        IMPLICIT NONE

!$OMP master
        DEALLOCATE(proza)
!$OMP end master

        RETURN
      END

!>    @brief function wrapper for [proza_lock] - avoid side effects
!>    @param[in] ProzA_lock processor lock variable
!>    @return index of processor
!>    @details
!>    prove entry (locking)\n
      INTEGER FUNCTION par_name(proza_lock)
        IMPLICIT NONE
        INTEGER proza_lock

!$OMP flush(ProzA_lock)
        par_name = proza_lock

        RETURN
      END

!>    @brief wrapper routine to set [proza_lock] - avoid side effects
!>    @param[out] ProzA_lock processor lock variable
!>    @details
!>    disable an entry (locking)\n
      SUBROUTINE par_disab(proza_lock)
        IMPLICIT NONE
        INTEGER proza_lock

        proza_lock = 1
!$OMP flush(ProzA_lock)

        RETURN
      END

!>    @brief initialize/reset [proza_lock] - clean marker of already computed blocks
!>    @param[out] ProzA_lock processor lock variable
      SUBROUTINE par_reset(proza_lock)
        use mod_blocking_size
        IMPLICIT NONE
        INTEGER proza_lock(bdim_i,bdim_j,bdim_k)
        INTEGER i, j, k

        DO k = 1, bdim_k
          DO j = 1, bdim_j
            DO i = 1, bdim_i
              proza_lock(i,j,k) = 0
            END DO
          END DO
        END DO

        RETURN
      END
