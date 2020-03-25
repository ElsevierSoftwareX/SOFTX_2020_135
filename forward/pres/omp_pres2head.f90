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

!>    @brief head setup
!>    @param[in] init flag: 0-init, 1-normal setup
!>    @param[in] ismpl local sample index
!>    @details
!> head setup currently disabled\n
      SUBROUTINE omp_pres2head(init,ismpl)
        use arrays
        use mod_flow
        use mod_genrl
        use mod_linfos
        use mod_OMP_TOOLS
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k
        INCLUDE 'OMP_TOOLS.inc'
        INTEGER init
        DOUBLE PRECISION psurf

!     presetting for dirichlet boundary conditions to avoid site effects
        IF (pres_active) THEN
           CALL set_dpbc(ismpl)
           CALL set_dtbc(ismpl)
           psurf = 1.0D5
!$OMP barrier
!
!$OMP do schedule(static) collapse(3)
           DO i = 1, i0
              DO j = 1, j0
                DO k = 1, k0
                  head(i,j,k,ismpl) = (pres(i,j,k,ismpl) - psurf)/(rref*grav) + delza(k)
                END DO
              END DO
           END DO
!$OMP end do nowait
        END IF
!
        RETURN
      END
