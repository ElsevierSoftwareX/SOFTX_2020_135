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

!>    @brief parallelisation wrapper for "omp_set_tsal"
!>    @param[in] ismpl local sample index
      SUBROUTINE set_tsal(ismpl)
        use mod_genrl
        use mod_OMP_TOOLS
        IMPLICIT NONE
        integer :: ismpl
        INCLUDE 'OMP_TOOLS.inc'

#ifdef fOMP
!$OMP parallel num_threads(Tlevel_1)
!$      call omp_binding(ismpl)
#endif
        CALL omp_set_tsal(ismpl)
#ifdef fOMP
!$OMP end parallel
#endif

        RETURN
      END

!>    @brief calculate total salinity
!>    @param[in] ismpl local sample index
!>    @details
!>calculate total salinity\n
      SUBROUTINE omp_set_tsal(ismpl)
        use arrays
        use mod_conc
        use mod_genrl
        use mod_genrlc
        use mod_linfos
        IMPLICIT NONE

        integer :: ismpl
        integer :: i, j, k
        INTEGER species
        DOUBLE PRECISION summ, fac

!$OMP master
        IF (linfos(3)>=2) WRITE(*,'(1A)') &
          '  ... setting total salinity'
!$OMP end master

!$OMP do schedule(static) collapse(3)
        DO k = 1, k0
          DO j = 1, j0
            DO i = 1, i0
              summ = 0.D0
              DO species = 1, ntrans
                fac = mmas_c(species)/mmas_nacl
                summ = summ + fac*conc(i,j,k,species,ismpl)
              END DO
              tsal(i,j,k,ismpl) = summ
            END DO
          END DO
        END DO
!$OMP end do nowait

        RETURN
      END
