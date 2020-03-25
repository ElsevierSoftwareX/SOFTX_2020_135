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

!>    @brief first initialisation for each realisation/sample/gradient/ensemble (original state variables)
!>    @param[in] ismpl local sample index
      SUBROUTINE single_init(ismpl)
        use arrays
        use mod_genrl
        use mod_conc
        use mod_time
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k
        INTEGER i_max, ilevel
        INTRINSIC max

! Here the parallelisation can be improved (second-inner level) !
! But beware the MASTER (inside) and BARRIER (outside) construct at the end !
!       save "OPTI" state (before)
        ilevel = cgen_opti
        CALL dcopy(i0*j0*k0,headold(1,ilevel,idx_master),1,headold(1,ilevel,ismpl),1)
        CALL dcopy(i0*j0*k0,tempold(1,ilevel,idx_master),1,tempold(1,ilevel,ismpl),1)
        CALL dcopy(i0*j0*k0*ntrans,concold(1,1,ilevel,idx_master),1,concold(1,1,ilevel,ismpl),1)
        CALL dcopy(i0*j0*k0,presold(1,ilevel,idx_master),1,presold(1,ilevel,ismpl),1)
!
!$OMP master
        i_max = max(maxunits,bc_maxunits)
        DO j = 1, nprop
          DO i = 1, i_max
            propunitold(i,j) = propunit(i,j,idx_master)
          END DO
        END DO
        DO k = 1, nbctp
          DO j = 1, 3
            DO i = 1, ngsmax
              bcperiodold(i,j,k) = bcperiod(i,j,k,idx_master)
            END DO
          END DO
        END DO
!$OMP end master
!
        RETURN
      END
