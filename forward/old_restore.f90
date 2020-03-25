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

!>    @brief parallelisation wrapper for "omp_old_restore"
!>    @param[in] level level number
!>    @param[in] ismpl local sample index
      SUBROUTINE old_restore(level,ismpl)
        use mod_genrl
        use mod_OMP_TOOLS
        IMPLICIT NONE
        integer :: ismpl

        INCLUDE 'OMP_TOOLS.inc'
        INTEGER level
        INTRINSIC abs

#ifdef fOMP
!$OMP parallel num_threads(Tlevel_1)
!$      call omp_binding(abs(ismpl))
#endif
        CALL omp_old_restore(level,ismpl)
#ifdef fOMP
!$OMP end parallel
#endif
!
        RETURN
      END

!>    @brief restores an old state/version
!>    @param[in] level level number (which old version)
!>    @param[in] ismpl local sample index
      SUBROUTINE omp_old_restore(level,ismpl)
        use arrays
        use mod_genrl
        use mod_conc
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k, l

        INTEGER level, tpos, tanz
        INTRINSIC abs, max

        CALL omp_part(i0*j0*k0,tpos,tanz)
        CALL ijk_m(tpos,i,j,k)
!     save state (before)
        CALL dcopy(tanz,headold(tpos,level,abs(ismpl)),1, head(i,j,k,max(1,ismpl)), 1)
        CALL dcopy(tanz,tempold(tpos,level,abs(ismpl)),1, temp(i,j,k,max(1,ismpl)), 1)
        CALL dcopy(tanz,presold(tpos,level,abs(ismpl)),1, pres(i,j,k,max(1,ismpl)), 1)
        DO l = 1, ntrans
          CALL dcopy(tanz,concold(tpos,l,level,abs(ismpl)),1, conc(i,j,k,l,max(1,ismpl)),1)
        END DO
!
        RETURN
      END
