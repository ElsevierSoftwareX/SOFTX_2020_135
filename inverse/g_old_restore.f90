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

!>    @brief parallelisation wrapper for "omp_g_old_restore"
!>    @param[in] level level number
!>    @param[in] ismpl local sample index
      SUBROUTINE g_old_restore(level,ismpl)
        use mod_genrl
        use mod_OMP_TOOLS
        IMPLICIT NONE
        integer :: ismpl
        INCLUDE 'OMP_TOOLS.inc'
        INTEGER level

#ifndef AD_RM
#ifdef fOMP
!$OMP parallel num_threads(Tlevel_1)
!$      call omp_binding(ismpl)
#endif
        CALL omp_g_old_restore(level,ismpl)
#ifdef fOMP
!$OMP end parallel
#endif
#endif
        RETURN
      END
      
#ifndef AD_RM
!>    @brief restore an old state
!>    @param[in] level level number
!>    @param[in] ismpl local sample index
      SUBROUTINE omp_g_old_restore(level,ismpl)
        use arrays
        use g_arrays
        use mod_genrl
        use mod_conc
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k, l

        INTEGER level, tpos, tanz


        IF (level>4 .OR. level<=0) THEN
          WRITE(*,'(1A)') &
            'error: "level" out of range in "g_old_restore" !'
          STOP
        END IF

        CALL omp_part(i0*j0*k0,tpos,tanz)
        CALL ijk_m(tpos,i,j,k)

!     save state of the gradients
        CALL dcopy(tanz,g_headold(tpos,level,ismpl),1, g_head(i,j,k,ismpl),1)
        CALL dcopy(tanz,g_tempold(tpos,level,ismpl),1, g_temp(i,j,k,ismpl),1)
        DO l = 1, ntrans
          CALL dcopy(tanz,g_concold(tpos,l,level,ismpl),1, g_conc(i,j,k,l,ismpl),1)
        END DO
        CALL dcopy(tanz,g_presold(tpos,level,ismpl),1, g_pres(i,j,k,ismpl),1)
        RETURN
      END
#endif
