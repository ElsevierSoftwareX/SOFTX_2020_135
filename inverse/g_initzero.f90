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

!>    @brief parallelisation wrapper for "omp_g_initzero"
!>    @param[in] ismpl local sample index
      SUBROUTINE g_initzero(ismpl)
        use mod_genrl
        use mod_OMP_TOOLS
        IMPLICIT NONE
        integer :: ismpl
        INCLUDE 'OMP_TOOLS.inc'

#ifdef AD 
#ifdef fOMP
!$OMP parallel num_threads(Tlevel_1)
!$      call omp_binding(ismpl)
#endif
        CALL omp_g_initzero(ismpl)
#ifdef fOMP
!$OMP end parallel
#endif
#endif

        RETURN
      END

#ifdef AD 
!>    @brief to initialize all derived-variables with zero
!>    @param[in] ismpl local sample index
!>    @details
!>    Initialize all derived-variables with zero needed before seeding\n
      SUBROUTINE omp_g_initzero(ismpl)
        use arrays
        use g_arrays
        use mod_genrl
        use mod_conc
        use mod_time
        use g_mod_time
        use mod_data
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k, l
!     OpenMP stuff
        INTEGER tpos, tanz
!     OpenMP stuff
        CALL omp_part(i0*j0*k0,tpos,tanz)
        CALL ijk_m(tpos,i,j,k)

        CALL set_dval(tanz,0.D0,g_a(i,j,k,ismpl))
        CALL set_dval(tanz,0.D0,g_b(i,j,k,ismpl))
        CALL set_dval(tanz,0.D0,g_c(i,j,k,ismpl))
        CALL set_dval(tanz,0.D0,g_d(i,j,k,ismpl))
        CALL set_dval(tanz,0.D0,g_e(i,j,k,ismpl))
        CALL set_dval(tanz,0.D0,g_f(i,j,k,ismpl))
        CALL set_dval(tanz,0.D0,g_g(i,j,k,ismpl))
        CALL set_dval(tanz,0.D0,g_w(i,j,k,ismpl))
        CALL set_dval(tanz,0.D0,g_x(i,j,k,ismpl))
! -----------------
        CALL set_dval(tanz,0.D0,g_head(i,j,k,ismpl))
        CALL set_dval(tanz,0.D0,g_temp(i,j,k,ismpl))
        CALL set_dval(tanz,0.D0,g_pres(i,j,k,ismpl))
        DO l = 1, ntrans
          CALL set_dval(tanz,0.D0,g_conc(i,j,k,l,ismpl))
        END DO
!     OpenMP call instead of the regular one !!!
        CALL omp_gp_old_save(cgen_fw,ismpl)
        CALL omp_gp_old_save(cgen_time,ismpl)
! only thread-0 for all others
!$OMP master
        g_simtime(ismpl) = 0.D0
        CALL set_dval(nunits*nprop,0.D0,g_propunit(1,1,ismpl))
        CALL set_dval(ngsmax*3*nbctp,0.D0,g_bcperiod(1,1,1,ismpl))
        IF (ndata>0) THEN
          CALL set_dval(ndata,0.D0,sdata(1,ismpl))
          CALL set_dval(ndata,0.D0,g_sdata(1,ismpl))
        END IF
        CALL set_dval(nbc_data*ndbc,0.D0,g_dbc_data(1,1,ismpl))
!$OMP end master
        RETURN
      END SUBROUTINE
#endif
