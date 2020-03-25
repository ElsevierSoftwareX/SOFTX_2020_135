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

!>    @brief computes the full Jacobi matrix (-> "jac")
!>    @param[in] idummy dummy holder for sample-index (ignored)
      SUBROUTINE jacobi_compute(idummy)
        use arrays
#ifdef AD
        use g_arrays
#endif
#ifdef AD_RM
        use arrays_ad
#endif
        use mod_genrl
        use mod_data
        use mod_time
        use mod_inverse
        use mod_OMP_TOOLS
        use mod_linfos
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j
        INCLUDE 'OMP_TOOLS.inc'
        INTEGER seed_komp, idummy,param_enabled
        DOUBLE PRECISION trun_oa, tend_oa
        external param_enabled

!
!         init master/backup copy - equal for all samples
          CALL DCOPY(mpara,main_input(1,idx_master),1,main_input_master,1)
#ifdef BENCH
          CALL sys_cputime(trun_oa)
#endif
          write_disable = .TRUE.
          write_iter_disable = .TRUE.

#ifndef AD_RM
!        Jacobian computation / seeding
!          sample computation, dynamic scheduling for OpenMP
#ifdef fOMP
!$OMP    parallel num_threads(Tlevel_0) &
!$OMP      default(none) private(ismpl, seed_komp) &
!$OMP      shared(iter_inv, update_dd, mpara, Tlevel_0) &
!$OMP      shared(main_input, g_main_input, main_output) &
!$OMP      shared(simtime_0, max_simtime, main_input_master, jac)
#endif
!         computing in parallel threads, sample index
          ismpl = omp_get_his_thread_num() + 1
!         sample init, each sample only one times
          CALL single_init(ismpl)
!$OMP     barrier
!
! !!! to avoid compiler bugs
          CALL omp_ordered_create(mpara)
!          write(*,*) "Mpara=",mpara
#ifdef fOMP
!$OMP    do schedule(dynamic,1) ordered
!----!$OMP    do schedule(dynamic,1)
#endif  

          DO seed_komp = 1, mpara
!       ### based on automatic differentiation (AD) ###
!           init all AD variables with zeros
            CALL g_initzero(ismpl)
#else
#ifdef fOMP
!$OMP    parallel num_threads(Tlevel_0) &
!$OMP      default(none) private(ismpl, seed_komp, i) &
!$OMP      shared(iter_inv, update_dd, mpara, ndata, Tlevel_0) &
!$OMP      shared(main_input, main_input_ad, main_output) &
!$OMP      shared(simtime_0, max_simtime, main_input_master, jac, jacT)
#endif
!         computing in parallel threads, sample index
          ismpl = omp_get_his_thread_num() + 1
!         sample init, each sample only one times
          CALL single_init(ismpl)
!$OMP     barrier
!
! !!! to avoid compiler bugs
!          CALL omp_ordered_create(ndata)!mpara)
#ifdef fOMP
!-----!$OMP    do schedule(dynamic,1) ordered
!$OMP    do schedule(static,1)
!schedule(dynamic,1)
#endif  

          DO seed_komp =  1, ndata
!            write(*,*) "seed_komp=",seed_komp
            CALL initzero_ad(ismpl)
#endif
!           sample init
            CALL prepare_realisation(ismpl)
!           init parameter
            CALL DCOPY(mpara,main_input_master,1,main_input(1,ismpl),1)
!
! --------
!           one derived function step
#ifndef AD_RM
!           get the seeding vector
            CALL seeding(seed_komp, g_main_input(1,ismpl), ismpl)
            CALL g_forward_compute(main_input(1,ismpl), g_main_input(1,ismpl), &
              main_output(1,ismpl), jac(1,seed_komp), &
              simtime_0,max_simtime,iter_inv,seed_komp,ismpl)

#else
!           get the seeding vector
            CALL seeding(seed_komp, main_input_ad(1,ismpl), ismpl)

            CALL forward_compute_ad(main_input, jacT(1,seed_komp), &
              main_output, main_input_ad(1,ismpl), &
              simtime_0,max_simtime,iter_inv,seed_komp,ismpl)
           do i=1,mpara
              if (param_enabled(i,ismpl)==2) then
                    jacT(i,seed_komp)=jacT(i,seed_komp)*main_input(i,ismpl)
             end if
          end do
#endif

! --------
!
            IF (update_dd) THEN
!           ### based on divided difference (DD) ###
!              save state of the gradients
              CALL gp_old_save(cgen_opti,ismpl)
!              compute DD-jacobian matrices
              CALL dd_iter(seed_komp,ismpl)
!              restore state of the gradients
              CALL g_old_restore(cgen_opti,ismpl)
            END IF            
!            write(*,*) "seed_komp=",seed_komp
          END DO
#ifdef fOMP
!$OMP    end do
#endif
! !!! to avoid compiler bugs
#ifndef AD_RM
          CALL omp_ordered_delete()
#endif

#ifdef fOMP
!$OMP    end parallel
#endif
          write_disable = .FALSE.
          write_iter_disable = .FALSE.

#ifdef AD_RM
        do i=1,mpara
            do j=1,ndata
                jac(j,i)=jacT(i,j)
            enddo
        enddo
#endif
!          WRITE(*,*) "jac=",jac
#ifdef BENCH
          CALL sys_cputime(tend_oa)
          WRITE(*,'(1A,F9.2,1A)') &
            '  [T] : full Jacobian matrix computation time =', &
            tend_oa - trun_oa, ' sec'
#endif
!         restore master state
          CALL DCOPY(mpara,main_input_master,1,main_input(1,idx_master),1)
!
        RETURN
      END
