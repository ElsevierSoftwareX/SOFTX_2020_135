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

!> @brief wrapper for ENKF/SIMUL initialisation and computation
!> @param[in] simtime_run simulation start time (begin intervall)
!> @param[in] simtime_end simulation end time (end intervall)
!> @param[in] iter_out iteration counter (ENKF, SIMUL)
!> @param[in] smode switch can contain "init" (for initialisation) and/or "run" (to run simulation); optional "mean" for special naming
!> @param[in] nrens number realisations/ensembles
!> @details
!> needed extra swaping, to solve parallel performance issue
      SUBROUTINE forward_multi_compute(simtime_run, simtime_end, iter_out, smode, nrens)
        use arrays
        use mod_genrl
        use mod_simul
        use mod_time
        use mod_OMP_TOOLS
        IMPLICIT NONE
        integer :: ismpl
        INCLUDE 'OMP_TOOLS.inc'
        DOUBLE PRECISION simtime_run, simtime_end
        integer  iter_out, nrens, irens
        character (len=*) :: smode
!       benchmark stuff
        DOUBLE PRECISION trun_oa, tend_oa
        INTRINSIC index

!
!       init master/backup copy - equal for all samples
        IF (index(smode,'init')>0) CALL DCOPY(mpara,main_input(1,idx_master),1,main_input_master,1)
#ifdef BENCH
        CALL sys_cputime(trun_oa)
#endif
!
!       sample computation, dynamic scheduling for OpenMP
#ifdef fOMP
!$OMP   parallel num_threads(Tlevel_0) &
!$OMP     default(none) private(ismpl, irens) &
!$OMP     shared(iter_out, nrens, mpara, nsmpl, Tlevel_0, smode) &
!$OMP     shared(main_input, main_output, main_input_master) &
!$OMP     shared(simtime_run, simtime_end, runmode) &
!$OMP     shared(vdefault)
#endif
!
#ifdef fOMP
!-- !$OMP do schedule(dynamic,1) ordered !!! to avoid compiler bugs !!!
          CALL omp_ordered_create(nrens)
#endif
!
          IF (index(smode,'init')>0) THEN
!           get memory one-times for each sample
#ifdef SIMUL_sgsim
            CALL alloc_sgsim()
#endif
#ifdef SIMUL_visim
            CALL alloc_visim()
#endif
!
#ifdef fOMP
!$OMP       do schedule(dynamic,1)
#endif

            DO irens = 1, nsmpl
!             first initialisation for each sample
              CALL single_init(irens)
            END DO
#ifdef fOMP
!$OMP       end do
!$OMP       barrier
#endif
          END IF
!
!         sample index - memory place for all computations of this thread
          ismpl = omp_get_his_thread_num() + 1
!         first work for this thread
          irens = ismpl
!
!       each thread-team computes one realisation (index [1 .. "Tlevel_0"])
!       -> to avoid state-swapping between currently parallel running threads
!         run simulation
          
          CALL omp_forward_multi_compute(main_input(1,ismpl), main_output(1,ismpl), &
            simtime_run, simtime_end, iter_out, smode, irens, nrens, ismpl)
!
#ifdef fOMP
!-- !$OMP do schedule(dynamic,1) ordered !!! to avoid compiler bugs !!!
!$OMP     do schedule(dynamic,1)
#endif
          DO irens = Tlevel_0 +1, nrens
!
!         now compute all other realisations (index ["Tlevel_0"+1 .. nrens])
!         -> state-swapping only between the last and the next realisation (not from a parallel thread)
!
!           swap state for realisations with index greater than "Tlevel_0"
!             ENKF needs lossless mode -> nrens==nsmpl!
            IF (runmode>=2) CALL swap_state(ismpl,irens,ismpl)
!
!           run simulation
            CALL omp_forward_multi_compute(main_input(1,ismpl), main_output(1,ismpl), &
              simtime_run, simtime_end, iter_out, smode, irens, nrens, ismpl)
!
!           swap back the state for realisations with index greater than "Tlevel_0"
            IF (runmode>=2) CALL swap_state(irens,ismpl,ismpl)
!
          END DO
#ifdef fOMP
!$OMP     end do
          CALL omp_ordered_delete()
!-- !$OMP end do !!! to avoid compiler bugs !!!
#endif
!
          IF (index(smode,'init')>0) THEN
!           free memory (for each sample)
#ifdef SIMUL_sgsim
            CALL dealloc_sgsim()
#endif
#ifdef SIMUL_visim
            CALL dealloc_visim()
#endif
          END IF
!
#ifdef fOMP
!$OMP   end parallel
#endif
!
#ifdef BENCH
        CALL sys_cputime(tend_oa)
        WRITE(*,'(1A,F9.2,1A)') '  [T] : over all **SIM computation time =', tend_oa - trun_oa, ' sec'
#endif
!
        RETURN
      END


!> @brief ENKF/SIMUL initialisation and computation
!> @param[in] dinput parameter vector for initial setup (dependency vector)
!> @param[out] doutput simulated/computed data/result values
!> @param[in] simtime_run simulation start time (begin intervall)
!> @param[in] simtime_end simulation end time (end intervall)
!> @param[in] iter_out iteration counter (ENKF, SIMUL)
!> @param[in] smode switch can contain "init" (for initialisation) and/or "run" (to run simulation)
!> @param[in] irens ensemble member/realization index
!> @param[in] nrens number realisations/ensembles
!> @param[in] ismpl local sample index
!> @details
      SUBROUTINE omp_forward_multi_compute(dinput, doutput, simtime_run, simtime_end, iter_out, smode, irens, nrens, ismpl)
        use arrays
        use mod_genrl
        use mod_genrlc
        use mod_simul
        use mod_time
        use mod_data
        use mod_OMP_TOOLS
        use mod_enkf, only:&
             enkf_log_out,&
             comp_real_out
        IMPLICIT NONE
        integer :: ismpl
        integer :: l
        INCLUDE 'OMP_TOOLS.inc'
!       input vector
        DOUBLE PRECISION dinput(mpara)
!       output vector
        DOUBLE PRECISION doutput(ndata)
!
        DOUBLE PRECISION simtime_run, simtime_end
        integer  iter_out, irens, nrens
        character (len=*) :: smode
        character (len=256) :: project_sfx_org
        INTRINSIC index, trim
!
! ---- stochastic initialisation
        IF (index(smode,'init')>0) THEN
          project_sfx_org = project_sfx(ismpl)
          project_sfx(ismpl) = trim(project_sfx_org)//'_init'
!         copy-in/initialise parameters for each realisation
          CALL prepare_realisation(ismpl)
!         init parameter
          CALL DCOPY(mpara,main_input_master,1,dinput,1)
!         generate parameter-initialisation/random-seed
          CALL simul_wrapper(dinput, ismpl, irens)
!         output of the physical values
          CALL forward_write(irens,ismpl)
!$OMP critical
          IF (index(smode,'enkf')>0 .and. enkf_log_out) &
               WRITE(37,*) '** ensemble created ** irens =', irens, 'of ', nrens
!$OMP end critical
          project_sfx(ismpl) = project_sfx_org
        END IF
! ---- end stochastic initialisation
!
! ---- progress realisation forward computation
        IF (index(smode,'run')>0) THEN
!         normal forward simulation for each realisation
!$OMP critical
          IF (index(smode,'enkf')>0 .and. enkf_log_out) &
               WRITE(37,*) '*** progress ensemble member ***', irens, &
               ' tstart  ', simtime_run/tunit, 'tend ', simtime_end/tunit
!$OMP end critical
          if(comp_real_out) then
             WRITE(*,'(1A,1I5,1A,1I4,1A,1I4)') &
                  '  [I] : compute realisation:', irens, ', sample index:', ismpl
          end if
          smon_idx(ismpl) = irens
!
!         forward computation
          CALL forward_compute(dinput, doutput, simtime_run, simtime_end, iter_out, irens, ismpl)
!
          l = irens
          IF (index(smode,'mean')>0 .AND. l==nsmpl) l=-4
!
!         keyword "verbose" disables the output
          IF (index(smode,'verb')<=0) THEN
!           output of the parameters
            CALL write_parameter(l,ismpl)
         IF (l>0) CALL write_hdfparameter2(l,ismpl)
!           output of the physical values           
              CALL forward_write(l,ismpl)
          END IF
        END IF
! ---- end forward computation
!
        RETURN
      END
