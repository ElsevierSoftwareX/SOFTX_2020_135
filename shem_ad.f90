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

!>    @brief main program deterministic inversion (based on AD)
!>    @details
!>
!> **SHEMAT-Suite (Simulator for HEat and MAss Transport)** is a
!> numerical code for computing flow, heat and species transport
!> equations in porous media. The governing equations of the code are
!> the groundwater flow equation, the heat transport equation and the
!> species transport equation. \n\n
!>
!> SHEMAT-Suite includes parameter estimation and data assimilation
!> approaches, both stochastic (Monte Carlo, ensemble Kalman filter)
!> and deterministic (Bayesian inversion using automatic
!> differentiation for calculating derivatives).\n\n
!>
!> SHEMAT-Suite is based on numerical modules developed by: \n
!> -  Jörn Bartels     (gtn neubrandenburg gmbh)\n
!> -  Michael Kühn     (tu hamburg-harburg)\n
!> -  Roland Wagner    (ag,rwth aachen)\n
!> -  Martin H. Bücker   (sc,rwth aachen)\n
!> -  Christoph Clauser (ag,rwth aachen)\n\n\n
!>
!> Note: To be able to use input file parsing with hdf5, the
!> hdf5-input-files have to be generated using the script:
!> `convert_to_hdf5.py`. This script can be found in the repository
!> `SHEMAT-Suite_Scripts` under
!> `python/preprocessing/convert_to_hdf5.py`.
      PROGRAM shem_ad
        use arrays
        use mod_genrl
        use mod_genrlc
        use mod_time
        use mod_data
        use mod_linfos
        use mod_OMP_TOOLS
        use mod_inverse
#ifndef noHDF
        use mod_input_file_parser_hdf5
#endif
        IMPLICIT NONE
        integer :: ismpl
        integer :: i
        double precision :: tsglobal, tfglobal
        double precision :: tflocal
        INCLUDE 'version.inc'
        INCLUDE 'OMP_TOOLS.inc'

        CHARACTER*80 filename
        DOUBLE PRECISION qval_old, rms_para_old, qval_new
        INTEGER start_iter_inv
        INTEGER lblank, max_it, criteria
        INTEGER startt,endt,rate,startloop,endloop,startfw,endfw
        integer startprepare,endprepare
        DOUBLE PRECISION qfunc
!     time stuff
!     if "restart" option is used
        LOGICAL test_option
        EXTERNAL qfunc, lblank, test_option
!     time of the last saving of the restart data
        DOUBLE PRECISION backuptime, depsilon


        CALL sys_cputime(tsglobal)
        backuptime = tsglobal
        write_disable = .FALSE.
        write_iter_disable = .FALSE.
        nested_build = .TRUE.
        ismpl = 1
        def_binary = 'inverse'
        lib_override = .FALSE.

        WRITE(*,*) ' '
        WRITE(*,*) '======================================'
        WRITE(*,*) version
        WRITE(*,*) '======================================'
        WRITE(*,*) ' '

        OPEN(66,file='shemade.job')

! -----------------------------------------

!     read new model
10      READ(66,'(A)',end=99999) filename
        CALL sys_cputime(tslocal)


        IF (filename(1:1)=='!') GO TO 10
        IF (filename(1:1)=='%') GO TO 10
        IF (filename(1:1)=='#') GO TO 10

        WRITE(*,*) ' '
        WRITE(*,*) ' '
        WRITE(*,*) ' *** NEW MODEL '

        project = filename
        runmode = 0
        transient = .FALSE.
        nperiod = 0
        is_init_flow_trafo_needed = .true.

!     no DD-mode
        update_dd = .FALSE.

#ifndef noHDF
        CALL h5parse_open_datafile(filename)
#endif

!     read forward model
        CALL read_model(filename,ismpl)
        CALL read_control(filename,ismpl)

!     read timestep parameter
        CALL read_time(filename,ismpl)

!     read data
        IF (runmode>0) CALL read_data(filename_data,ismpl)
        IF (runmode>1) CALL read_inverse(filename_inverse,ismpl)
!     split units
        CALL read_split_inv(filename,ismpl)

!     setup log file name
        status_log = filename(1:lblank(filename)) // '_status.log'
        status_log_inv = filename(1:lblank(filename)) // &
          '_status-inv.log'
        restart_name = filename(1:lblank(filename)) // '_restart'

!     initialisation for different solver (optimisation methods)
        CALL init_step()

        start_iter_inv = 1
        iter_inv = 1
        IF (test_option('restart')) THEN
!       read all necessary data from last save point
          CALL read_restartfw(restart_name,simtime(ismpl),itimestep_0, &
            ismpl)
          IF (runmode>=2) THEN
            CALL read_restartinv(restart_name,start_iter_inv,ismpl)
            iter_inv = start_iter_inv - 1
          END IF

!       add output for restarted optimization observation
          WRITE(*,'(3A)') '  [W] : "', status_log(1:lblank(status_log &
            )), '"'
          OPEN(76,file=status_log,status='unknown',position='append')
          WRITE(76,'(3A)') '%    restart Project: "', &
            filename(1:lblank(filename)), '"'
          WRITE(76,'(1A)') '%'
          IF (transient) THEN
            WRITE(76,'(2A)') '% <time step>', ' <cum time>'
          END IF
#ifdef head_base
          WRITE(76,'(4A)') '%    <iteration>',' <delta head>',' <delta temp>',' (<delta conc> ...)'
#endif
#ifdef pres_base
          WRITE(76,'(4A)') '%    <iteration>',' <delta pres>',' <delta temp>',' (<delta conc> ...)'
#endif
          CLOSE(76)
          IF ((ndata/=0) .AND. (mpara/=0) .AND. (runmode>=2)) THEN
            OPEN(76,file=status_log_inv,status='unknown', &
              position='append')
            WRITE(76,'(4A14,2A7,A11)') '% objec.func. ', &
              '         theta_obs  ', '    theta_par  ', &
              '       relax  ', ' inv_it', '    rel_it', ' time [min]'
            CLOSE(76)
          END IF
        ELSE
!       add output for new optimization observation
          WRITE(*,'(3A)') '  [W] : "', status_log(1:lblank(status_log &
            )), '"'
          OPEN(76,file=status_log)
          WRITE(76,'(2A)') '% Shemat-Suite version: ', version
          WRITE(76,'(2A)') '%          build: ', datum
          WRITE(76,'(2A)') '%   build command line: ', makecmd
          WRITE(76,'(3A)') '%        Project: "', &
            filename(1:lblank(filename)), '"'
          WRITE(76,'(1A)') '%'
          IF (transient) THEN
            WRITE(76,'(2A)') '% <time step>', ' <cum time>'
          END IF
#ifdef head_base
          WRITE(76,'(4A)') '%    <iteration>',' <delta head>',' <delta temp>',' (<delta conc> ...)'
#endif
#ifdef pres_base
          WRITE(76,'(4A)') '%    <iteration>',' <delta pres>',' <delta temp>',' (<delta conc> ...)'
#endif
          CLOSE(76)
          IF ((ndata/=0) .AND. (mpara/=0) .AND. (runmode>=2)) THEN
            OPEN(76,file=status_log_inv)
            WRITE(76,'(4A14,2A7,A11)') '% objec.func. ', &
              '         theta_obs  ', '    theta_par  ', &
              '       relax  ', ' inv_it', '    rel_it', ' time [min]'
            CLOSE(76)
          END IF
        END IF
#ifndef noHDF
        CALL h5parse_close_datafile()
#endif
!
!     forward modeling (initialization)
        CALL forward_init(ismpl)

!       forward modeling (nonlinear iteration)
        call system_clock(startfw,rate)
        CALL forward_iter(simtime_0,max_simtime,iter_inv,0,ismpl)
        call system_clock(endfw)
        write(*,*) "Time for forward_iter:",REAL(endfw-startfw)/REAL(rate)

!     --- RUNMODE ? ---
        IF (runmode>=2) THEN
!-------
          CALL prepare_jacobian(ismpl)
          qval_old = qfunc(linlog_input,apri_input,ismpl)
          WRITE(*,'(A,1e16.8)') '  objective function,     value=', qval_old
          CALL datapara_output(ndata,mpara,iter_inv,rms_data,rms_para,ismpl)
          CALL write_data(iter_inv,ismpl)
          CALL para_write(runmode,iter_inv,start_iter_inv,ismpl)
!--------

! ######################################################################
!   begin iteration loop
! ######################################################################
          DO iter_inv = start_iter_inv, maxiter_inv

            call system_clock(startloop,rate)
!           begin inversion
            WRITE(*,*) ''
            WRITE(*,*) ''
            WRITE(*,*) '*** Inverse iteration number ', iter_inv, &
              ' ***'

!           write out the last forward results
            IF (mod(iter_inv,iter_out)==0 .OR. iter_inv==maxiter_inv) THEN
              CALL forward_write(iter_inv,ismpl)
              IF (runmode>0) CALL write_data(iter_inv,ismpl)
            END IF
!           but this for each iteration
            CALL para_write(runmode,iter_inv,iter_inv,ismpl)

!           write backup of all necessary data after half an hour (1800se
!AW--           CALL sys_cputime(tfglobal)
!AW--           IF (tfglobal-backuptime.gt.1800.d0) THEN
!AW--             CALL resforward_write(restart_name, simtime(ismpl),itimestep<- geht nicht mehr,ismpl)
!AW--             CALL resinverse_write(restart_name, iter_inv,ismpl)
!AW--             backuptime = tfglobal
!AW--           END IF

!           transient   : use state from file init (necessary)
!           steady-state: use state from the last computation (optimisation)
            CALL old_restore(cgen_opti,ismpl)

!          ***********************************************
#ifndef JACOBI_FREE
            write(*,*) "full Jacobi compute"
            !           full Jacobi matrix computation - one times
            startprepare = 0
            endprepare = 0
            call system_clock(startprepare)
            CALL jacobi_compute(ismpl)
            call system_clock(endprepare)
#endif
!          ***********************************************

            write(*,*) "Time for full jacobi:",REAL(endprepare-startprepare)/REAL(rate)
            !           -- do normal iteration (optimization step's) --
            IF (mod(iter_inv,iter_out)==0 .OR. (iter_inv>=maxiter_inv)) THEN
              CALL write_jacw(ismpl)
            END IF

            CALL prepare_jacobian(ismpl)
!
!           get the quality function value for position p_k+1, (destroy all "tmp" vectors)
            IF (iter_inv==1) THEN
!             first iteration: compute the quality function
              qval_old = qfunc(linlog_input,apri_input,ismpl)
              rms_para_old = rms_para
            ELSE
!             later iterations: copy from the last iteration (before)
              qval_old = qval_new
              rms_para_old = rms_para
            END IF

!           only show values of the initialization
            IF ((iter_inv==1) .AND. (ndata/=0) .AND. (mpara/=0)) THEN
              WRITE(*,'(3A)') '  [W] : "', status_log_inv(1:lblank( &
                status_log_inv)), '"'
              OPEN(76,file=status_log_inv,status='unknown', &
                position='append')
              WRITE(76,'(4(e14.6),2(1X,I5,1X),F11.2)') qval_old, &
                rms_data,rms_para,1.D0, 0, 0, 0.D0
              CLOSE(76)
            END IF

call system_clock(startt,rate)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           Bayesian inversion step('s)
!             p^aprio == apri_input
!             p^k     == linlog_input

            CALL inv_step(qval_new,qval_old,linlog_input,apri_input,grad, &
              grad_sec,tmp_vec,tmp_mat,ismpl)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            call system_clock(endt,rate)

            CALL datapara_output(ndata,mpara,iter_inv,rms_data, &
              rms_para,ismpl)

            CALL data_output(ndata,iter_inv,maxiter_inv,tol_inv, &
              rms_data,ismpl)
            CALL prepare_jacobian(ismpl)

            IF ((covar/=0 .OR. resmat/=0) .AND. iter_inv/=0) THEN
              CALL calc_covares(mpara,tmp_vec,ismpl)
              CALL write_covariances(iter_inv,ismpl)
              CALL write_resolution(iter_inv,ismpl)
            END IF
            call system_clock(endloop)

!           inverse model converged?
            IF ((dsqrt(rms_data)<=tol_inv*dble(ndata)) .OR. (iter_inv>=maxiter_inv) .OR. abort_inv) THEN 
               GO TO 88888
            ELSE 
               write(*,*) "rms_data=",dsqrt(rms_data), " >",tol_inv*dble(ndata)
            endif
            write(*,*) "Time for inv_step: ",REAL(endt-startt)/REAL(rate)
            write(*,*) "Time for loop: ",REAL(endloop-startloop)/REAL(rate)

! ######################################################################
!   end inversion loop
! ######################################################################
          END DO


        END IF
!     --- end RUNMODE ---
88888   CONTINUE

!     write out the last iteration
        CALL forward_write(-1,ismpl)
        IF (runmode>0) CALL write_data(-1,ismpl)
        CALL para_write(runmode,-1,iter_inv,ismpl)
        CALL compare_run(ismpl)

        IF (runmode>=2) THEN
!         restore the state of PHI & TEMP of the derivation
          CALL old_restore(cgen_opti,ismpl)
#ifndef MF_STBAY
          CALL write_inv_hdf(ismpl)
#endif
        END IF

!     --- free memory
        CALL dealloc_arrays(ismpl)
        CALL props_end(ismpl)
        CALL dealloc_data(ismpl)

        IF (runmode>1) THEN
          CALL dealloc_inverse(ismpl)
#ifdef AD 
          CALL g_dealloc_arrays(ismpl)
#endif
#ifdef AD_RM
          CALL dealloc_arrays_ad(ismpl)
#endif
!AW-noexist-   call g_props_end
        END IF
!     ---

        WRITE(*,'(3A)') '  [W] : "', status_log(1:lblank(status_log)),'"'
        OPEN(76,file=status_log,status='unknown',position='append')
        CALL sys_cputime(tflocal)
        WRITE(76,'(A,F11.2,A)') '% job finished, total cpu time: ', &
          (tflocal-tslocal)/60.0D0, ' min'
        CLOSE(76)

!     Another file to load?
        GO TO 10

! -----------------------------------------

!     finis terrae
99999   CONTINUE
!
        CLOSE(66)
        CALL sys_cputime(tfglobal)
        tslocal = tfglobal - tsglobal
        i = tslocal/60.D0
        WRITE(*,'(1A,1I4,1A,1F5.2,1A)') ' total cpu time: ', i, ':', &
          tslocal - dble(i)*60.D0, ' min'
        WRITE(*,*) 'RUN O.K.'

      END
