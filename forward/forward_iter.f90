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

!> @brief time discretisation loop
!> @param[in] iter_out inverse iteration, SM realisation
!> @param[in] simtime_run start time of the simulation
!> @param[in] simtime_end finish time of the simulation
!> @param[in] iseed 0: FW simulation, 1 .. <mpara>: AD seeding index
!> @param[in] ismpl local sample index
!> @details
!> In-a-Nutshell description of this subroutine: \n
!> - Preprocessing before time step loop, initial variable values,
!>   monitoring output, extra steady-state initialisation, status_log\n
!> - Time loop: \n
!>   - before computations: time stepping, saving old variable arrays,
!>     output \n
!>   - computation: calling `forward_wrapper.f90`
!>   - after computations: save simulated data, update simtime,
!>     output, check divergence for variable step size \n
!> - Postprocessing: standard output
      subroutine forward_iter(simtime_run,simtime_end,iter_out,iseed,ismpl)

        use arrays
        use mod_genrl
        use mod_genrlc
        use mod_time
        use mod_linfos

        implicit none

        ! local sample index
        integer :: ismpl

        ! iter_out: inverse iteration, SM realisation
        INTEGER iter_out, iseed

        ! Time step index
        integer :: itimestep

        ! Size of a time period
        double precision :: deltt

        ! Start time of the simulation
        double precision, intent (in) :: simtime_run

        ! Finish time of the simulation
        double precision, intent (in) :: simtime_end

        double precision, external :: deltat
        integer, external :: lblank


        ! Preprocessing
        ! -------------

        ! initial values for some variables/arrays
        flag_1st_timestep(ismpl)=0
        itimestep = itimestep_0
        simtime(ismpl) = simtime_run
        deltt = deltat(simtime(ismpl),ismpl)
        tr_switch(ismpl) = .true.
        iter_nlold = maxiter_nl/2

        ! initial monitoring output
        if (transient .and. monitor .and. simtime_run == simtime_0) then
          call write_monitor(1,ismpl)
          call write_monitor_user(1,ismpl)
        end if

        ! runmode 2: extra steady state initialisation
        if (transient .and. runmode == 2) then

          tr_switch(ismpl) = .false.
          if (iseed == 0 .and. linfos(2) >= 0) write(*,'(1A)') '  [I] : extra steady state initialisation'

          call forward_wrapper(itimestep,iseed,ismpl)

          if (iseed == 0 .and. linfos(2) >= 0) write(*,'(1A)') '  [I] : normal transient process'
          tr_switch(ismpl) = .true.

        end if

        ! Write to status_log
        if (transient .and. iseed == 0 .and. (.not. write_iter_disable)) then

          open(76, file=status_log, status='unknown', position='append')
          write(76, fmt='(I8,1e14.6,1e14.6)') itimestep, deltt, simtime(ismpl)/tunit
          close(76)

        end if

        ! Time step loop for forward modeling
        ! -----------------------------------
1000    CONTINUE

        if (transient) then

          ! Advance time step
          itimestep = itimestep + 1

          ! Initialize flag for variable time step size
          flag_delt(ismpl) = 0

          ! Time stepping info to standard out
          if (linfos(1)>=1) then
            write(*,*) ' '
            write(*,'(1A,1I6)') '  >>>> new time step: ', itimestep
            write(*,'(1A,1e16.8,1A,1e16.8)') '  >>>>     cum. time= ', &
                (simtime(ismpl)+deltt)/tunit, '/', max_simtime/tunit
            write(*,'(1A,1e16.8)') '  >>>>     time step= ', &
                (deltt)/tunit

          end if

          ! Save old time level
          call old_save(cgen_time,ismpl)

        end if

! ######### Forward Iteration ######
        call forward_wrapper(itimestep,iseed,ismpl)
! ##################################

        ! save and collect the computed values for:
        ! - comparison with 'ddata(*,cid_pv)' (observed data)
        ! - data-output (write_data.f)
        call save_data(ismpl)

        if (transient) then

          ! Update simulation time
          simtime(ismpl) = simtime(ismpl) + deltt

          ! wrapper for output
          call write_outt(deltt,ismpl)

          ! monitoring output
          if (monitor .and. flag_delt(ismpl) /= -2) then
            call write_monitor(2,ismpl)
            call write_monitor_user(2,ismpl)
          end if

          ! Write to status_log
          if ( .not. write_iter_disable .and. iseed == 0) then

            ! Status log info to standard out
            if (linfos(1)>=1) then
              write(*,'(3A)') '  [W] : "', status_log(1:lblank(status_log)), '"'
            end if

            open(76, file=status_log, status='unknown', position='append')
            write(76, fmt='(I8,1e14.6,1e14.6)') itimestep, deltt, simtime(ismpl)/tunit
            close(76)

          end if

          ! Check for variable time stepping divergence flag
          if (flag_delt(ismpl) == -2) then
            ! Restore old values of variable arrays if time step was
            ! halfed (restoring simtime is handled in "deltat")
            call old_restore(cgen_time, ismpl)
          end if

          ! Set variable time stepping divergence flag to zero to
          ! avoid double calling of deltat.
          flag_delt(ismpl) = 0
          deltt = deltat(simtime(ismpl),ismpl)

        end if

        ! Important: This if statement (with the goto) needs to be
        ! outside of the other "transient" scopes, to generate
        ! reverse-mode code!
        if (transient .and. simtime(ismpl) < simtime_end) go to 1000
! --------------- 1000: return to next time step ----------------------

        ! Postprocessing
        ! --------------

        ! Write to status_log and status_log_inv
        if ((transient) .and. (runmode >= 2) .and. .not. write_iter_disable .and. iseed == 0) then
          write(*,'(3A)') '  [W] : "', status_log(1:lblank(status_log)), '"'

          open(76,file=status_log,status='unknown',position='append')
          write(76,'(1A)') key_char//' transient end'
          close(76)

          open(76,file=status_log_inv,status='unknown', position='append')
          write(76,'(1A,I8,1e14.6,1e14.6)') key_char//' transient: ', itimestep, deltt, simtime(ismpl)/tunit
          close(76)
        end if

        ! Standard output: steady state
        if (iter_out > 0 .and. linfos(1) >= 1 .and. .not. transient .and. iseed == 0) then
          write(*,'(29X,1A)') ' ===> leaving nonlinear iteration'
        end if

        ! Standard output: transient
        if (linfos(1) >= 0 .and. transient .and. iseed == 0) then
          write(*,'(1A,I8,1A,1e14.6)') '  [I] : final time step = ', &
              itimestep,', simulation time', simtime(ismpl)/tunit
        end if

        return

      end subroutine forward_iter
