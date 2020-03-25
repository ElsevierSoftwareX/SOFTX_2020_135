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

!> @brief get the size of the current time step
!> @param[in] sim_time current simulation time
!> @param[in] ismpl local sample index
!> @return deltat of the time step
!> @details
!> Return the deltat of the time step depending on the simulation time
!> [sim_time] and the time step table [delta_time]\n\n
!>
!> 1. Input time stepping: Add delta_time values until the cumulative
!> sum minus 1/100 times the last deltat is larger than the simulation
!> time input. Then, the final deltat is the output of the
!> function. \n\n
!>
!> 2. Variable time stepping: Set (a) the starting time step in the
!> beginning, (b) go back to the last simtime and half the previous
!> time step, when the nonlinear maxiter was reached, or (c) keep the
!> time step or double it if everything was alright, then look for the
!> maximum simulation time or upcoming output times and adjust the
!> time step accordingly.
      double precision function deltat(sim_time,ismpl)

        use arrays
        use mod_genrl
        use mod_time

        implicit none

        double precision, intent (in) :: sim_time

        integer :: ismpl

        ! Counter
        integer :: i

        ! Tolerance for sim_time smaller timesum
        double precision :: tentol

        ! Cumulative sum of delta_time values
        double precision :: timesum

        ! Time step counter
        integer :: it

        ! Output time id
        integer :: id

        if (.not. delt_vary) then
          ! Time stepping from `# time periods`
          ! -----------------------------------

          deltat = 0.0d0

          ! time step counter
          it = 1

          ! init time step size
          deltat = delta_time(it)

          ! 1% numerical tolerance
          tentol = 0.01d0*deltat

          ! init time sum
          timesum = simtime_0 + deltat
!
          do while (sim_time + tentol >= timesum .and. it < ntimestep)
            it = it + 1
            ! update time step size
            deltat = delta_time(it)
            ! 1% numerical tolerance
            tentol = 0.01d0*deltat
            ! update time sum
            timesum = timesum + deltat
          end do

        else
          ! Variable time stepping.
          ! -----------------------

          if (sim_time == simtime_0 .and. flag_1st_timestep(ismpl) == 0) then

            ! Set starting time step
            deltat = delt_start
            delt_old(ismpl) = deltat
            delt_count(ismpl) = 0

          else

            ! Set previous time step
            deltat = delt_old(ismpl)

          end if

          if (flag_delt(ismpl) == -2) then
            ! Iterative solver or Picard/Newton iteration reached maxiter

            !Restore old simulation time
            simtime(ismpl) = simtime(ismpl) - delt_old(ismpl)

            flag_1st_timestep(ismpl) = 1

            write(*,*) "Halfing time step!"
            write(*,*) "Old deltat: ", deltat

            ! Half the time step size
            deltat = 0.5*deltat

            write(*,*) "New deltat: ", deltat

            ! Check whether minimum deltat is reached
            if (deltat < delt_min) then
              write(*,*) "Minimum step size reached. Aborting."
              write(unit = *, fmt = *) "Current deltat = ", deltat
              write(unit = *, fmt = *) "Minimum deltat = ", delt_min
              stop
            end if

          else if (flag_delt(ismpl) == 1) then
            ! Iterative solver or Picard/Newton iteration fine

            ! Counter closer to doubling
            delt_count(ismpl) = delt_count(ismpl) + 1

            ! Doubling time step size if delt_count hits delt_double
            if (delt_count(ismpl) == delt_double) then

              ! Reset doubling counter
              delt_count(ismpl) = 0

              write(*,*) "Doubling time step size!"
              write(*,*) "Old deltat: ", deltat

              ! Double the time step size.
              deltat = 2*deltat

              write(*,*) "New deltat: ", deltat

              ! Check whether maximum delt is reached
              if (deltat > delt_max) then
                write(*,*) "Maximum deltat reached. No further timestep size increase possible."
                deltat = delt_max
                write(*,*) "deltat: ", deltat
                write(*,*) "delt_max: ", delt_max
              end if
            end if

            ! Is max_simtime reached?
            ! Info: add delt_old, because simtime update not yet realized
            if ((simtime(ismpl) + delt_old(ismpl) + deltat) > max_simtime) then
              write(*,*) "Changing step size to reach simulation end."

              ! Set deltat to reach max_simtime
              deltat = max_simtime - simtime(ismpl)-delt_old(ismpl)

              write(*,*) "deltat: ", deltat

            end if

            ! Find the upcoming output time id
            id = 0
            do i = 1, noutt +1

              if (outt(i) > simtime(ismpl) + delt_old(ismpl) ) then

                id = i
                exit

              end if

            end do

            ! Change deltat to match output time if output time is in
            ! the next time step.
            if ((simtime(ismpl) + delt_old(ismpl) + deltat) > outt(id)) then

              ! Match only if the new deltat is larger than 10 times
              ! the minimum step size
              if ((outt(id) - simtime(ismpl) - delt_old(ismpl))  > 10.0d0*delt_min) then

                deltat = outt(id) - simtime(ismpl) - delt_old(ismpl)

                write(*,*) "Changing deltat to match output time."
                write(*,*) "deltat: ", deltat

              end if

            end if

          end if

          delt_old(ismpl) = deltat

        end if

        return

      end function deltat
