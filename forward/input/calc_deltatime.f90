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

!> @brief compute time steps table delta_time
!> @param[in] ismpl local sample index
!> @details
!> Compute the time steps table delta time. \n\n
!>
!> Linear: dt(j) = (tend-tstart)/n \n\n
!>
!> Logarithmic: dt(j) = tunit * (10**(j/n*((tstart-tend)/tunit + 1)) -
!> 1) \n\n
!>
!> Varied logarithmic (tstart needs to be non-zero): dt(j) = tstart *
!> (10**(j/n * tend/tstart) - 10**((j-1)/n * tend/tstart)) \n\n
!>
!> Logarithmic decreasing: dt(j) = tunit *
!> (10**((n-j)/n*((tstart-tend)/tunit + 1)) - 1) \n\n
subroutine calc_deltatime(ismpl)

  use arrays, only: delta_time
  use mod_genrl, only: memory
  use mod_time, only: nperiod, ntimestep, tunit, iperiod, dperiod

  implicit none

  ! local sample index
  integer :: ismpl

  ! Counter: overall time step
  integer :: it

  ! Number of time steps in current time periods
  integer :: nstep

  ! Counters: Period, time step of current period
  integer :: iper, istep

  ! Logarithmic factor, used in logarithmic time stepping
  double precision :: facl

  ! End (t1) and start (t0) time of time step
  double precision :: t1, t0


  ! Re-allocate delta_time to total number of time periods
  deallocate(delta_time)
  memory = memory - 1
  allocate(delta_time(ntimestep))
  memory = memory + ntimestep

  ! Set intitial time period counter
  it = 0

  ! Loop over time periods
  do iper = 1, nperiod

    ! Number of time steps in current time periods
    nstep = iperiod(iper,1)

    ! Compute time steps according to type
    if (iperiod(iper,2) == 1) then
      ! type: linear

      do istep = 1, nstep
        it = it + 1
        delta_time(it) = (dperiod(iper,2) - dperiod(iper,1)) / dble(nstep)
      end do

    else if (iperiod(iper,2) == 2) then
      ! type: logarithmically increasing

      facl = log10((dperiod(iper,2) - dperiod(iper,1))/tunit + 1.D0) / dble(nstep)

      t0 = dperiod(iper,1)
      do istep = 1, nstep
        it = it + 1
        t1 = dperiod(iper,1) + tunit * 10.0d0**(dble(istep) * facl) - tunit
        delta_time(it) = t1 - t0
        t0 = t1
      end do

    else if (iperiod(iper,2)==3) then
      ! type: varied logarithmically increasing

      if (dperiod(iper,1) <= 0.0d0) then
        write(unit = *, fmt = *) 'error (read_time.f90): ',&
            'No varied Log with starting time 0.0d0'
        stop 1
      end if

      facl = log10(dperiod(iper,2) / dperiod(iper,1)) / dble(nstep)

      t0 = dperiod(iper,1)
      do istep = 1, nstep
        it = it + 1
        t1 = dperiod(iper,1) * 10.0d0**(dble(istep) * facl)
        delta_time(it) = t1 - t0
        t0 = t1
      end do

    else if (iperiod(iper,2)==-2) then
      ! type: logarithmically decreasing

      facl = log10((dperiod(iper,2) - dperiod(iper,1))/tunit + 1.0d0) / dble(nstep)

      t0 = dperiod(iper,2)
      do istep = nstep - 1, 0, -1
        it = it + 1
        t1 = dperiod(iper,1) + tunit * 10.0d0**(dble(istep) * facl) - tunit
        delta_time(it) = t0 - t1
        t0 = t1
      end do

    else
      write(*,'(1A,1I2,1A,1I2,1A)') &
          'error: wrong time step type ', iperiod(iper,2), ' at ', iper, &
          ' in "read_time.f" !'
      stop

    end if

  end do

  ! Sanity check: Total number of time steps
  if (it /= ntimestep) then
    write(*,'(1A)') 'error: lost some time steps in "read_time.f" !'
    stop
  end if

  return

end subroutine calc_deltatime
