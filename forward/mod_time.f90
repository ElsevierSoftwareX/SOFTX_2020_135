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

!>    @brief global variables for timestepping
module mod_time

      !> @brief Time unit factor.
      !> @details
      !> Time unit factor. \n
      !> Default: 1.0d0. \n\n
      !>
      !> Read in under keyword `# tunit`. \n\n
      !>
      !> All time inputs are multiplied by this factor. Since the
      !> general computation expects seconds, the following
      !> correspondences hold: \n
      !> - tunit = 1.0d0: unit [s] \n
      !> - tunit = 60.0d0: unit [min] \n
      !> - tunit = 3600.0d0: unit [h] \n
      !> - tunit = 86400.0d0: unit [day] \n
      !> - tunit = 31536000.0d0: unit [a] with 365 days per year \n
      double precision :: tunit

      !> @brief Default time unit factor.
      !> @details
      !> Default time unit factor. \n
      !> Corresponds to time unit [s], seconds. \n\n
      double precision, parameter :: tunit_const = 1.0d0

      !> @brief Transient computation switch.
      !> @details
      !> Transient computation switch. \n
      !> `.true.`: transient computation. \n
      !> `.false.`: steady state computation.
      logical transient

      !> @brief Maximum number of time periods.
      !> @details
      !> Maximum number of time periods. \n
      integer, parameter :: npermax = 1056

      !> @brief Maximum simulation time.
      !> @details
      !> Maximum simulation time. \n\n
      !>
      !> Read under `# time periods`. Set to the end time of the
      !> last time period.
      double precision max_simtime

      !> @brief Start of simulation time.
      !> @details
      !> Start of simulation time. \n\n
      !>
      !> Read under `# timestep control` (second line, fourth entry)
      !> or under `# tstart`.
      !>
      !> "simtime_0" is needed for inverse-restart
      double precision simtime_0

      !> @brief Total number of time steps.
      !> @details
      !> Total number of time steps. \n
      !> Default: 0. \n\n
      !>
      !> If `# time periods` is read, the number of time steps is
      !> recorded in this variable. For variable time step size, the
      !> value of ntimestep is zero.
      integer ntimestep

      !> @brief Initial time step index.
      !> @details
      !>  Initial time step index. \n
      !> The initial time step index. Default: 0. \n\n
      !>
      !> For forward runs, itimestep_0 will set an offset in the time
      !> index. A simulation of 100 time steps, will then run from
      !> `itimestep_0` to `100+itimestep_0`. Otherwise, the
      !> time-handling will be equivalent.\n\n
      !>
      !> Read in under `# titer`. \n\n
      !>
      !> "itimestep_0" is needed for inverse-restart.
      integer itimestep_0

      !> @brief Number of time periods.
      !> @details
      !> Number of time periods. \n
      !> The number of time periods is not used if variable time step
      !> size is defined.
      integer :: nperiod

      !> @brief Time period specification array.
      !> @details
      !> Time period specification array. \n
      !> - iperiod(ip,1): number of time steps of period ip \n
      !> - iperiod(ip,2): time step size of period ip \n\n
      !>
      !> time step type: 1=linear, 2=logarithmic, 3=jlo (logarithmic
      !> variation), -2= logarithmic decreasing
      integer, allocatable, dimension (:,:) :: iperiod

      !> @brief Time period start and end array.
      !> @details
      !> Time period start and end array. \n
      !> - dperiod(ip,1): starting time of period ip \n
      !> - dperiod(ip,2): ending time of period ip \n\n
      !>
      !> Start and end time of the time period are given in time unit
      !> `tunit`, and then converted into seconds for the rest of the
      !> computation. \n
      double precision, allocatable, dimension (:,:) :: dperiod

      !> @brief Static relaxation for flow.
      !> @details
      !> Static relaxation for flow. \n
      !>
      !> Read under `# timestep control`, first entry.
      double precision :: thetaf

      !> @brief Static relaxation for temperature.
      !> @details
      !> Static relaxation for temperature. \n
      !>
      !> Read under `# timestep control`, second entry.
      double precision :: thetat

      !> @brief Static relaxation for concentration.
      !> @details
      !> Static relaxation for concentration. \n
      !>
      !> Read under `# timestep control`, third entry.
      double precision :: thetac
!
      logical bctp
!     number of bc-tp entries
      integer ngsmax
      integer ngststep, nbctp
!

      !> @brief Switch for monitor output.
      !> @details
      !> Switch for monitor output. \n\n
      logical :: monitor

      integer nmon,imon(npermax,4)

      !> @brief Number of output times.
      !> @details
      !> Number of output times. \n
      !> Records entry under `# output times`.
      integer noutt

end module mod_time
