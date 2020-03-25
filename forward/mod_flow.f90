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

!>    @brief global variables for flow arrays
module mod_flow

      ! Flow: linear solver parameters
      ! ------------------------------

      !> @brief Flow linear solver error input.
      !> @details
      !> Flow linear solver error input. \n
      !> Read under `# lsolvef`, first entry. \n\n
      !> User in linear solver routines for break criteria.
      double precision :: errf

      !> @brief Flow linear solver explicit/implicit switch.
      !> @details
      !> Flow linear solver explicit/implicit switch. \n
      !> Currently set to default 1.0d0. \n\n
      !> Switching capability currently disabled.
      double precision :: aparf

      !> @brief Flow linear solver specification number.
      !> @details
      !> Flow linear solver specification number. \n
      !> Read under `# lsolvef`, second entry. \n\n
      !> Used in solve.f90 to specify solvername, solve criterion,
      !> preconditioner. \n\n
      !>
      !> [ ctrl = solvername + 16*criteria + 256*precondition ] \n
      !>     extract [ ctrl ] : \n
      !>      solvername = mod(ctrl,16) \n
      !>      ctrl = ctrl/16 \n
      !>      criteria = mod(ctrl,16) \n
      !>      ctrl = ctrl/16 \n
      !>      precondition = ctrl \n\n
      !>
      !> Examples:\n
      !> - 64: BiCGStab (parallel) solver, autodetected error, ILU
      !>        preconditioner
      !> - 67: PLU (serial) solver, autodetected error, ILU
      !>        preconditioner
      integer :: controlf

      !> @brief Flow linear solver maximum iteration number.
      !> @details
      !> Flow linear solver maximum iteration number. \n
      !> Read under `# lsolvef`, third entry. \n\n
      integer :: lmaxitf
!
      ! stopping criteria nonlinear outer loop, head or pressure
      ! --------------------------------------------------------

      !> @brief Flow nonlinear iteration tolerance.
      !> @details
      !> Flow nonlinear iteration tolerance. \n
      !> If the maximal difference between headold and head is smaller
      !> than this tolerance, the iteration converged.
      double precision :: nltolf

      !> @brief Flow nonlinear iteration relaxation.
      !> @details
      !> Flow nonlinear iteration relaxation. \n
      double precision :: nlrelaxf

      !> @brief Flow nonlinear iteration max.
      !> @details
      !> Flow nonlinear iteration max. \n
      !> Used in relaxation, `nl_relax`.
      double precision :: nlmaxf

!     stopping criteria nonlinear outer loop, saturation
      double precision nltols,nlrelaxs,nlmaxs
!     constant of gravitational force, compressibility of rock
      double precision grav,rref
!
!     Input and Output are in [MPa], internal computation in [Pa]
      double precision Pa_conv, Pa_conv1
        parameter (Pa_conv = 1.0d6)
        parameter (Pa_conv1 = 1.0d0 / Pa_conv)
!
end module mod_flow
