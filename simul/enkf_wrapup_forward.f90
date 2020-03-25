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

!> @brief Wrapup the forward simulation
!> @details
!> Rename project to original string, set the timestep, which may
!> be used by the inverse computation, increase the assimilation
!> counter.
subroutine enkf_wrapup_forward()

  use mod_genrlc, only:&
       project

  use mod_time, only:&
       itimestep_0,&
       ntimestep,&
       max_simtime,&
       simtime_0

  use mod_enkf, only:&
       project_org,&
       itimestep_ismpl,&
       ttend,&
       iter_out

  implicit none

  !----------------------------------------------------
  !--------------------- PROJECT ----------------------
  !----------------------------------------------------
  project = project_org

  !----------------------------------------------------
  !--------------------- TIMESTEP ---------------------
  !----------------------------------------------------
  itimestep_0 = max(itimestep_0 +1, &
       int( dble(ntimestep -itimestep_ismpl) *ttend /(max_simtime - simtime_0) ) )

  !----------------------------------------------------
  !--------------------- ITER_OUT ---------------------
  !----------------------------------------------------
  ! increase assimilation counter
  iter_out = iter_out +1

end subroutine enkf_wrapup_forward
