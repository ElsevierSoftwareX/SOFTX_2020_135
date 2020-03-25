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

!> @brief Prepare forward run
!> @details
!> First, the start time and end time of the forward simulation are
!> set. Then the subroutine that updates the parameters is
!> called. Finally, the project directory name is adapted.
subroutine enkf_prepare_forward()

  use arrays, only:&
       simtime

  use mod_genrlc, only:&
       project

  use mod_simul, only:&
       ssample_outdir

  use mod_time, only:&
       simtime_0,&
       tunit
       

  use mod_enkf, only:&
       nrens,&
       irobs,&
       nrobs_int,&
       obst,&
       project_org,&
       ttstart,&
       ttend,&
       iter_out,&
       sflags

  implicit none

  integer :: irens

  !--------------------------------------------------------
  !------------------------ TTSTART -----------------------
  !--------------------------------------------------------
  if( irobs == 1 ) then
     ! FIRST ttstart 
     ttstart = simtime_0

     ! INITIALIZE simtime, iter_out, obst(0)
     DO irens = 1, nrens
        simtime(irens) = simtime_0
     END DO
     iter_out = 0
     obst(0) = 0.D0
  else
     ! UPDATE ttstart
     ttstart = simtime(1)
  end if

  !--------------------------------------------------------
  !------------------------- TTEND -----------------------
  !--------------------------------------------------------
  ! ttend at observation time
  ttend = obst(irobs)*tunit

  !--------------------------------------------------------
  !------------------- UPDATE PARAMETERS-------------------
  !--------------------------------------------------------
  ! Update parameter space
  call enkf_update_parameters()

  !--------------------------------------------------------
  !------------------------- SFLAGS -----------------------
  !--------------------------------------------------------
  sflags = 'enkf+run+verbose'
  IF (irobs == nrobs_int) sflags = 'enkf+run'

  !--------------------------------------------------------
  !------------------------- PROJECT ----------------------
  !--------------------------------------------------------

  ! Use different project and store old one
  project_org = project
  project = trim(ssample_outdir)//trim(project_org)

end subroutine enkf_prepare_forward
