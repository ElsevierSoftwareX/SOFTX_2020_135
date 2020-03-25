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

!> @brief Wrapper routine for SGSim Intialization
!> @param[in] iter_enkf EnKF iteration counter
!> @details
!> First the number of realisations is set, then some
!> project-directories are renamed. Finally, the main subroutines
!> of the SGSim Intialization are called. Possibly a given variable
!> distribution is imposed.
subroutine enkf_simulate(iter_enkf)

  use arrays, only:&
       project_sfx,&
       smon_idx

  use mod_genrl, only:&
       nsmpl

  use mod_genrlc, only:&
       project

  use mod_simul, only:&
       sm_max,&
       ssample_outdir

  use mod_enkf,only:&
       ref_dist_switch,&
       project_org,&
       nrens

  implicit none

  integer, intent(in) :: iter_enkf

  integer :: irens
  
  
  !--- no. of realisations ---
  nrens = sm_max

  ! Specify project suffixes
  DO irens = 1, nsmpl
     project_sfx(irens) = '_E0'
  END DO

  ! Add output directory to project and save the old value
  ! in project_org
  project_org = project
  project = trim(ssample_outdir)//trim(project_org)


  !----------------------------------------------------
  !----------------------------------------------------
  ! SGSIM-INIT for EnKF
  ! Fill propunit
  CALL forward_multi_compute(0.0d0, 0.0d0, iter_enkf, 'enkf+init', nrens)
  !----------------------------------------------------
  !----------------------------------------------------

  ! initialize mean
  smon_idx(nsmpl) = 0
  CALL prepare_realisation(nsmpl)


  ! Put project back to the original value
  project = project_org

  !----------------------------------------------------
  !----------------------------------------------------
  !Fills head/propunit with samples

  if(ref_dist_switch) then
     call enkf_variable_dist()
  end if


  !----------------------------------------------------
  !----------------------------------------------------




end subroutine enkf_simulate
