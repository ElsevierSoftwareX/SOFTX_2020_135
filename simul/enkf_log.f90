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

!> @brief EnKF output written to enkf.log
!> @param[in] i_output index of different enkf.log outputs
!> @param[in] assim_id Observation variable index
!> @details
!> In `enkf.log` some output is saved.
subroutine enkf_log(i_output, assim_id)

  use mod_simul, only:&
       senkf_outdir

  use mod_enkf, only:&
       nrens,&
       obst,&
       irobs,&
       enkf_log_out,&
       enkf_variable_names,&
       normobs

  implicit none

  integer, intent(in) :: i_output
  integer, intent(in) :: assim_id

  integer :: irens

  !-------------------------- 0 -----------------------
  if (i_output == 0) then

     OPEN(37,file=senkf_outdir//'enkf.log')  

  !-------------------------- 1 -----------------------
  else if (i_output == 1 .and. enkf_log_out) then

     WRITE(37,*) 'ENKF - control file'
     WRITE(37,*) '==================='
     WRITE(37,*)     
     WRITE(37,*) &
          '           *** ensemble creation completed ***'
     WRITE(37,*)
     WRITE(37,*) 'number of ensemble members    nrens:', nrens
     WRITE(37,*)

  !-------------------------- 2 -----------------------
  else if (i_output == 2 .and. enkf_log_out) then

     WRITE(37,*)
     WRITE(37,*)
     WRITE(37,*)
     WRITE(37,*) '*** start assimilation at time  ', &
          obst(irobs)

  !-------------------------- 3 -----------------------
  else if (i_output == 3 .and. enkf_log_out) then
     WRITE(37,*) 'Assimilation type sequential'

  !-------------------------- 4 -----------------------
  else if (i_output == 4 .and. enkf_log_out) then
     WRITE(37,*) '*** enter assimilation for ', enkf_variable_names(assim_id), ' ***'

  !-------------------------- 5 -----------------------
  else if (i_output == 5 .and. enkf_log_out) then
     WRITE(37,*)'Assimilation type combined   scaling: ',normobs
     WRITE(37,*)
     WRITE(37,*) '*** enter assimilation ***'

  !-------------------------- 6 -----------------------
  else if (i_output == 6 .and. enkf_log_out) then
     WRITE(37,*) '           *** recompute ensembles ***'
     WRITE(37,*)

  !-------------------------- 7 -----------------------
  else if (i_output == 7) then

     close(unit = 37)

  else if (enkf_log_out) then
     write(unit = *, fmt = *) '[E] Error in enkf_log.f90'
     stop 1
  end if

end subroutine enkf_log
