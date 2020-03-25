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

!> @brief Read the maximum number of observations
!> @details
!> __INPUT__: The number of observations in the observations file and
!> the name of the observation file.
!>
!> _Remark_: Use nrobs_int_pure, to go trough the
!> observation file once. For iterative/dual the nrobs_int_pure
!> observations are distrubted later to the nrobs_int observation times
!> of the scheme.
!>
!> __OUTPUT__: The maximum number of observation points at any
!> observation time.
subroutine enkf_read_obs_len()


  use mod_enkf, only: &
       nrobs_int_pure,&
       obs_name,&
       dual_enkf_switch,&
       max_obs_loc

  implicit none

  integer :: iobs, iloc
  
  integer :: obs_loc

  integer :: idu
  double precision    :: adu


  !open observation data file
  obs_name = trim(adjustl(obs_name))
  open(unit = 36, file = obs_name)

  ! Skip header
  read(unit = 36, fmt = *)

  ! Loop over observation intervals
  max_obs_loc = 0
  do iobs = 1, nrobs_int_pure
     read(unit = 36,fmt = *) idu, adu, obs_loc, idu, idu, idu, idu, idu, idu
     ! Save maximum number of observation locations
     max_obs_loc = max(max_obs_loc,obs_loc)
     ! Skip entries
     do iloc = 1, obs_loc
        read(unit = 36,fmt = *) idu, idu, idu, adu, adu, adu, adu, adu, adu
     end do
  end do

  close(unit = 36)



end subroutine enkf_read_obs_len
