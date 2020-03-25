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

!> @brief Set observations and observation variances for NS-EnKF
!> @details
!> Call the subroutine that sets the measurment values to the
!> regime of the normalized Gaussian. Then normalize obsvar_ns to
!> the observation variance of the original state vector variables
!> at the measurement locations.
subroutine enkf_normal_score_observations_obsvars()

  use mod_enkf, only:&
       obsvar_ns,&
       irobs,&
       nrobs_loc,&
       num_enkf_vars,&
       obs_ns_orig_vars,&
       act_o,&
       obsvar
  
  implicit none

  integer :: ienkfvar, iloc
  integer :: nrobs

  nrobs = nrobs_loc(irobs)

  
  allocate(obsvar_ns(nrobs,num_enkf_vars))

  do ienkfvar = 1, num_enkf_vars
     if(act_o(irobs,ienkfvar) == 1) then

        ! Normal score transform of measurement value
        call enkf_normalize_obs_mean_obsvars(ienkfvar)
        
        ! Normalization of obsvars_ns
        do iloc = 1, nrobs
           obsvar_ns(iloc,ienkfvar) = obsvar(ienkfvar)/obs_ns_orig_vars(iloc,ienkfvar)
        end do
        
     end if
  end do

end subroutine enkf_normal_score_observations_obsvars
