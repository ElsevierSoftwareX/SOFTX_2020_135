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

!> @brief Set observation values in regime of normalized Gaussian
!> @param[in] ienkfvar EnKF variable index
!> @details
!> SET: Measurement values in var_obs in the regime fo the
!> normalized Gaussian.
!>
!> var_obs is set at each measurement location iloc (outer
!> loop).
!>
!> The inner loop runs over the ensemble members and the sorted
!> arrays obs_ns_original and obs_ns_gauss are used. In
!> obs_ns_original, we look for the first entry that is larger than
!> the measurement value. Depending on whether this is the first,
!> the last, or one of the intermediate values in obs_ns_orginial,
!> an interpolation for the measurement value is carried out. We
!> explain the procedure for intermediate values:
!>
!> In the double-precision variable fac, we normalize the distance
!> between the measurement value (var_obs-entry) and the largest
!> state vector value that is smaller than the measurement value
!> (obs_ns_original at iens-1) with respect to the difference
!> between the two state vector values surrounding the measurement
!> value (obs_ns_original at iens-1 and obs_ns_original at
!> iens).
!>
!> With this ratio, we set the new measurement value (var_obs) to
!> the corresponding value at the same ration of distances between
!> the normalized gaussian state vector points at the same rank as
!> the original state vector points (obs_ns_gauss at iens-1 and
!> obs_ns_gauss at iens).
subroutine enkf_normalize_obs_mean_obsvars(ienkfvar)

  use mod_enkf, only:&
       nrobs_loc,&
       irobs,&
       var_obs,&
       nrens,&
       obs_ns_original,&
       obs_ns_gauss
  
  implicit none

  integer, intent(in) :: ienkfvar

  integer :: iloc, iens

  double precision :: fac

  fac = 0.0d0
  
  do iloc = 1, nrobs_loc(irobs)
     do iens = 1, nrens
        
        if(obs_ns_original(iens,iloc,ienkfvar) > var_obs(irobs,iloc,ienkfvar)) then
           ! First value in obs_ns_original
           if(iens == 1) then
              fac = (obs_ns_original(iens,iloc,ienkfvar)-var_obs(irobs,iloc,ienkfvar))&
                   /(obs_ns_original(nrens,iloc,ienkfvar)-obs_ns_original(iens,iloc,ienkfvar))
              var_obs(irobs,iloc,ienkfvar) = obs_ns_gauss(iens,iloc,ienkfvar)&
                   - fac * (obs_ns_gauss(nrens,iloc,ienkfvar)-obs_ns_gauss(iens,iloc,ienkfvar))
           ! Intermediate value in obs_ns_original
           else
              fac = (var_obs(irobs,iloc,ienkfvar)-obs_ns_original(iens-1,iloc,ienkfvar))&
                   /(obs_ns_original(iens,iloc,ienkfvar)-obs_ns_original(iens-1,iloc,ienkfvar))
              var_obs(irobs,iloc,ienkfvar) = fac * obs_ns_gauss(iens,iloc,ienkfvar) &
                   + (1-fac) * obs_ns_gauss(iens-1,iloc,ienkfvar)
           end if
           exit
        ! Last value in obs_ns_original
        else if (iens == nrens) then
           fac = (var_obs(irobs,iloc,ienkfvar)-obs_ns_original(nrens,iloc,ienkfvar))&
                /(obs_ns_original(nrens,iloc,ienkfvar)-obs_ns_original(1,iloc,ienkfvar))
           var_obs(irobs,iloc,ienkfvar) = obs_ns_gauss(nrens,iloc,ienkfvar)&
                + fac * (obs_ns_gauss(nrens,iloc,ienkfvar)-obs_ns_gauss(1,iloc,ienkfvar))
           exit
        end if

     end do
  end do

end subroutine enkf_normalize_obs_mean_obsvars
