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

!> @brief Allocating arrays for EnKF update
!> @details
!> Allocating the state vector, average vector, variances vector,
!> system variances vector and, if needed, the correlation matrix
!> and the rhos-matrix for covariance localization.
subroutine enkf_alloc_state_vector()

  use mod_enkf, only: &
       mem,&
       ave,&
       var,&
       corr_matrix,&
       nstate,&
       nrens,&
       sysvarmem,&
       num_cov_ref,&
       nrobs_loc,&
       irobs,&
       rhos,&
       num_enkf_vars,&
       switch_cov_ref,&
       cov_loc_switch

  implicit none

  ALLOCATE(mem(nstate,nrens))
  ALLOCATE(ave(nstate))
  ALLOCATE(var(nstate))
  allocate(sysvarmem(nstate))

  if(switch_cov_ref) then
     allocate(corr_matrix(nstate,max(num_cov_ref,1)))
  end if

  if(cov_loc_switch) then
     allocate(rhos(nstate,nrobs_loc(irobs),num_enkf_vars))
  end if

end subroutine enkf_alloc_state_vector
