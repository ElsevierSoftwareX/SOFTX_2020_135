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

!> @brief Calculate submatrices of hybrid covariance matrix
!> @param[in] ivar Variable index
!> @param[in] nrobs Number of observation locations
!> @details
!> Calculate submatrices of hybrid covariance matrix. Theses are
!> the parts of the full covariance matrix used in the EnKF update.
subroutine enkf_hybrid_subcovs(ivar,nrobs)

  use mod_enkf, only:&
       pb,&
       pb_r,&
       pb_reps,&
       nstate,&
       act_s,&
       lstate,&
       i_obs,&
       j_obs,&
       k_obs,&
       irobs
  
  implicit none

  integer, intent(in) :: ivar
  integer, intent(in) :: nrobs

  integer :: imem, iobs, jobs, ihyb, jhyb

  integer :: ii,ij,ik,iivar, l
  
  integer, external :: index_loc_to_lin
  
  allocate(pb_r(nrobs,nrobs))
  allocate(pb_reps(nstate,nrobs))
  pb_r(:,:) = 0.0d0
  pb_reps(:,:) = 0.0d0
  
  do jobs = 1, nrobs
     ! Find linear index with i_obs,j_obs,k_obs
     l = index_loc_to_lin(i_obs(irobs,jobs),j_obs(irobs,jobs),k_obs(irobs,jobs))
     ! Calculate ihyb (full variable index)
     jhyb = (ivar-1)*lstate + l
     do iobs = 1, nrobs
        l = index_loc_to_lin(i_obs(irobs,iobs),j_obs(irobs,iobs),k_obs(irobs,iobs))
        ihyb = (ivar-1)*lstate + l
        pb_r(iobs,jobs) = pb(ihyb,jhyb)
     end do
  end do
  
  do jobs = 1, nrobs
     l = index_loc_to_lin(i_obs(irobs,jobs),j_obs(irobs,jobs),k_obs(irobs,jobs))
     jhyb = (ivar-1)*lstate + l
     do imem = 1, nstate
        ! sysindx not taken into account!
        call enkf_index_mem_to_loc(imem,ii,ij,ik,iivar)
        l = index_loc_to_lin(ii,ij,ik)
        ihyb = (iivar-1)*lstate + l
        pb_reps(imem,jobs) = pb(ihyb,jhyb)
     end do
  end do
  
end subroutine enkf_hybrid_subcovs
