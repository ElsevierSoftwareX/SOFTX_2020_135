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

!> @brief Set the hybrid covariance matrix
!> @details
!> Set the hybrid covariance matrix
subroutine enkf_hybrid_cov()

  use mod_enkf, only:&
       pb,&
       nstate,&
       lstate,&
       mem,&
       nrens,&
       ave,&
       hybrid_cov_kind,&
       num_enkf_vars
  
  implicit none

  integer :: nhyb, ihyb, jhyb, iloc, jloc
  integer :: ii,ij,ik,ivar
  integer :: ji,jj,jk,jvar
  double precision :: hybrid_var
  double precision, allocatable :: A(:,:)
  integer :: i_ens

  ! nhyb: All Variables are needed, because the observations are not
  ! necessarily in the EnKF-state vector
  nhyb = num_enkf_vars*lstate
  
  allocate(pb(nhyb,nhyb))

  hybrid_cov_kind = 1
  if(hybrid_cov_kind==1) then
     ! Given variances for the different variables and no covariances
     ! as the static covariance in hybrid EnKF.
     do ihyb = 1, nhyb
        do jhyb = 1, nhyb
           if (ihyb==jhyb) then
              ! Get the variances of the variable/parameter
              ! corresponding to the hybrid index
              call enkf_hybrid_var_covs(ihyb,hybrid_var)
              pb(ihyb,jhyb) = hybrid_var
           else
              pb(ihyb,jhyb) = 0.0d0
           end if
        end do
     end do

  else if (hybrid_cov_kind==2) then
     ! Incoming state vector covariance matrix at first assimilation
     ! step as static covariance in hybrid EnKF.
     !
     ! Deprecated!!!!!!!!!!!!!!!!!!! 
     allocate(A(nstate,nrens))
     do i_ens = 1, nrens
        A(:,i_ens) = mem(:,i_ens) - ave(:)
     end do

     call dgemm('n','t',nstate,nstate,nrens,1.0d0/dble(nrens-1),A,nstate,A,nstate,0.0d0,pb,nstate)
     deallocate(A)
  else
     write(unit = *, fmt = *) "[E1] Error in enkf_hybrid_cov, wrong hybrid_cov_kind: ", hybrid_cov_kind
  end if
  
end subroutine enkf_hybrid_cov

