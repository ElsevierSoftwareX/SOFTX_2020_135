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

!> @brief Compute correlation matrix with respect to one location
!> @details
!> This subroutine takes the location and variable information from
!> the inputfile and calculates the correlation of all state-vector
!> variables with respect to the specified state-vector variable in
!> corr_matrix.
!>
!> The reference point variable is taken from the state vector mem
!> if ivar_cov_ref in {1,2,3,4,5,6}, or it is calculated from the
!> underlying arrays if ivar_cov_ref in {11,12,13,14,15,16}.
subroutine enkf_compute_covs()
  
  use mod_enkf, only:&
       mem,&
       ave,&
       var,&
       nstate,&
       nrens,&
       corr_matrix,&
       mat_cov_ref,&
       num_cov_ref

  implicit none
  
  integer :: istate, irens, iref

  integer :: istate_cov_ref

  double precision :: average_temp
  double precision :: average_ref
  double precision :: stddev_temp
  double precision :: stddev_ref

  double precision :: norm_const

  double precision, dimension(nrens) :: mem_temp
  double precision, dimension(nrens) :: mem_ref

  integer :: i_cov_ref
  integer :: j_cov_ref
  integer :: k_cov_ref
  integer :: ivar_cov_ref

  double precision, external :: varpar

  integer, external :: index_loc_to_mem

  norm_const = real(nrens,8) - 1.0d0
  
  ! Loop over number of covariance reference locations
  do iref = 1, num_cov_ref

     ! Get the location and variable information
     i_cov_ref = mat_cov_ref(iref,1)
     j_cov_ref = mat_cov_ref(iref,2)
     k_cov_ref = mat_cov_ref(iref,3)
     ivar_cov_ref = mat_cov_ref(iref,4)

     ! Reference from mem
     if(ivar_cov_ref < 10) then

        istate_cov_ref = index_loc_to_mem(i_cov_ref,j_cov_ref,k_cov_ref,&
             ivar_cov_ref)
        
        do irens = 1, nrens
           mem_ref(irens) = mem(istate_cov_ref,irens)
        end do
        average_ref = ave(istate_cov_ref)
        stddev_ref = dsqrt(var(istate_cov_ref))
        
     ! Reference from varpar
     else if( ivar_cov_ref > 10) then

        ivar_cov_ref = ivar_cov_ref - 10

        do irens = 1,nrens
           mem_ref(irens) = varpar(i_cov_ref,j_cov_ref,k_cov_ref,irens,ivar_cov_ref)
        end do
        call enkf_vector_average(mem_ref,nrens,average_ref)
        call enkf_vector_variance(mem_ref,nrens,average_ref,stddev_ref)
        stddev_ref = dsqrt(stddev_ref)

     ! Reference Error
     else
        write(unit = *, fmt = *) '[E] Error in enkf_compute covs'
        stop
     end if

     do istate = 1, nstate

        ! State vector
        average_temp = ave(istate)
        stddev_temp = dsqrt(var(istate))

        do irens = 1, nrens
           mem_temp(irens) = mem(istate,irens)
        end do

        !Covariance matrix
        corr_matrix(istate,iref) = 0.0d0
        do irens = 1, nrens
           corr_matrix(istate,iref) = corr_matrix(istate,iref) &
                + (mem_temp(irens)-average_temp)&
                *(mem_ref(irens)-average_ref)
        end do
        corr_matrix(istate,iref) = corr_matrix(istate,iref)/norm_const

        !Correlation matrix
        corr_matrix(istate,iref) = corr_matrix(istate,iref)&
             /(stddev_ref*stddev_temp)

     end do

  end do

end subroutine enkf_compute_covs
