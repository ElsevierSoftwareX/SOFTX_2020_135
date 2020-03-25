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

!> @brief Generate Normal Score arrays
!> @details
!> Set `mem_ns_original`, `mem_ns_gauss`, `mem` (Gets the Gaussian
!> zero mean, unit variance, but not sorted). More information in
!> mod_enkf.f90.
subroutine enkf_normal_score()

  use m_random
  
  use mod_enkf, only: &
       nstate,&
       nrens,&
       mem,&
       mem_ns_original,&
       mem_ns_gauss

  implicit none

  integer, dimension(nrens) :: mem_ind
  integer :: istate, i

  allocate(mem_ns_original(nrens,nstate))
  allocate(mem_ns_gauss(nrens,nstate))
  mem_ns_original(:,:) = 0.0d0
  mem_ns_gauss(:,:) = 0.0d0

  do istate = 1, nstate
     
     mem_ns_original(:,istate) = mem(istate,:)

     ! Sort mem_ns_original and track changes in mem_ind
     call enkf_sort_vector_ind(mem_ns_original(:,istate),nrens,mem_ind)

     ! Ordered gaussian ensemble, zero mean, variance one.
     call enkf_generate_gaussians('cumu',mem(istate,:),nrens)
     mem_ns_gauss(:,istate) = mem(istate,:)

     ! Original order for mem
     call enkf_sort_vector_back(mem(istate,:),nrens,mem_ind)

  end do

end subroutine enkf_normal_score
