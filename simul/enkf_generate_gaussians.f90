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

!> @brief Generate sample vectr of Gaussian distribution
!> @param[in] choice Switch for cumulative distribution or random generation('cumu'/'rand')
!> @param[inout] vector Vector that will contain the Gaussian distribution
!> @param[in] length Length of the vector
!> @details
!> A Gaussian distriubution is generated, either according to a
!> cumulative distribution (choice:`cumu`) oder randomly, and then
!> sorted (choice:`rand`).
subroutine enkf_generate_gaussians(choice,vector,length)

  use m_random

  implicit none

  character (len=4), intent (in) :: choice
  double precision, intent(inout), dimension(length) :: vector
  integer, intent(in) :: length
  
  double precision :: total_weight
  double precision :: weight
  double precision :: cum_prob
  double precision :: cum_prob_last
  double precision :: gaussian_value

  integer :: i, ierr

  double precision, allocatable, dimension(:) :: gausswork
  
  select case (choice)
  case ('cumu')

     total_weight = dble(length)
     cum_prob = 0.0d0
     cum_prob_last = 0.0d0
     weight = 0.0d0

     do i = 1, length

        cum_prob = cum_prob + 1.0d0/total_weight
        weight = (cum_prob + cum_prob_last)*0.5d0

        call gauinv(weight,gaussian_value,ierr)

        if(ierr == 1) then
           write(unit = *, fmt = *) "[E1] Error in enkf_generate_gaussians!"
        end if

        cum_prob_last = cum_prob

        vector(i) = gaussian_value

     end do
  case ('rand')

     allocate(gausswork(length))

     call random(gausswork,length)

     call enkf_sort_vector(gausswork,length)

     do i = 1, length
        vector(i) = gausswork(i)
     end do
     
     deallocate(gausswork)

  case default

     write(unit = *, fmt = *) "[E2] Error in enkf_generate_gaussians."
  end select


end subroutine enkf_generate_gaussians
