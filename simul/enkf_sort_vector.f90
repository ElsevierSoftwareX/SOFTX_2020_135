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

!> @brief Sorting a vector and saving the index manipulation
!> @param[inout] vector Vector to be sorted
!> @param[in] length Length of the vector
!> @param[inout] indices Updated index vector
!> @details
!> The double precision input vector is sorted and 
!> the original order is saved in the index vector
subroutine enkf_sort_vector_ind(vector,length,indices)

  implicit none

  integer, intent(in) :: length
  double precision, intent(inout), dimension(length) :: vector
  integer, intent(out), dimension(length) :: indices

  double precision :: r_dummy
  integer :: i_dummy
  integer :: i, j

  r_dummy = 0.0d0

  do i = 1, length
     indices(i) = i
  end do
  
  do i = 1, length-1
     do j = i+1, length

        if (vector(j) < vector(i)) then
           !Exchange values at places i,j
           r_dummy = vector(j)
           vector(j) = vector(i)
           vector(i) = r_dummy

           !Exchange indices i,j
           i_dummy = indices(j)
           indices(j) = indices(i)
           indices(i) = i_dummy
        end if

     end do
  end do

end subroutine enkf_sort_vector_ind

!> @param[inout] vector Vector to be sorted
!> @param[in] length Length of the vector
!> @brief Sort a vector
!> @details
!> Sort a double-precision vector.
subroutine enkf_sort_vector(vector,length)

  implicit none

  integer, intent(in) :: length
  double precision, intent(inout), dimension(length) :: vector

  double precision :: r_dummy
  integer :: i, j

  r_dummy = 0.0d0

  do i = 1, length-1
     do j = i+1, length

        if (vector(j) < vector(i)) then
           !Exchange values at places i,j
           r_dummy = vector(j)
           vector(j) = vector(i)
           vector(i) = r_dummy
        end if

     end do
  end do

end subroutine enkf_sort_vector
