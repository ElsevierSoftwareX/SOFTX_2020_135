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

!> @brief Variance of vector
!> @param[in] vec Vector for variance calculation
!> @param[in] len Length of vector
!> @param[in] ave Mean of vector
!> @param[out] var Variance
subroutine enkf_vector_variance(vec,len,ave,var)

  implicit none

  integer, intent(in) :: len
  double precision, intent(in) :: ave
  double precision, intent(in) :: vec(len)
  double precision, intent(out) :: var

  integer :: i
  
  var = 0.0d0
  do i = 1, len
     var = var + (vec(i)-ave)*(vec(i)-ave)
  end do
  var = (1.0d0/real(len-1,8))*var
     
end subroutine enkf_vector_variance
