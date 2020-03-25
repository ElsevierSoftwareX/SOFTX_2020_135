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

!> @brief Set global random seed according to two seed-seeds
!> @param[in] seed_1 Seed-seed number one
!> @param[in] seed_2 Seed-seed number two
!> @details
!> The seed-seeds are used as the initial two numbers of a
!> Fibonacci series modulo 97. The length of the final random seed
!> is determined by the compiler subroutine random_seed; with the
!> same subroutine, the random seed is put.
!> Through this subroutine, the global random seed can be specified
!> which allows for reproducible runs.
subroutine enkf_set_random_seed(seed_1,seed_2)

  implicit none

  integer, intent(in) :: seed_1
  integer, intent(in) :: seed_2

  integer :: seed_size, i
  integer, allocatable, dimension(:) :: seed_put

  !Setting the seed for the random number production
  call random_seed(size = seed_size)
  allocate(seed_put(seed_size))
  
  do i = 1, seed_size
     if(i == 1) then
        !First two seed numbers were input
        seed_put(i) = seed_1
     else if (i == 2) then
        seed_put(i) = seed_2
     else
        !Fibonacci modulo 97 as seed-sequence
        seed_put(i) = mod(seed_put(i-2) + seed_put(i-1),97)
     end if
  end do
  call random_seed(put = seed_put)

  deallocate(seed_put)

end subroutine enkf_set_random_seed
