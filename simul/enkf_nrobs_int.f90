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

!> @brief Change `nrobs_int` for Dual EnKF and Iterative EnKF
subroutine enkf_nrobs_int()

  use mod_enkf, only: &
       nrobs_int,&
       nrobs_int_pure,&
       dual_enkf_switch,&
       iterative_switch,&
       iterative_nrobs_int

  implicit none

  ! Save input value
  nrobs_int_pure = nrobs_int

  ! Dual EnKF (2*n)
  if(dual_enkf_switch) then
     nrobs_int = 2*nrobs_int_pure
  end if

  ! Iterative EnKF (1/2*n'*(n'+1) + n)
  if(iterative_switch) then
     nrobs_int = ( iterative_nrobs_int*((nrobs_int_pure-1)/iterative_nrobs_int)&
          *((nrobs_int_pure-1)/iterative_nrobs_int + 1) )/2  +  nrobs_int_pure
  end if
  
end subroutine enkf_nrobs_int
