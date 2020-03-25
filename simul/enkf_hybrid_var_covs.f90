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

!> @brief Set standard hybrid variance 
!> @param[in] ihyb Hybrid index running over all possible state vector variables
!> @param[out] hybrid_cov Standard hybrid variance
!> @details
!> Set standard hybrid variance.
subroutine enkf_hybrid_var_covs(ihyb,hybrid_cov)

  use mod_enkf, only:&
       lstate

  implicit none

  integer,intent(in) :: ihyb

  double precision, intent(out) :: hybrid_cov

  integer :: ivar

  ! Get the variable index ivar corresponding to the hybrid index ihyb
  ivar = ((ihyb-1)/lstate) + 1
  
  select case (ivar)
  case(1)                       !head
     hybrid_cov = 0.5
  case(2)                       !temp
     hybrid_cov = 0.5
  case(3)                       !conc
     hybrid_cov = 1.0d-9
  case(4)                       !kz
     hybrid_cov = 0.5
  case(5)                       !lz
     hybrid_cov = 0.5
  case(6)                       !por
     hybrid_cov = 0.05
  case default
     write(unit = *, fmt = *) '[E1] Error in enkf_hybrid_var_covs'
  end select
  
end subroutine enkf_hybrid_var_covs
