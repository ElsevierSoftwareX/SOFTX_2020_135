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

!> @brief Deallocate EnKF update arrays
!> @details
!> These arrays are allocated and deallocated at every update.
subroutine enkf_dealloc_state_vector()

  use mod_enkf, only: &
       mem,&
       ave,&
       var,&
       sysvarmem,&
       corr_matrix,&
       rhos,&
       pp_state_inds,&
       pp_lstate_inds

  implicit none

  DEALLOCATE(mem)
  DEALLOCATE(ave)
  DEALLOCATE(var)
  deallocate(sysvarmem)
  if (allocated(corr_matrix)) deallocate(corr_matrix)
  if (allocated(rhos)) deallocate(rhos)
  if (allocated(pp_state_inds)) deallocate(pp_state_inds)
  if (allocated(pp_lstate_inds)) deallocate(pp_lstate_inds)

end subroutine enkf_dealloc_state_vector
