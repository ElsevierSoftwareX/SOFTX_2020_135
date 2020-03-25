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

!> @brief Deallocate Normal-Score EnKF arrays
!> @details
!> The normal score arrays allocated in the subroutines
!> `enkf_normal_score`, `enkf_normal_score_obsvars` and
!> `enkf_normal_score_observations_obsvars` are deallocated.
subroutine enkf_normal_score_dealloc()
  
  use mod_enkf, only: &
       mem_ns_original,&
       mem_ns_gauss,&
       obs_ns,&
       obs_ns_original,&
       obs_ns_gauss,&
       obs_ns_orig_vars,&
       obsvar_ns

  implicit none

  ! From enkf_normal_score
  if(allocated(mem_ns_original)) deallocate(mem_ns_original)
  if(allocated(mem_ns_gauss)) deallocate(mem_ns_gauss)
  
  ! From enkf_normal_score_obsvars
  if(allocated(obs_ns)) deallocate(obs_ns)
  if(allocated(obs_ns_original)) deallocate(obs_ns_original)
  if(allocated(obs_ns_gauss)) deallocate(obs_ns_gauss)
  if(allocated(obs_ns_orig_vars)) deallocate(obs_ns_orig_vars)

  ! From enkf_normal_score_observations_obsvars
  if(allocated(obsvar_ns)) deallocate(obsvar_ns)
  
end subroutine enkf_normal_score_dealloc
