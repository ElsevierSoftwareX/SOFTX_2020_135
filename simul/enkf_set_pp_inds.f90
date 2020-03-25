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

!> @brief Set state vector indices of pilot point locations
!> @details
!> Fill the array `pp_state_inds` with the state vector indices of the
!> Pilot Point locations . Then, fill the array `pp_lstate_inds` with
!> the linear (one-variable) index of the pilot points.
!>
!> Find the state vector indices of the Pilot-Point variable. Save the
!> first state vector index of the Pilot-Point variable in
!> `ipp_start`. Save the last index of the Pilot Point variable
!> `ipp_end`.
subroutine enkf_set_pp_inds ()

  use mod_genrl, only: &
       i0,&
       j0,&
       k0
  
  use mod_enkf, only: &
       pp_state_inds,&
       pp_lstate_inds,&
       lstate,&
       num_pp,&
       pp_inds,&
       sysindx,&
       pp_ivar,&
       ipp_start,&
       ipp_end

  implicit none
  
  integer :: ipp, i_pp, j_pp, k_pp, irens

  integer, external :: index_loc_to_mem

  ! State vector indices of Pilot-Point locations
  allocate(pp_state_inds(num_pp))
  allocate(pp_lstate_inds(num_pp))
  do ipp = 1, num_pp
     i_pp = pp_inds(ipp,1)
     j_pp = pp_inds(ipp,2)
     k_pp = pp_inds(ipp,3)

     pp_lstate_inds(ipp) = sysindx(i_pp + (j_pp-1)*i0 + (k_pp-1)*i0*j0)
     pp_state_inds(ipp) = index_loc_to_mem(i_pp,j_pp,k_pp,pp_ivar)
  end do

  ! Save first and last index of full Pilot-Point variable in mem
  ipp_start = index_loc_to_mem(1,1,1,pp_ivar)
  ipp_end = index_loc_to_mem(i0,j0,k0,pp_ivar)
  if(.not. ipp_end-ipp_start+1 == lstate) then
     write (unit = *, fmt = *) '[E1] Error with indices in enkf_set_pp_inds.'
  end if

end subroutine enkf_set_pp_inds
