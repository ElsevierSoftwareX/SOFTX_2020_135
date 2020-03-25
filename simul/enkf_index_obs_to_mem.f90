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

!> @brief Calculate mem index from observation location index and variable index
!> @param[in] iobs Index of observation location between 1 and nrobs_loc(irobs)
!> @param[in] ivar Variable index
!> @details
!> Calculate `mem` index from observation location index and variable
!> index.
integer function index_obs_to_mem(iobs,ivar)

  use mod_enkf, only:&
       irobs,&
       i_obs,&
       j_obs,&
       k_obs,&
       lstate,&
       nstate,&
       rank_mem

  implicit none


  integer, intent(in) :: iobs
  integer, intent(in) :: ivar
  
  integer :: imem

  integer :: l

  integer, external :: index_loc_to_lin

  !Linear index inside the grid
  l = index_loc_to_lin(i_obs(irobs,iobs),j_obs(irobs,iobs),k_obs(irobs,iobs))

  if(l< 1) then

     imem = 0

  else if(rank_mem(ivar) < 0) then

     imem = -1

  else

     !Index inside the state vector
     imem = rank_mem(ivar)*lstate + l

     if(imem < 1 .or. imem > nstate) then
        write(unit = *, fmt = *) '[E2] Error in enkf_index_obs_to_mem:'
        stop 1
     end if

  end if

  index_obs_to_mem = imem
  
end function index_obs_to_mem
