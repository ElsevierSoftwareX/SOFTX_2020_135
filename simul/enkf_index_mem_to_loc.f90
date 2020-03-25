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

!> @brief Calculate location and variable indices from mem index
!> @param[in] imem Index of the state vector
!> @param[out] i grid index in x direction
!> @param[out] j grid index in y direction
!> @param[out] k grid index in z direction
!> @param[out] ivar Index of the variable between 1 and num_enkf_vars
!> @details
!> `sysindx` has to be trivial (i.e., no excluded units), otherwise
!> this simple approach does not work
subroutine enkf_index_mem_to_loc(imem,i,j,k,ivar)

  use mod_genrl, only:&
       i0,&
       j0,&
       k0

  use mod_enkf, only:&
       lstate,&
       lstate0,&
       rank_mem,&
       num_enkf_vars

  implicit none

  integer, intent(in) :: imem

  integer, intent(out) :: i
  integer, intent(out) :: j
  integer, intent(out) :: k
  integer, intent(out) :: ivar

  integer :: l, rmem, jvar

  if(lstate == lstate0) then
     
     ! Find ivar
     rmem = (imem-1)/lstate
     if(rmem>=0 .and. rmem<= num_enkf_vars) then
        do jvar = 1, num_enkf_vars
           if(rank_mem(jvar) == rmem) then
              ivar = jvar
           end if
        end do
     else
        write(unit = *, fmt = *) "[E5] Error in enkf_index_mem_to_loc"
        stop
     end if
     
     ! Find i, j, k
     l = mod(imem-1,lstate)+1
     if(l<1 .or. l>lstate) then
        write(unit = *, fmt = *) "[E2] Error in enkf_index_mem_to_loc"
        stop
     end if

     i = mod(l-1,i0)+1
     l = l-i
     if(mod(l,i0)==0) then
        l = l/i0 + 1
     else
        write(unit = *, fmt = *) "[E3] Error in enkf_index_mem_to_loc"
        stop
     end if

     j = mod(l-1,j0)+1
     l = l-j
     if(0<= l .and. l<=k0 .and. mod(l,j0)==0) then
        k = l/j0 + 1
     else
        write(unit = *, fmt = *) "[E4] Error in enkf_index_mem_to_loc"
        stop
     end if

  else
     ! non-trivial sysindx does not allow easy inverting
     write(unit = *, fmt = *) "[E1] Error in enkf_index_mem_to_loc"
     stop
  end if



end subroutine enkf_index_mem_to_loc
