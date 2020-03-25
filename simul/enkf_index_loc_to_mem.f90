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

!> @brief Calculate mem index from location and variable indices
!> @param[in] i grid index in x direction
!> @param[in] j grid index in y direction
!> @param[in] k grid index in z direction
!> @param[in] ivar Index of the variable between 1 and num_enkf_vars
!> @details
!> The location and variable index are checked against rank_mem and
!> sysindx and then, the __linear state vector index__ is
!> calculated.
!>
!> __Warnings and errors__:
!> 1. Wrong `sysindx` (location not in state vector): `imem = 0`
!> 2. Wrong `rank_mem` (variable not in state vector): `imem = -1`
!> 3. Wrong `imem` (`imem` outside `1, nstate`): Error code.
integer function index_loc_to_mem(i,j,k,ivar)

  use mod_enkf, only:&
       lstate,&
       rank_mem,&
       nstate

  implicit none
  
  integer, intent(in) :: i
  integer, intent(in) :: j
  integer, intent(in) :: k
  integer, intent(in) :: ivar

  integer :: imem

  integer :: l

  integer, external :: index_loc_to_lin

  !Linear grid index
  l =  index_loc_to_lin(i,j,k)

  if(l < 1) then

     imem = 0

  else if(rank_mem(ivar) < 0) then

     imem = -1

  else

     !Index in state vector mem
     imem = rank_mem(ivar)*lstate + l

     if(imem < 1 .or. imem > nstate) then

        write(unit = *, fmt = *) '[E] Error in enkf_index_loc_to_mem'
        stop 1

     end if

  end if

  index_loc_to_mem = imem

end function index_loc_to_mem

!> @brief Calculate linear index from location indices
!> @param[in] i grid index in x direction
!> @param[in] j grid index in y direction
!> @param[in] k grid index in z direction
!> @details
!> A linear index of state-active locations is calculated from the
!> location indices and from sysindx.
!>
!> At this point, the variable/parameter does not play a role.
integer function index_loc_to_lin (i,j,k)

  use mod_genrl, only:&
       i0,&
       j0

  use mod_enkf, only:&
       sysindx

  implicit none

  integer, intent(in) :: i
  integer, intent(in) :: j
  integer, intent(in) :: k

  integer :: l

  ! Linear grid index
  l =  i + (j-1)*i0 + (k-1)*i0*j0

  ! Check if location in state vector
  index_loc_to_lin = sysindx(l)

end function index_loc_to_lin
