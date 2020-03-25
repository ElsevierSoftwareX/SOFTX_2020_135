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

!> @brief Save overall residual between mean and true model
!> @param[in] a_before_after bef/aft assimilation
!> @details
!> USE: Variables describing the length of the state vector, number of
!> active species, number of observation times and, most
!> importantly, the true values and average values for each state
!> vector variable.
!>
!> SET: The overall residual between the true values and the ensemble
!> averages of the state vector variables in resvar. This is
!> calculated as sum of squares of the differences - no
!> normalization is calculated yet!
subroutine enkf_resid(a_before_after)

  use mod_genrl, only: &
       i0,&
       j0,&
       k0

  use mod_enkf, only: &
       var_true,&
       ave,&
       act_s,&
       irobs,&
       num_enkf_vars,&
       resvar

  implicit none

  character (len=3), intent(in) :: a_before_after

  integer :: i, j, k, imem, lmem, ivar
  integer :: ires

  integer, external :: index_loc_to_lin
  integer, external :: index_loc_to_mem

  ! Set ires
  if(a_before_after == 'bef') then
     ires = 2*irobs -1
  else if(a_before_after == 'aft') then
     ires = 2*irobs
  else
     write(unit = *, fmt = *) '[E1] Error in enkf_resid()'
     stop 1
  end if

  ! Set resvar
  do ivar = 1, num_enkf_vars
     if(act_s(ivar)==1) then

        do k = 1, k0
           do j = 1, j0
              do i = 1, i0
                 imem = index_loc_to_mem(i,j,k,ivar)
                 lmem = index_loc_to_lin(i,j,k)
                 if (imem > 0) then

                    resvar(ires,ivar) = resvar(ires,ivar) &
                         + (var_true(lmem,ivar)-ave(imem))**2

                 end if
              end do
           end do
        end do

     end if
  end do

end subroutine enkf_resid
