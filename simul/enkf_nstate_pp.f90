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

!> @brief Set nstate for Pilot-Point EnKF
!> @param[in] befaft 'bef'/'aft', switch before or after assimilation
!> @details
!> Before the update, `nstate` is reduced, the number of state vector
!> entries of the Pilot-Point variable are reduced from `lstate` to
!> `num_pp`.
!>
!> After the update, the old `nstate` is restored.
subroutine enkf_nstate_pp(befaft)

  use mod_enkf, only: &
       nstate, &
       nstate_pp_temp,&
       lstate,&
       num_pp

  implicit none

  character (len=3), intent (in) :: befaft

  if (befaft == 'bef') then

     nstate_pp_temp = nstate
     nstate = nstate_pp_temp-lstate+num_pp

  else if (befaft == 'aft') then

     nstate = nstate_pp_temp

  else

     write(unit = *, fmt = *) "[E1] Error in enkf_nstate_pp"
     stop 1

  end if

end subroutine enkf_nstate_pp
