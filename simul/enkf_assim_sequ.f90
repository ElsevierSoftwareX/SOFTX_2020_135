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

!> @brief Invoke assimilation step sequentially for each variable
!> @details
!> The assimilation wrapper for update according to observations
!> of a single variable is called for each observed
!> variable. Additionally: Calling functions for debugging output.
subroutine enkf_assim_sequ()

  use mod_enkf, only: &
       act_o,&
       irobs,&
       num_enkf_vars

  implicit none
  
  integer :: ivar

  integer :: ifile2
  character (len=3) :: afile2

  call enkf_output_assimstp(0)
  
  call enkf_log(3,0)

  !-------------------------------------------------------------------------------
  ! Sequentially assimilate HEAD(ivar = 1), TEMP(2), CONC(3), KZ(4), LZ(5), POR(6)
  !-------------------------------------------------------------------------------
  do ivar = 1, num_enkf_vars
     if (act_o(irobs,ivar)==1) then
        call enkf_assim_single(ivar)
     end if
  end do

  call enkf_output_assimstp(1)
  
end subroutine enkf_assim_sequ
