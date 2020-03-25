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

!> @brief EnKF assimilation step output wrapper
!> @param[in] i_assimstp index of different assimilation output calls
!> @details
!> Assimilation step output wrapper.
!> - 0: open assimstp file
!> - 1: close assimstp file 
subroutine enkf_output_assimstp(i_assimstp)

  use arrays, only:&
       project_sfx

  use mod_simul, only:&
       senkf_outdir

  use mod_enkf, only:&
       irobs,&
       assimstp_switch

  implicit none

  integer, intent(in) :: i_assimstp

  integer :: ifile2
  character (len=3) :: afile2

  !-------------------------- 0 -----------------------
  if (i_assimstp == 0 .and. assimstp_switch) then

     ifile2 = irobs + 100
     WRITE(afile2,'(I3)') ifile2
     afile2 = trim(afile2)
     OPEN(unit=12,file=senkf_outdir//'assimstp'//trim(project_sfx(1))//'_'//afile2)
     
  !-------------------------- 1 -----------------------
  else if (i_assimstp == 1 .and. assimstp_switch) then

     close(12)
        
  else if (assimstp_switch) then
     write(unit = *, fmt = *) '[E] Error in enkf_output_assimstp.f90'
     stop 1
  end if

end subroutine enkf_output_assimstp
