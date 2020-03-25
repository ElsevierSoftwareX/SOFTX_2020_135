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

!> @brief Output single updated thermal conductivities
!> @param[in] a_before_after switch before or after assimilation ('bef'/'aft')
!> @details
!> One line of output contains all the realizations of the
!> thermal conductivity in that order order: tcon, tcon_2, ...,
!> tcon_nrens
subroutine enkf_output_tcon(a_before_after)

  use mod_simul, only: &
       senkf_outdir

  use mod_enkf, only: &
       irobs,&
       num_pres_vel,&
       nrens,&
       nstate,&
       mem,&
       num_tcon
  
  implicit none

  character (len=3), intent(in) :: a_before_after
  character (len=5) :: c_num_tcon

  integer :: imem, irens

  ! Open/Reopen output file
  if (irobs == 1) then
     open(unit = 20, file = senkf_outdir//'tcon'//&
          '_'//trim(a_before_after)//'.txt')
  else
     open(unit = 20, file = senkf_outdir//'tcon'//&
          '_'//trim(a_before_after)//'.txt', position="append")
  end if

  ! First line
  if (irobs == 1) then
     write(unit = 20, fmt = '(2a,i5)') &
          'Thermal conductivities '//a_before_after//' assimilations (num_tcon,nrens):', &
          num_tcon, nrens
  end if

  ! String with number of realizations times number of tcons
  write(c_num_tcon, fmt = '(i5)') num_tcon*nrens
  
  ! Output thermal conductivities
  write(unit = 20, fmt = '('//trim(c_num_tcon)//'e16.8)') &
       ((mem(nstate-num_pres_vel-num_tcon+imem,irens),irens=1,nrens),imem=1,num_tcon)

  ! Close output file
  close(20)

end subroutine enkf_output_tcon
