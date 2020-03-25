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

!> @brief Perturb the prescribed velocity
!> @details
!> Add Gaussian noise to the prescribed velocity vdefault and call
!> the subroutine `enkf_velocity_dbc()` if head is active.
subroutine enkf_pres_vel()

  use m_random
  
  use arrays, only:&
       vdefault
       
  use mod_genrl, only:&
       nsmpl, head_active
  
  use mod_enkf, only:&
       num_pres_vel,&
       vdefault_stddevs
  
  implicit none

  integer :: ismpl, i

  double precision, allocatable, dimension(:) :: normedgauss

  allocate(normedgauss(nsmpl*num_pres_vel))

  call random(normedgauss,nsmpl*num_pres_vel)
  
  do ismpl = 1, nsmpl
     do i = 1, num_pres_vel
        vdefault(i,ismpl) = vdefault(i,ismpl) + vdefault_stddevs(i)*normedgauss((i-1)*nsmpl+ismpl)
     end do
  end do

  deallocate(normedgauss)

  if(head_active) then
     do ismpl = 1, nsmpl
        call enkf_velocity_dbc(ismpl)
     end do
  end if
  
end subroutine enkf_pres_vel
