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

!> @brief Initial uncertainty for thermal conductivity
!> @details
!> When only one thermal conductivity is supposed to be estimated
!> for a whole unit, this routine adds a random component to the
!> value at the beginning of the simulation.
subroutine enkf_tcon()

  use m_random
  
  use arrays, only: &
       idx_lz,&
       propunit

  use mod_genrl, only: &
       nsmpl

  use mod_enkf, only: &
       tcon_inds,&
       tcon_stddev,&
       num_tcon
  
  implicit none

  integer :: ismpl, itc
  
  double precision, allocatable, dimension(:) :: perturb_1

  allocate(perturb_1(num_tcon*nsmpl))
  
  call random(perturb_1,num_tcon*nsmpl)

  do ismpl = 1, nsmpl
     do itc = 1, num_tcon
        propunit(tcon_inds(itc),idx_lz,ismpl) = &
             propunit(tcon_inds(itc),idx_lz,ismpl) + tcon_stddev * perturb_1((itc-1)*num_tcon+ismpl)
     end do
  end do
  
  deallocate(perturb_1)
  
end subroutine enkf_tcon
