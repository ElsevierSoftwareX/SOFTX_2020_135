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

!> @brief __NOT RELIABLE__: Perturb boundary conditions.
!> @details
!> Head boundary conditions can be perturbed. This needs to be
!> developed further.
subroutine enkf_stoch_bc()

  use m_random

  use arrays, only: &
       propunit,&
       nbc_data,&
       ibc_data,&
       dbc_data,&
       cbc_i,&
       cbc_j,&
       cbc_pv,&
       pv_head

  use mod_genrl, only: &
       nsmpl

  use mod_enkf, only: &
       stoch_bc_stddevs,&
       stoch_bc_nbc,&
       stoch_bc_ibc,&
       num_stoch_bc

  implicit none

  integer :: i, ismpl

  double precision, allocatable, dimension(:) :: perturb_1, perturb_2

  allocate(stoch_bc_nbc(num_stoch_bc))
  allocate(stoch_bc_ibc(nbc_data))
  stoch_bc_nbc(:) = 0
  stoch_bc_ibc(:) = 0

  allocate(perturb_1(nsmpl))
  allocate(perturb_2(nsmpl))

  call random(perturb_1,nsmpl)
  call random(perturb_2,nsmpl)

  !First unit, head bc
  do i = 1, nsmpl
     propunit(1,18,i)=propunit(1,18,i) + stoch_bc_stddevs(1)*perturb_1(i)
  end do
  do ismpl = 1, nsmpl
     do i = 1, nbc_data
        ! i=1
        if(ibc_data(i,cbc_i) == 1)then
           ! Head
           if(ibc_data(i,cbc_pv) == pv_head) then
              dbc_data(i,1,ismpl)=dbc_data(i,1,ismpl) + stoch_bc_stddevs(1)*perturb_1(ismpl)
              stoch_bc_nbc(1) = i
              stoch_bc_ibc(i) = 1
           end if
        ! j=1
        else if(ibc_data(i,cbc_j)== 1) then
           ! Head
           if(ibc_data(i,cbc_pv) == pv_head) then
              dbc_data(i,1,ismpl)=dbc_data(i,1,ismpl) + stoch_bc_stddevs(1)*perturb_2(ismpl)
              stoch_bc_nbc(2) = i
              stoch_bc_ibc(i) = 2
           end if
        end if
     end do
  end do

  deallocate(perturb_1)
  deallocate(perturb_2)

end subroutine enkf_stoch_bc

