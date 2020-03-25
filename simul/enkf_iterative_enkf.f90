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

!> @brief Preparations for Iterative EnKF update
!> @details
!> First update: `head`, `temp`, `conc` and the simulation time are
!> saved. The number of starts is set to one. Every update:
!> Calculate the number of update before the next restart. If this
!> number is reached, reset head, temp, conc and the simulation
!> time to the start values.
subroutine enkf_iterative_enkf()

  use arrays, only: &
       head,&
       temp,&
       conc,&
       simtime
       
  use mod_genrl, only:&
       i0,&
       j0,&
       k0,&
       nsmpl
  
  use mod_conc, only:&
       ntrans

  use mod_enkf, only: &
       irobs,&
       nrens,&
       head_tmp,&
       temp_tmp,&
       conc_tmp,&
       simtime_tmp,&
       nrobs_int,&
       iterative_nrobs_int,&
       iterative_irobs

  implicit none

  integer :: i, j, k , iens, iconc, ismpl
  integer :: irobs_next_restart
  
  ! Before first update: Save variable values and simtime
  if (irobs == 1) then
     allocate(head_tmp(i0,j0,k0,nrens))
     allocate(temp_tmp(i0,j0,k0,nrens))
     allocate(conc_tmp(i0,j0,k0,max(ntrans,1),nrens))
     allocate(simtime_tmp(nsmpl))

     do iens = 1, nrens
        do k = 1, k0
           do j = 1, j0
              do i = 1, i0
                 head_tmp(i,j,k,iens) = head(i,j,k,iens)
                 temp_tmp(i,j,k,iens) = temp(i,j,k,iens)
                 do iconc = 1, max(ntrans,1)
                    conc_tmp(i,j,k,iconc,iens) = conc(i,j,k,iconc,iens)
                 end do
              end do
           end do
        end do
     end do
     do ismpl = 1, nsmpl
        simtime_tmp(ismpl) = simtime(ismpl)
     end do

     iterative_irobs = 1
  end if

  ! Number of updates before next restart
  irobs_next_restart = iterative_nrobs_int*(iterative_irobs-1)*iterative_irobs/2

  if(irobs == irobs_next_restart+1)then
     do iens = 1, nrens
        do k = 1, k0
           do j = 1, j0
              do i = 1, i0
                 head(i,j,k,iens) = head_tmp(i,j,k,iens)
                 temp(i,j,k,iens) = temp_tmp(i,j,k,iens)
                 do iconc = 1, max(ntrans,1)
                    conc(i,j,k,iconc,iens) = conc_tmp(i,j,k,iconc,iens)
                 end do
              end do
           end do
        end do
     end do
     do ismpl = 1, nsmpl
        simtime(ismpl) = simtime_tmp(ismpl)
     end do

     iterative_irobs = iterative_irobs + 1
  end if

  ! Last step: Deallocate temp-arrays
  if (irobs == nrobs_int) then
     if(allocated(head_tmp)) deallocate(head_tmp)
     if(allocated(temp_tmp)) deallocate(temp_tmp)
     if(allocated(conc_tmp)) deallocate(conc_tmp)
     if(allocated(simtime_tmp)) deallocate(simtime_tmp)
  end if

end subroutine enkf_iterative_enkf
