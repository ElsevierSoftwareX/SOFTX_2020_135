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

!> @brief Preparations for Dual EnKF update
!> @details
!> For an odd update, set the activity of the states to zero,
!> introduce the new number of active state variables and activity
!> vector and save the dynamic state variable values and save the
!> simulation time. For an even update, return to the save values
!> of the dynamic variables and the simulation time and set the
!> activity of the paramters to zero, introduce the new number of
!> active state variables and the activity vector.
subroutine enkf_dual_enkf()

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
       n_act_s_tmp,&
       nrens,&
       n_act_s,&
       n_act_s_param,&
       n_act_s_state,&
       nstate,&
       nstate_param,&
       nstate_state,&
       head_tmp,&
       temp_tmp,&
       conc_tmp,&
       simtime_tmp,&
       nrobs_int,&
       act_s_tmp,&
       act_s

  implicit none

  integer :: remainder_dual_enkf
  integer :: i, j, k, iens, ismpl, iconc, ivar

  ! Allocate temporary arrays
  if (irobs == 1) then
     allocate(head_tmp(i0,j0,k0,nrens))
     allocate(temp_tmp(i0,j0,k0,nrens))
     allocate(conc_tmp(i0,j0,k0,max(ntrans,1),nrens))
     allocate(simtime_tmp(nsmpl))
     n_act_s_tmp = n_act_s
     act_s_tmp = act_s
  end if

  ! Remainder: Parameter or state update
  remainder_dual_enkf = mod(irobs,2)

  select case (remainder_dual_enkf)
  case (1)                      !Parameter update
     n_act_s = n_act_s_param
     nstate = nstate_param
     act_s(1) = 0
     act_s(2) = 0
     act_s(3) = 0
     act_s(4) = act_s_tmp(4)
     act_s(5) = act_s_tmp(5)
     act_s(6) = act_s_tmp(6)
     ! Save values of state variables
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
     ! Save simulation time
     do ismpl = 1, nsmpl
        simtime_tmp(ismpl) = simtime(ismpl)
     end do
  case(0)                       !State update
     n_act_s = n_act_s_state
     nstate = nstate_state
     act_s(1) = act_s_tmp(1)
     act_s(2) = act_s_tmp(2)
     act_s(3) = act_s_tmp(3)
     act_s(4) = 0
     act_s(5) = 0
     act_s(6) = 0
     ! Return to saved values of state variables
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
     ! Return to save simulation time
     do ismpl = 1, nsmpl
        simtime(ismpl) = simtime_tmp(ismpl)
     end do
  end select

  ! Deallocate temporary arrays
  if (irobs == nrobs_int) then
     if(allocated(head_tmp)) deallocate(head_tmp)
     if(allocated(temp_tmp)) deallocate(temp_tmp)
     if(allocated(conc_tmp)) deallocate(conc_tmp)
     if(allocated(simtime_tmp)) deallocate(simtime_tmp)
  end if

end subroutine enkf_dual_enkf
