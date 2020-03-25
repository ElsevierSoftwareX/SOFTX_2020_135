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

!> @brief Specify state vector size
!> @details
!> At first the number of active variables and
!> parameters in the state vector is determined.
!> A variable/parameter is not the same as a (EnKF) state variable,
!> which can the value of both a variable/ parameter at a single grid
!> cell.
!>
!> Then `sysindx` is determined using the EnKF-Input `ex_unit`.  All
!> units less than or equal ex_unit will be excluded from the state
!> vector.  According to `sysindx` the number of active grid cells
!> (`lstate0`) and the number of active state variables (`nstate`) are
!> determined and saved.
subroutine enkf_state_len()

  use arrays, only: &
       uindex

  use mod_genrl, only: &
       i0,&
       j0,&
       k0

  use mod_enkf, only: &
       ex_unit,&
       lstate0,&
       act_s,&
       num_enkf_vars,&
       n_act_s_param,&
       n_act_s_state,&
       nstate_param,&
       nstate_state,&
       enkf_log_out,&
       stoch_bc_switch,&
       num_stoch_bc,&
       pres_vel_switch,&
       num_pres_vel,&
       tcon_switch,&
       num_tcon,&
       sysindx,&
       n_act_s,&
       lstate,&
       nstate

  implicit none

  integer :: i, j, k, l, l0

  ! Define lstate0
  lstate0 = i0*j0*k0

  ! Allocate sysindx
  allocate(sysindx(lstate0))

  ! Number of active variables/parameters in state vector
  n_act_s = sum(act_s)
  n_act_s_param = act_s(4) + act_s(5) + act_s(6)
  n_act_s_state = act_s(1)  + act_s(2)  + act_s(3)
  !Check
  if( (n_act_s < 0) .or. (n_act_s > num_enkf_vars) ) then
     write(unit = *, fmt = *) '[E1] Error in enkf_state_len.f90'
     stop
  end if

  if( (.not. stoch_bc_switch) .and. (.not. tcon_switch) .and. (.not. pres_vel_switch)) then
     if (n_act_s .le. 0) then
        write(unit = *, fmt = *) '[E2] Error in enkf_state_len.f90'
        stop
     end if
  end if
  
  !----------------------------------------------------------

  ! Determine sysindx
  if(enkf_log_out) write(37,*)'sysindx'
  ! Some units excluded
  IF (ex_unit.ne.0) then
     l0=0
     DO k=1,k0
        DO j=1,j0
           DO i=1,i0
              l = (k-1)*i0*j0 + (j-1)*i0 + i
              sysindx(l)=0
              IF(uindex(i,j,k).gt.ex_unit) then
                 l0=l0+1
                 sysindx(l)=l0
              END IF
              if(enkf_log_out) write(37,*) i,j,k,l,sysindx(l)
           END DO
        END DO
     END DO
     lstate = l0
  ! No units excluded
  ELSE
     DO k=1,k0
        DO j=1,j0
           DO i=1,i0
              l = (k-1)*i0*j0 + (j-1)*i0 + i
              sysindx(l)=l
              if(enkf_log_out) write(37,*)i,j,k,l,sysindx(l)
           END DO
        END DO
     END DO
     lstate = lstate0
  END IF

  ! Length of the state vector
  nstate = lstate*n_act_s
  ! Dual: Length of the parameters in the state vector
  nstate_param = lstate*n_act_s_param
  ! Dual: Length of the dynamic variables inthe state vector
  nstate_state = lstate*n_act_s_state
  
  
  ! Boundary conditions: Add num_stoch_bc entries to the state vector
  if (stoch_bc_switch) then
     nstate = nstate + num_stoch_bc
  end if

  ! Thermal conductivity
  if (tcon_switch) then
     nstate = nstate + num_tcon
  end if
  
  ! Prescribed velocity: Add two/three entries to state vector
  if (pres_vel_switch) then
     nstate = nstate + num_pres_vel
  end if
  
end subroutine enkf_state_len
