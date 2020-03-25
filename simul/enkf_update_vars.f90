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

!> @brief Set update variables after ENKF
!> @details
!> - Outer loop over variables (ivar)
!> - Loop over the samples (irens) and the location (i,j,k).
!> - Calculate index in state vector mem (imem from
!> index_loc_to_mem).
!> - Call enkf_set_varpar
!> If they are included in the EnKF, update boundary conditions,
!> constant thermal conductivity or prescribed velocity.
subroutine enkf_update_vars()

  use arrays, only: &
       propunit,&
       nbc_data,&
       dbc_data,&
       vdefault,&
       idx_lz
  
  use mod_genrl, only: &
       i0,&
       j0,&
       k0,&
       head_active

  use mod_enkf, only: &
       nrens,&
       lstate,&
       sysindx,&
       act_s,&
       num_enkf_vars,&
       mem,&
       nstate,&
       num_stoch_bc,&
       stoch_bc_switch,&
       stoch_bc_ibc,&
       pres_vel_switch,&
       num_pres_vel,&
       tcon_switch,&
       num_tcon,&
       tcon_inds

  implicit none

  integer :: irens, i, j, k, ivar
  integer :: imem
  
  integer, external :: index_loc_to_mem

  do ivar = 1, num_enkf_vars
     if(act_s(ivar) == 1) then

        DO irens = 1, nrens
           DO k = 1, k0
              DO j = 1, j0
                 DO i = 1, i0

                    imem = index_loc_to_mem(i,j,k,ivar)

                    if (imem > 0) then
                       call enkf_set_varpar(i,j,k,irens,imem,ivar)
                    end if

                 END DO
              END DO
           END DO
        END DO

     end if
  end do

  ! Update boundary conditions
  if(stoch_bc_switch) then
     
     do irens = 1, nrens
        do i = 1, nbc_data

           ! Use index-array
           if(stoch_bc_ibc(i)>0) then
              dbc_data(i,1,irens)=mem(nstate-num_stoch_bc+stoch_bc_ibc(i),irens)
           end if
           
        end do
     end do
     
  end if


  if (tcon_switch) then
     do irens = 1, nrens
        do i = 1, num_tcon
           propunit(tcon_inds(i),idx_lz,irens) = mem(nstate-num_pres_vel-num_tcon+i,irens)
        end do
     end do
  end if
  
  ! Update prescribed flow
  if(pres_vel_switch) then
     
     do irens = 1, nrens
        do i = 1, num_pres_vel
           vdefault(i,irens) = mem(nstate-num_pres_vel+i,irens)
        end do
     end do

     if(head_active) then
        do irens = 1, nrens
           call enkf_velocity_dbc(irens)
        end do
     end if
     
  end if

end subroutine enkf_update_vars
