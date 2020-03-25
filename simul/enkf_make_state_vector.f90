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

!> @brief Put state vectors/variances into mem/sysvarmem
!> @details
!> USE: Size of the state-vector/variance matrices, activity and
!> variances of the parameters. Import function: varpar, which
!> sets the variables/parameters.
!>
!> SET: rank_mem, mem, sysvarmem
subroutine enkf_make_state_vector()

  use arrays, only: &
       dbc_data,&
       vdefault,&
       propunit,&
       idx_lz
  
  use mod_genrl, only: &
       i0,&
       j0,&
       k0
 
  use mod_enkf, only: &
       nstate,&
       nrens,&
       lstate,&
       sysindx,&
       act_s,&
       num_enkf_vars,&
       sysvar,&
       stoch_bc_switch,&
       num_stoch_bc,&
       stoch_bc_nbc,&
       stoch_bc_stddevs,&
       pres_vel_switch,&
       num_pres_vel,&
       vdefault_sysvarmem,&
       tcon_switch,&
       tcon_inds,&
       num_tcon,&
       tcon_sysvarmem,&
       mem,&
       sysvarmem,&
       rank_mem


  implicit none

  double precision :: kz, lz, por
  external kz, lz, por

  integer, external :: index_loc_to_mem

  double precision, external :: varpar

  integer :: irens, i, j, k, ivar, iv
  integer :: rmem, imem
  
  ! Set rank_mem
  rmem = 0
  do ivar = 1, num_enkf_vars
     if(act_s(ivar) == 1) then
        rank_mem(ivar) = rmem
        rmem = rmem +1
     else
        rank_mem(ivar) = -1
     end if
  end do

  ! Set mem and sysvarmem
  do ivar = 1, num_enkf_vars
     if(act_s(ivar) == 1) then

        DO irens = 1, nrens
           DO k = 1, k0
              DO j = 1, j0
                 DO i = 1, i0

                    imem = index_loc_to_mem(i,j,k,ivar)

                    if (imem > 0) then
                       mem(imem,irens) = varpar(i,j,k,irens,ivar)
                       sysvarmem(imem) = sysvar(ivar)
                    end if

                 END DO
              END DO
           END DO
        END DO

     end if
  end do

  ! Boundary conditions
  if (stoch_bc_switch) then
     do irens = 1, nrens
        do imem = nstate-num_stoch_bc+1, nstate
           mem(imem,irens) = dbc_data(stoch_bc_nbc(imem-nstate+num_stoch_bc),1,irens)
           ! Set system variance according to perturbation
           sysvarmem(imem) = stoch_bc_stddevs(1)/100.0d0
        end do
     end do
  end if
  
  ! Thermal conductivity
  if (tcon_switch) then
     
     do irens = 1, nrens
        do imem = nstate-num_pres_vel-num_tcon+1,nstate-num_pres_vel
           mem(imem,irens) = propunit(tcon_inds(imem-nstate+num_pres_vel+num_tcon),idx_lz,irens)
           sysvarmem(imem) = tcon_sysvarmem**2
        end do
     end do
     
  end if
  
  ! Prescribed Flow
  if (pres_vel_switch) then

     do irens = 1, nrens
        do imem = nstate-num_pres_vel+1, nstate
           iv = imem-nstate+num_pres_vel
           mem(imem,irens) = vdefault(iv,irens)
           sysvarmem(imem) = vdefault_sysvarmem(iv)**2
        end do
     end do
     
  end if

end subroutine enkf_make_state_vector
