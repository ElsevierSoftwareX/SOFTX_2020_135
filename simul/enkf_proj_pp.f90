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

!> @brief Calculate Pilot Point EnKF projection operator and interpolation
!> @details
!> Calculate pilot point projection operator and interpolation using
!> the covariance matrix `gss`.
!> 1. Fill matrix `gpp`, the part of the covariance matrix of the
!> state vector containing covariances between parameters at Pilot
!> Points.
!> 2. Fill matrix `grp`, the part of the covariance matrix of the
!> state vector containing covariances between parameter at non-pilot
!> points and parameters at Pilot Points.
!> 3. Calculate LU-decomposition of `gpp`.
!> 4. Calculate inverse of `gpp`. The inverse of `gpp` is saved in
!> `gpp`.
!> 5. Calculate `prp` or `Grp*Gpp^{-1}`, the part of the projection
!> operator needed for the interpolation.
!> 6. Fill `mem_up_pp` with parameter updates at Pilot Points.
!> 7. Calculate the interpolation. From `prp` and `mem_up_pp`, calculate
!> `mem_up_rp`, the interpolated parameter updates between Pilot Points.
!> 8. Generate full state vector `mem` by merging `mem_temp`, the
!> values of dynamic states and parameters at Pilot Points with
!> `mem_up_rp` the interpolated parameters updates between Pilot Points.
!> 9. Deallocate Pilot Point arrays.
subroutine enkf_proj_pp()

  use mod_enkf, only: &
       mem,&
       nstate,&
       lstate,&
       nrens,&
       gss,&
       pp_state_inds,&
       pp_lstate_inds,&
       num_pp,&
       mem_pp,&
       mem_rp,&
       ipp_start,&
       ipp_end

  implicit none

  integer :: irens, imem, ipp, jpp, irp
  double precision, allocatable :: gpp(:,:)
  double precision, allocatable :: grp(:,:)
  double precision, allocatable :: prp(:,:)
  double precision, allocatable :: mem_up_pp(:,:)
  double precision, allocatable :: mem_up_rp(:,:)
  double precision, allocatable :: mem_temp(:,:)
  integer, allocatable :: gpp_ipiv(:)
  integer :: ierr, lwork
  double precision, allocatable :: work(:)

  ! 1. Allocate and fill gpp
  allocate(gpp(num_pp,num_pp))
  do jpp = 1, num_pp
     do ipp = 1, num_pp
        gpp(ipp,jpp) = gss(pp_lstate_inds(ipp),pp_lstate_inds(jpp))
     end do
  end do

  ! 2. Allocate and fill grp
  allocate(grp(lstate-num_pp,num_pp))
  do jpp = 1, num_pp
     irp = 0
     do imem = 1, lstate
        if( .not. any(pp_lstate_inds==imem)) then
           irp = irp+1
           grp(irp,jpp) = gss(imem,pp_lstate_inds(jpp))
        end if
     end do
  end do

  ! 3. Calculate LU decomposition of gpp
  allocate(gpp_ipiv(num_pp))
  call dgetrf(num_pp,num_pp,gpp,num_pp,gpp_ipiv,ierr)
  if (ierr/=0) then
     write (unit = *, fmt = *) '[E] enkf_proj_pp: Error in dgetrf.f'
  end if

  ! 4. Invert LU decomposition of gpp, inverse matrix saved in gpp
  lwork = num_pp
  allocate(work(lwork))
  call dgetri(num_pp,gpp,num_pp,gpp_ipiv,work,lwork,ierr)
  if (ierr/=0) then
     write (unit = *, fmt = *) '[E] enkf_proj_pp: Error in dgetri.f'
  end if

  ! 5. Calculate prp, part of the projection operator
  allocate(prp(lstate-num_pp,num_pp))
  CALL dgemm('n','n',lstate-num_pp,num_pp,num_pp,1.0d0,grp,lstate-num_pp,gpp,num_pp, &
       0.0d0,prp,lstate-num_pp)

  ! 6. Calculate mem_up_pp, Pilot Point updates
  allocate(mem_up_pp(num_pp,nrens))
  do irens = 1, nrens
     ipp = 1
     do imem = ipp_start, ipp_start+num_pp-1
        mem_up_pp(ipp,irens) = mem(imem,irens) - mem_pp(imem,irens)
        ipp = ipp + 1
     end do
  end do

  ! 7. Calculate interpolation of updates and save in mem_up_rp
  allocate(mem_up_rp(lstate-num_pp,nrens))
  call dgemm('n','n',lstate-num_pp,nrens,num_pp,1.0d0,prp,lstate-num_pp,mem_up_pp,num_pp,0.0d0,mem_up_rp,lstate-num_pp)

  ! 8. Merge (mem_rp + mem_up_rp) and mem_temp to get full state vector
  allocate(mem_temp(nstate-lstate+num_pp,nrens))
  do irens = 1,nrens
     do ipp = 1, nstate-lstate+num_pp
        mem_temp(ipp,irens) = mem(ipp,irens)
     end do
  end do

  deallocate(mem)
  allocate(mem(nstate,nrens))
  do irens = 1, nrens
     irp = 1
     ipp = 1
     jpp = 1
     do imem = 1, nstate
        if (imem<ipp_start) then
           mem(imem,irens) = mem_temp(ipp,irens)
           ipp = ipp + 1
        else if (imem > ipp_end) then
           mem(imem,irens) = mem_temp(ipp,irens)
           ipp = ipp + 1
        else
           if ((jpp <= num_pp) .and. (jpp >= 1)) then
              if (imem == pp_state_inds(jpp)) then
                 mem(imem,irens) = mem_temp(ipp,irens)
                 ipp = ipp + 1
                 jpp = jpp + 1
              else
                 mem(imem,irens) = mem_rp(irp,irens) + mem_up_rp(irp,irens)
                 irp = irp + 1
              end if
           else
              mem(imem,irens) = mem_rp(irp,irens) + mem_up_rp(irp,irens)
              irp = irp + 1
           end if
        end if
     end do
  end do

  ! 9. Deallocate Pilot-Point arrays
  deallocate(gpp)
  deallocate(grp)
  deallocate(gpp_ipiv)
  deallocate(work)
  deallocate(prp)
  deallocate(mem_up_pp)
  deallocate(mem_up_rp)
  deallocate(mem_temp)
  deallocate(mem_pp)
  deallocate(mem_rp)

end subroutine enkf_proj_pp
