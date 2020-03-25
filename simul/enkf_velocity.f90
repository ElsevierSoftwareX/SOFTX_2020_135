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

!> @brief Set dbc_data according to prescribed velocity
!> @param[in] irens ensemble member/realization index
!> @details
!> Stepwise documentation
!> - 1. Parameter a: Factor turning physical distance into head difference
!> - 2. Set boundary condition data: Front, Left, Back, Right
!>  + a. Front: Start with given value at the front left corner from
!>              input file and add according to distance and
!>              velocity.
!>  + b. Left: Start with given value at the fron left corner and add.
!>  + c. Back: Start with final value of Left and add.
!>  + d. Right: Start with final value of Front and add.
subroutine enkf_velocity_dbc(irens)

  use arrays, only: &
       dbc_data, vdefault, delx, dely,&
       first_flow

  use mod_genrl, only: &
       i0, j0, k0, nsmpl
  
  use mod_flow, only: &
       grav

  implicit none

  integer, intent(in) :: irens

  double precision, allocatable, dimension(:,:) :: a
  integer :: i, ib
  
  DOUBLE PRECISION rhof
  EXTERNAL rhof
  DOUBLE PRECISION visf
  EXTERNAL visf
  DOUBLE PRECISION kz
  EXTERNAL kz

  allocate(a(3,nsmpl))

  !  Parameter a: dhx=a*dx, dhy=a*dy, dhz=a*dz
  do i=1,3
     a(i,irens) = (-visf(1,1,1,irens)/(kz(1,1,1,irens)*rhof(1,1,1,irens)*grav))*vdefault(i,irens)
  end do

  ! Front:
  call check_bc_index(first_flow,1,1,1)
  do i = 2, i0
     ib = first_flow-1+i
     dbc_data(ib,1,irens) = dbc_data(ib-1,1,irens)+a(1,irens)*0.5d0*(delx(i)+delx(i-1))
  end do
  call check_bc_index(first_flow-1+i0,i0,1,1)
  
  ! Left:
  call check_bc_index(first_flow+i0,1,2,1)
  dbc_data(first_flow+i0,1,irens) = dbc_data(first_flow,1,irens)+a(2,irens)*0.5d0*(dely(1)+dely(2))
  do i = 3, j0
     ib = first_flow-1+i0+i-1
     dbc_data(ib,1,irens) = dbc_data(ib-1,1,irens)+a(2,irens)*0.5d0*(dely(i)+dely(i-1))
  end do
  call check_bc_index(first_flow+i0+j0-2,1,j0,1)
  
  ! Back:
  call check_bc_index(first_flow+i0+j0-1,2,j0,1)
  do i = 2, i0
     ib = first_flow-1+i0+j0+i-2
     dbc_data(ib,1,irens) = dbc_data(ib-1,1,irens)+a(1,irens)*0.5d0*(delx(i)+delx(i-1))
  end do
  call check_bc_index(first_flow+2*i0+j0-3,i0,j0,1)

  ! Right:
  call check_bc_index(first_flow+2*i0+j0-2,i0,2,1)
  dbc_data(first_flow+2*i0+j0-2,1,irens) = dbc_data(first_flow-1+i0,1,irens) + a(2,irens)*0.5d0*(dely(1)+dely(2))
  do i = 3, j0-1
     ib = first_flow-1+2*i0+j0+i-3
     dbc_data(ib,1,irens) = dbc_data(ib-1,1,irens)+a(2,irens)*0.5d0*(dely(i)+dely(i-1))
  end do
  call check_bc_index(first_flow+2*i0+2*j0-5,i0,j0-1,1)
  
  if(allocated(a)) deallocate(a)
  
end subroutine enkf_velocity_dbc

!> @brief Set dbc_dataold array according to prescribed velocity
!> @details
!> Stepwise documentation
!> - 1. Parameter a: Factor turning physical distance into head difference
!> - 2. Set boundary condition data: Front, Left, Back, Right
!>  + a. Front: Start with given value at the front left corner from
!>              input file and add according to distance and
!>              velocity.
!>  + b. Left: Start with given value at the fron left corner and add.
!>  + c. Back: Start with final value of Left and add.
!>  + d. Right: Start with final value of Front and add.
subroutine enkf_velocity_dbcold()
  use arrays, only: &
       dbc_dataold, vdefault, delx, dely, &
       first_flow

  use mod_genrl, only: &
       i0, j0, k0
  
  use mod_flow, only: &
       grav

  implicit none

  double precision, allocatable, dimension(:) :: a
  integer :: i, ib
  
  DOUBLE PRECISION rhof
  EXTERNAL rhof
  DOUBLE PRECISION visf
  EXTERNAL visf
  DOUBLE PRECISION kz
  EXTERNAL kz

  allocate(a(3))

  !  Parameter: dhx=a*dx, dhy=a*dy, dhz=a*dz 
  do i=1,3
     a(i) = (-visf(1,1,1,1)/(kz(1,1,1,1)*rhof(1,1,1,1)*grav))*vdefault(i,1)
  end do

  ! Front:
  call check_bc_index(first_flow,1,1,1)
  do i = 2, i0
     ib = first_flow-1+i
     dbc_dataold(ib) = dbc_dataold(ib-1)+a(1)*0.5d0*(delx(i)+delx(i-1))
  end do
  call check_bc_index(first_flow-1+i0,i0,1,1)
  
  ! Left:
  call check_bc_index(first_flow+i0,1,2,1)
  dbc_dataold(first_flow+i0) = dbc_dataold(first_flow)+a(2)*0.5d0*(dely(1)+dely(2))
  do i = 3, j0
     ib = first_flow-1+i0+i-1
     dbc_dataold(ib) = dbc_dataold(ib-1)+a(2)*0.5d0*(dely(i)+dely(i-1))
  end do
  call check_bc_index(first_flow+i0+j0-2,1,j0,1)
  
  ! Back:
  call check_bc_index(first_flow+i0+j0-1,2,j0,1)
  do i = 2, i0
     ib = first_flow-1+i0+j0+i-2
     dbc_dataold(ib) = dbc_dataold(ib-1)+a(1)*0.5d0*(delx(i)+delx(i-1))
  end do
  call check_bc_index(first_flow+2*i0+j0-3,i0,j0,1)

  ! Right:
  call check_bc_index(first_flow+2*i0+j0-2,i0,2,1)
  dbc_dataold(first_flow+2*i0+j0-2) = dbc_dataold(first_flow-1+i0) + a(2)*0.5d0*(dely(1)+dely(2))
  do i = 3, j0-1
     ib = first_flow-1+2*i0+j0+i-3
     dbc_dataold(ib) = dbc_dataold(ib-1)+a(2)*0.5d0*(dely(i)+dely(i-1))
  end do
  call check_bc_index(first_flow+2*i0+2*j0-5,i0,j0-1,1)

  if(allocated(a)) deallocate(a)

end subroutine enkf_velocity_dbcold

!> @brief Check consistency of boundary conditions index and location
!> @param[in] ib Boundary index
!> @param[in] i Checked x-axis index
!> @param[in] j Checked y-axis index
!> @param[in] k Checked z-axis index
!> @details
!> Check that boundary condition index ib is associated with the
!> cell at location i, j, k.
subroutine check_bc_index(ib,i,j,k)

  use arrays, only:&
       ibc_data, cbc_i, cbc_j, cbc_k
  
  implicit none

  integer, intent(in) :: ib
  integer, intent(in) :: i
  integer, intent(in) :: j
  integer, intent(in) :: k

  if (ibc_data(ib,cbc_i) .ne. i &
       .or. ibc_data(ib,cbc_j) .ne. j &
       .or. ibc_data(ib,cbc_k) .ne. k) then
     write(unit = *, fmt = *) "[E] enkf_velocity.f90"
     write(unit = *, fmt = *) "  Needed (ib, i, j, k):"
     write(unit = *, fmt = *) ib, i, j, k
     write(unit = *, fmt = *) "  Boundary condition index (ib, i, j, k):"
     write(unit = *, fmt = *) ib, ibc_data(ib,cbc_i), ibc_data(ib,cbc_j), ibc_data(ib,cbc_k)
     stop 1
  end if
  
end subroutine check_bc_index
