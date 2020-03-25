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

!> @brief Set covariance matrix for Pilot Point EnKF
!> @details
!> Calculate the covariance matrix `gss`, which will be used in
!> Pilot-Point EnKF to define both, the projection operator and
!> interpolation operator.
!> 1. Extract the Pilot-Point part of the state vector and subtract
!> the mean (saved in array `memp`)
!> 2. Calculate the ensemble covariance matrix using dgemm (saved in
!> `gss`)
!> 3. Read existing `gss` from file or write `gss` to file.
!>
!> __ATTENTION__: If `gss` is written, a large number of ensemble
!> members should be used (minimum 10000) in order to keep sampling
!> errors at a minimum.
subroutine enkf_cov_pp()

  use mod_enkf, only: &
       mem,&
       ave,&
       lstate,&
       nrens,&
       gss,&
       pp_get_out,&
       ipp_start,&
       ipp_end

  implicit none

  integer :: i
  integer :: j

  double precision, allocatable :: memp(:,:)

  integer, external :: index_loc_to_mem

  allocate(gss(lstate,lstate))

  ! Calculate gss
  select case (pp_get_out)
  case ('def','out')
     allocate(memp(lstate,nrens))

     ! Subtract mean from mem at right indices
     do i = ipp_start, ipp_end
        do j = 1, nrens
           memp(1+i-ipp_start,j) = mem(i,j)-ave(i)
        end do
     end do

     ! Calculate ensemble covariance: 1/(N-1)*XX^T, where X=memp
     CALL dgemm('n','t',lstate,lstate,nrens,1.0d0/(real(nrens)-1.0d0),memp,lstate,memp,lstate, &
          0.0d0,gss,lstate)

     deallocate(memp)
  end select

  ! Read or write gss
  select case (pp_get_out)
  case ('get')
     open(unit=63, file='gss.dat', status='old', action='read')
     read(unit=63, fmt = '(e20.10)') gss
     close(unit=63)
  case ('out')
     open(unit=63, file='gss.dat', status='replace', action='write')
     write (unit = 63, fmt = '(e20.10)') gss
     close(unit=63)
     write (unit = *, fmt = *) "[OK] Written gss to gss.dat."
     stop
  end select

end subroutine enkf_cov_pp
