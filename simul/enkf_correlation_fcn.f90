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

!> @brief Correlation as function of Euclidean distance
!> @param[in] lx x-component of covariance localization length scale
!> @param[in] ly y-component of covariance localization length scale
!> @param[in] lz z-component of covariance localization length scale
!> @param[in] bx x-component of Euclidean distance
!> @param[in] by y-component of Euclidean distance
!> @param[in] bz z-component of Euclidean distance
!> @param[out] corr Correlation as function of Euclidean distance
!> @details
!> The correlation is calculated according to __Gaspari,Cohn 1999__
!> "Construction of Correlation Functions in Two and Three
!> Dimensions": Fifth order correlation function
subroutine enkf_correlation_fcn(lx,ly,lz,bx,by,bz,corr)

  implicit none

  double precision, intent(in) :: lx
  double precision, intent(in) :: ly
  double precision, intent(in) :: lz

  double precision, intent(in) :: bx
  double precision, intent(in) :: by
  double precision, intent(in) :: bz

  double precision, intent(out) :: corr

  double precision :: q
  double precision :: ax, ay, az

  !----------------------------------------------------
  
  ! Euclidean distance relative to length scale
  ! Including additional factor
  ax = dsqrt(10.0d0/3.0d0)*lx
  ay = dsqrt(10.0d0/3.0d0)*ly
  az = dsqrt(10.0d0/3.0d0)*lz

  q = dsqrt((bx/ax)**2 + (by/ay)**2 + (bz/az)**2)

  ! Function evaluation
  if( bx < 0 .or. by < 0 .or. bz < 0 .or. q < 0) then
     write(unit = *, fmt = *) '[E] Error in enkf_correlation_fcn.f90.'
  else if ( (0 .le. q) .and. (q .le. 1) ) then
     corr = -0.25d0 *         q*q*q*q*q&
          + 0.5d0 *           q*q*q*q&
          + 0.625d0 *         q*q*q&
          - 1.66666666666d0 * q*q&
          + 1.0d0
  else if ( (1 < q) .and. (q .le. 2.0d0) ) then
     corr = 0.08333333333d0 * q*q*q*q*q&
          - 0.5d0 *           q*q*q*q&
          + 0.625d0 *         q*q*q&
          + 1.66666666666d0 * q*q&
          - 5.0d0 *           q&
          + 4.0d0&
          -0.66666666666666d0/q
  else if (2.0d0 < q) then
     corr = 0.0d0
  end if

end subroutine enkf_correlation_fcn
