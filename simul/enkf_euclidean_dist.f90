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

!> @brief Calculate Euclidean distance between two locations
!> @param[in] i x-index of first location
!> @param[in] j y-index of first location
!> @param[in] k z-index of first location
!> @param[in] l x-index of second location
!> @param[in] m y-index of second location
!> @param[in] n z-index of second location
!> @param[out] dist Euclidean distance
!> @param[out] distx x-component of Euclidean distance
!> @param[out] disty y-component of Euclidean distance
!> @param[out] distz z-component of Euclidean distance
!> @details
!> Euclidean distance and the components of the distance in x, y, z
!> direction are calculated.
subroutine enkf_euclidean_dist(i,j,k,l,m,n,dist,distx,disty,distz)

  use arrays, only:&
       delxa,&
       delya,&
       delza

  implicit none

  !First point location
  integer, intent(in) :: i
  integer, intent(in) :: j
  integer, intent(in) :: k

  !Second point location
  integer, intent(in) :: l
  integer, intent(in) :: m
  integer, intent(in) :: n

  double precision, intent(out) :: dist
  double precision, intent(out) :: distx
  double precision, intent(out) :: disty
  double precision, intent(out) :: distz

  ! Euclidean distance
  dist = (delxa(i)-delxa(l))*(delxa(i)-delxa(l))&
       + (delya(j)-delya(m))*(delya(j)-delya(m))&
       + (delza(k)-delza(n))*(delza(k)-delza(n))
  dist = dsqrt(dist)

  ! Distance x, y, z components
  distx = abs(delxa(i)-delxa(l))
  disty = abs(delya(j)-delya(m))
  distz = abs(delza(k)-delza(n))

end subroutine enkf_euclidean_dist
