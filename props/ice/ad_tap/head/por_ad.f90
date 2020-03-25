!        Generated by TAPENADE     (INRIA, Ecuador team)
!  tapenade 3.x
!
!  Differentiation of por in reverse (adjoint) mode (with options noISIZE i8):
!   gradient     of useful results: *propunit por
!   with respect to varying inputs: *propunit
!   Plus diff mem management of: propunit:in
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
!>    @brief assign porosity to cell
!>    @param[in] i grid indices
!>    @param[in] j grid indices
!>    @param[in] k grid indices
!>    @param[in] ismpl local sample index
!>    @return  porosity                            porlocal [-]
SUBROUTINE POR_AD0(i, j, k, ismpl, por_ad)
  use arrays

  USE ARRAYS_AD

  IMPLICIT NONE
  INTEGER :: i, j, k, ismpl
  DOUBLE PRECISION :: por_ad
  DOUBLE PRECISION :: por
  propunit_ad(uindex(i, j, k), idx_por, ismpl) = propunit_ad(uindex(i, j&
&   , k), idx_por, ismpl) + por_ad
END SUBROUTINE POR_AD0

