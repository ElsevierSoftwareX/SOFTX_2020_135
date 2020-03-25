!        Generated by TAPENADE     (INRIA, Ecuador team)
!  tapenade 3.x
!
!  Differentiation of rhocm in reverse (adjoint) mode (with options noISIZE i8):
!   gradient     of useful results: *temp *propunit rhocm
!   with respect to varying inputs: *temp *propunit
!   Plus diff mem management of: temp:in propunit:in
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
!>    @brief calculates product of heat capacity and density of rock.
!>    @param[in] i cell index, direction I0
!>    @param[in] j cell index, direction J0
!>    @param[in] k cell index, direction K0
!>    @param[in] ismpl local sample index
!>    @return  rhocm                                rc   [J / (m^3 K)]
!>    @details
!>    input:\n
!>      pressure                            plocal [Mpa]\n
!>      temperature                         tlocal in [C]\n
SUBROUTINE RHOCM_AD0(i, j, k, ismpl, rhocm_ad)
  use arrays

  USE ARRAYS_AD

  USE MOD_TEMP
  IMPLICIT NONE
  INTEGER :: i, j, k, ismpl
  DOUBLE PRECISION :: plocal, tlocal
  DOUBLE PRECISION :: tlocal_ad
  DOUBLE PRECISION :: temporary_ad
  DOUBLE PRECISION :: rhocm
  DOUBLE PRECISION :: rhocm_ad
  tlocal = temp(i, j, k, ismpl)
  propunit_ad(uindex(i, j, k), idx_rc, ismpl) = propunit_ad(uindex(i, j&
&   , k), idx_rc, ismpl) + (cma1+cma2*tlocal+cma3*tlocal**2)*rhocm_ad
  temporary_ad = propunit(uindex(i, j, k), idx_rc, ismpl)*rhocm_ad
  tlocal_ad = (cma2+2*tlocal*cma3)*temporary_ad
  temp_ad(i, j, k, ismpl) = temp_ad(i, j, k, ismpl) + tlocal_ad
END SUBROUTINE RHOCM_AD0
