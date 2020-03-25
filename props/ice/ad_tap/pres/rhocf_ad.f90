!        Generated by TAPENADE     (INRIA, Ecuador team)
!  tapenade 3.x
!
!  Differentiation of rhocf in reverse (adjoint) mode (with options noISIZE i8):
!   gradient     of useful results: *temp *pres rhocf
!   with respect to varying inputs: *temp *pres
!   Plus diff mem management of: temp:in pres:in
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
!>    @brief calculates heat capacity times density of water.
!>    @param[in] i cell index, direction I0
!>    @param[in] j cell index, direction J0
!>    @param[in] k cell index, direction K0
!>    @param[in] ismpl local sample index
!>    @return  rhocf                                rhocf[J / (m^3 * K)]
!>    @details
!>    input:\n
!>      heat capacity of water              cpf [J / kg / K]\n
!>      density                             rhof [kg/m^3]\n
!>      pressure                            plocal [Pa]\n
!>      temperature                         tlocal [C]\n
SUBROUTINE RHOCF_AD(i, j, k, ismpl, rhocf_adv)
  use arrays

  USE ARRAYS_AD

  IMPLICIT NONE
  double precision :: rhocf_adv
  INTEGER :: i, j, k, ismpl
  EXTERNAL CPF, RHOF
  EXTERNAL CPF_AD, RHOF_AD
  DOUBLE PRECISION :: plocal, tlocal, rfluid, cfluid, CPF, &
& RHOF
  DOUBLE PRECISION :: rfluid_ad, cfluid_ad
  DOUBLE PRECISION :: rhocf
!     water density [kg/m**3]
  rfluid = RHOF(i, j, k, ismpl)
!     water isobariC heat CapaCity [J/(kg*K)]
  cfluid = CPF(i, j, k, ismpl)
  rfluid_ad = cfluid*rhocf_adv
  cfluid_ad = rfluid*rhocf_adv
  CALL CPF_AD(i, j, k, ismpl, cfluid_ad)
  CALL RHOF_AD(i, j, k, ismpl, rfluid_ad)
END SUBROUTINE RHOCF_AD
