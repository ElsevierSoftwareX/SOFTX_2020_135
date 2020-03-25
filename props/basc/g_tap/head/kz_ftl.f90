!        Generated by TAPENADE     (INRIA, Ecuador team)
!  tapenade 3.x
!
!  Differentiation of kz in forward (tangent) mode:
!   variations   of useful results: kz
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
!> @brief assign permeability in z direction to cell
!> @param[in] i grid indices
!> @param[in] j grid indices
!> @param[in] k grid indices
!> @param[in] ismpl local sample index
!> @return  permeability                        (m^2)
!> @details
!> kz returns the permeability in z-direction[m2] at node(i,j,k) from
!> the input file.\n
DOUBLE PRECISION FUNCTION g_KZ(i, j, k, ismpl, kz)
  USE ARRAYS

  USE g_ARRAYS

  IMPLICIT NONE
! Location indices
  INTEGER, INTENT(IN) :: i
  INTEGER, INTENT(IN) :: j
  INTEGER, INTENT(IN) :: k
! Sample index
  INTEGER :: ismpl
  DOUBLE PRECISION :: kz
  g_kz = g_propunit(uindex(i, j, k), idx_kz, ismpl)
  kz = propunit(uindex(i, j, k), idx_kz, ismpl)
  RETURN
END FUNCTION g_KZ
