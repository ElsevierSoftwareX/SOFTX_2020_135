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

!> @brief compute "i,j,k"-index of "m", depending on "i0,j0,k0"
!> @param[in] m continuous memory cell-index
!> @param[out] i cell index, direction I0
!> @param[out] j cell index, direction J0
!> @param[out] k cell index, direction K0
!> @details
!> Compute k cell index: \n\n
!>
!>    k = (m-1) / (i0*j0) + 1 \n\n
!>
!> Explanations:\n
!> 1) m-1: The integer steps between k-indices are at elements m =
!>    n*(i0*j0 + 1) for some positive integer n [f.e.: m=i0*j0 is
!>    still in the lowest layer, m=i0*j0+1 is in the second lowest
!>    layer]. The steps of the division are at m_div = n*(i0*j0). \n
!> 2) +1 at end: The division gives k = 0 for the lowest layer, k = 1
!>    for the second-lowest, and so on. The layer-counting starts at
!>    k=1, so 1 is added. \n\n
!>
!> Compute j cell index: \n\n
!>
!>    j = (m-1 - (k-1)*i0*j0) / i0 + 1 \n\n
!>
!> Explanations:\n
!> 1) (m-1 - (k-1)*i0*j0): The index is projected onto its
!>    partner-index on the lowest layer by subtracting the number of
!>    z-layers on top, then one is subtracted, since the steps between
!>    j-indices should be at m=n*(i0+1), but the steps of the division
!>    by i0 are at m=n*i0. \n
!> 2) +1 at end: The division by i0 gives j = 0 for the first line in
!>    x-direction, j = 1 for the second, and so on. The line-counting
!>    should start at j=1, so 1 is added. \n\n
!>
!> Compute i cell index: \n\n
!>
!>    i = m - (k-1)*i0*j0 - (j-1)*i0\n\n
!>
!> Explanations:\n
!> 1) (m - (k-1)*i0*j0): The index is projected onto its
!>    partner-index on the lowest layer by subtracting the number of
!>    z-layers on top. \n
!> 2) - (j-1)*i0: The index is projected onto its
!>    partner-index on the first line by subtracting the number of
!>    y-lines . \n\n
      subroutine ijk_m(m,i,j,k)

        use mod_genrl

        implicit none

        ! Continuous linear index
        integer, intent (in) :: m

        ! i-cell/j-cell/k-cell index
        integer, intent (out) :: i, j, k


        ! Compute k-cell index
        k = (m-1) / (i0*j0) + 1

        ! Compute j-cell index
        j = (m-1 - (k-1)*i0*j0) / i0 + 1

        ! Compute i-cell index
        i = m - (k-1)*i0*j0 - (j-1)*i0

        return

      end subroutine ijk_m
