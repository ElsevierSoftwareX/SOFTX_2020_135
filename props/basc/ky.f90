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

!> @brief assign permeability in y direction to cell
!> @param[in] i grid indices
!> @param[in] j grid indices
!> @param[in] k grid indices
!> @param[in] ismpl local sample index
!> @return  permeability                        (m^2)
!> @details
!> ky returns the permeability in y-direction [m2] at node(i,j,k) from
!> the input file.\n\n
!>
!> The permeability in y-direction is the product of the permeability
!> in z-direction and the anisotropy factor for the y-direction.
      double precision function ky(i,j,k,ismpl)
        use arrays

        implicit none

        ! Location indices
        integer, intent (in) :: i
        integer, intent (in) :: j
        integer, intent (in) :: k

        ! Sample index
        integer :: ismpl


        ky = propunit(uindex(i,j,k),idx_kz,ismpl)* &
          propunit(uindex(i,j,k),idx_an_ky,ismpl)

        return

      end function ky
