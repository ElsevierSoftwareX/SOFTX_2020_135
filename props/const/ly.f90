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

!> @brief calculates effective thermal conductivity of the cell
!> @param[in] i cell index, direction I0
!> @param[in] j cell index, direction J0
!> @param[in] k cell index, direction K0
!> @param[in] ismpl local sample index
!> @return  thermal conductivity ly[W/(m*K)]
!> @details
!> calculates effective thermal conductivity of the two phase system
!> matrix-porosity .\n\n
!>
!> input:\n
!> porosity                            porlocal [-]\n
!> thermal conductivity of fluid           lamf [W/(m*K)]\n
!> thermal conductivity of matrix       lammref [W/(m*K)]\n
      double precision function ly(i,j,k,ismpl)
        use arrays

        implicit none

        ! Location indices
        integer, intent (in) :: i
        integer, intent (in) :: j
        integer, intent (in) :: k

        ! Sample index
        integer :: ismpl

        ! Local fluid thermal conductivity [W/(m*K)]
        double precision :: lamfluid
        double precision, external :: lamf

        ! Matrix thermal conductivity [W/(m*K)]
        double precision :: lammref

        ! Local porosity [-]
        double precision, external :: por


        ! Local fluid thermal conductivity [W/(m*K)]
        lamfluid = lamf(i,j,k,ismpl)

        ! Reference matrix thermal conductivity [W/(m*K)]
        lammref = propunit(uindex(i,j,k),idx_lz,ismpl)* &
          propunit(uindex(i,j,k),idx_an_ly,ismpl)

        ly = lammref**(1.0d0-por(i,j,k,ismpl))*lamfluid**por(i,j,k,ismpl)

        return

      end function ly
