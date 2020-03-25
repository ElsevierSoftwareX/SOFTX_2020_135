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

!> @brief calculates volumetric heat capacity of the cell
!> @param[in] i cell index, direction I0
!> @param[in] j cell index, direction J0
!> @param[in] k cell index, direction K0
!> @param[in] ismpl local sample index
!> @return volumetric heat capacity
!> @details
!> calculates volumetric heat capacity of the system
!> matrix-porosity [J/(K*m3)].\n
      double precision function rhoceff(i,j,k,ismpl)

        use arrays
        ! use mod_temp

        implicit none

        ! Location indices
        integer, intent (in) :: i
        integer, intent (in) :: j
        integer, intent (in) :: k

        ! Sample index
        integer :: ismpl

        ! Local temperature [degC]
        double precision :: tlocal

        ! Local porosity [-]
        double precision :: porlocal
        double precision, external :: por

        ! Matrix fraction in cell
        double precision :: fm

        ! Fluid fraction in cell
        double precision :: ff

        ! Heat capacity of the matrix
        double precision, external :: rhocm

        ! Heat capacity of the fluid
        double precision, external :: rhocf


        ! Local Temperature in degC
        tlocal = temp(i,j,k,ismpl)

        ! Local porosity
        porlocal = por(i,j,k,ismpl)

        ! Matrix fraction
        fm = 1.D0 - porlocal

        ! Fluid fraction
        ff = porlocal

        ! Heat capacity in cell, arithmetic mean
        rhoceff = ff*rhocf(i,j,k,ismpl) + fm*rhocm(i,j,k,ismpl)

        return

      end function rhoceff
