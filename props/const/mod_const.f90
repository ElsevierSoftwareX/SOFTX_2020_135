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

!> @brief global variables and constants for PROPS=const.
module mod_const
  
  !> @brief Dimension of fluid property index array.
  !> @details
  !> Dimension of fluid property index array. \n
  !> Dimension of fluid property index array fprops.
  integer, parameter ::  npropsf = 5

  !> @brief Index of rhof in fluid property index array.
  !> @details
  !> Index of rhof in fluid property index array. \n
  integer, parameter :: pconst_rhof = 1

  !> @brief Index of compf in fluid property index array.
  !> @details
  !> Index of compf in fluid property index array. \n
  integer, parameter :: pconst_compf = 2

  !> @brief Index of cpf in fluid property index array.
  !> @details
  !> Index of cpf in fluid property index array. \n
  integer, parameter :: pconst_cpf = 3

  !> @brief Index of lamf in fluid property index array.
  !> @details
  !> Index of lamf in fluid property index array. \n
  integer, parameter :: pconst_lamf = 4

  !> @brief Index of visf in fluid property index array.
  !> @details
  !> Index of visf in fluid property index array. \n
  integer, parameter :: pconst_visf = 5

  !> @brief Fluid property index array.
  !> @details
  !> Fluid property index array. \n
  !> The array contains indices for parameters rhof, compf, cpf, lamf,
  !> visf.
  double precision, dimension (npropsf) :: fprops

end module mod_const
