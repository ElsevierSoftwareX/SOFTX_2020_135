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

!> @brief Read the synthetic true model
!> @details
!> __USED__: Variables defining the length of the state vector and
!> which nodes/variables are taken into account or not.
!> Names of the observation input files (normal and
!> chemical)
!>
!> __SET__: Values of the true model for: Permeability, hydraulic
!> head, temperature, concentration, thermal conductivity
!> and porosity in the array `var_true`.
subroutine enkf_read_true()

  use mod_enkf, only: &
       true_name,&
       true_chem_name,&
       lstate0,&
       ex_unit,&
       act_s,&
       var_true

  implicit none

  double precision :: adu
  integer :: idu
  integer :: i, j, k
  double precision :: kz_t, h_t, t_t, c_t, lz_t, por_t
  integer :: l, iunit, lmem

  integer, external :: index_loc_to_lin

  !-----------------------------------------------------------
  ! NO CELL-CENTERED allowed for True input
  !-----------------------------------------------------------
  open(unit = 19, file = true_name, status = 'old')

  read(unit = 19,fmt = *)
  read(unit = 19,fmt = *)
  read(unit = 19,fmt = *)

  do l = 1, lstate0
     read(unit = 19, fmt = *) i, j, k, adu, adu, adu, iunit, &
          idu, idu, idu, h_t, t_t, adu, adu, adu, adu, adu, adu, por_t, adu, &
          adu, kz_t, adu, adu, lz_t, adu, adu, adu, adu, adu, &
          adu, adu, adu, adu, adu, adu

     if (iunit .gt. ex_unit) then
        lmem=index_loc_to_lin(i,j,k)
        var_true(lmem,1) = h_t
        var_true(lmem,2) = t_t
        var_true(lmem,4) = kz_t
        var_true(lmem,5) = lz_t
        var_true(lmem,6) = por_t
     end if
  end do

  close(unit = 19)

  if(act_s(3) == 1) then
     !true_chem_name: Set in the *.enkf input file
     open(unit = 19, file = true_chem_name, status = 'old')

     read(unit = 19, fmt = *)
     read(unit = 19, fmt = *)
     read(unit = 19, fmt = *)

     do l = 1, lstate0
        read(unit = 19, fmt = *) adu, adu, adu, iunit, i, j, k, &
             idu, idu, idu, h_t, t_t, adu, adu, adu, adu, adu, c_t
        if (iunit .gt. ex_unit) then
           lmem=index_loc_to_lin(i,j,k)
           var_true(lmem,3) = c_t
        end if
     end do

     close(unit = 19)

  end if

end subroutine enkf_read_true
