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

!> @brief Return a variable/parameter value
!> @param[in] i grid index in x direction
!> @param[in] j grid index in y direction
!> @param[in] k grid index in z direction
!> @param[in] irens ensemble member/realization index
!> @param[in] ienkfvar EnKF variable/parameter index
!> @details
!> Return a variable/parameter value according to the location, the
!> realization and the variable-index of the state vector
!> candidates head, temp, conc, kz, lz, por.
!>
!> head, temp, conc are read from the corresponding arrays from the
!> Fortran-module arrays.f90, whereas kz, lz, por are read from the
!> corresponding functions from the SHEMAT-property module. These
!> functions themselves read the values from the array propunit and
!> sometimes manipulate them.
!>
!> At the moment, lz is written directly from the
!> propunit-array. The function lz will return a mixture of matrix
!> thermal conductivity and fluid thermal conductivity. The value
!> read from propunit is the pure matrix conductivity, which is
!> more suitable for EnKF-updates.
double precision function varpar(i,j,k,irens,ienkfvar)

  use arrays, only:&
       head,&
       temp,&
       conc,&
       propunit,&
       uindex,&
       idx_lz
       
  implicit none

  integer, intent(in) :: i, j, k
  integer, intent(in) :: irens
  integer, intent(in) :: ienkfvar

  double precision :: kz, lz, por
  external kz, lz, por

  select case (ienkfvar)
  case(1)
     varpar = head(i,j,k,irens)
  case(2)
     varpar = temp(i,j,k,irens)
  case(3)
     varpar = conc(i,j,k,1,irens)
  case(4)
     varpar = dlog10(kz(i,j,k,irens)+1.0d-30)
  case(5)
     varpar = propunit(uindex(i,j,k),idx_lz,irens)!!lz(i,j,k,irens)
  case(6)
     varpar = por(i,j,k,irens)
  case default
     write(unit = *, fmt = *) '[E1] Error in enkf_varpar()'
     stop 1
  end select
  
  return
     
end function varpar

!> @brief Setting variables/parameters from state vector
!> @param{in] i x-index of variable/parameter location
!> @param{in] j y-index of variable/parameter location
!> @param{in] k z-index of variable/parameter location
!> @param[in] irens ensemble member/realization index
!> @param[in] imem Index of variable/parameter location in state vector
!> @param[in] ivar Index of variable/parameter in state vector
!> @details
!> This subroutine sets an entry of a dynamic variable (head, temp,
!> conc) or a parameter (propunit, kz, lz, por) specified by the
!> index ivar. The location is specified by the location-indices i,
!> j, k and the number of the ensemble member by the ensemble index
!> irens. The index in the EnKF-state vector is given by imem.
!>
!> __ATTENTION__: For conc, only one kind of tracer is supported. For
!> the parameters saved in propunit, the subroutine stab_param is
!> called to ensure that the assimilated values are inside the
!> physical regime.
subroutine enkf_set_varpar(i,j,k,irens,imem,ivar)

  use arrays, only:&
       head, temp, conc,&
       propunit,&
       uindex,&
       idx_kz, idx_lz, idx_por,&
       nunits, nprop

  use mod_enkf, only:&
       mem
  
  implicit none

  integer, intent(in) :: i,j,k
  integer, intent(in) :: imem
  integer, intent(in) :: irens
  integer, intent(in) :: ivar

  select case (ivar)
  case(1)
     head(i,j,k,irens) = mem(imem,irens)
  case(2)
     temp(i,j,k,irens) = mem(imem,irens)
  case(3)
     ! concentration field may contain a number of species, here only
     ! the first is considered
     conc(i,j,k,1,irens) = mem(imem,irens)
  case(4)
     propunit(uindex(i,j,k),idx_kz,irens) = 10.0d0**mem(imem,irens)
     call stab_param(propunit(uindex(i,j,k),idx_kz,irens),idx_kz,uindex(i,j,k))
  case(5)
     propunit(uindex(i,j,k),idx_lz,irens) = mem(imem,irens)
     call stab_param(propunit(uindex(i,j,k),idx_lz,irens),idx_lz,uindex(i,j,k))
  case(6)
     propunit(uindex(i,j,k),idx_por,irens) = mem(imem,irens)
     call stab_param(propunit(uindex(i,j,k),idx_por,irens),idx_por,uindex(i,j,k))
  end select
  
end subroutine enkf_set_varpar
