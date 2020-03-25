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

!>    @brief Compute static relaxation for flow and temperature
!>    @param[in] ijk number of cells
!>    @param[inout] ismpl local sample index
!>    @details
!> Static relaxation is computed for variable arrays head, pres, temp:
!> \n\n
!>
!> var = (1-theta)*varold + theta*var \n
subroutine static_relaxation(ijk,ismpl)

  use arrays
  use mod_genrl
  use mod_time

  implicit none

  ! local sample index
  integer :: ismpl
  
  ! Number of cells
  integer, intent (in) :: ijk


  ! flow
  if (thetaf /= 1.0d0)  then

#ifdef head_base
    CALL dscal(ijk,thetaf,head(1,1,1,ismpl),1)
    CALL daxpy(ijk,1.0D0-thetaf,headold(1,cgen_time,ismpl),1,head(1,1,1,ismpl),1)
#endif
#ifdef pres_base
    CALL dscal(ijk,thetaf,pres(1,1,1,ismpl),1)
    CALL daxpy(ijk,1.0D0-thetaf,presold(1,cgen_time,ismpl),1,pres(1,1,1,ismpl),1)
#endif

  end if

  ! temperature
  if (thetat /= 1.0d0) then

    CALL dscal(ijk,thetat,temp(1,1,1,ismpl),1)
    CALL daxpy(ijk,1.0D0-thetat,tempold(1,cgen_time,ismpl),1,temp(1,1,1,ismpl),1)

  end if

  return
  
end subroutine static_relaxation
