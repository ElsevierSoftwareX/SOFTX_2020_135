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

!>    @brief Set variable time step size variables delt_count and flag_delt
!>    @param[in] iter_nl nonlinear iteration counter
!>    @param[in] ismpl local sample index
!>    @details
!> Set variable time step size variables delt_count and flag_delt
subroutine set_var_deltat(iter_nl, ismpl)

  use arrays
  use mod_genrl
  use mod_linfos

  implicit none

  ! local sample index
  integer :: ismpl

  ! nonlinear iteration counter
  integer, intent (in) :: iter_nl

  ! Set delt_count
  ! ------------------

  ! Test whether to hasten or prevent time-step doubling
  IF (iter_nl+iter_nlold.LT.maxiter_nl/2) THEN
    delt_count(ismpl) = MIN(delt_count(ismpl)+2,delt_double-1)
  ELSE IF (iter_nl+iter_nlold.LT.maxiter_nl) THEN
    delt_count(ismpl) = MIN(delt_count(ismpl)+1,delt_double-1)
  ELSE IF (iter_nl+iter_nlold.GT.maxiter_nl) THEN
    delt_count(ismpl) = MAX(delt_count(ismpl)-1,0)
  END IF

  ! Standard output
  IF (linfos(3)>=2) WRITE (*,*) 'delt_count', delt_count(ismpl), &
      'delt_double', delt_double, 'iter_nlold', iter_nlold

  ! Set counter for next iteration
  iter_nlold = iter_nl/4 + (3*iter_nlold)/4

  ! Set flag_delt
  ! ------------------

  ! Test if nonlinear iteration reached maximum iteration count
  IF (((iter_nl.EQ.maxiter_nl).OR.(flag_delt(ismpl).EQ.(-1)))) THEN

    IF (iter_nl.EQ.maxiter_nl .AND. (nlconverge .eq. 0)) THEN
      WRITE(*,*) "Nonlinear iteration reached maximum iteration ", iter_nl, ". Set flag_delt to -2."
    END IF

    IF (flag_delt(ismpl).EQ.(-1)) THEN
      WRITE(*,*) "Iterative Solver reached maximum iteration. Set flag_delt to -2."
    END IF

    flag_delt(ismpl) = -2

  ELSE IF (flag_delt(ismpl).EQ.0) THEN

    flag_delt(ismpl) = flag_delt(ismpl) + 1

  END IF

end subroutine set_var_deltat
