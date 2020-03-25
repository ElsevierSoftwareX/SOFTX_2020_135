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

!> @brief Copy variables/parameters from previous EnKF-iteration
!> @details
!> old_restore is called, which copies formerly save values into
!> the current variable arrays. Then the boundary condition data
!> is copied.
subroutine enkf_old_restore()

  use arrays, only:&
       nbc_data,&
       dbc_data,&
       dbc_dataold

  use mod_genrl, only:&
       nsmpl,&
       cgen_opti

  implicit none

  integer :: irens, i

  ! INIT to recompute ensembles
  DO irens = 1, nsmpl
     ! init state-variables state (before)
     CALL old_restore(cgen_opti,irens)
     ! init "DBC" state (before)
     DO i = 1, nbc_data
        dbc_data(i,1,irens) = dbc_dataold(i)
     END DO
  END DO

end subroutine enkf_old_restore
