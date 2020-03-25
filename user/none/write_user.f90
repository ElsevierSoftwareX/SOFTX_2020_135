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

!> @brief special output defined by the developer of this working directory (example routine)
!> @param[in] ident iteration number
!> @param[in] ismpl local sample index
!> @details
!> This routine can be taken as a starting point for writing
!> additional user-defined output routines.
      subroutine write_user(ident,ismpl)

        implicit none

        ! Sample index
        integer :: ismpl

        ! iteration number
        integer, intent (in) ::ident

        ! common character variables
        character (len=80) :: filename
        character (len=8) :: snumber


!     transform the typical suffix into a character variable for the filename
        if (ident>=0) then
          write(snumber,'(1I7)') ident
        else if (ident==-1) then
          write(snumber,'(A8)') 'final'
        else if (ident==-2) then
          write(snumber,'(A8)') 'debug'
        else if (ident==-3) then
          write(snumber,'(A8)') 'ens_mean'
        else if (ident==-4) then
          write(snumber,'(A8)') 'mean'
        else if (ident==-5) then
          write(snumber,'(A8)') 'ens_mean'
        end if
!
! ----- begin main body -----
!     ...
! ------ end main body ------
!
        return

      end subroutine write_user
