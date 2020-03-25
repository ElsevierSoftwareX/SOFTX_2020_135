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

!>    @brief wrapper for the forward simulation call
!>    @param[in] iseed 0: FW simulation, 1 .. <mpara>: AD seeding index
!>    @param[in] iter time iteration counter
!>    @param[in] ismpl local sample index
!>    @details
!> seperates the equation system solving from the\n
!> specific nonlinear solver (Picard based)\n
      SUBROUTINE forward_wrapper(iter,iseed,ismpl)
        use arrays
        use mod_genrl
        use mod_genrlc
        use mod_time
        IMPLICIT NONE
        integer :: ismpl

        INTEGER iter, iseed

        IF (lread_joutt) THEN
          WRITE(project_sfx(ismpl),'(1A,1I7.7)') '_all', iter
!         read data of each time step
          CALL read_outt_hdf(ismpl)
          project_sfx(ismpl) = ' '
        ELSE
! ######### Forward Iteration ######
          CALL forward_picard(iter,ismpl)
! ##################################
        END IF
        RETURN
      END
