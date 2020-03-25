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

!>    @brief initialisation routine, no "reinjection" functionality
!>    @param[in] ismpl local sample index
      SUBROUTINE user_init(ismpl)
        use arrays
        use mod_genrl
        use mod_linfos
        IMPLICIT NONE
        integer :: ismpl
        integer :: i

        IF (linfos(3)>=2) WRITE(*,*) ' ... user_init (REINJECTION)'
!     avoid side effects
        DO i = 1, nbc_data
!        init from the first set
          dbc_dataold(i) = dbc_data(i,1,1)
        END DO
        RETURN
      END

!>    @brief reinjection "dummy" routine, no "reinjection" functionality
!>    @param[in] ismpl local sample index
      SUBROUTINE calc_user(ismpl)
        use arrays
        use mod_genrl
        use mod_linfos
        IMPLICIT NONE
        integer :: ismpl
        integer :: i
        INTRINSIC int

        IF (linfos(3)>=2) WRITE(*,*) ' ... calc_user'
!     Dummy body
        i = int(head(1,1,1,ismpl))
        RETURN
      END
