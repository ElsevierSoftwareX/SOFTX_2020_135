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

!>    @brief convert 4 discrete numbers c1,c2,c3,c4 (each 4bit) into the number code ctrl0 (16bit)
!>    @param[out] ctrl number code
!>    @param[in] c1 first discrete value
!>    @param[in] c2 second discrete value
!>    @param[in] c3 third discrete value
!>    @param[in] c4 fourth discrete value
      SUBROUTINE encntrl4(ctrl,c1,c2,c3,c4)
        IMPLICIT NONE
        INTEGER ctrl, c1, c2, c3, c4
!
        ctrl = c1 + 16*c2 + 256*c3 + 4096*c4
!
        RETURN
      END
