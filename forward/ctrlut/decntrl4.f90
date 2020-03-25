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

!>    @brief convert the number code ctrl0 (16bit) into 4 discrete numbers c1,c2,c3,c4 (each 4bit)
!>    @param[in] ctrl0 number code
!>    @param[out] c1 first discrete value
!>    @param[out] c2 second discrete value
!>    @param[out] c3 third discrete value
!>    @param[out] c4 fourth discrete value
!>    @details
!>    decode [ctrl]\n
!>    ctrl = c1 + 16*c2 + 256*c3 + 4096*c4\n
      SUBROUTINE decntrl4(ctrl0,c1,c2,c3,c4)
        IMPLICIT NONE
        INTEGER ctrl, ctrl0, c1, c2, c3, c4
!
        ctrl = ctrl0
        c1 = mod(ctrl,16)
!
        ctrl = ctrl/16
        c2 = mod(ctrl,16)
!
        ctrl = ctrl/16
        c3 = mod(ctrl,16)
!
        ctrl = ctrl/16
        c4 = mod(ctrl,16)
!
        RETURN
      END
