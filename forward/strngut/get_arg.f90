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

!>    @brief returns the begin and the end of the argument of the key in the line string
!>    @param[in] key search string
!>    @param[in] line current character line (find here the key)
!>    @param[out] first begin index of the argument
!>    @param[out] last end index of the argument
      SUBROUTINE get_arg(key,line,first,last)
        IMPLICIT NONE
        character (len=*) :: key, line
        INTEGER first, last, locstr, cfirst
        EXTERNAL locstr, cfirst
        INTRINSIC achar

        last = 0
        first = locstr(line,key)
        IF (first==0) RETURN
!
        first = first + len(key)
        first = first - 1 + cfirst(line(first:len(line)))
!
        DO last = first, len(line) - 1
          IF (line(last+1:last+1)==' ' .OR. line(last+1:last+1)==',' &
            .OR. line(last+1:last+1)==';' .OR. &
            line(last+1:last+1)==':' .OR. line(last+1:last+1)=='=' &
            .OR. line(last+1:last+1)==achar(9)) RETURN
        END DO
        last = len(line)
!
        RETURN
      END
