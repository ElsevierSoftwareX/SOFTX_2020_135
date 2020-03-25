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

!>    @brief returns the position of the last non-blank character in string
!>    @param[in] string search string
!>    @return position of the last non-blank character
!>    @details
!>    (actually there is a standard unix function of the same name\n
!>    that does the exactly same thing. this is provided for\n
!>    compatibility on non-unix systems.)\n
      INTEGER FUNCTION lblank(string)
        IMPLICIT NONE
        character (len=*) :: string
        INTRINSIC achar

        DO lblank = len(string), 1, -1
          IF (string(lblank:lblank)/=' ' .AND. &
            string(lblank:lblank)/=achar(9)) RETURN
        END DO
        lblank = 0
        RETURN
      END
