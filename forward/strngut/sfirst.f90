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

!>    @brief returns position before the first separator character in string (special version)
!>    @param[in] string search string
!>    @return position before the first separator character
!>    @details
!>    special version for reading extended parameters (used in "no_ext_link")\n
      INTEGER FUNCTION sfirst(string)
        IMPLICIT NONE
        character (len=*) :: string
        INTEGER i
        INTRINSIC achar

        sfirst = 0
        DO i = 1, len(string)
          IF (string(i:i)==',' .OR. string(i:i)==';' .OR. &
            string(i:i)==':' .OR. string(i:i)=='=') RETURN
          IF (string(i:i)/=' ' .AND. string(i:i)/=achar(9)) sfirst = i
        END DO
        RETURN
      END
