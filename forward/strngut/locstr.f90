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

!>    @brief returns the position of str2 in str1 ignores case
!>    @param[in] str1 string where to find str2
!>    @param[in] str2 search string (key)
!>    @return position of str2 in str1
!>    @details
!>    returns 0 if str2 not found in str1\n
      INTEGER FUNCTION locstr(str1,str2)
        IMPLICIT NONE
        character (len=*) :: str1, str2
        INTEGER i, j, capdif
        LOGICAL same

        locstr = 0
        capdif = ichar('a') - ichar('A')
!
        DO i = 1, len(str1) - len(str2) + 1
          same = .TRUE.
          DO j = 1, len(str2)
            same = same .AND. (str1(i+j-1:i+j-1)==str2(j:j) .OR. 'A'<= &
              str2(j:j) .AND. str2(j:j)<='Z' .AND. ichar(str1(i+j-1:i+ &
              j-1))==ichar(str2(j:j))+capdif .OR. 'a'<=str2(j:j) .AND. &
              str2(j:j)<='z' .AND. ichar(str1(i+j-1:i+ &
              j-1))==ichar(str2(j:j))-capdif)
          END DO
          IF (same) THEN
            locstr = i
            RETURN
          END IF
        END DO
!
        RETURN
      END
