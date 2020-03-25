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

!>    @brief returns positions of first and last non-blank character
!>    @param[in] text string line
!>    @param[out] ianf first non-blank character
!>    @param[out] iend last non-blank character
!>    @details
!> returns positions of first and last non-blank character in string\n
      SUBROUTINE chln(text,ianf,iend)
        IMPLICIT NONE
        INTEGER i, j, l, ianf, iend
        character (len=*) :: text
        LOGICAL lead, trail

        lead = .TRUE.
        trail = .TRUE.
        l = len(text)
!
        DO i = 1, l
          IF (lead .AND. text(i:i)/=' ') THEN
!           first non-blank CharaCter: store ianf
            ianf = i
            lead = .FALSE.
          END IF
!
          j = l - i + 1
          IF (trail .AND. text(j:j)/=' ') THEN
!           last non-blank CharaCter: store iend
            iend = j
            trail = .FALSE.
          END IF
          IF (( .NOT. lead) .AND. ( .NOT. trail)) RETURN
        END DO
!
!     error return, all string blank
        ianf = 0
        iend = 0
!
        RETURN
      END
