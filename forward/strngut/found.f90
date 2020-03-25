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

!>    @brief locate keyword in file 'ifil', begin at position 1
!>    @param[in] ifil file handler of the opened file
!>    @param[in] key section name
!>    @param[out] block_s current/last readed line from file
!>    @param[in] need_it "true": stop when not found (generate an error for important sections)
!>    @return "true" when section (key) found in file
      LOGICAL FUNCTION found(ifil,key,block_s,need_it)
        USE mod_genrl
        IMPLICIT NONE
        INTEGER ifil, iwrd, kb, ke
        character (len=*) :: key, block_s
        LOGICAL need_it

        INTEGER lblank, locstr
        EXTERNAL lblank, locstr

        iwrd = 1
        found = .FALSE.

!     search dataset
        CALL chln(key,kb,ke)

10      CONTINUE
        READ(ifil,'(a)',end=15,err=15) block_s
        found = locstr(block_s,key(kb:ke)) == 1
!     no matching, try next (10)
        IF ( .NOT. found .OR. lblank(block_s)==0 .OR. &
          block_s(1:1)/=key_char) GO TO 10
        RETURN

!     end of file -> reading file again
15      CONTINUE
        REWIND ifil

20      CONTINUE
        READ(ifil,'(a)',end=25,err=25) block_s
        found = locstr(block_s,key(kb:ke)) == 1
!     no matching, try next (20)
        IF ( .NOT. found .OR. lblank(block_s)==0 .OR. &
          block_s(1:1)/=key_char) GO TO 20
        RETURN

!     end of file -> array not in file
25      REWIND ifil
        IF (need_it) THEN
          WRITE(*,*) 'error: reading field "', key(kb:ke), &
            '", not found!'
          STOP
        END IF

        RETURN
      END
