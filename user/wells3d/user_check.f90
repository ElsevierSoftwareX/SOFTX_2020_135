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

!>    @brief check current user directory choice
!>    @param[in] ismpl local sample index
      SUBROUTINE user_check(ismpl)
        use mod_genrlc
        IMPLICIT NONE
        INTEGER ismpl
        character (len=20) :: ldef_user
        PARAMETER (ldef_user='wells3d')
        LOGICAL test_option
        EXTERNAL test_option
        INTRINSIC trim

#ifndef USER_wells3d
        WRITE(*,'(3A)') 'error: this source was written for USER=', &
          ldef_user, ', please correct this check in "user_check.f"!'
        STOP
#endif
        IF ( .NOT. test_option('USER='//trim(def_user))) THEN
          IF (ldef_user/=def_user) THEN
            WRITE(*,'(7A)') 'error: model file needs an executable', &
              ' build from USER=', trim(def_user), &
              ', but the current', ' consist of USER=', &
              trim(ldef_user), '!'
            STOP
          END IF
        END IF
        RETURN
      END
