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

!>    @brief test for a specific runtime option
!>    @param[in] soption expected runtime option (string)
!>    @return existence of the expected runtime option
      LOGICAL FUNCTION test_option(soption)
        IMPLICIT NONE
        INTEGER i
        character (len=80) :: param
        character (len=*) :: soption
        INTEGER p

        test_option = .FALSE.
        i = command_argument_count()
        DO p = 1, i
          CALL get_command_argument(p,param)
          IF (param==soption) THEN
            test_option = .TRUE.
            RETURN
          END IF
        END DO
        RETURN
      END

!>    @brief integer value for a specific runtime option
!>    @param[in] soption runtime option (string)
!>    @return integer value for the defined runtime option
      INTEGER FUNCTION get_ioptval(soption)
        IMPLICIT NONE
        INTEGER i
        character (len=80) :: param
        character (len=*) :: soption
        INTEGER p, t

        i = command_argument_count()
        DO p = 1, i
          CALL get_command_argument(p,param)
          IF (param==soption) THEN
            IF (i>p) THEN
              t = p + 1
              CALL get_command_argument(t,param)
              READ(param,*) get_ioptval
              RETURN
            ELSE
              WRITE(*,'(3A)') 'error: no option value for "', &
                soption, '"!'
              STOP
            END IF
          END IF
        END DO
        RETURN
      END

!>    @brief returns the string value for a specific runtime option
!>    @param[in] soption runtime option (string)
!>    @param[out] param string value
      SUBROUTINE get_coptval(soption,param)
        IMPLICIT NONE
        INTEGER i
        character (len=*) :: param
        character (len=*) :: soption
        INTEGER p, t

        i = command_argument_count()
        DO p = 1, i
          CALL get_command_argument(p,param)
          IF (param==soption) THEN
            IF (i>p) THEN
              t = p + 1
              CALL get_command_argument(t,param)
              RETURN
            ELSE
              WRITE(*,'(3A)') 'error: no option value for "', &
                soption, '"!'
              STOP
            END IF
          END IF
        END DO
        RETURN
      END
