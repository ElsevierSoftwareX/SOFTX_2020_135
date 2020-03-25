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


!>    @brief check current PROPS choice
!>    @param[in] ismpl local sample index
      SUBROUTINE props_check(ismpl)


        use mod_genrlc
        IMPLICIT NONE
        INTEGER ismpl
        character (len=10) :: ldef_props
        PARAMETER (ldef_props='ice')
        LOGICAL test_option
        EXTERNAL test_option
        INTRINSIC trim


#ifndef PROPS_ice
        WRITE(*,'(3A)') 'error: this source was written for PROPS=', &
          ldef_props, &
          ', please correct this check in "props_check.f"!'
        STOP
#endif
        IF ( .NOT. test_option('PROPS='//trim(def_props))) THEN
          IF (ldef_props/=def_props) THEN
            WRITE(*,'(7A)') 'error: model file needs an executable', &
              ' build from PROPS=', trim(def_props), &
              ', but the current', ' consist of PROPS=', &
              trim(ldef_props), '!'
            STOP
          END IF
        END IF

        RETURN
      END
