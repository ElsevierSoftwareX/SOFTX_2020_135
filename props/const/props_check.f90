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

!> @brief check current PROPS choice
!> @param[in] ismpl local sample index
!> @details
!> Check the local/current PROPS ldef_props against the PROPS choice
!> in the input file (def_props).
      subroutine props_check(ismpl)
        use mod_genrlc

        implicit none

        ! Sample index
        integer :: ismpl

        ! Local PROPS
        character (len=10), parameter :: ldef_props = 'const'

        ! Test options of command line input
        logical, external :: test_option

        intrinsic trim


#ifndef PROPS_const
        write(*,'(3A)') 'error: this source was written for PROPS=', &
          ldef_props, &
          ', please correct this check in "props_check.f"!'
        stop
#endif
        if ( .not. test_option('PROPS='//trim(def_props))) then
          if (ldef_props/=def_props) then
            write(*,'(7A)') 'Error: model file needs an executable', &
              ' build from PROPS=', trim(def_props), &
              ', but the current', ' consist of PROPS=', &
              trim(ldef_props), '!'
            stop
          end if
        end if

        return

      end subroutine props_check
