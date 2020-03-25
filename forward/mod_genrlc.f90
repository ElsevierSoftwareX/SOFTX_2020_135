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

!>    @brief commonly used character global variables
module mod_genrlc
      character (len=80) :: title
      character (len=256) :: project
      character (len=256) :: restart_name, status_log_inv

      !> @brief Filename of status log output file.
      !> @details
      !> Filename of status log output file. \n
      character (len=256) :: status_log
!
!     defines the PROPS and USER choice
      character (len=80) :: def_props, def_user
!
!     binary can be one of: 'forward', 'simul' or 'inverse'
      character (len=7) :: def_binary
!
!     external file control for: simul, enkf, inverse, data
      character (len=256) :: filename_simul, filename_enkf, filename_inverse
      character (len=256) :: filename_data
end module mod_genrlc
