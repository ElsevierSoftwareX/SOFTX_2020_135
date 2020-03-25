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

!>    @brief info, verbosity level
module mod_linfos

  !> @brief Array of verbosity level switches.
  !> @details
  !> Array of verbosity level switches. \n\n
  !> linfos(1)=\n
  !> - `0`, only necessary outputs\n
  !> - `1`, time steps information\n
  !> - `2`, with initialisation and file reading outputs\n
  !>
  !> linfos(2)=\n
  !> - `0`, 'inversion/simulation' with only necessary outputs\n
  !> - `1`, 'inversion/simulation' with some more outputs\n
  !> - `2`, 'inversion/simulation' with all outputs (optimisation outputs)\n
  !> - `3`, 'inversion/simulation' with all debug outputs\n
  !>
  !> linfos(3)=\n
  !> - `-1`, 'non linear solver' disables some error outputs\n
  !> - `0`, 'non linear solver' with only necessary outputs\n
  !> - `1`, 'non linear solver' with some stage outputs\n
  !> - `2`, 'non linear solver' with all outputs\n
  !>
  !> linfos(4)=\n
  !> - `0`, 'system solver' without any outputs\n
  !> - `1`, 'system solver' with some outputs\n
  !> - `2`, 'system solver' with some debug outputs\n
  !> - `3`, 'system solver' with all debug outputs\n
  integer, dimension (4) :: linfos

end module mod_linfos
