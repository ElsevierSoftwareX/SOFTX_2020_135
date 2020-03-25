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

!>    @brief defines the reinjection parameters, user directory switch
module mod_wells3d
!     number of IN sources
      integer num_in
        parameter (num_in = 1)
!     number of OUT productions
      integer num_pro
        parameter (num_pro = 1)
!
!     finish time
      double precision stoptime
! 1d      parameter (stoptime = 86400.d0)
! 19 h!!!
        parameter (stoptime = 68400.d0)
!
      integer iin(num_in),   jin(num_in),   kin(num_in)
      integer ipro(num_pro), jpro(num_pro), kpro(num_pro)
!
        data iin /13/
        data jin /8/
        data kin /13/
        data ipro /13/
        data jpro /18/
        data kpro /13/
end module mod_wells3d
