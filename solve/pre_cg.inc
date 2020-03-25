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

!
!     Constant : at the beginning of a new system
      integer START
!
!     Constant : Matrix-Vector-Product [v]:=[A]x[s]
      integer MVP
!
!     Constant : two Matrix-Vector-Product [v]:=[A]x[s] & [z]:=[A]x[x]
      integer MVPx
!
!     Constant : nothing to do anymore ...
      integer FINE
!
!     Constant : must break
      integer ABORT
!
      parameter (START = 0)
      parameter (MVP   = 1)
      parameter (MVPx  = 2)
      parameter (FINE  = 3)
      parameter (ABORT = 4)

!
!     offsets for temporary vectors saved in [locTMP]
!     ld_lv    : leading dimensions (lokal number of elements)
!     locTMP(1,1) ^= [b_hat]
!     locTMP(1,2) ^= [r]
!     locTMP(1,3) ^= [z]
!     locTMP(1,4) ^= [s]
!     locTMP(1,5) ^= [p]
!     locTMP(1,6) ^= [q]
!     locTMP(1,7) ^= [v]
!     locTMP(1,8) ^= [mt] (only for EST)
!     locTMP(1,9) ^= [best] X value
      integer b_hat,r,z,s,p,q,v,mt,best
      parameter(b_hat=1)
      parameter(r=2)
      parameter(z=3)
      parameter(s=4)
      parameter(p=5)
      parameter(q=6)
      parameter(v=7)
      parameter(mt=8)
      parameter(best=9)
