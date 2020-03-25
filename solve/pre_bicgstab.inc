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
!     K=L_K*R_K : preconditioner
!
!     Constant : aborting, ([r0_hat]^T*[r])=0
      integer ABORT
!
!     Constant : at the beginning of a new system
      integer START
!
!     Constant : preconditioner [y]:=[K^-1]x[p],
!                matrix-vector product [v]:=[A]x[y]
      integer Do_y_p_v
!
!     Constant : preconditioner [z]:=[K^-1]x[s],
!                preconditioner [s_pc]:=[L_K^-1]x[s],
!                matrix-vector product [t]:=[A]x[z],
!                preconditioner [t_pc]:=[L_K^-1]x[t]
      integer Do_z_s_t
!
!     Constant : preconditioner [y]:=[K^-1]x[p],
!                matrix-vector product [v]:=[A]x[y],
!                matrix-vector product [z]:=[A]x[x] for advanced precision
      integer more_y_p_v
!
!     Constant : nothing to do anymore ...
      integer FINE
!
      parameter (ABORT      =-1)
      parameter (START      = 0)
      parameter (Do_y_p_v   = 1)
      parameter (Do_z_s_t   = 2)
      parameter (more_y_p_v = 3)
      parameter (FINE       = 4)

!
!     offsets for temporary vectors saved in [locTMP]
!     ld_lv    : leading dimension (local number of elements)
!     locTMP(1,1) ^= [b_hat] (only for EST)
!     locTMP(1,2) ^= [r0_tilde] (only for EST)
!     locTMP(1,3) ^= [r]
!     locTMP(1,4) ^= [z]
!     locTMP(1,5) ^= [s]
!     locTMP(1,6) ^= [t]
!     locTMP(1,7) ^= [v]
!     locTMP(1,8) ^= [p]
!     locTMP(1,9) ^= [y]
!     locTMP(1,10) ^= [t_pc]
!     locTMP(1,11) ^= [s_pc]
!     locTMP(1,12)^= [mt] (only for EST)
!     locTMP(1,13)^= [best] X value
      integer b_hat,r0_tilde,r,z,s,t,v,p,y,t_pc,s_pc,mt,best
      parameter(b_hat=1)
      parameter(r0_tilde=2)
      parameter(r=3)
      parameter(z=4)
      parameter(s=5)
      parameter(t=6)
      parameter(v=7)
      parameter(p=8)
      parameter(y=9)
      parameter(t_pc=10)
      parameter(s_pc=11)
      parameter(mt=12)
      parameter(best=13)
