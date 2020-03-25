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

!>    @brief solver constants and prozessor cache & block size information
module mod_blocking_size
!     processor-cache size, cache-block size
      integer cache_size, bl_size
!
      integer block_i, block_j, block_k
      integer bdim_i, bdim_j, bdim_k
      integer max_blocks
!
!     chunk size optimisation (for different cases)
      double precision bldiv_cg
      double precision bldiv_bicg(2)
      double precision bldiv_mvp
      double precision bldiv_dot(3)
!
!     maximal number of temporary vectors needed in a linear system solver
!     1..9: "pre_CG" system solver
!     1..13: "pre_BiCGStab" system solver
      integer max_locTMP
        parameter (max_locTMP = 13)
end module mod_blocking_size
