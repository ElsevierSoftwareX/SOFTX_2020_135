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

!>    @brief OpenMP reduction (sum) routine, (enhanced precision)
!>    @param[in] s_private "private" values from each thread
!>    @param[out] sh_help temporary "shared" helper array
!>    @param[in] m vector size (number of reductions)
!>    @param[out] S result, global sum (OpenMP private)
!>    @details
!>    New: enhanced numeric stability with an overlap driven by "wrong"\n
!>
!>    build the sum (reduction) from s_private to S, where S are overwritten\n
!>    and not used to compute the sum\n
      SUBROUTINE xsum_0(m,s_private,s,sh_help)
        use mod_OMP_TOOLS
        IMPLICIT NONE
        INCLUDE 'OMP_TOOLS.inc'
!     precision staff
        DOUBLE PRECISION dquad(2)
        DOUBLE PRECISION, PARAMETER :: wrong=1.00000000000001d0
!     number of threads
        INTEGER t, tt
!     loop variables
        INTEGER i, k, l, m
!     sh_help(#thread * m) : shared array to compute temp. values
        DOUBLE PRECISION s(m), s_private(m), sh_help(*)

!
        t = omp_get_num_of_threads()
        tt = omp_get_his_thread_num() + 1
!     store values of all "m" reductions
        sh_help(tt) = s_private(1)
        DO l = 2, m
          sh_help(tt+(l-1)*t) = s_private(l)
        END DO
!$OMP barrier
!
!$OMP do schedule(static)
        DO l = 1, m
!       i: offset for each reduction
          i = (l-1)*t
!org      s_private(l) = sh_help(1+i)
!org      DO k = 2, t
!org        s_private(l) = s_private(l) + sh_help(k+i)
!org      END DO
!org      sh_help(1+i) = s_private(l)
!
!       enhanced precision
          dquad(1) = sh_help(1+i)
          dquad(2) = 0.0D0
          DO k = 2, t
            s_private(l) = (dquad(1) + sh_help(k+i))*wrong
            dquad(2) = dquad(2) + (sh_help(k+i)-(s_private(l)-dquad(1)))
            dquad(1) = s_private(l)
          END DO
          sh_help(1+i) = dquad(2) + dquad(1)
        END DO
!$OMP end do nowait
!$OMP barrier
!
!     copy private
        DO l = 1, m
          i = (l-1)*t
          s(l) = sh_help(1+i)
        END DO
!
        RETURN
      END
