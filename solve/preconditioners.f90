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

!>    @brief preconditioner "dummy" routine, no preconditioning (copy instead)
!>    @param[in] N length of all vectors
!>    @param[in] r vector [r]
!>    @param[out] z vector [z]
!>    @details
!>    compute [z], the solution of [K]x[z]=[r]\n
!>    [K] is substitude from [A] and work as a preconditioner\n
!>    here simple example : only copy [z]:=[r]\n
      SUBROUTINE myprco(n,r,z)
        use mod_OMP_TOOLS
        IMPLICIT NONE
        INCLUDE 'OMP_TOOLS.inc'
!     N : length of all vectors r,z
        INTEGER n
!     vectors [r], [z]
        DOUBLE PRECISION r(*), z(*)
!     thread-stuff
        INTEGER tpos, tanz

        CALL omp_part(n,tpos,tanz)
!     very simple Pre-Cond.
        CALL dcopy(tanz,r(tpos),1,z(tpos),1)
        RETURN
      END

!>    @brief simple diagonal preconditioner
!>    @param[in] N length of all vectors
!>    @param[in] r vector [r]
!>    @param[in] D diagonal vector for preconditioning
!>    @param[out] z vector [z]
!>    @details
!>    compute [z], the solution of [K]x[z]=[r]\n
!>    [K] is substitude from [A] and work as a preconditioner\n
!>    here simple example : [z]:=[r]/[d]\n
      SUBROUTINE diagprco(n,d,r,z)
        use mod_OMP_TOOLS
        IMPLICIT NONE
        INCLUDE 'OMP_TOOLS.inc'
!     N : length of all vectors r,z
        INTEGER n, i
!     vectors [r], [z]
        DOUBLE PRECISION r(n), z(n), d(n)
!     thread-stuff
        INTEGER tpos, tanz

        CALL omp_part(n,tpos,tanz)
!     works if there are not NULL-elements in [D]
!     very simple Pre-Cond.
        DO i = tpos, tpos + tanz - 1
!AW        z(i)=r(i)/D(i)
          IF (d(i)/=0.0D0) z(i) = r(i)/d(i)
        END DO
        RETURN
      END
