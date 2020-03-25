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

!>    @brief return normalisation-vector for initial normalisation of the linear-system
!>    @param[in] mw right side, solution, diagonal (system matrix)
!>    @param[in] mbc_mask boundary mask
!>    @param[in] N number of elements, vector length
!>    @param[in] md diagonal of the system matrx
!>    @param[in] mx currentr solution vector (used for step by step updates)
!>    @param[out] dnormalise updated (new) normalisation vector
!>    @param[in] soddmx second offset (dominant) diagonal*mx
      SUBROUTINE get_norm2(n,mw,mbc_mask,md,mx,dnormalise,soddmx)
        IMPLICIT NONE
        INTEGER n, i
        DOUBLE PRECISION dnormalise(n), dcrit_nrm
!     right side, solution, diagonal (system matrix)
        DOUBLE PRECISION mw(n), mx(n), md(n), soddmx(n)
!     boundary mask
        CHARACTER mbc_mask(n)
        INTRINSIC dsqrt, dabs, dble
        LOGICAL test_null
        EXTERNAL test_null

!     correction for 2-Norm
        dcrit_nrm = dsqrt(dble(n))

        DO i = 1, n
          dnormalise(i) = 1.D0
!        add, since it is no boundary condition element
          IF (mbc_mask(i)=='+') THEN
!           important: diagonal dominance
            IF ( .NOT. test_null(mx(i))) THEN
              dnormalise(i) = 1.D0 /(max(dabs(md(i)*mx(i)), soddmx(i))*dcrit_nrm)
            ELSE
              dnormalise(i) = 1.D0 /(max(dabs(md(i)), dabs(mw(i)), soddmx(i))*dcrit_nrm)
            END IF
          END IF
        END DO
!
        RETURN
      END
