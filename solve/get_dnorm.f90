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

!>    @brief return modified diagonal prepared for normalisation
!>    @param[in] md diagonal (system matrix)
!>    @param[in] mbc_mask boundary mask
!>    @param[in] N number of elements, vector length
!>    @param[out] dnrm normalisation vector
      SUBROUTINE get_dnorm(n,mbc_mask,md,dnrm)
        IMPLICIT NONE
        INTEGER n, i
        DOUBLE PRECISION dnrm(n), dcrit_nrm
!     diagonal (system matrix)
        DOUBLE PRECISION md(n)
!     boundary mask
        CHARACTER mbc_mask(n)
        INTRINSIC dsqrt, dabs, dble
        LOGICAL test_null
        EXTERNAL test_null

!     correction for 2-Norm
        dcrit_nrm = dsqrt(dble(n))

        DO i = 1, n
          dnrm(i) = 0.D0
!        add, since it is no boundary condition element
          IF (mbc_mask(i)=='+') THEN
            IF (test_null(md(i))) THEN
              WRITE(*,*) 'zero at ', i, md(i), mbc_mask(i)
            END IF
!           important: diagonal dominance
!             used for adaptive scaling the residuum
            dnrm(i) = 1.D0/dabs(md(i)*dcrit_nrm)
          END IF
        END DO
!
        RETURN
      END
