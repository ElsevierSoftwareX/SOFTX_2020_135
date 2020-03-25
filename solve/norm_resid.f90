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

!>    @brief adaptive update normalisation of the residuum
!>    @param[in] N vector length
!>    @param[in] dnrm main part of the normalisation vector
!>    @param[in] xval current result vector, needed for normalisation
!>    @param[out] rval update vector
      SUBROUTINE norm_resid(n,dnrm,xval,rval)
        IMPLICIT NONE
        INTEGER n, i
        DOUBLE PRECISION dtmp, dnrm(n), xval(n), rval(n)
        LOGICAL test_null
        EXTERNAL test_null
        INTRINSIC dabs

        DO i = 1, n
          IF ( .NOT. test_null(dnrm(i))) THEN
            IF ( .NOT. test_null(xval(i))) THEN
              dtmp = dnrm(i)/dabs(xval(i))
              rval(i) = rval(i)*dtmp
            END IF
          END IF
        END DO
!
        RETURN
      END
