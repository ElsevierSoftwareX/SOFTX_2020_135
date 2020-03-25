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

!>    @brief normalise the linear system, matrix and the right side
!>    @param[in] N size of all given vectors
!>    @param[in] mbc_mask boundary mask
!>    @param[in,out] mw right side, solution
!>    @param[in] mx [x]-vector
!>    @param[in,out] MA 1. diagonal of the system matrix
!>    @param[in,out] MB 2. diagonal of the system matrix
!>    @param[in,out] MC 3. diagonal of the system matrix
!>    @param[in,out] MD 4. diagonal of the system matrix
!>    @param[in,out] ME 5. diagonal of the system matrix
!>    @param[in,out] MF 6. diagonal of the system matrix
!>    @param[in,out] MG 7. diagonal of the system matrix
!>    @param[in,out] dnormalise normalising vector for adaptive resiuum (disabled function - working on matnormalise), as input initialised with second offset (dominant) diagonal*x
!>    @param[out] matnormalise normalising vector for matrix system
      SUBROUTINE norm_linsys2(n,mbc_mask,mw,mx,ma,mb,mc,md,me,mf,mg,dnormalise,matnormalise)
        IMPLICIT NONE
        INTEGER n, i
!     normalising value
        DOUBLE PRECISION dnormalise(n), matnormalise(n)
!     matrix coefficients
        DOUBLE PRECISION ma(n), mb(n), mc(n), md(n), me(n), mf(n), &
          mg(n)
!     right side, solution
        DOUBLE PRECISION mw(n), mx(n)
        CHARACTER mbc_mask(n)


!     get the norm.-vector for the normalisation of the system (in "matnormalise")
        CALL get_norm2(n,mw,mbc_mask,md,mx,matnormalise,dnormalise)

!     normalise the system matrix
        DO i = 1, n
          ma(i) = ma(i)*matnormalise(i)
        END DO
        DO i = 1, n
          mb(i) = mb(i)*matnormalise(i)
        END DO
        DO i = 1, n
          mc(i) = mc(i)*matnormalise(i)
        END DO
        DO i = 1, n
          md(i) = md(i)*matnormalise(i)
        END DO
        DO i = 1, n
          me(i) = me(i)*matnormalise(i)
        END DO
        DO i = 1, n
          mf(i) = mf(i)*matnormalise(i)
        END DO
        DO i = 1, n
          mg(i) = mg(i)*matnormalise(i)
        END DO
!     normalise the right side
        DO i = 1, n
          mw(i) = mw(i)*matnormalise(i)
        END DO

!     disable the adaptive normalisation (normalisation vector = 0)
        CALL set_dval(n,0.0d0,dnormalise)

        RETURN
      END
