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

!>    @brief calculate the mean of a 3D array
!>    @param[in] values 3D array
!>    @param[out] meanv average value
!>    @param[in] ismpl local sample index (ignored here)
      SUBROUTINE mean3d(values,meanv,ismpl)
        use arrays
        use mod_genrl
        use mod_linfos
        IMPLICIT NONE
        DOUBLE PRECISION meanv, values(i0,j0,k0)
        integer :: i, j, k
        integer :: ismpl
        INTRINSIC dble

        IF (linfos(3)>=2) WRITE(*,*) ' ... mean 3D working'
!
        meanv = 0.0d0
        DO k = 1, K0
          DO j = 1, J0
            DO i = 1, I0
              meanv = meanv + values(i,j,k)
            END DO
          END DO
        END DO
!
        meanv = meanv/dble(i0*j0*k0)
!
        RETURN
      END
