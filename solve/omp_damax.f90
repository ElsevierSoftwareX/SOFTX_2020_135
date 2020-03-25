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

!>    @brief computes the maximum of vector [x] and broadcast to all, (OpenMP version)
!>    @param[in] N length of vector [x]
!>    @param[in] X vector [x]
!>    @param[in] sh_max openmp-shared help variable
!>    @param[out] MAX_X maximum of [x]
      SUBROUTINE omp_damax(n,x,max_x,sh_max)
        use mod_OMP_TOOLS
        IMPLICIT NONE
        INCLUDE 'OMP_TOOLS.inc'
!     N    : length of vector [x]
!     ind  : index of maximal element
        INTEGER n, i, ind
!     vector [x]
        DOUBLE PRECISION x(n), sh_max, max_x
        INTEGER idamax
        EXTERNAL idamax
        INTRINSIC dabs

!      orphaning feature needed
!$OMP  master
!      clear before compute
        sh_max = 0.0D0
!$OMP  end master
!$OMP  barrier
!
!      very simple (slow) variation
!$OMP  do schedule(static) reduction(max:sh_max)
        DO i = 1, omp_get_num_of_threads()
!        compute local maximum
          ind = idamax(n,x,1)
          IF (ind>=1 .AND. n>=1) sh_max = max(sh_max,dabs(x(ind)))
        END DO
!$OMP  end do
!      barrier here ...
!
        max_x = sh_max
        RETURN
      END

!>    @brief computes the maximum of vector [x], serial (no OpenMP) implementation, see above
!>    @param[in] N length of vector [x]
!>    @param[in] X vector [x]
!>    @param[out] MAX_X maximum of [x]
      SUBROUTINE s_damax(n,x,max_x)
        IMPLICIT NONE
!     N    : length of vector [x]
!     ind  : index of maximal element
        INTEGER n, ind
!     vector [x]
        DOUBLE PRECISION x(n), max_x
        INTEGER idamax
        EXTERNAL idamax
        INTRINSIC dabs

!     very simple (slow) variation
!       compute local maximum
        max_x = 0.0D0
        ind = idamax(n,x,1)
        IF (ind>=1 .AND. n>=1) max_x = dabs(x(ind))
        RETURN
      END
