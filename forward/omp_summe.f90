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

!>    @brief OpenMP "reduction" collection used in Courant/Peclet/Neumann computations
!>    @param[in,out] dval_maxx max X value
!>    @param[in,out] dval_minx min X value
!>    @param[in,out] dval_avgx sum, X value
!>    @param[in,out] dval_maxy max Y value
!>    @param[in,out] dval_miny min Y value
!>    @param[in,out] dval_avgy sum, Y value
!>    @param[in,out] dval_maxz max Z value
!>    @param[in,out] dval_minz min Z value
!>    @param[in,out] dval_avgz sum, Z value
!>    @param[in,out] c1 number in X
!>    @param[in,out] c2 number in Y
!>    @param[in,out] c3 number in Z
!>    @param[in] ismpl local sample index
!>    @details
!> computes global (inter) thread maximum, minimum and sum (float and integer)\n
      SUBROUTINE omp_summe(dval_maxx,dval_minx,dval_avgx,dval_maxy, &
          dval_miny,dval_avgy,dval_maxz,dval_minz,dval_avgz,c1,c2,c3, &
          ismpl)
        use arrays
        use mod_genrl
        use mod_OMP_TOOLS
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j
        INCLUDE 'OMP_TOOLS.inc'
        DOUBLE PRECISION dval_maxx, dval_minx, dval_avgx, dval_maxy, &
          dval_miny, dval_avgy, dval_maxz, dval_minz, dval_avgz
        INTEGER c1, c2, c3, myid

!       each thread needs to set its private values in the global arrays
        myid = omp_get_his_thread_num() + 1
!
        omp_dglobal(myid,1,ismpl) = dval_maxx
        omp_dglobal(myid,2,ismpl) = dval_minx
        omp_dglobal(myid,3,ismpl) = dval_avgx
!
        omp_dglobal(myid,4,ismpl) = dval_maxy
        omp_dglobal(myid,5,ismpl) = dval_miny
        omp_dglobal(myid,6,ismpl) = dval_avgy
!
        omp_dglobal(myid,7,ismpl) = dval_maxz
        omp_dglobal(myid,8,ismpl) = dval_minz
        omp_dglobal(myid,9,ismpl) = dval_avgz
!
        omp_iglobal(myid,1,ismpl) = c1
        omp_iglobal(myid,2,ismpl) = c2
        omp_iglobal(myid,3,ismpl) = c3
!
!     global OpenMP reduction
!$OMP barrier
!$OMP master
!     fast OpenMP reduction MAX
        DO i = 2, tlevel_1
          omp_dglobal(1,1,ismpl) = max(omp_dglobal(1,1,ismpl), &
            omp_dglobal(i,1,ismpl))
        END DO
        DO i = 2, tlevel_1
          omp_dglobal(1,4,ismpl) = max(omp_dglobal(1,4,ismpl), &
            omp_dglobal(i,4,ismpl))
        END DO
        DO i = 2, tlevel_1
          omp_dglobal(1,7,ismpl) = max(omp_dglobal(1,7,ismpl), &
            omp_dglobal(i,7,ismpl))
        END DO
!     fast OpenMP reduction MIN
        DO i = 2, tlevel_1
          omp_dglobal(1,2,ismpl) = min(omp_dglobal(1,2,ismpl), &
            omp_dglobal(i,2,ismpl))
        END DO
        DO i = 2, tlevel_1
          omp_dglobal(1,5,ismpl) = min(omp_dglobal(1,5,ismpl), &
            omp_dglobal(i,5,ismpl))
        END DO
        DO i = 2, tlevel_1
          omp_dglobal(1,8,ismpl) = min(omp_dglobal(1,8,ismpl), &
            omp_dglobal(i,8,ismpl))
        END DO
!     fast OpenMP reduction SUM
        DO i = 2, tlevel_1
          omp_dglobal(1,3,ismpl) = omp_dglobal(1,3,ismpl) + &
            omp_dglobal(i,3,ismpl)
        END DO
        DO i = 2, tlevel_1
          omp_dglobal(1,6,ismpl) = omp_dglobal(1,6,ismpl) + &
            omp_dglobal(i,6,ismpl)
        END DO
        DO i = 2, tlevel_1
          omp_dglobal(1,9,ismpl) = omp_dglobal(1,9,ismpl) + &
            omp_dglobal(i,9,ismpl)
        END DO
!     fast OpenMP reduction SUM (integer)
        DO j = 1, 3
          DO i = 2, tlevel_1
            omp_iglobal(1,j,ismpl) = omp_iglobal(1,j,ismpl) + &
              omp_iglobal(i,j,ismpl)
          END DO
        END DO
!$OMP end master
!$OMP barrier
!
!       each thread gets its own copy of the global reduction
        dval_maxx = omp_dglobal(1,1,ismpl)
        dval_minx = omp_dglobal(1,2,ismpl)
        dval_avgx = omp_dglobal(1,3,ismpl)
        dval_maxy = omp_dglobal(1,4,ismpl)
        dval_miny = omp_dglobal(1,5,ismpl)
        dval_avgy = omp_dglobal(1,6,ismpl)
        dval_maxz = omp_dglobal(1,7,ismpl)
        dval_minz = omp_dglobal(1,8,ismpl)
        dval_avgz = omp_dglobal(1,9,ismpl)
        c1 = omp_iglobal(1,1,ismpl)
        c2 = omp_iglobal(1,2,ismpl)
        c3 = omp_iglobal(1,3,ismpl)
!
        RETURN
      END
