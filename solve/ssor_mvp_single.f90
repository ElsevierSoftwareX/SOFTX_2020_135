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

!>    @brief apply 7point-star matrix multiply [as]:=[M]x[s]
!>    @param[in] N_I lengths of I dimension of local matrix [M]
!>    @param[in] N_J lengths of J dimension of local matrix [M]
!>    @param[in] N_K lengths of K dimension of local matrix [M]
!>    @param[in] s vector [s]
!>    @param[out] as vector [as]
!>    @param[out] t temporary vector, used for computation
!>    @param[in] MA 1. diagonal of the system matrix [M]
!>    @param[in] MB 2. diagonal of the system matrix [M]
!>    @param[in] MC 3. diagonal of the system matrix [M]
!>    @param[in] MD 4. (main) diagonal of the system matrix [M]
!>    @param[in] ME 5. diagonal of the system matrix [M]
!>    @param[in] MF 6. diagonal of the system matrix [M]
!>    @param[in] MG 7. diagonal of the system matrix [M]
!>    @details
!>    OpenMP parallelised, general version - no special blocking\n
!>    apply 7point-star matrix multiply\n
!>    compute [as]:=[M]x[s], [s],[as],[M] given in 3-D-structure\n
!>    Data-Cube :\n
!> @image html cube.png
!       k     * * * *
!     /     *     * *
!    0 -j * * * *   *
!    |    *     *   *
!    i    *     * *
!         * * * *
      SUBROUTINE ssor_mvp_single(n_i,n_j,n_k,s,as,t,ma,mb,mc,md,me,mf, &
          mg)
        use mod_OMP_TOOLS
        IMPLICIT NONE
        INCLUDE 'OMP_TOOLS.inc'
        INTEGER n_i, n_j, n_k
!
        DOUBLE PRECISION s(n_i,n_j,n_k), as(n_i,n_j,n_k)
        DOUBLE PRECISION ma(n_i,n_j,n_k), mb(n_i,n_j,n_k), &
          mc(n_i,n_j,n_k)
        DOUBLE PRECISION md(n_i,n_j,n_k), me(n_i,n_j,n_k), &
          mf(n_i,n_j,n_k)
        DOUBLE PRECISION mg(n_i,n_j,n_k)
!     t : tmp. vector
        DOUBLE PRECISION t(n_i,n_j,n_k)
!     loop variable
        INTEGER i, j, k
!
        INTEGER panz, ppos
        EXTERNAL panz, ppos
!     THREAD-stuff
        ! INTEGER tpos, tanz
!
        INTEGER ijk


!       compute boundary size
        ijk = n_i*n_j*n_k


! #####################################################################
! Main
! given : A=L+D+U

! t:=(D+U)^(-1)*s

!$OMP master
        CALL dcopy(ijk,s,1,t,1)


        DO k = n_k, 1, -1
          DO j = n_j, 1, -1
            DO i = n_i, 1, -1
              IF (i<n_i) t(i,j,k) = t(i,j,k) - me(i,j,k)*t(i+1,j,k)
              IF (j<n_j) t(i,j,k) = t(i,j,k) - mf(i,j,k)*t(i,j+1,k)
              IF (k<n_k) t(i,j,k) = t(i,j,k) - mg(i,j,k)*t(i,j,k+1)
!AW         ugly hack ...
              IF (md(i,j,k)/=0.0D0) THEN
                t(i,j,k) = t(i,j,k)/md(i,j,k)
              END IF
            END DO
          END DO
        END DO



!CC t2 : tot, t3 = as
! t2 =(s - D*t)
! t3 =(D+L)^(-1) * t2

        CALL dcopy(ijk,s,1,as,1)

!AW   call daxpy(ijk,-1.0d0,t,1,as,1)
        DO k = 1, n_k
          DO j = 1, n_j
            DO i = 1, n_i
              as(i,j,k) = as(i,j,k) - t(i,j,k)*md(i,j,k)
            END DO
          END DO
        END DO


        DO k = 1, n_k
          DO j = 1, n_j
            DO i = 1, n_i
              IF (i>1) as(i,j,k) = as(i,j,k) - mc(i,j,k)*as(i-1,j,k)
              IF (j>1) as(i,j,k) = as(i,j,k) - mb(i,j,k)*as(i,j-1,k)
              IF (k>1) as(i,j,k) = as(i,j,k) - ma(i,j,k)*as(i,j,k-1)
!AW         ugly hack ...
              IF (md(i,j,k)/=0.0D0) THEN
                as(i,j,k) = as(i,j,k)/md(i,j,k)
              END IF
            END DO
          END DO
        END DO


!CC t3 = as, t4 : tot
! t4 = t + t3
! A^T*s = D * t4

        CALL daxpy(ijk,1.0D0,t,1,as,1)

        DO k = 1, n_k
          DO j = 1, n_j
            DO i = 1, n_i
              as(i,j,k) = as(i,j,k)*md(i,j,k)
            END DO
          END DO
        END DO
!$OMP end master
!$OMP barrier


        RETURN
      END
