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

!>    @brief solve of : D*[D+L]^(-1) x [b] = [b^]
!>    @param[in] N_I lengths of I-dimension of local matrix [M]
!>    @param[in] N_J lengths of J-dimension of local matrix [M]
!>    @param[in] N_K lengths of K-dimension of local matrix [M]
!>    @param[in] b vector [b]
!>    @param[in] MA diagonals of matrix [M]
!>    @param[in] MB diagonals of matrix [M]
!>    @param[in] MC diagonals of matrix [M]
!>    @param[in] MD diagonals of matrix [M]
!>    @param[out] b_hat the solution vector [b^]
!>    @details
!>    solve of : D*[D+L]^(-1) x [b] = [b^]\n
!>               with [M]=[L]+[D]+[R]\n
      SUBROUTINE ddl(n_i,n_j,n_k,b,b_hat,ma,mb,mc,md)
        use mod_OMP_TOOLS
        IMPLICIT NONE
        INCLUDE 'OMP_TOOLS.inc'
        INTEGER n_i, n_j, n_k, i, j, k
        DOUBLE PRECISION b(n_i,n_j,n_k), b_hat(n_i,n_j,n_k)
        DOUBLE PRECISION ma(n_i,n_j,n_k), mb(n_i,n_j,n_k)
        DOUBLE PRECISION mc(n_i,n_j,n_k), md(n_i,n_j,n_k)

!     ###################### not parallel Code !!! #########################
!     prepare [b^]
!     [b^] = D*(D+L)^(-1) * [b]
!$OMP master
        DO k = 1, n_k
          DO j = 1, n_j
            DO i = 1, n_i
              b_hat(i,j,k) = b(i,j,k)
              IF (i>1) b_hat(i,j,k) = b_hat(i,j,k) - &
                mc(i,j,k)*b_hat(i-1,j,k)
              IF (j>1) b_hat(i,j,k) = b_hat(i,j,k) - &
                mb(i,j,k)*b_hat(i,j-1,k)
              IF (k>1) b_hat(i,j,k) = b_hat(i,j,k) - &
                ma(i,j,k)*b_hat(i,j,k-1)
!AW           ugly hack ...
              IF (md(i,j,k)/=0.0D0) THEN
                b_hat(i,j,k) = b_hat(i,j,k)/md(i,j,k)
              END IF
            END DO
          END DO
        END DO
        DO k = 1, n_k
          DO j = 1, n_j
            DO i = 1, n_i
              b_hat(i,j,k) = b_hat(i,j,k)*md(i,j,k)
            END DO
          END DO
        END DO
!$OMP end master
!     ################### above not parallel Code !!! ######################
        RETURN
      END

!>    @brief solve of : [D+U]^(-1) x [x^] = [x]
!>    @param[in] N_I lengths of I-dimension of local matrix [M]
!>    @param[in] N_J lengths of J-dimension of local matrix [M]
!>    @param[in] N_K lengths of K-dimension of local matrix [M]
!>    @param[in] x_hat vector [x^]
!>    @param[in] MD diagonals of matrix [M]
!>    @param[in] ME diagonals of matrix [M]
!>    @param[in] MF diagonals of matrix [M]
!>    @param[in] MG diagonals of matrix [M]
!>    @param[out] x the solution vector [x]
      SUBROUTINE du(n_i,n_j,n_k,x_hat,x,md,me,mf,mg)
        IMPLICIT NONE
        INTEGER n_i, n_j, n_k, i, j, k
        DOUBLE PRECISION x(n_i,n_j,n_k), x_hat(n_i,n_j,n_k)
        DOUBLE PRECISION md(n_i,n_j,n_k), me(n_i,n_j,n_k)
        DOUBLE PRECISION mf(n_i,n_j,n_k), mg(n_i,n_j,n_k)

!     ###################### not parallel Code !!! #########################
!     compute [x]
!     [x] = (D+U)^(-1) * [x^]
        DO k = n_k, 1, -1
          DO j = n_j, 1, -1
            DO i = n_i, 1, -1
              x(i,j,k) = x_hat(i,j,k)
              IF (i<n_i) x(i,j,k) = x(i,j,k) - me(i,j,k)*x(i+1,j,k)
              IF (j<n_j) x(i,j,k) = x(i,j,k) - mf(i,j,k)*x(i,j+1,k)
              IF (k<n_k) x(i,j,k) = x(i,j,k) - mg(i,j,k)*x(i,j,k+1)
!AW           ugly hack ...
              IF (md(i,j,k)/=0.0D0) THEN
                x(i,j,k) = x(i,j,k)/md(i,j,k)
              END IF
            END DO
          END DO
        END DO
!     ################### above not parallel Code !!! ######################
        RETURN
      END
