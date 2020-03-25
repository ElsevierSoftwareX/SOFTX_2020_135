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

!>    @brief calculate the covariances a posteriori
!>    @param[in] tmp temporary space, normal [g_p_,ndata] (now enough space for [g_p_,g_p_] too)
!>    @param[in] g_p_ number of independent parameters
!>    @param[in] ismpl local sample index
!>    @details
!>    calculate the covariances a posteriori\n
!>    and the resolution martrix\n
      SUBROUTINE calc_covares(g_p_,tmp,ismpl)
        use arrays
        use mod_genrl
        use mod_genrlc
        use mod_inverse
        use mod_data
        use mod_linfos
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k

        INTEGER ii, g_p_
!     "g_p_" should be "mpara"
        DOUBLE PRECISION wp2_ps, wd2_ps
        EXTERNAL wp2_ps, wd2_ps

!     normal [g_p_,ndata] : now enough space for [g_p_,g_p_] too
        DOUBLE PRECISION tmp(g_p_,max(ndata,g_p_),3)


        IF (covar==0 .AND. resmat==0) RETURN

        IF (linfos(2)>=1) WRITE(*,*) ' .... calculate covariances'

!     setup normal equation lhs to temporary loCation
! i:  tmp2 = Stwd2_ps = jac^T * wd2_ps
!            [g_p__ x ndata]
        DO j = 1, ndata
          DO i = 1, g_p_
            tmp(i,j,2) = 0.0D0
            DO k = 1, ndata
              tmp(i,j,2) = tmp(i,j,2) + jac(k,i)*wd2_ps(k,j &
                ,ismpl)
            END DO
          END DO
        END DO

! i:  tmp3 = tmp2 *jac
!            [g_p_^2]
        DO j = 1, g_p_
          DO i = 1, g_p_
            tmp(i,j,3) = 0.0D0
            DO k = 1, ndata
              tmp(i,j,3) = tmp(i,j,3) + tmp(i,k,2)*jac(k,j)
            END DO
          END DO
        END DO

!-------------------------------------------------
!     need for covariance & resolution matrix

!     i:  tmp_mat = tmp3 +wp2_ps
!     M = [jac^T * wd2_ps *jac]
!     [g_p_^2]
        DO j = 1, g_p_
          DO i = 1, g_p_
            tmp_mat(i,j) = tmp(i,j,3) + wp2_ps(i,j,ismpl)
          END DO
        END DO

!     prepare right side
        DO j = 1, g_p_
          DO i = 1, g_p_
            tmp(i,j,1) = 0.0D0
            covar_p(i,j) = 0.0D0
          END DO
        END DO

!     Warning: need "g_p_" =< "ndata" !!!
!     [tmp1] = [I]
        DO j = 1, g_p_
          tmp(j,j,1) = 1.0D0
        END DO

!     solve for parameter CovarianCe a posteriori
!     [covar_p] = [J*] = [M^-1] : [g_p_]x[g_p_]
        CALL dense_solve(g_p_,g_p_,tmp_mat,tmp(1,1,1),covar_p)

!     now covar_p is the inverse of the matrix M
!-------------------------------------------------

!     CalCulate errors
        DO i = 1, g_p_
          IF (covar_p(i,i)<0.0D0) WRITE(*,'(A,I6)') &
            'warning: Covar-P(i,i) < zero at i=', i
          CALL set_optie(i,dsqrt(covar_p(i,i)),ismpl)
        END DO

        IF (resmat/=0) THEN
          DO j = 1, g_p_
            DO i = 1, g_p_
              resmat_p(i,j) = 0.0D0
              DO k = 1, g_p_
!                 [resmat_p] = [J*]x[wp2_ps] : [g_p_]x[g_p_]
                resmat_p(i,j) = resmat_p(i,j) + &
                  covar_p(i,k)*wp2_ps(k,j,ismpl)
              END DO
            END DO
!           [resmat_p] = [I] -[J*]x[wp2_ps] = [I] -[M^-1]x[wp2_ps]
!           matlab: CT= eye(size(U))-U*CI; U=A^-1=M^-1=J*; CI=CIp=wp2_ps
            resmat_p(j,j) = 1.0D0 - resmat_p(j,j)
          END DO
        END IF

        RETURN
      END
