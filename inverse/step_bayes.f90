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

!>    @brief calculate one bayes step from p_k to p_k+1
!>    @param[in] p_k parameter (iteration k)
!>    @param[in] p_apr apriori parameter value
!>    @param[out] p_delta delta to the estimated parameters
!>    @param[in,out] tmp temporary vectors (unspecified values)
!>    @param[in,out] tmpM temporary martix (unspecified values)
!>    @param[in] ismpl local sample index
!>    @details
!>calculate one bayes step from p_k to p_k+1, \n
!>given jacobian and parameter set p_k. \n
!>corresponding to the definition of Tarantola (2004) p.79\n
!> -   p: parameter array\n
!> -   c: common block array\n
!> -   f: external function\n
!> -   i: compute inside this function (temporary variable)\n
!> -   p_k+1 = p_k +(S^T *C_dd^-1 *S + reg_func() *C_pp^-1)^-1 \n
!>      *(S^T *C_dd^-1 *(d -g(p_k)) - reg_func() *C_pp^-1 *(p_k -p_apr))\n
!> -   p_k+1 = p_k +\n
!>        (S^T *C_dd^-1 *S + reg_func() *C_pp^-1)^-1 *\n
!>        (S^T *C_dd^-1 *(d -g(p_k)) -\n
!>        reg_func() *C_pp^-1 *(p_k -p_apr))\n
!> -   [mpara] = [mpara] +\n
!>         ([mpara x ndata] *[ndata^2] *[ndata x mpara] +[1]*[mpara^2]) *\n
!>         ([mpara x ndata] *[ndata^2] *([ndata] -[ndata]) -\n
!>         [1]*[mpara^2] *([mpara] -[mpara]))\n
!>
!> ==>\n
!>   StWd2 = S^T *C_dd^-1\n
!>
!>   p_k+1 = p_k +(StWd2 *S + reg_func() *C_pp^-1)^-1 \n
!>        *(StWd2 *(d -g(p_k)) - reg_func() *C_pp^-1 *(p_k -p_apr))\n
!> ==>\n
!> - p:  forw == g(p_k)\n
!>     [ndata]\n
!> - p:  jacobi == S\n
!>     [ndata x mpara]\n
!> - c:  d == data(:,1)\n
!>     [ndata]\n
!> - f:  C_pp == I .* d?unit(:); C_pp^-1 == I ./ d?unit(:) = Wp2\n
!>                            [mpara^2]\n
!> - f:  C_dd == I .* data(:,2); C_dd^-1 == I ./ data(:,2) = Wd2\n
!>                            [ndata^2]\n
!> - i:  tmp1 = StWd2 = jacobi^T *Wd2\n
!>       [mpara x ndata]\n
!>
!>p_k+1 = p_k +(tmp1 *jacobi + reg_func() *Wp2)^-1 \n
!>           *(tmp1 *(d -forw) - reg_func() *Wp2 *(p_k -p_apr))\n
!> ==>\n
!> - i:  tmp2 = p_k -p_apr\n
!>       [mpara]\n
!> - i:  tmp3 = d -forw\n
!>       [ndata]\n
!> - i:  tmp4 = tmp1 *jacobi\n
!>       [mpara^2]\n
!>
!>    p_k+1 = p_k +(tmp4 + reg_func() *Wp2)^-1 *(tmp1 *tmp3 - reg_func() *Wp2 *tmp2)\n
!> ==>\n
!> - i:  tmp5 = tmp1 *tmp3 - reg_func() *Wp2 *tmp2\n
!>       [mpara]\n
!> - i:  tmpM = tmp4 + reg_func() *Wp2\n
!>       [mpara^2]\n
!>
!>    p_k+1 = p_k +tmpM^-1 *tmp5\n
!> ==>\n
!> - i:  tmp6 = tmpM^-1 *tmp5\n
!>       [mpara]\n
!> - i:  p_k+1 = p_k +tmp6\n
!>        [mpara]\n
      SUBROUTINE step_bayes(p_k,p_apr,tmp,tmpm,p_delta,ismpl)
        use arrays
        use mod_genrl
        use mod_inverse
        use mod_data
        use mod_linfos
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k
!
        DOUBLE PRECISION p_k(mpara), p_apr(mpara), p_delta(mpara)
!     not dense !!!
        DOUBLE PRECISION tmp(mpara,max(mpara,ndata),5), &
          tmpm(mpara,mpara)
!
        DOUBLE PRECISION wp2_ps, wd2_ps, ddot, reg_func
!                                [mpara^2], [ndata^2]
        EXTERNAL wp2_ps, wd2_ps, ddot, reg_func

!
        IF (linfos(2)>=1) WRITE(*,*) ' ... calculate bayes step'
!
! i:  tmp1 = StWd2 = jacobi^T *Wd2
!        [mpara x ndata]
        DO j = 1, ndata
          DO i = 1, mpara
            tmp(i,j,1) = 0.0D0
            DO k = 1, ndata
              tmp(i,j,1) = tmp(i,j,1) + jac(k,i)*wd2_ps(k,j &
                ,ismpl)
            END DO
          END DO
        END DO

! i:  tmp2 = p_k -p_apr
!        [mpara]
        DO i = 1, mpara
          tmp(i,1,2) = p_k(i) - p_apr(i)
        END DO

! i:  tmp3 = d -forw
!        [ndata]
        DO j = 1, ndata
          tmp(1,j,3) = ddata(j,cdd_pv) - sdata(j,ismpl)
        END DO

! i:  tmp4 = tmp1 *jacobi
!        [mpara^2]
        DO j = 1, mpara
          DO i = 1, mpara
            tmp(i,j,4) = 0.0D0
            DO k = 1, ndata
              tmp(i,j,4) = tmp(i,j,4) + tmp(i,k,1)*jac(k,j)
            END DO
          END DO
        END DO

! i:  tmp5 = tmp1 *tmp3 -reg_func() *wp2 *tmp2
!        [mpara]
        DO i = 1, mpara
          tmp(i,1,5) = 0.0D0
          DO k = 1, ndata
            tmp(i,1,5) = tmp(i,1,5) + tmp(i,k,1)*tmp(1,k,3)
          END DO
        END DO

        DO i = 1, mpara
          DO k = 1, mpara
            tmp(i,1,5) = tmp(i,1,5) - reg_func()*wp2_ps(i,k,ismpl)*tmp(k &
              ,1,2)
          END DO
        END DO

! i:  tmpM = tmp4 +reg_func() *Wp2
!        [mpara^2]
        DO j = 1, mpara
          DO i = 1, mpara
            tmpm(i,j) = tmp(i,j,4) + reg_func()*wp2_ps(i,j,ismpl)
          END DO
        END DO

! i:  tmp6 = tmpM^-1 *tmp5
!        [mpara]
!     full_solve(A,B,x)
        CALL set_dval(mpara,0.D0,p_delta)
        CALL dense_solve(mpara,1,tmpm(1,1),tmp(1,1,5),p_delta)
!
        IF (linfos(2)>=1) THEN
          WRITE(*,*) ' [I] : nrm2(p_k-p_(k+1)):', dsqrt(ddot(mpara,p_delta,1,p_delta,1))
          WRITE(*,*) ' [I] : end bayes step'
        END IF
!
        RETURN
      END
