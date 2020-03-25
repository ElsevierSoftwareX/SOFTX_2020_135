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

!>    @brief calculate one d-space step from p_k to p_k+1
!>    @param[in] p_k parameter (iteration k)
!>    @param[in] p_apr apriori parameter value
!>    @param[out] p_delta delta to the estimated parameters
!>    @param[in,out] tmp temporary vectors (unspecified values)
!>    @param[in,out] tmpM temporary martix (unspecified values)
!>    @param[in] ismpl local sample index
!>    @details
!>calculate one d-space step from p_k to p_k+1,\n
!>given jacobian and parameter set p_k\n
!> -   p: parameter array\n
!> -   c: common block array\n
!> -   f: external function\n
!> -   i: compute inside this function (temporary variable)\n
!> -    regularisation parameter (reg_func()) skipped !
!> -   p_k+1  = p_apr +C_pp *S^T *(C_dd +S *C_pp *S^T)^-1 *\n
!>             (d -g(p_k) +S *(p_k -p_apr))\n
!> -   p_k+1  = p_apr +\n
!>             C_pp *S^T *\n
!>             (C_dd +S *C_pp *S^T)^-1 *\n
!>             (d -g(p_k) +S *(p_k -p_apr))\n
!> -   [mpara] = [mpara] +\n
!>             [mpara^2] *[mpara x ndata] *\n
!>             ([ndata^2] +[ndata x mpara] *[mpara^2] *[mpara x ndata]) *\n
!>             ([ndata] -[ndata] +[ndata x mpara] *([mpara] -[mpara]))\n
!>
!> ==>\n
!>   Wp2St = C_pp *S^T\n
!>
!>   p_k+1  = p_apr +Wp2St *(C_dd +S *Wp2St)^-1 *\n
!>             (d -g(p_k) +S *(p_k -p_apr))\n
!> ==>\n
!> - p:  forw == g(p_k)\n
!>         [ndata]\n
!> - p:  jacobi == S\n
!>         [ndata x mpara]\n
!> - c:  d == data(:,1)\n
!>         [ndata]\n
!> - f:  C_pp == I .* d?unit(:) = wp2_ds\n
!>            [mpara^2]\n
!> - f:  C_dd == I .* data(:,2) = wd2_ds\n
!>            [ndata^2]\n
!> - i:  tmp1 = Wp2St =  wp2_ds *jacobi^T\n
!>           [mpara x ndata]\n
!>
!>    p_k+1  = p_apr +tmp1 *(wd2_ds +jacobi *tmp1)^-1 *\n
!>             (d -forw +jacobi *(p_k -p_apr))\n
!> ==>\n
!> - i:  tmp2 = jacobi *(p_k -p_apr)\n
!>           [ndata]\n
!> - i:  tmp3 = d -forw +tmp2\n
!>           [ndata]\n
!> - i:  tmpM = wd2_ds +jacobi *tmp1\n
!>           [ndata^2]\n
!>
!>    p_k+1 = p_apr +tmp1 *(tmpM^-1 *tmp3)\n
!> ==>ps\n
!> - i:  tmp4 = tmpM^-1 *tmp3
!>           [ndata]\n
!> - i:  p_k+1 = p_apr +tmp1 *tmp4
!>            [mpara]\n
      SUBROUTINE step_d_space(p_k,p_apr,tmp,tmpm,p_delta,ismpl)
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
        DOUBLE PRECISION tmp(mpara,max(mpara,ndata),4)
!     high memory requirement
        DOUBLE PRECISION tmpm(ndata,ndata)

        DOUBLE PRECISION wp2_ds, wd2_ds, ddot
!                      [mpara^2], [ndata^2]
        EXTERNAL wp2_ds, wd2_ds, ddot

        IF (linfos(2)>=1) WRITE(*,*) ' ... calculate d-space step'
!
        IF (ndata<1) THEN
          WRITE(*,'(A)') &
            'error: in "step_d_space", "ndata" equals to zero !'
          STOP
        END IF
!
! i:  tmp1 = Wp2 =  wp *jacobi^T
!            [mpara x ndata]
        DO j = 1, ndata
          DO i = 1, mpara
            tmp(i,j,1) = 0.0D0
            DO k = 1, mpara
              tmp(i,j,1) = tmp(i,j,1) + jac(j,k)*wp2_ds(i,k,ismpl)
            END DO
          END DO
        END DO

! i:  tmp2 = jacobi *(p_k -p_apr) <- (tmp2 = p_k -p_apr)
!            [ndata]           <- ([mpara])
        DO j = 1, ndata
          tmp(1,j,2) = 0.0D0
          DO i = 1, mpara
            tmp(1,j,2) = tmp(1,j,2) + jac(j,i)*(p_k(i)-p_apr(i))
          END DO
        END DO

! i:  tmp3 = d -forw +tmp2 <- (tmp3 = d -forw)
!            [ndata]
        DO j = 1, ndata
          tmp(j,1,3) = ddata(j,cdd_pv) - sdata(j,ismpl) + tmp(1,j,2)
        END DO

! i:  tmpM = wd2_ds +jacobi *tmp1
!            [ndata^2]
        DO j = 1, ndata
          DO i = 1, ndata
            tmpm(i,j) = wd2_ds(i,j,ismpl)
            DO k = 1, mpara
              tmpm(i,j) = tmpm(i,j) + jac(i,k)*tmp(k,j,1)
            END DO
          END DO
        END DO

! i:  tmp4 = tmpM^-1 *tmp3
!            [ndata]
!     full_solve(A,B,x)
        CALL dense_solve(ndata,1,tmpm(1,1),tmp(1,1,3),tmp(1,1,4))

! i:  p_delta = tmp1 *tmp4
!             [mpara]
        DO i = 1, mpara
          p_delta(i) = 0.0D0
          DO k = 1, ndata
            p_delta(i) = p_delta(i) + tmp(i,k,1)*tmp(k,1,4)
          END DO
        END DO

! i:  p_k+1 = p_apr +p_delta
!     -> p_delta! = p_apr +p_delta -p_k
!     -> p_k+1 = p_k +p_delta!
!             [mpara]
        DO i = 1, mpara
          p_delta(i) = p_delta(i) + p_apr(i) - p_k(i)
        END DO
!
        IF (linfos(2)>=1) THEN
          WRITE(*,*) ' [I] : nrm2(p_k-p_(k+1)):', dsqrt(ddot(mpara,p_delta,1,p_delta,1))
          WRITE(*,*) ' [I] : end d-space step'
        END IF
!
        RETURN
      END
