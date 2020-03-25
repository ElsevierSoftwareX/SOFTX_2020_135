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

!>    @brief calculate parameter weighting matrix wp^2
!>    @param[in] i matrix index, first dimension
!>    @param[in] j matrix index, second dimension
!>    @param[in] ismpl local sample index
!>    @return weighting value
!>    @details
!>in this case wp^2  = Cpp^-1 is assumed to\n
!>be diagonal and the inverse parameter covariance \n
!>apriori in bayesian estimation.\n
!>    C_pp == I .* d?unit(:)\n
!>    C_pp^-1 == I ./ d?unit(:) = wp2_ps =wp2_ps\n
!>               [g_p_^2]\n
      DOUBLE PRECISION FUNCTION wp2_ps(i,j,ismpl)
        use arrays
        use mod_inverse
        use mod_genrl
        use mod_data
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j
        INTEGER s_k, s_u
        DOUBLE PRECISION wgt, get_optid
        EXTERNAL get_optid

        IF (para_weight==1) THEN
          wp2_ps = covar_prior_p(i,j)
        ELSE IF (i==j) THEN
          s_k = seed_para(1,i)
          s_u = seed_para(2,i)
!         wgt = 0.0d0
!         do k = 1, ndata
!            wgt = wgt +jac(k,i)**2
!         enddo
          wgt = 1.0D0
          IF ((s_k<=lastidx) .AND. (s_k>=1)) THEN
            wgt = propwgt(s_k-firstidx+1)
          ELSE IF ((s_k<=bc_lastidx) .AND. (s_k>=1)) THEN
            wgt = bcwgt(s_k-bc_firstidx+1)
          ELSE IF (s_k==-1) THEN
            wgt = tpwgt(opti_tp(1,s_u),opti_tp(2,s_u),opti_tp(3,s_u))
          END IF
          wp2_ps = wgt/(get_optid(i,ismpl)*get_optid(i,ismpl))
        ELSE
          wp2_ps = 0.0D0
        END IF
!
        RETURN
      END

!>    @brief calculate data weighing matrix wd^2
!>    @param[in] i matrix index, first dimension
!>    @param[in] j matrix index, second dimension
!>    @param[in] ismpl local sample index
!>    @return weighting value
!>    @details
!>in this case wd^2  = Cdd^-1 is assumed to\n
!>be diagonal and the inverse data covariance \n
!>apriori in bayesian estimation.\n
!>    C_dd == I .* data(:,2)\n
!>    C_dd^-1 == I ./ data(:,2) = Cddm1 = wd2_ps\n
!>               [ndata^2]\n
      DOUBLE PRECISION FUNCTION wd2_ps(i,j,ismpl)
        use arrays
        use mod_inverse
        use mod_data
        IMPLICIT NONE
        INTEGER i, j, k, ismpl
        DOUBLE PRECISION wgt

        IF (data_weight==1) THEN
          wd2_ps = covar_prior_d(i,j)
        ELSE IF (i==j) THEN
!        wgt=0.0d0
!        do k=1,mpara
!          wgt = wgt +jac(i,k)**2
!        enddo
          wgt = 1.0D0
!    wgt = wdt(type(i))
          wd2_ps = 1.0D0/(wgt*ddata(i,cdd_w)*ddata(i,cdd_w))
        ELSE
          wd2_ps = 0.0D0
        END IF
!
        RETURN
      END
