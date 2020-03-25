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
!>
      SUBROUTINE step_bayes_mf(p_k,p_apr,tmp,p_delta,ismpl)
        use arrays
        use mod_genrl
        use mod_inverse
        use mod_data
        use mod_flow
        use mod_temp
        use mod_conc
        use mod_OMP_TOOLS
        use mod_linfos
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, k
        INTRINSIC max, min, dabs, dsqrt
        DOUBLE PRECISION p_k(mpara), p_apr(mpara), p_delta(mpara)
        DOUBLE PRECISION tmp(max(mpara,ndata),4)
        DOUBLE PRECISION wp2_ps, ddot, reg_func
        DOUBLE PRECISION depsilon, dmaxval
        EXTERNAL wp2_ps, ddot, reg_func

!
        IF (linfos(2)>=1) WRITE(*,*) ' ... calculate (matrix free) bayes step'
!
!$OMP parallel shared(ismpl) num_threads(Tlevel_0)
!         sample init, each sample only one times
          CALL single_init(ismpl)
!$OMP end parallel


!      (S^T *C_dd^-1 *(d -g(p_k)) - reg_func() *C_pp^-1 *(p_k -p_apr))

! i:  tmp3 = d -forw
!        [ndata]
        CALL dcopy(ndata,ddata(1,cdd_pv),1,tmp(1,3),1)
        CALL daxpy(ndata,-1.0d0,sdata(1,ismpl),1,tmp(1,3),1)
!
! i:  tmp4 = StWd2 = jacobi^T *Wd2 *tmp3
!        [mpara x ndata]
        CALL step_bayes_tmp1(tmp(1,3),tmp(1,1),tmp(1,4),ismpl)
!
! i:  tmp2 = p_k -p_apr
!        [mpara]
        CALL dcopy(mpara,p_k,1,tmp(1,2),1)
        CALL daxpy(mpara,-1.0d0,p_apr,1,tmp(1,2),1)
!
! i:  tmp4 = tmp4 -reg_func() *wp2 *tmp2
!        [mpara]
        DO i = 1, mpara
          DO k = 1, mpara
            tmp(i,4) = tmp(i,4) - reg_func()*wp2_ps(i,k,ismpl)*tmp(k,2)
          END DO
        END DO

! i:  p_delta = tmpM^-1 *tmp4
!        [mpara]
        CALL set_dval(mpara,0.D0,p_delta)
!       --
!       use exstimated relative precision and computes an absolute estimation
!         best precision of all physical values / state variables
          depsilon = 1.0d-16
          IF (head_active) THEN
            CALL s_damax(I0*J0*K0,head,dmaxval)
            depsilon = max(depsilon,nltolf/min(1.0d99,dmaxval))
          ENDIF
          IF (pres_active) THEN
            CALL s_damax(I0*J0*K0,pres,dmaxval)
            depsilon = max(depsilon,nltolf/min(1.0d99,dmaxval))
          ENDIF
          IF (temp_active) THEN
            CALL s_damax(I0*J0*K0,temp,dmaxval)
            depsilon = max(depsilon,nltolt/min(1.0d99,dmaxval))
          ENDIF
          IF (trans_active) THEN
            CALL s_damax(I0*J0*K0,conc,dmaxval)
            depsilon = max(depsilon,nltolc/min(1.0d99,dmaxval))
          ENDIF
!         highest parameter value
          dmaxval = p_k(1)
          DO i = 2, mpara
            dmaxval = max(dmaxval, p_k(i))
          END DO
!         estimated precision feasibility
          depsilon = depsilon *dmaxval
#ifdef DEBUG
          WRITE(*,*) ' [I] : -> BiCGstab solver runs until residuum <',depsilon
#endif
!       --
!       for "dnrm" dummy vector
        CALL set_dval(mpara,0.D0,tmp(1,3))
!       compute "p_delta'
        WRITE(*,*) ' '
        CALL omp_bayes_solve(mpara,p_delta,tmp(1,4),r, depsilon,max(mpara,ndata),0, &
          lss_lloctmp(1,1,1,ismpl),tmp(1,1),tmp(1,3),ismpl)
        WRITE(*,*) ' [I] : finish linear system for one bayes step'
        WRITE(*,*) ' '
!
! ----------------
#ifdef DEBUG
!       debug test: S^T *C_dd^-1 *S + reg_func() *C_pp^-1) * p_delta
!       - matrix free system -
        call step_bayes_tmpm(p_delta,tmp(1,1),tmp(1,3),ismpl)
!       test result
        do i = 1, min(mpara,10)
          WRITE(*,*) ' [d] : res =',dabs((tmp(i,3)-tmp(i,4))/tmp(i,4))
        enddo
        IF (mpara>10) WRITE(*,*) ' [d] : res = ...'
#endif
! ----------------
!
        IF (linfos(2)>=1) THEN
          WRITE(*,*) ' [I] : nrm2(p_k-p_(k+1)):', dsqrt(ddot(mpara,p_delta,1,p_delta,1))
          WRITE(*,*) ' [I] : end bayes step'
        END IF
!
        RETURN
      END


!>    @brief calculate "Jacobi^T vector product"; tmpm = <tmp4> *tmpv + (reg_func() *C_pp^-1) *tmpv
!>    @param[in] tmpv input vector
!>    @param[in,out] tmp temporary helper vectors (space)
!>    @param[out] tmpm result vector
!>    @param[in] ismpl local sample index
!>    @details
!> tmpm = (<tmp4> + reg_func() *C_pp^-1) *tmpv\n
!> [mpara]\n
!>
!> - tmpm = <tmp4> *tmpv\n
!> - tmpm = tmpm + (reg_func() *C_pp^-1) *tmpv\n
      SUBROUTINE step_bayes_tmpm(tmpv,tmp,tmpm,ismpl)
        use arrays
        use mod_data
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j
        DOUBLE PRECISION tmpm(mpara), tmpv(mpara)
        DOUBLE PRECISION tmp(ndata,2)
        DOUBLE PRECISION reg_func, wp2_ps, d_max
        LOGICAL test_null
        EXTERNAL reg_func, wp2_ps, test_null

!
!       test for trivial solution, to avoid useless computation
        CALL s_damax(mpara,tmpv,d_max)
        IF (test_null(d_max)) THEN
          WRITE(*,*) ' [I] : -> Bayes matrix vector product - trivial solution, skip computation'
          CALL set_dval(mpara,0.0d0,tmpm)
          RETURN
        END IF
!
        write(*,*) "[I] : Bayes matrix vector product - computing J^T *(J*x)"
        CALL step_bayes_tmp4(tmpv,tmp,tmpm,ismpl)
        DO j = 1, mpara
          DO i = 1, mpara
            tmpm(j) = tmpm(j) + reg_func()*wp2_ps(j,i,ismpl)*tmpv(i)
          END DO
        END DO
!
        RETURN
      END


!>    @brief calculate tmp4 = <tmp1> *jacobi *tmpv
!>    @param[in] tmpv input vector
!>    @param[in,out] tmp temporary helper vectors (space)
!>    @param[out] tmp4 result vector
!>    @param[in] ismpl local sample index
!>    @details
!> tmp4 = <tmp1> *jacobi *tmpv\n
!> [mpara]\n
!>
!> - tmp(:,1) = jacobi *tmpv\n
!> - tmp4 = <tmp1> *tmp(:,1)\n
      SUBROUTINE step_bayes_tmp4(tmpv,tmp,tmp4,ismpl)
        use arrays
        use mod_data
        IMPLICIT NONE
        integer :: ismpl
        DOUBLE PRECISION tmp4(mpara), tmpv(mpara)
        DOUBLE PRECISION tmp(ndata,2)
        integer start_jacmv,end_jacmv,rate

!
        call system_clock(start_jacmv,rate)
        CALL jac_mv(tmpv,tmp(1,1),ismpl)
        call system_clock(end_jacmv)
        write(*,*) "J*x timing: ",REAL(end_jacmv-start_jacmv)/REAL(rate)
        CALL step_bayes_tmp1(tmp(1,1),tmp(1,2),tmp4,ismpl)
!
        RETURN
      END


!>    @brief calculate tmp1 = jacobi^T *Wd2 *tmpv
!>    @param[in] tmpv input vector
!>    @param[in,out] tmp temporary helper vector (space)
!>    @param[out] tmp1 result vector
!>    @param[in] ismpl local sample index
!>    @details
!> tmp1 = jacobi^T *Wd2 *tmpv\n
!> [mpara]\n
!>
!> - tmp = Wd2 *tmpv\n
!> - tmp1 = jacobi^T *tmp\n
      SUBROUTINE step_bayes_tmp1(tmpv,tmp,tmp1,ismpl)
        use arrays
        use mod_data
        IMPLICIT NONE
        integer :: ismpl
        integer :: j, k
        DOUBLE PRECISION tmp1(mpara)
        DOUBLE PRECISION tmp(ndata), tmpv(ndata)
        DOUBLE PRECISION wd2_ps
        integer start_jactmv,end_jactmv,rate
        EXTERNAL wd2_ps

!
        DO j = 1, ndata
          tmp(j) = 0.0d0
          DO k = 1, ndata
            tmp(j) = tmp(j) + wd2_ps(j,k,ismpl)*tmpv(k)
          END DO
        END DO
        call system_clock(start_jactmv,rate)
        CALL jact_mv(tmp,tmp1,ismpl)
        call system_clock(end_jactmv)

        write(*,*) "J^T*x timing: ",REAL(end_jactmv-start_jactmv)/REAL(rate)
!
        RETURN
      END
