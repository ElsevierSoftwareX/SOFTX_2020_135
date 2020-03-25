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

!>    @brief simple line search: 1./2^relax-max < mu < 1
!>    @param[in] qval_old old quality function value
!>    @param[out] qval_new quality function value
!>    @param[in,out] p_k new parameter vector (iteration k)
!>    @param[in] p_apr apriori parameter vector
!>    @param[in,out] p_delta delta to the estimated parameters
!>    @param[in] ismpl_dummy local sample index (ignored)
!>    @details
!> line search by computing different relaxations of the parameter update\n
      SUBROUTINE update_model(qval_new,qval_old,p_k,p_apr,p_delta, &
          ismpl_dummy)
        use arrays
#ifdef AD
        use g_arrays
#endif
#ifdef AD_RM
        use arrays_ad
#endif
        use mod_genrl
        use mod_genrlc
        use mod_flow
        use mod_temp
        use mod_time
        use mod_inverse
        use mod_data
        use mod_OMP_TOOLS
        use mod_linfos
        IMPLICIT NONE
        integer :: ismpl
        double precision :: tflocal
        INCLUDE 'OMP_TOOLS.inc'
!
        DOUBLE PRECISION p_k(mpara), p_delta(mpara), p_apr(mpara)
        DOUBLE PRECISION, ALLOCATABLE :: tmp_p_k(:,:), qval_cur(:)
        DOUBLE PRECISION, ALLOCATABLE :: rms_data_ismpl(:)
        DOUBLE PRECISION, ALLOCATABLE :: rms_para_ismpl(:)
        DOUBLE PRECISION, ALLOCATABLE :: sdata_best(:)
        DOUBLE PRECISION, ALLOCATABLE :: src_best(:)
!     relaxation
        DOUBLE PRECISION qval_new, qval_old, relax
        INTEGER lblank, ismpl_dummy, relax_best, relax_von, relax_bis, relax_iter
        DOUBLE PRECISION qfunc
        EXTERNAL qfunc, lblank
        INTRINSIC min, dble


!     init
        ALLOCATE(tmp_p_k(mpara, relax_max +1))
        ALLOCATE(qval_cur(relax_max +1))
        ALLOCATE(rms_data_ismpl(relax_max +1))
        ALLOCATE(rms_para_ismpl(relax_max +1))
        IF (ndata > 0) ALLOCATE(sdata_best(ndata))
        ALLOCATE(src_best(I0*J0*K0))
        qval_new = qval_old

!     init for relaxation iteration
        relax_bis = 0
        write_disable = .TRUE.
        write_iter_disable = .TRUE.
!       -> none of them (out of range)
        relax_best = relax_max +2
! **********************************************************************
100     CONTINUE
!       start-intervall
        relax_von = relax_bis +1
!       end-intervall
        relax_bis = min(relax_bis + Tlevel_0, relax_max +1)

#ifdef fOMP
!$OMP   parallel num_threads(Tlevel_0) &
!$OMP     default(none) private(ismpl, relax_iter, relax) &
!$OMP     shared(relax_von, relax_bis, relax_max, mpara, Tlevel_0) &
!$OMP     shared(p_k, tmp_p_k, p_delta, simtime_0, max_simtime) &
!$OMP     shared(rms_data_ismpl, qval_old, tol_inv, ndata, relax_best) &
!$OMP     shared(rms_para_ismpl, runmode, transient) &
!$OMP     shared(main_input, main_output) &
!$OMP     shared(sdata, sdata_best, w, src_best, I0, J0, K0) &
!$OMP     shared(iter_inv, qval_cur, p_apr, rms_data, rms_para)
#endif
! !!! to avoid compiler bugs
        CALL omp_ordered_create(relax_max +1)
#ifdef fOMP
!--- !$OMP    do schedule(dynamic,1) ordered
!$OMP   do schedule(dynamic,1)
#endif
        DO relax_iter = relax_von, relax_bis
          ismpl = omp_get_his_thread_num() + 1
          CALL prepare_realisation(ismpl)
!         compute relaxation
          relax = 0.5D0**(relax_iter -1)
          IF (relax_iter>1) WRITE(*,'(A,F6.4,A)') &
            ' -> try relaxation=', relax, ' !'
!         init p_k
          CALL dcopy(mpara,p_k,1,tmp_p_k(1, relax_iter),1)
!         apply the update : linlog_input(:) == p^k+1 := p^k + relax*p_delta
          CALL update_param(tmp_p_k(1, relax_iter),p_delta,relax,ismpl)
!         update parameter in physical space
          CALL restore_jacobian(tmp_p_k(1, relax_iter),ismpl)
!         forward modeling
          CALL forward_compute(main_input(1,ismpl),main_output(1,ismpl),simtime_0,max_simtime,iter_inv,0,ismpl)
!$OMP critical
!         get the function value for position p_k+1
          qval_cur(relax_iter) = qfunc(tmp_p_k(1, relax_iter),p_apr,ismpl)
          rms_data_ismpl(relax_iter) = rms_data
          rms_para_ismpl(relax_iter) = rms_para

!         save sdata/src-init state, when no better found
          IF (relax_iter==1 .AND. relax_best > relax_max) THEN
            IF (ndata > 0) CALL dcopy(ndata,sdata(1,ismpl),1,sdata_best,1)
            CALL dcopy(I0*J0*K0,w(1,1,1,ismpl),1,src_best,1)
          END IF

!         when better and it is the first
          IF ((qval_cur(relax_iter)<qval_old .OR. &
              rms_data_ismpl(relax_iter)<=tol_inv*dble(ndata)) .AND. &
              relax_iter<relax_best) THEN
            relax_best = relax_iter
!           save best sdata/src-state
            IF (ndata > 0) CALL dcopy(ndata,sdata(1,ismpl),1,sdata_best,1)
            CALL dcopy(I0*J0*K0,w(1,1,1,ismpl),1,src_best,1)
!           steady-state: it is usefull to restore the computed state
!           transient: no need to re-use the computed state
!             runmode=2: init from steady-state, runmode=3: init from initial reading
!             re-use not needed time-state: negative ismpl indicates the master-copy (1)
            IF (.NOT. transient .AND. runmode==2) CALL old_save(cgen_time,-ismpl)
          END IF
!$OMP end critical
        END DO
#ifdef fOMP
!$OMP    end do
#endif
! !!! to avoid compiler bugs
          CALL omp_ordered_delete()
#ifdef fOMP
!$OMP    end parallel
#endif

!       next try
        IF (relax_bis<relax_max +1) THEN
          IF (relax_best>relax_max +1) GO TO 100
        END IF
! **********************************************************************

        ismpl = idx_master
        write_disable = .FALSE.
        write_iter_disable = .FALSE.
        IF (relax_best>relax_max +1) relax_best = 1
!       setup best relaxiation
        qval_new = qval_cur(relax_best)
        rms_data = rms_data_ismpl(relax_best)
        rms_para = rms_para_ismpl(relax_best)
        relax = 0.5D0**(relax_best -1)
        CALL dcopy(mpara,tmp_p_k(1, relax_best),1,p_k,1)
        IF (ndata > 0) CALL dcopy(ndata,sdata_best,1,sdata(1,ismpl),1)
        CALL dcopy(I0*J0*K0,src_best,1,w(1,1,1,ismpl),1)
!       update parameter in physical space
        CALL restore_jacobian(p_k,ismpl)
        IF (.NOT. transient .AND. runmode==2) THEN
!         steady-state: it is usefull to restore the computed state
!         transient: no need to use the computed state
!           runmode=2: init from steady-state, runmode=3: init from initial reading
!           re-use, not needed time-state: negative ismpl indicates the master-copy (1)
!         -> get master-copy
          CALL old_restore(cgen_time,-ismpl)
!         -> save master-copy for all
          CALL old_save(cgen_opti,ismpl)
        ENDIF

        IF (linfos(2)>=1) THEN
          IF (relax_max==0) THEN
            WRITE(*,'(2(A,1e16.8))') '  objective function:', &
              qval_old, '->', qval_new
          ELSE
            WRITE(*,'(3(A,1e16.8),A,I2)') '  relaxation:', &
              relax, ', objective function:', qval_old, '->', &
              qval_new, ', iteration:', relax_best
          END IF
        END IF

        IF ((ndata/=0) .AND. (mpara/=0) .AND. (relax_out/=0)) THEN
          WRITE(*,'(3A)') '  [W] : "', status_log_inv(1:lblank(status_log_inv)), '"'
          OPEN(76,file=status_log_inv,status='unknown', position='append')
          CALL sys_cputime(tflocal)
          WRITE(76,'(4(e14.6),2(1X,I5,1X),F11.2)') qval_new, &
            rms_data, rms_para, relax, iter_inv, relax_best, &
            (tflocal-tslocal)/60.0D0
          CLOSE(76)
        END IF

        DEALLOCATE(src_best)
        IF (ndata > 0) DEALLOCATE(sdata_best)
        DEALLOCATE(rms_para_ismpl)
        DEALLOCATE(rms_data_ismpl)
        DEALLOCATE(qval_cur)
        DEALLOCATE(tmp_p_k)

        RETURN
      END
