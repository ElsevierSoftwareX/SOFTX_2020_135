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

!>    @brief compute gradients by divided differences
!>    @param[in] seed_komp seeding component/index, needed from parent subroutine
!>    @param[in] ismpl local sample index
      SUBROUTINE dd_iter(seed_komp,ismpl)
        use arrays
#ifdef AD_RM 
        use arrays_ad
#else
        use g_arrays
#endif
        use mod_genrl
        use mod_genrlc
        use mod_data
        use mod_conc
        use mod_time
        use mod_inverse
        IMPLICIT NONE
        INTEGER i, j, k, ismpl
!     needed parameters from parent subroutine
        INTEGER seed_komp, s_u, s_k
        CHARACTER*16 pname
        DOUBLE PRECISION, ALLOCATABLE :: hh_head(:), hh_temp(:), &
          hh_conc(:)
        DOUBLE PRECISION, ALLOCATABLE :: hh_pres(:)
        DOUBLE PRECISION, ALLOCATABLE :: loctmp(:), lsdata(:)

        INTEGER lblank
        EXTERNAL lblank

!     for DD use
        LOGICAL test_null, hh_next, hh_first
        DOUBLE PRECISION ddot
        EXTERNAL test_null, ddot
        INTRINSIC dsqrt
        DOUBLE PRECISION hh, hh_bak, s_diff
        DOUBLE PRECISION ad_head, ad_temp, ad_conc, ad_pres
        DOUBLE PRECISION dd_head, dd_temp, dd_conc, dd_pres
        DOUBLE PRECISION dif_head, dif_temp, dif_conc, dif_pres
        DOUBLE PRECISION hh_difn1, hh_difn2, hh_difn3, hh_difn4, hh_difn5
        LOGICAL hh_difc1, hh_difc2, hh_difc3, hh_difc4, hh_difc5
!     backup of 'hh_difn*'
        DOUBLE PRECISION hh_difn1_old, hh_difn2_old, hh_difn3_old
        DOUBLE PRECISION hh_difn4_old, hh_difn5_old
!     save the right 'hh' value for output
        DOUBLE PRECISION hh_1, hh_2, hh_3, hh_4, hh_5, hh_step
        INTEGER how_often, how_omax
        PARAMETER (how_omax=2)
!     sqrt(0.1)
        PARAMETER (hh_step=0.3162277660168379332D0)


        ALLOCATE(hh_head(i0*j0*k0))
        ALLOCATE(hh_temp(i0*j0*k0))
        ALLOCATE(hh_conc(i0*j0*k0*max(ntrans,1)))
        ALLOCATE(hh_pres(i0*j0*k0))
        ALLOCATE(loctmp(i0*j0*k0*max(ntrans,1)))
        IF (ndata>0) ALLOCATE(lsdata(ndata))
        hh = 1.0D-2
        how_often = 0
        hh_1 = hh
        hh_2 = hh
        hh_3 = hh
        hh_4 = hh
        hh_5 = hh
        hh_difn1_old = 1.0D100
        hh_difn2_old = 1.0D100
        hh_difn3_old = 1.0D100
        hh_difn4_old = 1.0D100
        hh_difn5_old = 1.0D100
#ifndef AD_RM 
        CALL s_damax(i0*j0*k0,g_head(1,1,1,ismpl),ad_head)
        CALL s_damax(i0*j0*k0,g_temp(1,1,1,ismpl),ad_temp)
        CALL s_damax(i0*j0*k0*ntrans,g_conc(1,1,1,1,ismpl),ad_conc)
        CALL s_damax(i0*j0*k0,g_pres(1,1,1,ismpl),ad_pres)
#else
        CALL s_damax(i0*j0*k0,head_ad(1,1,1,ismpl),ad_head)
        CALL s_damax(i0*j0*k0,temp_ad(1,1,1,ismpl),ad_temp)
        CALL s_damax(i0*j0*k0*ntrans,conc_ad(1,1,1,1,ismpl),ad_conc)
        CALL s_damax(i0*j0*k0,pres_ad(1,1,1,ismpl),ad_pres)
#endif
        dd_head = 0.D0
        dd_temp = 0.D0
        dd_conc = 0.D0
        dd_pres = 0.D0
        dif_head = 0.D0
        dif_temp = 0.D0
        dif_conc = 0.D0
        dif_pres = 0.D0
        hh_first = .TRUE.
!     save orig. 'sdata'
        IF (ndata>0) CALL dcopy(ndata,sdata(1,ismpl),1,lsdata,1)
! ################################
1000    CONTINUE

!     restore state (before)
        CALL old_restore(cgen_opti,ismpl)

!     -- just only apply the step size '+h' (x_=x+h) --
        CALL dd_seeding(seed_komp,hh,hh_bak,ismpl)
!     one function step '+h', but may change AD variables (in case of none picard iters)
        CALL forward_iter(simtime_0,max_simtime,iter_inv,0,ismpl)
!     restore the org. values (x=x_-h), compute DD values in the pv-arrays
        CALL dd_jacobian(seed_komp,hh,hh_bak,ismpl)

!     save DD jacobian (x+h)
        CALL dcopy(i0*j0*k0,head(1,1,1,ismpl),1,hh_head,1)
        CALL dcopy(i0*j0*k0,temp(1,1,1,ismpl),1,hh_temp,1)
        CALL dcopy(i0*j0*k0*ntrans,conc(1,1,1,1,ismpl),1,hh_conc,1)
        CALL dcopy(i0*j0*k0,pres(1,1,1,ismpl),1,hh_pres,1)

!     restore state (before)
        CALL old_restore(cgen_opti,ismpl)

        hh = -hh
!     -- just only apply the step size '-h' (x_=x-h) --
        CALL dd_seeding(seed_komp,hh,hh_bak,ismpl)
!     one function step '-h', but may change AD variables (in case of none picard iters)
        CALL forward_iter(simtime_0,max_simtime,iter_inv,0,ismpl)
!     restore the org. values (x=x_+h), compute DD values in the pv-arrays
        CALL dd_jacobian(seed_komp,hh,hh_bak,ismpl)
!     restore '+hh'
        hh = -hh

        hh_next = .FALSE.

!     test the quality of (x+h) with (x-h)
!     # head
        CALL s_damax(i0*j0*k0,hh_head,s_diff)
        CALL dcopy(i0*j0*k0,hh_head,1,loctmp,1)
        CALL daxpy(i0*j0*k0,-1.0D0,head(1,1,1,ismpl),1,loctmp,1)
        IF ( .NOT. test_null(s_diff)) CALL dscal(i0*j0*k0, &
          1.0D0/s_diff,loctmp,1)
        hh_difn1 = dsqrt(ddot(i0*j0*k0,loctmp,1,loctmp,1))
!     # temp
        CALL s_damax(i0*j0*k0,hh_temp,s_diff)
        CALL dcopy(i0*j0*k0,hh_temp,1,loctmp,1)
        CALL daxpy(i0*j0*k0,-1.0D0,temp(1,1,1,ismpl),1,loctmp,1)
        IF ( .NOT. test_null(s_diff)) CALL dscal(i0*j0*k0, &
          1.0D0/s_diff,loctmp,1)
        hh_difn2 = dsqrt(ddot(i0*j0*k0,loctmp,1,loctmp,1))
!     # conc
        CALL s_damax(i0*j0*k0*ntrans,hh_conc,s_diff)
        CALL dcopy(i0*j0*k0*ntrans,hh_conc,1,loctmp,1)
        CALL daxpy(i0*j0*k0*ntrans,-1.0D0,conc(1,1,1,1,ismpl),1, &
          loctmp,1)
        IF ( .NOT. test_null(s_diff)) CALL dscal(i0*j0*k0*ntrans, &
          1.0D0/s_diff,loctmp,1)
        hh_difn3 = dsqrt(ddot(i0*j0*k0*ntrans,loctmp,1,loctmp,1))
!     # pres
        CALL s_damax(i0*j0*k0,hh_pres,s_diff)
        CALL dcopy(i0*j0*k0,hh_pres,1,loctmp,1)
        CALL daxpy(i0*j0*k0,-1.0D0,pres(1,1,1,ismpl),1,loctmp,1)
        IF ( .NOT. test_null(s_diff)) CALL dscal(i0*j0*k0, &
          1.0D0/s_diff,loctmp,1)
        hh_difn5 = dsqrt(ddot(i0*j0*k0,loctmp,1,loctmp,1))

!     save best-state
!     # head
        IF (hh_difn1<hh_difn1_old .AND. .NOT. test_null(hh_difn1)) &
            THEN
!        m = m + [x-h]
          CALL daxpy(i0*j0*k0,1.0D0,head(1,1,1,ismpl),1,hh_head,1)
!        m = m/2
          CALL dscal(i0*j0*k0,0.5D0,hh_head,1)
          hh_difn1_old = hh_difn1
          IF ( .NOT. test_null(hh_bak)) THEN
            hh_1 = hh_bak*hh
          ELSE
            hh_1 = hh
          END IF
!        try improvment (another iteration)
          hh_next = .TRUE.
          hh_difc1 = .TRUE.
        ELSE
          hh_difc1 = .FALSE.
        END IF
!     # temp
        IF (hh_difn2<hh_difn2_old .AND. .NOT. test_null(hh_difn2)) &
            THEN
!        m = m + [x-h]
          CALL daxpy(i0*j0*k0,1.0D0,temp(1,1,1,ismpl),1,hh_temp,1)
!        m = m/2
          CALL dscal(i0*j0*k0,0.5D0,hh_temp,1)
          hh_difn2_old = hh_difn2
          IF ( .NOT. test_null(hh_bak)) THEN
            hh_2 = hh_bak*hh
          ELSE
            hh_2 = hh
          END IF
!        try improvment (another iteration)
          hh_next = .TRUE.
          hh_difc2 = .TRUE.
        ELSE
          hh_difc2 = .FALSE.
        END IF
!     # conc
        IF (hh_difn3<hh_difn3_old .AND. .NOT. test_null(hh_difn3)) &
            THEN
!        m = m + [x-h]
          CALL daxpy(i0*j0*k0*ntrans,1.0D0,conc(1,1,1,1,ismpl),1, &
            hh_conc,1)
!        m = m/2
          CALL dscal(i0*j0*k0*ntrans,0.5D0,hh_conc,1)
          hh_difn3_old = hh_difn3
          IF ( .NOT. test_null(hh_bak)) THEN
            hh_3 = hh_bak*hh
          ELSE
            hh_3 = hh
          END IF
!        try improvment (another iteration)
          hh_next = .TRUE.
          hh_difc3 = .TRUE.
        ELSE
          hh_difc3 = .FALSE.
        END IF
!     # pres
        IF (hh_difn5<hh_difn5_old .AND. .NOT. test_null(hh_difn5)) &
            THEN
!        m = m + [x-h]
          CALL daxpy(i0*j0*k0,1.0D0,pres(1,1,1,ismpl),1,hh_pres,1)
!        m = m/2
          CALL dscal(i0*j0*k0,0.5D0,hh_pres,1)
          hh_difn5_old = hh_difn5
          IF ( .NOT. test_null(hh_bak)) THEN
            hh_5 = hh_bak*hh
          ELSE
            hh_5 = hh
          END IF
!        try improvment (another iteration)
          hh_next = .TRUE.
          hh_difc5 = .TRUE.
        ELSE
          hh_difc5 = .FALSE.
        END IF

!--- debug ------------------
!       head
        IF (hh_difc1 .OR. hh_first) THEN
          CALL s_damax(i0*j0*k0,hh_head,dd_head)
          CALL dcopy(i0*j0*k0,hh_head,1,loctmp,1)
#ifndef AD_RM
          CALL daxpy(i0*j0*k0,-1.0D0,g_head(1,1,1,ismpl),1,loctmp,1)
#else
           CALL daxpy(i0*j0*k0,-1.0D0,head_ad(1,1,1,ismpl),1,loctmp,1)
#endif
          dif_head = dsqrt(ddot(i0*j0*k0,loctmp,1,loctmp,1))
          WRITE(*,'(4(A,1e16.8))') '   up-head : nrm2(DD-AD)=', &
          dif_head, ' max(AD)=', ad_head, ' max(DD)=', dd_head, ' h=', &
          hh_1
        ELSE
          IF (.NOT. test_null(hh_difn1)) &
            WRITE(*,'(4(A,1e16.8))') '      head : nrm2(DD-AD)=', &
            dif_head, ' max(AD)=', ad_head, ' max(DD)=', dd_head, ' h=', &
            hh_1
        END IF
!       temp
        IF (hh_difc2 .OR. hh_first) THEN
          CALL s_damax(i0*j0*k0,hh_temp,dd_temp)
          CALL dcopy(i0*j0*k0,hh_temp,1,loctmp,1)
#ifndef AD_RM
          CALL daxpy(i0*j0*k0,-1.0D0,g_temp(1,1,1,ismpl),1,loctmp,1)
#else
          CALL daxpy(i0*j0*k0,-1.0D0,temp_ad(1,1,1,ismpl),1,loctmp,1)
#endif
          dif_temp = dsqrt(ddot(i0*j0*k0,loctmp,1,loctmp,1))
          WRITE(*,'(4(A,1e16.8))') '   up-temp: nrm2(DD-AD)=', &
          dif_temp, ' max(AD)=', ad_temp, ' max(DD)=', dd_temp, ' h=', &
          hh_2
        ELSE
          IF (.NOT. test_null(hh_difn2)) &
            WRITE(*,'(4(A,1e16.8))') '      temp: nrm2(DD-AD)=', &
            dif_temp, ' max(AD)=', ad_temp, ' max(DD)=', dd_temp, ' h=', &
            hh_2
        END IF
!       conc
        IF (hh_difc3 .OR. hh_first) THEN
          CALL s_damax(i0*j0*k0*ntrans,hh_conc,dd_conc)
          CALL dcopy(i0*j0*k0*ntrans,hh_conc,1,loctmp,1)
#ifndef AD_RM
          CALL daxpy(i0*j0*k0*ntrans,-1.0D0,g_conc(1,1,1,1,ismpl),1, &
            loctmp,1)
#else
           CALL daxpy(i0*j0*k0*ntrans,-1.0D0,conc_ad(1,1,1,1,ismpl),1, &
            loctmp,1)
#endif
          dif_conc = dsqrt(ddot(i0*j0*k0*ntrans,loctmp,1,loctmp,1))
          WRITE(*,'(4(A,1e16.8))') '   up-conc: nrm2(DD-AD)=', &
          dif_conc, ' max(AD)=', ad_conc, ' max(DD)=', dd_conc, ' h=', &
          hh_3
        ELSE
          IF (.NOT. test_null(hh_difn3)) &
            WRITE(*,'(4(A,1e16.8))') '      conc: nrm2(DD-AD)=', &
            dif_conc, ' max(AD)=', ad_conc, ' max(DD)=', dd_conc, ' h=', &
            hh_3
        END IF
!       pres
        IF (hh_difc5 .OR. hh_first) THEN
          CALL s_damax(i0*j0*k0,hh_pres,dd_pres)
          CALL dcopy(i0*j0*k0,hh_pres,1,loctmp,1)
#ifndef AD_RM 
          CALL daxpy(i0*j0*k0,-1.0D0,g_pres(1,1,1,ismpl),1,loctmp,1)
#else
          CALL daxpy(i0*j0*k0,-1.0D0,pres_ad(1,1,1,ismpl),1,loctmp,1)
#endif
          dif_pres = dsqrt(ddot(i0*j0*k0,loctmp,1,loctmp,1))
          WRITE(*,'(4(A,1e16.8))') '   up-pres: nrm2(DD-AD)=', &
          dif_pres, ' max(AD)=', ad_pres, ' max(DD)=', dd_pres, ' h=', &
          hh_5
        ELSE
          IF (.NOT. test_null(hh_difn5)) &
            WRITE(*,'(4(A,1e16.8))') '      pres: nrm2(DD-AD)=', &
            dif_pres, ' max(AD)=', ad_pres, ' max(DD)=', dd_pres, ' h=', &
            hh_5
        END IF
!--- debug ------------------

        IF ( .NOT. hh_next) THEN
!       how often is it more inexactly
          how_often = how_often + 1
        END IF
        hh_first = .FALSE.

!     decrease about a half number
        hh = hh*hh_step
        WRITE(*,'(1A,1e24.16,1A,1e24.16)') '   --> rel.h=', hh, &
          ', base=', hh_bak
        IF (hh>=1.0D-8 .AND. how_often<how_omax) GO TO 1000
        hh = hh/hh_step
! ################################

! ----------------------------------------------------------------------
        OPEN(88,file=project(1:lblank(project))//'_dd.log', &
          status='unknown',position='append')
        s_k = seed_para(1,seed_komp)
        s_u = seed_para(2,seed_komp)
        CALL param_name(seed_komp,pname,ismpl)
        IF ((s_k<=lastidx) .AND. (s_k>=1)) THEN
!        parameter units
          WRITE(88,'(3A,I2,A,I2,A,I5,A,I5)') key_char//' seeding: ', pname, &
            ', component= ', s_k - firstidx + 1, '/', &
            lastidx - firstidx + 1, ', unit= ', s_u, '/', maxunits
        ELSE IF ((s_k<=bc_lastidx) .AND. (s_k>=1)) THEN
!        bc units
          WRITE(88,'(3A,I5,A,I5)') key_char//' seeding: ', pname, ', unit= ', &
            s_u, '/', bc_maxunits
        ELSE IF (s_k==-1) THEN
!        tp units
          WRITE(88,'(3A,I5,A,I5)') key_char//' seeding: ', pname, ', unit= ', &
            s_u, '/', mpara_tp
        END IF

!       head
        WRITE(*,'(4(A,1e16.8))') '   head : nrm2(DD-AD)=', dif_head, &
          ' max(AD)=', ad_head, ' max(DD)=', dd_head, ' h=', hh_1
        WRITE(88,'(4(A,1e16.8))') 'head : nrm2(DD-AD)=', dif_head, &
          ' max(AD)=', ad_head, ' max(DD)=', dd_head, ' h=', hh_1
!       temp
        WRITE(*,'(4(A,1e16.8))') '   temp: nrm2(DD-AD)=', dif_temp, &
          ' max(AD)=', ad_temp, ' max(DD)=', dd_temp, ' h=', hh_2
        WRITE(88,'(4(A,1e16.8))') 'temp: nrm2(DD-AD)=', dif_temp, &
          ' max(AD)=', ad_temp, ' max(DD)=', dd_temp, ' h=', hh_2
!       conc
        WRITE(*,'(4(A,1e16.8))') '   conc: nrm2(DD-AD)=', dif_conc, &
          ' max(AD)=', ad_conc, ' max(DD)=', dd_conc, ' h=', hh_3
        WRITE(88,'(4(A,1e16.8))') 'conc: nrm2(DD-AD)=', dif_conc, &
          ' max(AD)=', ad_conc, ' max(DD)=', dd_conc, ' h=', hh_3
!       pres
        WRITE(*,'(4(A,1e16.8))') '   pres: nrm2(DD-AD)=', dif_pres, &
          ' max(AD)=', ad_pres, ' max(DD)=', dd_pres, ' h=', hh_5
        WRITE(88,'(4(A,1e16.8))') 'pres: nrm2(DD-AD)=', dif_pres, &
          ' max(AD)=', ad_pres, ' max(DD)=', dd_pres, ' h=', hh_5

        CLOSE(88)
! ----------------------------------------------------------------------
!     restore orig. 'sdata'
        IF (ndata>0) CALL dcopy(ndata,lsdata,1,sdata(1,ismpl),1)

        IF (ndata>0) DEALLOCATE(lsdata)
        DEALLOCATE(loctmp)
        DEALLOCATE(hh_pres)
        DEALLOCATE(hh_conc)
        DEALLOCATE(hh_temp)
        DEALLOCATE(hh_head)

        RETURN
      END
