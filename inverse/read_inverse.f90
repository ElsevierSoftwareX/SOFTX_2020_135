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

!>    @brief read inverse parameters and allocate fields
!>    @param[in] filename
!>    @param[in] ismpl local sample index
      SUBROUTINE read_inverse(filename,ismpl)
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
        use mod_conc
        use mod_inverse
        use mod_time
        use mod_data
        use mod_OMP_TOOLS
        use mod_linfos
        IMPLICIT NONE
        integer i, j, k, ismpl 
        INCLUDE 'OMP_TOOLS.inc'
        CHARACTER filename*80, line*5000, ctmp*4
        CHARACTER*4 stmp(max(nprop,nbc))
        INTEGER locstr, lblank, itmp, i2tmp, c1, c2, c3, c4, tmplen
        INTEGER mpara_p, mpara_bc
        LOGICAL found, no_ext_link, no_ext_link_int
        DOUBLE PRECISION dtmp
        DOUBLE PRECISION, ALLOCATABLE :: datmp(:,:)
        EXTERNAL locstr, found, no_ext_link, no_ext_link_int, lblank
!     not global anymore (better use opti_props or opti_bc) !
        INTEGER proplevel(nprop)
        INTEGER, ALLOCATABLE :: unitlevel(:)
        CHARACTER*4 sll(0:2)
        DATA sll/'    ', ' lin', ' log'/
        INTRINSIC trim


!       open file
        OPEN(79,file=filename,status='old')
        WRITE(*,*)
        WRITE(*,*) '  reading inverse parameter'
        WRITE(*,*)
!
!       init HDF5 support, when available
        CALL open_hdf5(' ')

        ALLOCATE(unitlevel(nunits))
!       default init
        reg_para = 1.0D0
        reg_type = 1
        maxiter_inv = 1
        tol_inv = 1.0D0
        iter_out = 1

        IF (found(79,key_char//' inverse',line,.FALSE.)) THEN
          WRITE(*,*) ' [R] : inversion parameters'
          READ(79,*) reg_type
          READ(79,*) maxiter_inv, tol_inv
          READ(79,*) resmat, covar
          READ(79,*) iter_out
        ELSE
          resmat = 0
          covar = 0
        END IF
#ifdef JACOBI_FREE
        WRITE(*,*) ' [I] : The use of the AD reverse mode implies (Jacobi) matrix free computation!'
        IF (resmat/=0.OR.covar/=0) THEN
          WRITE(*,*) ' [I] : Computation of Covariances and Resolution matrix DISABLED!'
          WRITE(*,*) '       Otherwise it can create too big Jacobi matrix!'
        END IF
        resmat = 0
        covar = 0
#endif

        IF (found(79,'# regular',line,.FALSE.)) THEN
          READ(79,*) reg_para
          WRITE(*,*) ' [R] : regularisation parameter: ', reg_para
        ELSE
          WRITE(*,*) ' <D> : no regularisation parameter !'
        END IF

!------------------------------------------------------------------------

#ifdef head_base
            nltol_g(pv_head) = nltolf
#endif
#ifdef pres_base
            nltol_g(pv_pres) = nltolf
#endif
            IF (head_active) THEN
              IF (found(79,key_char//' grad nliterf',line,.FALSE.)) THEN
#ifdef head_base
                READ(79,*,err=111,end=111) nltol_g(pv_head)
#endif
#ifdef pres_base
                READ(79,*,err=111,end=111) nltol_g(pv_pres)
#endif
                WRITE(*,*) ' [R] : gradient flow nonlinear iteration tolerance (rel.)'
111             CONTINUE
              END IF
            END IF
            nltol_g(pv_temp) = nltolt
            IF (temp_active) THEN
              IF (found(79,key_char//' grad nlitert',line,.FALSE.)) THEN
                READ(79,*,err=112,end=112) nltol_g(pv_temp)
                WRITE(*,*) ' [R] : gradient temperature &
                  &nonlinear iteration tolerance (rel.)'
112             CONTINUE
              END IF
            END IF
            nltol_g(pv_conc) = nltolc
            IF (trans_active) THEN
              IF (found(79,key_char//' grad nliterc',line,.FALSE.)) THEN
                READ(79,*,err=113,end=113) nltol_g(pv_conc)
                WRITE(*,*) ' [R] : gradient transport nonlinear &
                  &iteration tolerance (rel.)'
113             CONTINUE
              END IF
            END IF

!------------------------------------------------------------------------

        IF (found(79,key_char//' errors',line,.FALSE.)) THEN
          ALLOCATE(datmp(nprop_load,maxunits))
          IF (no_ext_link(nprop_load,maxunits,1,datmp,'errors',line)) &
              THEN
            CALL read_array(79,nprop_load,maxunits,datmp,'# errors', &
              ismpl)
          END IF
          IF (nprop_load==lastidx-firstidx+1) THEN
!          load all, dense entries (no specific index)
            DO i = 1, maxunits
              DO j = 1, nprop_load
                d_propunit(i,firstidx-1+j) = datmp(j,i)
              END DO
            END DO
          ELSE
!          needs to handle manual reading with specific index
            WRITE(*,'(1A)') &
              'error in "read_inverse.f": reading errors: specific index handling not implemented, use dense entries instead'
            STOP
          END IF

          DO i = 1, maxunits
!          because of logarithimc scale, suppress zeros
            DO j = firstidx, lastidx
              d_propunit(i,j) = max(1.0D-10,d_propunit(i,j))
            END DO
            IF (linfos(2)>=2) WRITE(*,'('//c_npropunit//'e12.4,1I8)') &
              (d_propunit(i,j),j=firstidx,lastidx), i

          END DO
          WRITE(*,*) ' [R] : variances apriori'
          DEALLOCATE(datmp)
        ELSE
          WRITE(*,*) ' <D> : variances apriori !'
          DO i = 1, maxunits
            DO j = firstidx, lastidx
              d_propunit(i,j) = 1.0D0
            END DO
            IF (linfos(2)>=2) WRITE(*,'('//c_npropunit//'e12.4,1I8)') &
              (d_propunit(i,j),j=firstidx,lastidx), i
          END DO
        END IF

        IF (found(79,key_char//' bcerrors',line,.FALSE.)) THEN
          DO i = 1, bc_maxunits
            DO j = bc_firstidx, bc_lastidx
              d_propunit(i,j) = 1.0D0
            END DO
          END DO
          CALL get_arg('records',line,i,j)
          IF (i<1 .OR. j<i) THEN
            READ(79,*) tmplen
          ELSE
            READ(line(i:j),*) tmplen
          END IF

          DO i = 1, tmplen
            READ(79,*) k, dtmp, ctmp
            IF ((k>bc_maxunits) .OR. (k<1)) THEN
              WRITE(*,'(A,I7,A,I7,A)') 'error: "bcerrors": &
                &bc unit number out of range or not used, &
                &(', k, ') at line ', i, '!'
              STOP
            END IF
            IF ((ctmp/='head') .AND. (ctmp/='temp') .AND. &
                (ctmp/='pres') .AND. (ctmp/='conc')) THEN
              WRITE(*,'(A,A1,A,I3,A)') &
                'error: "bcerrors": bc unit type not allowed, "', &
                ctmp, '" at line ', i, '!'
              STOP
            END IF

            IF (ctmp=='head') d_propunit(k,idx_hbc) = max(1.D-10,dtmp)
            IF (ctmp=='temp') d_propunit(k,idx_tbc) = max(1.D-10,dtmp)
            IF (ctmp=='conc') d_propunit(k,idx_cbc) = max(1.D-10,dtmp)
            IF (ctmp=='pres') d_propunit(k,idx_hbc) = max(1.D-10,dtmp)
          END DO
          WRITE(*,'(A,I3)') '  [R] : bc variances apriori, records=', &
            tmplen
        ELSE
          WRITE(*,*) ' <D> : bc variances apriori !'
          DO i = 1, bc_maxunits
            DO j = bc_firstidx, bc_lastidx
              d_propunit(i,j) = 1.0D0
            END DO
          END DO
        END IF

!------------------------------------------------------------------------

        IF (found(79,key_char//' apriori',line,.FALSE.)) THEN
          ALLOCATE(datmp(nprop_load,maxunits))
          IF (no_ext_link(nprop_load,maxunits,1,datmp,'apriori',line)) &
              THEN
            CALL read_array(79,nprop_load,maxunits,datmp,'# apriori', &
              ismpl)
          END IF
          IF (nprop_load==lastidx-firstidx+1) THEN
!          load all, dense entries (no specific index)
            DO i = 1, maxunits
              DO j = 1, nprop_load
                a_propunit(i,firstidx-1+j) = datmp(j,i)
              END DO
            END DO
          ELSE
!          needs to handle manual reading with specific index
            WRITE(*,'(1A)') &
             'error in "read_inverse.f": reading apriori: specific index handling not implemented, use dense entries instead'
            STOP
          END IF

          DO i = 1, maxunits
!          because of logarithimc scale, suppress zeros
            a_propunit(i,idx_por) = max(prop_min(idx_por), &
              a_propunit(i,idx_por))
            a_propunit(i,idx_an_kx) = max(prop_min(idx_an_kx), &
              a_propunit(i,idx_an_kx))
            a_propunit(i,idx_an_ky) = max(prop_min(idx_an_ky), &
              a_propunit(i,idx_an_ky))
            a_propunit(i,idx_kz) = max(prop_min(idx_kz), &
              a_propunit(i,idx_kz))
!?           a_propunit(i,idx_Comp)    = max(prop_min(idx_Comp),
!?     &       a_propunit(i,idx_Comp))
            a_propunit(i,idx_an_lx) = max(prop_min(idx_an_lx), &
              a_propunit(i,idx_an_lx))
            a_propunit(i,idx_an_ly) = max(prop_min(idx_an_ly), &
              a_propunit(i,idx_an_ly))
            a_propunit(i,idx_lz) = max(prop_min(idx_lz), &
              a_propunit(i,idx_lz))
            a_propunit(i,idx_q) = max(prop_min(idx_q), &
              a_propunit(i,idx_q))
            a_propunit(i,idx_rc) = max(prop_min(idx_rc), &
              a_propunit(i,idx_rc))
            a_propunit(i,idx_df) = max(prop_min(idx_df), &
              a_propunit(i,idx_df))
            a_propunit(i,idx_ec) = max(prop_min(idx_ec), &
              a_propunit(i,idx_ec))
!?           a_propunit(i,idx_lC)    = max(prop_min(idx_lC),
!?     &       a_propunit(i,idx_lC))
            IF (linfos(2)>=2) WRITE(*,'('//c_npropunit//'e12.4,1I8)') &
              (a_propunit(i,j),j=firstidx,lastidx), i
          END DO
          WRITE(*,*) ' [R] : unit properties apriori'
          DEALLOCATE(datmp)
        ELSE
          WRITE(*,*) ' <D> : unit properties apriori !'
          DO i = 1, maxunits
            a_propunit(i,idx_por) = max(prop_min(idx_por), &
              propunit(i,idx_por,ismpl))
            a_propunit(i,idx_an_kx) = max(prop_min(idx_an_kx), &
              propunit(i,idx_an_kx,ismpl))
            a_propunit(i,idx_an_ky) = max(prop_min(idx_an_ky), &
              propunit(i,idx_an_ky,ismpl))
            a_propunit(i,idx_kz) = max(prop_min(idx_kz), &
              propunit(i,idx_kz,ismpl))
            a_propunit(i,idx_comp) = propunit(i,idx_comp,ismpl)
            a_propunit(i,idx_an_lx) = max(prop_min(idx_an_lx), &
              propunit(i,idx_an_lx,ismpl))
            a_propunit(i,idx_an_ly) = max(prop_min(idx_an_ly), &
              propunit(i,idx_an_ly,ismpl))
            a_propunit(i,idx_lz) = max(prop_min(idx_lz), &
              propunit(i,idx_lz,ismpl))
            a_propunit(i,idx_q) = max(prop_min(idx_q), &
              propunit(i,idx_q,ismpl))
            a_propunit(i,idx_rc) = max(prop_min(idx_rc), &
              propunit(i,idx_rc,ismpl))
            a_propunit(i,idx_df) = max(prop_min(idx_df), &
              propunit(i,idx_df,ismpl))
            a_propunit(i,idx_ec) = max(prop_min(idx_ec), &
              propunit(i,idx_ec,ismpl))
            a_propunit(i,idx_lc) = propunit(i,idx_lc,ismpl)
            IF (linfos(2)>=2) WRITE(*,'('//c_npropunit//'e12.4,1I8)') &
              (a_propunit(i,j),j=firstidx,lastidx), i
          END DO
        END IF

        IF (found(79,key_char//' bcapriori',line,.FALSE.)) THEN
          DO i = 1, bc_maxunits
            a_propunit(i,idx_hbc) = max(0.D0,propunit(i,idx_hbc,ismpl))
            a_propunit(i,idx_tbc) = propunit(i,idx_tbc,ismpl)
            a_propunit(i,idx_cbc) = propunit(i,idx_cbc,ismpl)
          END DO
          CALL get_arg('records',line,i,j)
          IF (i<1 .OR. j<i) THEN
            READ(79,*) tmplen
          ELSE
            READ(line(i:j),*) tmplen
          END IF

          DO i = 1, tmplen
            READ(79,*) k, dtmp, ctmp
            IF ((k>bc_maxunits) .OR. (k<1)) THEN
              WRITE(*,'(A,I7,A,I7,A)') 'error: "bcapriori": &
                &bc unit number out of range or not used, &
                &(', k, ') at line ', i, ' !!!'
              STOP
            END IF
            IF ((ctmp/='head') .AND. (ctmp/='temp') .AND. &
                (ctmp/='pres') .AND. (ctmp/='conc')) THEN
              WRITE(*,'(A,A1,A,I3,A)') &
                'error: "bcapriori": bc unit type not allowed, "', &
                ctmp, '" at line ', i, ' !!!'
              STOP
            END IF

            IF (ctmp=='head') a_propunit(k,idx_hbc) = max(0.D0,dtmp)
            IF (ctmp=='temp') a_propunit(k,idx_tbc) = dtmp
            IF (ctmp=='conc') a_propunit(k,idx_cbc) = dtmp
            IF (ctmp=='pres') a_propunit(k,idx_hbc) = dtmp
          END DO
          WRITE(*,'(A,I3)') &
            '  [R] : bc unit properties apriori, records=', tmplen
        ELSE
          WRITE(*,*) ' <D> : bc unit properties apriori !'
          DO i = 1, bc_maxunits
            a_propunit(i,idx_hbc) = max(0.D0,propunit(i,idx_hbc,ismpl))
            a_propunit(i,idx_tbc) = propunit(i,idx_tbc,ismpl)
            a_propunit(i,idx_cbc) = propunit(i,idx_cbc,ismpl)
          END DO
        END IF

! -----------------------------------------------------------------------

!     read opt. switChes ---
        IF (found(79,key_char//' enable property',line,.FALSE.)) THEN
          WRITE(*,*) ' [R] : property map (on/off)'
!         default - disable the property
          DO i = 1, nprop_load
            proplevel(i) = 0
          END DO
!         read level
          READ(79,*,err=241,end=241) (proplevel(i),i=1,nprop_load)
          GOTO 200
241       IF (i-1<nprop_load) WRITE(*,'(A,I2,A,I2,A)') &
            '  <D> : WARNING, to few properties specified ( read ', i-1, ' of ', nprop_load, ') !'
200       CONTINUE
        ELSE
          WRITE(*,*) ' <D> :  transformation types for properties !'
          DO i = 1, nprop_load
            proplevel(i) = 0
          END DO
          WRITE(*,'(A,13I2)') '        #', (proplevel(i),i=1, &
            nprop_load)
        END IF

!  --------------------
!     init
        DO i = 1, maxunits
          unitlevel(i) = 1
        END DO
!     set
        IF (found(79,key_char//' enable unit',line,.FALSE.)) THEN
          WRITE(*,*) ' [R] : unit map (on/off)'
          READ(79,*) (unitlevel(i),i=1,maxunits)
        ELSE
          WRITE(*,*) ' <D> : unit map (on) !'
        END IF

!  --------------------
!     init
        DO j = 1, maxunits
          DO i = 1, nprop_load
            IF ((proplevel(i)>0) .AND. (unitlevel(j)>0)) THEN
              opti_props(i,j) = proplevel(i)
            ELSE
              opti_props(i,j) = 0
            END IF
          END DO
        END DO
!     set
        IF (found(79,key_char//' optimize property',line,.FALSE.)) THEN
          CALL get_arg('records',line,i,j)
          IF (i>0 .AND. j>=i) THEN
            READ(line(i:j),*) k
          ELSE
            READ(79,*) k
          END IF
          WRITE(*,'(A,I3)') &
            '  [R] : property optimization table, records=', k

          DO i = 1, k
            READ(79,'(A)') line
            READ(line,*) itmp, i2tmp
            IF (itmp>maxunits) THEN
              WRITE(*,'(A)') 'error: "optimize property" &
                &unit number out of range !!!'
              STOP
            END IF
            IF (locstr(line,'por')>=1) opti_props(idx_por,itmp) &
              = i2tmp
            IF (locstr(line,'a_kx')>=1) opti_props(idx_an_kx,itmp) &
              = i2tmp
            IF (locstr(line,'a_ky')>=1) opti_props(idx_an_ky,itmp) &
              = i2tmp
            IF (locstr(line,'kz')>=1) opti_props(idx_kz,itmp) = i2tmp
            IF (locstr(line,'comp')>=1) opti_props(idx_comp,itmp) &
              = i2tmp
            IF (locstr(line,'a_lx')>=1) opti_props(idx_an_lx,itmp) &
              = i2tmp
            IF (locstr(line,'a_ly')>=1) opti_props(idx_an_ly,itmp) &
              = i2tmp
            IF (locstr(line,'lz')>=1) opti_props(idx_lz,itmp) = i2tmp
            IF (locstr(line,'q')>=1) opti_props(idx_q,itmp) = i2tmp
            IF (locstr(line,'rc')>=1) opti_props(idx_rc,itmp) = i2tmp
            IF (locstr(line,'df')>=1) opti_props(idx_df,itmp) = i2tmp
            IF (locstr(line,'ec')>=1) opti_props(idx_ec,itmp) = i2tmp
            IF (locstr(line,'lc')>=1) opti_props(idx_lc,itmp) = i2tmp
          END DO
        END IF
!     final count
        mpara_p = 0
        DO j = 1, maxunits
          DO i = 1, lastidx - firstidx + 1
            IF (opti_props(i,j)>0) mpara_p = mpara_p + 1
          END DO
        END DO

!  --------------------
!     init
        DO j = 1, bc_maxunits
          DO i = 1, nbc
            opti_bc(i,j) = 0
          END DO
        END DO
!     set
        IF (found(79,key_char//' optimize bc',line,.FALSE.)) THEN
          CALL get_arg('records',line,i,j)
          IF (i>0 .AND. j>=i) THEN
            READ(line(i:j),*) k
          ELSE
            READ(79,*) k
          END IF
          WRITE(*,'(A,I3)') '  [R] : bc optimization table, records=', k
          IF (k>0 .AND. transient) THEN
            WRITE(*,'(1A)') 'error: bc optimization not allowed for time dependend models!'
            WRITE(*,'(2A)') '  -> please, update this as tp optimization with the ', &
              'simulation begin as the start time.'
            STOP
          END IF
!
          DO i = 1, k
!          avaiting: [bc-unit   (1|2)  (head|temp|conc)]
!                 == [unit-id  lin/log      pv-type    ]
            READ(79,'(A)') line
            READ(line,*) itmp, i2tmp
            IF (itmp>bc_maxunits) THEN
              WRITE(*,'(A)') &
                'error: "optimize bc" unit number out of range !!!'
              STOP
            END IF
            IF (locstr(line,'head')>=1) opti_bc(1,itmp) = i2tmp
            IF (locstr(line,'temp')>=1) opti_bc(2,itmp) = i2tmp
            IF (locstr(line,'conc')>=1) opti_bc(3,itmp) = i2tmp
            IF (locstr(line,'pres')>=1) opti_bc(3,itmp) = i2tmp
          END DO
        END IF
!     final count
        mpara_bc = 0
        DO j = 1, bc_maxunits
          DO i = 1, bc_lastidx - bc_firstidx + 1
            IF (opti_bc(i,j)>0) mpara_bc = mpara_bc + 1
          END DO
        END DO

!  --------------------
!     init
        CALL set_ival(3*ngsmax*mopti_tp*nbctp,0,opti_tp)
        mpara_tp = 0
!     set
        IF (transient) THEN
          IF (found(79,key_char//' optimize tp',line,.FALSE.)) THEN
            CALL get_arg('records',line,i,j)
            IF (i>0 .AND. j>=i) THEN
              READ(line(i:j),*) c1
            ELSE
              READ(79,*) c1
            END IF
            WRITE(*,'(A,I3)') '  [R] : tp optimization table, records=', c1
            DO i = 1, c1
!             [tp entry index, bc index number, alpha-beta keyword]
              READ(79,'(A)') line
!             Here the bc-index needs to be assoziated with a unit number and
!               defines this the dependency to the tp period table (where to find the tp-entry)
              READ(line,*) itmp, i2tmp
!             sanity check
              IF (itmp>ngsmax) THEN
                WRITE(*,'(A)') 'error: "optimize tp" period number out of range !!!'
                STOP
              END IF
              IF (i2tmp>nbctp) THEN
                WRITE(*,'(A)') 'error: "optimize tp" bc index out of range !!!'
                STOP
              END IF
!           default: linear, log = disabled
!aw              c2 = 1
!aw              IF (locstr(line,'lin')>=1) c2 = 1
!aw              IF (locstr(line,'log')>=1) c2 = 2
!           al:alpha, be:beta (2 types -> mopti_tp=2)
              IF (locstr(line,'al')>=1) THEN
                mpara_tp = mpara_tp + 1
                opti_tp(1,mpara_tp) = itmp
                opti_tp(2,mpara_tp) = 1
                opti_tp(3,mpara_tp) = i2tmp
!aw              opti_tp(4,mpara_tp) = C2
              ELSE IF (locstr(line,'be')>=1) THEN
                mpara_tp = mpara_tp + 1
                opti_tp(1,mpara_tp) = itmp
                opti_tp(2,mpara_tp) = 2
                opti_tp(3,mpara_tp) = i2tmp
!aw              opti_tp(4,mpara_tp) = C2
              ELSE
                WRITE(*,'(1A,1I4,3A)') '  [!] : ignore ',i,'. line = "',trim(line),'"'
              END IF
            END DO
!           sanity check
            IF (mpara_tp>ngsmax*mopti_tp*max(nbctp,1)) THEN
              WRITE(*,'(A)') 'error: number of time period optimisation entries greater than expected !!!'
              STOP
            END IF
          END IF
        END IF

!     set "mpara" dynamic, must be done before any "inverse" allocation
        mpara = mpara_p + mpara_bc + mpara_tp
        IF (mpara==0) THEN
          WRITE(*,'(1A)') &
            'error: no optimization parameters specified!'
          STOP
        END IF
!     --- end opt. switches ---

        IF (found(79,key_char//' read output',line,.FALSE.)) THEN
          lread_joutt = .TRUE.
          resmat = 1
          covar = 1
          ALLOCATE(seed_index(nsmpl+1))
          WRITE(*,*) ' [I] : READ-MODE for state variables'
!          IF (mpara /= Tlevel_0) THEN
!!           avoid "barrier in a loop" errors in "read_joutt_hdf"
!            Tlevel_0 = 1
!            WRITE(*,*) ' [I] : cutting sample parallelisation to 1, to avoid errors for READ-MODE'
!          END IF
        ELSE IF (found(79,key_char//' write output',line,.FALSE.)) THEN
          lwrite_joutt = .TRUE.
          resmat = 1
          covar = 1
          WRITE(*,*) ' [I] : WRITE-MODE for state variables'
        END IF

!  --------------------

        IF ((transient) .AND. (mpara_tp>0)) THEN
!         todo: read values from file
          CALL set_dval(ngsmax*2*nbctp,1.D0,tpwgt)
          CALL set_dval(ngsmax*2*nbctp,0.D0,e_bcperiod)
!         default, when not all readed
          CALL set_dval(ngsmax*2*nbctp,1.0D0,d_bcperiod)
          IF (found(79,key_char//' tperrors',line,.FALSE.)) THEN
            CALL get_arg('records',line,i,j)
            IF (i>0 .AND. j>=i) THEN
              READ(line(i:j),*) c1
            ELSE
              READ(79,*) c1
            END IF
            WRITE(*,'(A,I3)') '  [R] : tp variances apriori, records=', c1
!
            DO i = 1, c1
              READ(79,'(A)') line
!             p-idx, bc-idx, value
              READ(line,*) itmp, i2tmp, dtmp
              IF (itmp>ngsmax) THEN
                WRITE(*,'(A)') 'error: "tperrors" period number out of range !!!'
                STOP
              END IF
              IF (i2tmp>nbctp) THEN
                WRITE(*,'(A)') 'error: "tperrors" bc index out of range !!!'
                STOP
              END IF
!             al:alpha, be:beta
              IF (locstr(line,'al')>=1) THEN
                d_bcperiod(itmp,1,i2tmp) = dtmp
              ELSE IF (locstr(line,'be')>=1) THEN
                d_bcperiod(itmp,2,i2tmp) = dtmp
              ELSE
                WRITE(*,'(A)') 'error: "tperrors", "alpha" or "beta" type must be specified !!!'
              END IF
            END DO
          ELSE
            WRITE(*,*) ' <D> : tp variances apriori !'
          END IF
!
!         default, when not all readed
          DO k = 1, nbctp
            DO j = 1, 2
              DO i = 1, ngsmax
                a_bcperiod(i,j,k) = bcperiod(i,j+1,k,ismpl)
              END DO
            END DO
          END DO
          IF (found(79,key_char//' tpapriori',line,.FALSE.)) THEN
            CALL get_arg('records',line,i,j)
            IF (i>0 .AND. j>=i) THEN
              READ(line(i:j),*) c1
            ELSE
              READ(79,*) c1
            END IF
            WRITE(*,'(A,I3)') '  [R] : time periods apriori, records=', c1
!
            DO i = 1, c1
              READ(79,'(A)') line
!             p-idx, bc-idx, value
              READ(line,*) itmp, i2tmp, dtmp
              IF (itmp>ngsmax) THEN
                WRITE(*,'(A)') 'error: "tpapriori" period number out of range !!!'
                STOP
              END IF
              IF (i2tmp>nbctp) THEN
                WRITE(*,'(A)') 'error: "tpapriori" bc index out of range !!!'
                STOP
              END IF
!             al:alpha, be:beta
              IF (locstr(line,'al')>=1) THEN
                a_bcperiod(itmp,1,i2tmp) = dtmp
              ELSE IF (locstr(line,'be')>=1) THEN
                a_bcperiod(itmp,2,i2tmp) = dtmp
              ELSE
                WRITE(*,'(A)') 'error: "tpapriori", "alpha" or "beta" type must be specified !!!'
              END IF
            END DO
          ELSE
            WRITE(*,*) ' <D> : time periods apriori !'
          END IF
        END IF

!     read weights
!     init
        DO i = 1, nprop_load
          propwgt(i) = 1.D0
        END DO
!     set
        IF (found(79,key_char//' weight property',line,.FALSE.)) THEN
          WRITE(*,*) ' [R] : weights for properties'
          READ(79,*,err=271) (propwgt(i),i=1,nprop_load)
          GO TO 272
271       WRITE(*,'(A,I2,A)') ' warning: to few weights specified (', &
            nprop_load, ') ! 1.0d0 assumed'
272       CONTINUE
        ELSE
          WRITE(*,*) ' <D> : uniform weights for properties used !'
        END IF

!     init
        DO i = 1, nbc
          bcwgt(i) = 1.D0
        END DO
!     set
        IF (found(79,key_char//' weight bc',line,.FALSE.)) THEN
          WRITE(*,*) ' [R] : weights for boundary conditions'
          READ(79,*,err=281) (bcwgt(i),i=1,nbc)
          GO TO 282
281       WRITE(*,'(A,I2,A)') ' warning: to few weights specified (', &
            nbc, ') ! 1.0d0 assumed'
282       CONTINUE
        ELSE
          WRITE(*,*) &
            ' <D> : uniform weights for boundary conditions !'
        END IF

!     init
        para_weight = 0
!     set
        IF (.NOT. lib_override) THEN
          IF (found(79,key_char//' covar prior para',line,.FALSE.)) THEN
            para_weight = 1
            ALLOCATE(covar_prior_p(mpara,mpara))
            memory = memory + mpara*mpara
!
            CALL get_arg('inverse',line,i,j)
            IF (i>0 .AND. j>=i) THEN
              icovarp = 1
              IF (no_ext_link(mpara,mpara,1,covar_prior_p, &
                  'covar_prior_p',line)) THEN
                DO j = 1, mpara
                  READ(79,*) (covar_prior_p(i,j),i=1,mpara)
                END DO
              END IF
              WRITE(*,*) &
                ' [R] : prior parameter covariance matrix (inverse)'
            ELSE
              icovarp = 0
              WRITE(*,*) ' [R] : prior parameter covariance matrix'
              IF (no_ext_link(mpara,mpara,1,covar_prior_p, &
                  'covar_prior_p',line)) THEN
                DO j = 1, mpara
                  READ(79,*) (covar_prior_p(i,j),i=1,mpara)
                END DO
              END IF
            END IF
          END IF
        ELSE
          WRITE(*,*) ' <D> : ignore section "# covar prior para" when in LIBRARY mode'
        END IF

!     init
        data_weight = 0
!     set
        IF (.NOT. lib_override) THEN
          IF (found(79,key_char//' covar prior data',line,.FALSE.)) THEN
            data_weight = 1
            ALLOCATE(covar_prior_d(ndata,ndata))
            memory = memory + ndata*ndata
!
            CALL get_arg('inverse',line,i,j)
            IF (i>0 .AND. j>=i) THEN
              icovard = 1
              IF (no_ext_link(ndata,ndata,1,covar_prior_d, &
                  'covar_prior_p',line)) THEN
                DO j = 1, ndata
                  READ(79,*) (covar_prior_d(i,j),i=1,ndata)
                END DO
              END IF
              WRITE(*,*) &
                ' [R] : prior data covariance matrix (inverse)'
            ELSE
              icovard = 0
              WRITE(*,*) ' [R] : prior data covariance matrix'
              IF (no_ext_link(ndata,ndata,1,covar_prior_d, &
                  'covar_prior_p',line)) THEN
                DO j = 1, ndata
                  READ(79,*) (covar_prior_d(i,j),i=1,ndata)
                END DO
              END IF
            END IF
          END IF
        ELSE
          WRITE(*,*) ' <D> : ignore section "'//key_char//' covar prior data" when in LIBRARY mode'
        END IF

!     print optimizations tables
        IF (linfos(2)>=1) THEN
!        parameter units
          IF (maxunits>=1) THEN
            WRITE(*,*) ' '
            WRITE(*,'(6X,A)') 'optimization matrix [properties]:'
            WRITE(*,'(8X,2A,'//c_npropunit//'(A4,1X),A1)') '| ', 'unit ', &
              (properties(i),i=firstidx,lastidx), '|'
            DO j = 1, maxunits
              DO k = 1, nprop_load
                stmp(k) = sll(opti_props(k,j))
              END DO
              WRITE(*,'(8X,A2,I4,1X,'//c_npropunit//'(A4,1X),A1)') '| ', j, &
                (stmp(i),i=1,nprop_load), '|'
            END DO
          END IF

!        bc units
          IF (bc_maxunits>=1) THEN
            WRITE(*,*) ' '
            WRITE(*,'(6X,A)') &
              'optimization matrix [boundary condition]:'
            WRITE(*,'(8X,1A2,1A6,1X,'//c_nbcunit//'(A4,1X),1A1)') '| ', 'BCunit', &
              'Flow', 'Temp', 'Conc', '|'
            DO j = 1, bc_maxunits
              DO i = 1, nbc
                stmp(i) = sll(opti_bc(i,j))
              END DO
              WRITE(*,'(8X,1A2,1I6,1X,'//c_nbcunit//'(A4,1X),1A1)') '| ', j, &
                (stmp(i),i=1,nbc), '|'
            END DO
          END IF

!        bctp
          IF (mpara_tp>=1) THEN
            WRITE(*,*) ' '
            WRITE(*,'(6X,A)') 'optimization matrix [time &
              &depended boundary condition]:'
            c1 = ibcperiod(1)
            DO i = 2, nbctp
              c1 = max(c1,ibcperiod(i))
            END DO
            WRITE(line,'(1000(I4,1X))') (i,i=1,c1)
            WRITE(*,'(8X,A12,1A,A1)') '| bctp | tp:', line(1:c1*5), &
              '|'
            DO i = 1, nbctp
              line = ' '

              DO j = 1, mpara_tp
                k = (opti_tp(1,j)-1)*5 + 1
                c2 = opti_tp(2,j)
                IF (opti_tp(3,j)==i .AND. line(k:k)==' ' .AND. c2==1) &
                  line(k:k+4) = 'alfa '
                IF (opti_tp(3,j)==i .AND. line(k:k)==' ' .AND. c2==2) &
                  line(k:k+4) = 'beta '
                IF (opti_tp(3,j)==i .AND. line(k:k)=='b' .AND. c2==1) &
                  line(k:k+4) = ' a+b '
                IF (opti_tp(3,j)==i .AND. line(k:k)=='a' .AND. c2==2) &
                  line(k:k+4) = ' a+b '
              END DO
              WRITE(*,'(8X,A2,1I4,1A2,4X,1A,A1)') '| ', i, ' |', &
                line(1:c1*5), '|'
            END DO
          END IF
!        print out dependenCy between bC and bCtp
          CALL show_bcdep(ismpl)
        END IF

!  --------------------

        DEALLOCATE(unitlevel)
        WRITE(*,*)
        CALL alloc_inverse(ismpl)

        IF (runmode>0) THEN
!        memory managment (2)
          WRITE(*,*) "Allocating ad arrays"
#ifdef AD
          CALL g_alloc_arrays(ismpl)
#endif
#ifdef AD_RM
          CALL alloc_arrays_ad(ismpl)
#endif
        END IF

!     finish HDF5 support, when available
        CALL close_hdf5()
        CLOSE(79)
        RETURN

      END
