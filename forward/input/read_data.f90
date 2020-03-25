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

!>    @brief read observed data and allocate fields
!>    @param[in] filename file name
!>    @param[in] ismpl local sample index
!>    @details
!> details see in documentation "read_observed_data.pdf"\n
      SUBROUTINE read_data(filename,ismpl)
        use arrays
        use mod_genrl
        use mod_genrlc
        use mod_flow
        use mod_temp
        use mod_time
        use mod_conc
        use mod_data
        use mod_linfos
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k, l
        character (len=80) :: filename
        character (len=80) :: line
        INTEGER locstr, lblank, type, ozone, i2, i3, ll, &
          ndata_sections, tmplen
        ! INTEGER c_i, c_j, c_k
        INTEGER level_timer, i_si, i_obs
        DOUBLE PRECISION d_timer
        LOGICAL found, no_ext_link, no_ext_link_int
        ! LOGICAL incomp_bc
        LOGICAL read_species, read_obs, read_absolute
        INTEGER, ALLOCATABLE :: tmp_idata(:,:)
        DOUBLE PRECISION, ALLOCATABLE :: tmp_data(:,:)
        EXTERNAL locstr, found, no_ext_link, no_ext_link_int, lblank


!     read file
        OPEN(79,file=filename,status='old')

!     init HDF5 support, when available
        CALL open_hdf5(' ')

        WRITE(*,*)
        WRITE(*,*) '  reading observed data'
        WRITE(*,*) '    from file "', filename(:lblank(filename)), &
          '"'
        WRITE(*,*)

!     init counter
        ndata_h = 0
        ndata_t = 0
        ndata_c = 0
        ndata_p = 0
        ndata_s = 0
        ndata_b = 0

!     count the number of observed data entries (sections and entries)
        REWIND 79
        ndata_sections = 0
        ndata = 0

10      READ(79,'(1A)',end=11) line
        i = locstr(line,key_char//' data')
!       data line?
        IF (i==1) THEN
!         increase number of sections
          ndata_sections = ndata_sections + 1
          CALL get_arg('records',line,i,j)
          IF (i<1 .OR. j<i) THEN
            READ(79,*) tmplen
          ELSE
            READ(line(i:j),*) tmplen
          END IF
!         increase number of entries
          ndata = ndata + tmplen
        END IF
!     read next line, up to the end of file
        GO TO 10

!     restart file
11      REWIND 79

        ALLOCATE(tmp_idata(max(ndata,1),n_idata))
        ALLOCATE(tmp_data(max(ndata,1),n_ddata))
        ll = 0

        DO i3 = 1, ndata_sections
          IF (found(79,key_char//' data',line,.TRUE.)) THEN
            CALL get_arg('records',line,i,j)
            IF (i<1 .OR. j<i) THEN
              READ(79,*) tmplen
            ELSE
              READ(line(i:j),*) tmplen
            END IF

!        time dependent data, if switch used, read one column more
            d_timer = max_simtime/tunit
!        0: do not read timers
!        1: read timers
!        2: use the given timer for all records
            level_timer = 0
            CALL get_arg('timer',line,i,j)
            IF (i>=1 .AND. j>=i) THEN
!           prove "r" instead of "read"
              IF (line(i:i)=='r') THEN
                level_timer = 1
              ELSE
                level_timer = 2
                READ(line(i:j),*,err=901) d_timer
!             write(*,'(1A,1e16.8)') '        timer=',d_timer
              END IF
            END IF

!        species for concentration
            i_si = 0
            read_species = .TRUE.
            CALL get_arg('species',line,i,j)
            IF (i>=1 .AND. j>=i) THEN
              READ(line(i:j),*) i_si
              read_species = .FALSE.
!           write(*,'(1A,1I3)') '      speCies=',i_si
            END IF

!        position type
            read_absolute = .FALSE.
            CALL get_arg('pos',line,i,j)
            IF (i>=1 .AND. j>=i) THEN
!           if (line(i:i).eq.'index') read_absolute = .false.
!           instead of prove "abs" or "index"
              IF (line(i:i)=='a') read_absolute = .TRUE.
!           write(*,*) '         position absolute=',read_absolute
            END IF

!        observation point index
            i_obs = 0
            read_obs = .TRUE.
            CALL get_arg('obs',line,i,j)
            IF (i>=1 .AND. j>=i) THEN
              READ(line(i:j),*) i_obs
              read_obs = .FALSE.
!           write(*,'(1A,1I3)') '          obs=',i_obs
            END IF

            WRITE(*,'(A,I6)') '  [R] : observed data nodes, records=' &
              , tmplen
!        sanity check
            IF (level_timer>=1 .AND. .NOT. transient) THEN
              WRITE(*,'(2A)') &
                'error: "timer" is specified for the data', &
                ' input, but the model is steady state !'
              STOP
            END IF

            IF (no_ext_link_int(tmplen,n_idata,1,tmp_idata,'idata', &
                line) .OR. no_ext_link(tmplen,n_ddata,1,tmp_data, &
                'data',line)) THEN
              DO j = 1, tmplen
!              record counter
                ll = ll + 1
                READ(79,'(1A)') line
                tmp_data(ll,cdd_time) = d_timer
                tmp_idata(ll,cid_si) = i_si
                tmp_idata(ll,cid_obs) = i_obs

!           ----------
                IF ( .NOT. read_absolute) THEN
!              index position
                  IF (level_timer==1 .AND. read_obs) THEN
                    READ(line,*,err=1001,end=1001) tmp_data(ll, &
                      cdd_pv), tmp_data(ll,cdd_w), &
                      tmp_data(ll,cdd_time), (tmp_idata(ll,k),k=cid_i, &
                      cid_k), tmp_idata(ll,cid_pv), &
                      tmp_idata(ll,cid_obs)
!                 re-read with sub-index (species)
                    IF (tmp_idata(ll,cid_pv)==pv_conc .AND. &
                      read_species) READ(line,*,err=1002,end=1002) &
                      tmp_data(ll,cdd_pv), tmp_data(ll,cdd_w), &
                      tmp_data(ll,cdd_time), (tmp_idata(ll,k),k=cid_i, &
                      cid_k), tmp_idata(ll,cid_pv), &
                      tmp_idata(ll,cid_si), tmp_idata(ll,cid_obs)
                  ELSE IF (level_timer==1 .AND. .NOT. read_obs) THEN
                    READ(line,*,err=1003,end=1003) tmp_data(ll, &
                      cdd_pv), tmp_data(ll,cdd_w), &
                      tmp_data(ll,cdd_time), (tmp_idata(ll,k),k=cid_i, &
                      cid_k), tmp_idata(ll,cid_pv)
!                 re-read with sub-index (species)
                    IF (tmp_idata(ll,cid_pv)==pv_conc .AND. &
                      read_species) READ(line,*,err=1004,end=1004) &
                      tmp_data(ll,cdd_pv), tmp_data(ll,cdd_w), &
                      tmp_data(ll,cdd_time), (tmp_idata(ll,k),k=cid_i, &
                      cid_k), tmp_idata(ll,cid_pv), &
                      tmp_idata(ll,cid_si)
                  ELSE IF (level_timer/=1 .AND. read_obs) THEN
                    READ(line,*,err=1005,end=1005) tmp_data(ll, &
                      cdd_pv), tmp_data(ll,cdd_w), &
                      (tmp_idata(ll,k),k=cid_i,cid_k), &
                      tmp_idata(ll,cid_pv), tmp_idata(ll,cid_obs)
!                 re-read with sub-index (species)
                    IF (tmp_idata(ll,cid_pv)==pv_conc .AND. &
                      read_species) READ(line,*,err=1006,end=1006) &
                      tmp_data(ll,cdd_pv), tmp_data(ll,cdd_w), &
                      (tmp_idata(ll,k),k=cid_i,cid_k), &
                      tmp_idata(ll,cid_pv), tmp_idata(ll,cid_si), &
                      tmp_idata(ll,cid_obs)
                  ELSE IF (level_timer/=1 .AND. .NOT. read_obs) THEN
                    READ(line,*,err=1007,end=1007) tmp_data(ll, &
                      cdd_pv), tmp_data(ll,cdd_w), &
                      (tmp_idata(ll,k),k=cid_i,cid_k), &
                      tmp_idata(ll,cid_pv)
!                 re-read with sub-index (species)
                    IF (tmp_idata(ll,cid_pv)==pv_conc .AND. &
                      read_species) READ(line,*,err=1008,end=1008) &
                      tmp_data(ll,cdd_pv), tmp_data(ll,cdd_w), &
                      (tmp_idata(ll,k),k=cid_i,cid_k), &
                      tmp_idata(ll,cid_pv), tmp_idata(ll,cid_si)
                  END IF
!           ----------
                ELSE
!           ----------
!             absolute position
                  IF (level_timer==1 .AND. read_obs) THEN
                    READ(line,*,err=2001,end=2001) tmp_data(ll, &
                      cdd_pv), tmp_data(ll,cdd_w), &
                      tmp_data(ll,cdd_time), (tmp_data(ll,k),k=cdd_i, &
                      cdd_k), tmp_idata(ll,cid_pv), &
                      tmp_idata(ll,cid_obs)
!                 re-read with sub-index (species)
                    IF (tmp_idata(ll,cid_pv)==pv_conc .AND. &
                      read_species) READ(line,*,err=2002,end=2002) &
                      tmp_data(ll,cdd_pv), tmp_data(ll,cdd_w), &
                      tmp_data(ll,cdd_time), (tmp_data(ll,k),k=cdd_i, &
                      cdd_k), tmp_idata(ll,cid_pv), &
                      tmp_idata(ll,cid_si), tmp_idata(ll,cid_obs)
                  ELSE IF (level_timer==1 .AND. .NOT. read_obs) THEN
                    READ(line,*,err=2003,end=2003) tmp_data(ll, &
                      cdd_pv), tmp_data(ll,cdd_w), &
                      tmp_data(ll,cdd_time), (tmp_data(ll,k),k=cdd_i, &
                      cdd_k), tmp_idata(ll,cid_pv)
!                 re-read with sub-index (species)
                    IF (tmp_idata(ll,cid_pv)==pv_conc .AND. &
                      read_species) READ(line,*,err=2004,end=2004) &
                      tmp_data(ll,cdd_pv), tmp_data(ll,cdd_w), &
                      tmp_data(ll,cdd_time), (tmp_data(ll,k),k=cdd_i, &
                      cdd_k), tmp_idata(ll,cid_pv), &
                      tmp_idata(ll,cid_si)
                  ELSE IF (level_timer/=1 .AND. read_obs) THEN
                    READ(line,*,err=2005,end=2005) tmp_data(ll, &
                      cdd_pv), tmp_data(ll,cdd_w), &
                      (tmp_data(ll,k),k=cdd_i,cdd_k), &
                      tmp_idata(ll,cid_pv), tmp_idata(ll,cid_obs)
!                 re-read with sub-index (species)
                    IF (tmp_idata(ll,cid_pv)==pv_conc .AND. &
                      read_species) READ(line,*,err=2006,end=2006) &
                      tmp_data(ll,cdd_pv), tmp_data(ll,cdd_w), &
                      (tmp_data(ll,k),k=cdd_i,cdd_k), &
                      tmp_idata(ll,cid_pv), tmp_idata(ll,cid_si), &
                      tmp_idata(ll,cid_obs)
                  ELSE IF (level_timer/=1 .AND. .NOT. read_obs) THEN
                    READ(line,*,err=2007,end=2007) tmp_data(ll, &
                      cdd_pv), tmp_data(ll,cdd_w), &
                      (tmp_data(ll,k),k=cdd_i,cdd_k), &
                      tmp_idata(ll,cid_pv)
!                 re-read with sub-index (species)
                    IF (tmp_idata(ll,cid_pv)==pv_conc .AND. &
                      read_species) READ(line,*,err=2008,end=2008) &
                      tmp_data(ll,cdd_pv), tmp_data(ll,cdd_w), &
                      (tmp_data(ll,k),k=cdd_i,cdd_k), &
                      tmp_idata(ll,cid_pv), tmp_idata(ll,cid_si)
                  END IF
                END IF
!           ----------

                IF (read_absolute) THEN
!                 search index values for x,y,z
                  tmp_idata(ll,cid_i) = 0
                  DO k = 1, i0
                    IF (delxa(k)-0.5D0*delx(k)<=tmp_data(ll,cdd_i)) &
                      tmp_idata(ll,cid_i) = k
                  END DO
                  tmp_idata(ll,cid_j) = 0
                  DO k = 1, j0
                    IF (delya(k)-0.5D0*dely(k)<=tmp_data(ll,cdd_j)) &
                      tmp_idata(ll,cid_j) = k
                  END DO
                  tmp_idata(ll,cid_k) = 0
                  DO k = 1, k0
                    IF (delza(k)-0.5D0*delz(k)<=tmp_data(ll,cdd_k)) &
                      tmp_idata(ll,cid_k) = k
                  END DO
                ELSE
!                 convert index i,j,k into absolute postion
                  tmp_data(ll,cdd_i) = delxa(tmp_idata(ll,cid_i))
                  tmp_data(ll,cdd_j) = delya(tmp_idata(ll,cid_j))
                  tmp_data(ll,cdd_k) = delza(tmp_idata(ll,cid_k))
                END IF

!              modify with 'tunit'
                tmp_data(ll,cdd_time) = tmp_data(ll,cdd_time)*tunit
!              sanity checks
                IF (tmp_data(ll,cdd_time)<simtime_0 .OR. &
                    tmp_data(ll,cdd_time)>max_simtime) THEN
                  WRITE(*,'(1A,1I6,1A)') &
                    'error: timer value out of range, at line ', j, '!'
                  WRITE(*,'(3A)') '  given:"', line, '"'
                  STOP
                END IF
                IF (tmp_data(ll,cdd_w)<dabs(1.0d-10*tmp_data(ll,cdd_pv))) THEN
                  tmp_data(ll,cdd_w)=max(dabs(1.0d-10*tmp_data(ll,cdd_pv)),1.0d-99)
                  WRITE(*,'(1A,1I6,1A,1G12.4,1A)') &
                    'warning: given error seems to be to small, at line ', j, &
                    ', cutting them to ',tmp_data(ll,cdd_w),'!'
                END IF
                IF (tmp_idata(ll,cid_i)<1 .OR. &
                    tmp_idata(ll,cid_i)>i0 .OR. &
                    tmp_idata(ll,cid_j)<1 .OR. &
                    tmp_idata(ll,cid_j)>j0 .OR. &
                    tmp_idata(ll,cid_k)<1 .OR. tmp_idata(ll,cid_k)>k0) &
                    THEN
                  WRITE(*,'(1A,1I6,1A)') &
                    'error: index out of range, at line ', j, '!'
                  WRITE(*,'(3A)') '  given:"', line, '"'
                  STOP
                END IF
                IF (tmp_data(ll,cdd_i)<0.0D0 .OR. &
                    tmp_data(ll,cdd_i)>delxa(i0)+0.5D0*delx(i0) .OR. &
                    tmp_data(ll,cdd_j)<0.0D0 .OR. &
                    tmp_data(ll,cdd_j)>delya(j0)+0.5D0*dely(j0) .OR. &
                    tmp_data(ll,cdd_k)<0.0D0 .OR. &
                    tmp_data(ll,cdd_k)>delza(k0)+0.5D0*delz(k0)) THEN
                  WRITE(*,'(1A,1I6,1A)') &
                    'error: absolute position out of range, at line ', &
                    j, '!'
                  WRITE(*,'(3A)') '  given:"', line, '"'
                  STOP
                END IF
                IF (tmp_idata(ll,cid_pv)<1 .OR. &
                    tmp_idata(ll,cid_pv)>npv) THEN
                  WRITE(*,'(1A,1I1,1A,1I6,1A)') &
                    'error: physical value index out of range [1..', &
                    npv, '], at line ', j, '!'
                  WRITE(*,'(3A)') '  given:"', line, '"'
                  STOP
                END IF
                IF (tmp_idata(ll,cid_si)<0 .OR. &
                    tmp_idata(ll,cid_si)>ntrans) THEN
                  WRITE(*,'(1A,1I6,1A)') &
                    'error: species index out of range, at line ', j, &
                    '!'
                  WRITE(*,'(3A)') '  given:"', line, '"'
                  STOP
                END IF
              END DO
            END IF
          END IF
        END DO
!     sanity check
        IF (ll/=ndata) THEN
          WRITE(*,'(1A)') 'error: lost some data in "read_data.f" !'
          WRITE(*,*) ' searching:', ndata, ', but reading:', ll
          STOP
        END IF

        IF (ndata>=1) THEN
!        compare data and boundaries, eleminate one of them when on the same position
          k = 0
          DO j = 1, ndata
            DO i = 1, nbc_data
!            only dirichlet boundaries
              IF ((tmp_idata(j,cid_i)==ibc_data(i, &
                  cbc_i)) .AND. (tmp_idata(j,cid_j)==ibc_data(i, &
                  cbc_j)) .AND. (tmp_idata(j,cid_k)==ibc_data(i, &
                  cbc_k)) .AND. (bt_diri==ibc_data(i, &
                  cbc_bt)) .AND. (tmp_idata(j,cid_pv)==ibc_data(i, &
                  cbc_pv))) THEN
!              disable them, later delete/ignore
#ifndef DbB
!              BbD, Boundaries before Data
                tmp_idata(j,cid_pv) = 0
#else
!              DbB, Data before Boundaries
                ibc_data(i,cbc_pv) = -ibc_data(i,cbc_pv)
#endif
                k = k + 1
              ELSE
!              dont change, leave them
              END IF
            END DO
          END DO

          ll = 0
!        compare data with data, eleminate one of them when on the same position
          DO j = 1, ndata
            DO i = 1, ndata
!            only dirichlet boundaries
              IF ((tmp_idata(i,cid_i)==tmp_idata(j, &
                  cid_i)) .AND. (tmp_idata(i,cid_j)==tmp_idata(j, &
                  cid_j)) .AND. (tmp_idata(i,cid_k)==tmp_idata(j, &
                  cid_k)) .AND. (tmp_idata(i,cid_pv)==tmp_idata(j, &
                  cid_pv)) .AND. (tmp_idata(i,cid_si)==tmp_idata(j, &
                  cid_si)) .AND. (tmp_data(i,cdd_time)==tmp_data(j, &
                  cdd_time)) .AND. (i/=j) .AND. &
                  tmp_idata(j,cid_pv)/=0) THEN
!              disable them, later delete/ignore
                tmp_idata(j,cid_pv) = 0
                ll = ll + 1
!              skip this i-loop
                GO TO 30
              END IF
            END DO
30          CONTINUE
          END DO

#ifdef DbB
          IF (k>0) WRITE(*,'(A,I6,2A)') 'warning: ', k, ' &
            &boundary points have been deleted, because of &
            &conflicting data'
!        rearange full boundary arrays
!          j - "copy from" index
!          i - "copy to" index
          j = 1
          DO i = 1, nbc_data - k
100         IF (j<=nbc_data) THEN
              IF (ibc_data(j,cbc_pv)<0) THEN
!              delete boundary property "j"
!              and increase "j" to the next "to be copy" position
                j = j + 1
                GO TO 100
              ELSE
!              fill entry for boundary point "i", copy "j" to "i" (now dense)
                DO i2 = 1, nibc
                  ibc_data(i,i2) = ibc_data(j,i2)
                END DO
                DO i2 = 1, ndbc
                  dbc_data(i,i2,ismpl) = dbc_data(j,i2,ismpl)
                END DO
                j = j + 1
              END IF
            END IF
          END DO
          nbc_data = nbc_data - k
          k = 0
!        reinit first/last_<pv>
          CALL sort_bc(ismpl)
#endif
          IF (k>0) WRITE(*,'(A,I6,2A)') 'warning: ', k, ' &
            &data have been deleted, because of conflicting &
            &boundaries'
          IF (ll>0) WRITE(*,'(A,I6,2A)') 'warning: ', ll, &
            ' data have been deleted, because of conflicting data'
          ll = ll + k
          ndata = ndata - ll

          CALL alloc_data(ismpl)

          j = 0
!        ***  copy ordered !!!  ***
          DO l = 1, npv
!          l=1: copy-in HEAD elements
!            2: copy-in TEMP elements
!            3: copy-in CONC elements
!            5: copy-in PRES elements
!            6: copy-in BHPR elements
            DO i = 1, ndata + ll
              IF (tmp_idata(i,cid_pv)==l) THEN
                j = j + 1
                DO i2 = 1, n_ddata
                  ddata(j,i2) = tmp_data(i,i2)
                END DO
                DO i2 = 1, n_idata
                  idata(j,i2) = tmp_idata(i,i2)
                END DO
              END IF
            END DO
          END DO
!        sanity check
          IF (j/=ndata) THEN
            WRITE(*,'(1A,1I5,1A,1I5,1A)') 'error: lost some DATA,', &
              j, '/', ndata, ', in "read_data"!!!'
            STOP
          END IF

!        need "ndata_p/t/c/e" before "alloc_inverse"
          ndata_h = 0
          ndata_t = 0
          ndata_c = 0
          ndata_p = 0
          ndata_s = 0
          ndata_b = 0
          DO l = 1, ndata
            i = idata(l,cid_i)
            j = idata(l,cid_j)
            k = idata(l,cid_k)
            type = idata(l,cid_pv)
            ozone = idata(l,cid_obs)
            IF (type==pv_head) ndata_h = ndata_h + 1
            IF (type==pv_pres) ndata_p = ndata_p + 1
            IF (type==pv_temp) ndata_t = ndata_t + 1
            IF (type==pv_conc) ndata_c = ndata_c + 1
            IF (type==pv_bhpr) ndata_b = ndata_b + 1
          END DO
!        sanity check
          IF (ndata_h+ndata_t+ndata_c+ndata_p+ndata_s+ndata_b/=ndata) THEN
            WRITE(*,'(1A,1I5,1A,1I5,1A)') 'error: lost some DATA,', &
              ndata_h + ndata_t + ndata_c + ndata_p + ndata_s + ndata_b, '/', &
              ndata, ', in "read_data"(2)!!!'
            STOP
          END IF
          WRITE(*,'(1A,1I7)') &
            '  [I] : data and data node adresses, records=', ndata
        ELSE
          WRITE(*,*) ' <D> : no observed data found, records=0'
        END IF

        DEALLOCATE(tmp_idata)
        DEALLOCATE(tmp_data)

!     finish HDF5 support, when available
        CALL close_hdf5()
        CLOSE(79)

        RETURN

!     error handler
901     WRITE(*,'(2A)') 'error: expecting "read" or a timer value', &
          ' behind the keyword "timer=" !!!'
        STOP
1001    WRITE(*,'(2A,1I6,1A)') &
          'error: expecting [value error timer i j k value-type obs', &
          '] at line ', j, ' !!!'
        WRITE(*,'(1A)') &
          '  cell position index: i,j,k are integer values'
        WRITE(*,'(3A)') '  given:"', line, '"'
        STOP
1002    WRITE(*,'(2A,1I6,1A)') 'error: expecting [value error &
          &timer i j k value-type species', ' obs] at line ', j, &
          ' !!!'
        WRITE(*,'(1A)') &
          '  cell position index: i,j,k are integer values'
        WRITE(*,'(3A)') '  given:"', line, '"'
        STOP
1003    WRITE(*,'(2A,1I6,1A)') &
          'error: expecting [value error timer i j k value-type', &
          '] at line ', j, ' !!!'
        WRITE(*,'(1A)') &
          '  cell position index: i,j,k are integer values'
        WRITE(*,'(3A)') '  given:"', line, '"'
        STOP
1004    WRITE(*,'(2A,1I6,1A)') 'error: expecting [value error &
          &timer i j k value-type species', '] at line ', j, ' !!!'
        WRITE(*,'(1A)') &
          '  cell position index: i,j,k are integer values'
        WRITE(*,'(3A)') '  given:"', line, '"'
        STOP
1005    WRITE(*,'(2A,1I6,1A)') &
          'error: expecting [value error i j k value-type obs', &
          '] at line ', j, ' !!!'
        WRITE(*,'(1A)') &
          '  cell position index: i,j,k are integer values'
        WRITE(*,'(3A)') '  given:"', line, '"'
        STOP
1006    WRITE(*,'(2A,1I6,1A)') &
          'error: expecting [value error i j k value-type species obs', &
          '] at line ', j, ' !!!'
        WRITE(*,'(1A)') &
          '  cell position index: i,j,k are integer values'
        WRITE(*,'(3A)') '  given:"', line, '"'
        STOP
1007    WRITE(*,'(2A,1I6,1A)') &
          'error: expecting [value error i j k value-type', &
          '] at line ', j, ' !!!'
        WRITE(*,'(1A)') &
          '  cell position index: i,j,k are integer values'
        WRITE(*,'(3A)') '  given:"', line, '"'
        STOP
1008    WRITE(*,'(2A,1I6,1A)') &
          'error: expecting [value error i j k value-type species', &
          '] at line ', j, ' !!!'
        WRITE(*,'(1A)') &
          '  cell position index: i,j,k are integer values'
        WRITE(*,'(3A)') '  given:"', line, '"'
        STOP
2001    WRITE(*,'(2A,1I6,1A)') &
          'error: expecting [value error timer x y z value-type obs', &
          '] at line ', j, ' !!!'
        WRITE(*,'(1A)') &
          '  absolute cell position: x,y,z are floating point values'
        WRITE(*,'(3A)') '  given:"', line, '"'
        STOP
2002    WRITE(*,'(2A,1I6,1A)') 'error: expecting [value error &
          &timer x y z value-type species', ' obs] at line ', j, &
          ' !!!'
        WRITE(*,'(1A)') &
          '  absolute cell position: x,y,z are floating point values'
        WRITE(*,'(3A)') '  given:"', line, '"'
        STOP
2003    WRITE(*,'(2A,1I6,1A)') &
          'error: expecting [value error timer x y z value-type', &
          '] at line ', j, ' !!!'
        WRITE(*,'(1A)') &
          '  absolute cell position: x,y,z are floating point values'
        WRITE(*,'(3A)') '  given:"', line, '"'
        STOP
2004    WRITE(*,'(2A,1I6,1A)') 'error: expecting [value error &
          &timer x y z value-type species', '] at line ', j, ' !!!'
        WRITE(*,'(1A)') &
          '  absolute cell position: x,y,z are floating point values'
        WRITE(*,'(3A)') '  given:"', line, '"'
        STOP
2005    WRITE(*,'(2A,1I6,1A)') &
          'error: expecting [value error x y z value-type obs', &
          '] at line ', j, ' !!!'
        WRITE(*,'(1A)') &
          '  absolute cell position: x,y,z are floating point values'
        WRITE(*,'(3A)') '  given:"', line, '"'
        STOP
2006    WRITE(*,'(2A,1I6,1A)') &
          'error: expecting [value error x y z value-type species obs', &
          '] at line ', j, ' !!!'
        WRITE(*,'(1A)') &
          '  absolute cell position: x,y,z are floating point values'
        WRITE(*,'(3A)') '  given:"', line, '"'
        STOP
2007    WRITE(*,'(2A,1I6,1A)') &
          'error: expecting [value error x y z value-type', &
          '] at line ', j, ' !!!'
        WRITE(*,'(1A)') &
          '  absolute cell position: x,y,z are floating point values'
        WRITE(*,'(3A)') '  given:"', line, '"'
        STOP
2008    WRITE(*,'(2A,1I6,1A)') &
          'error: expecting [value error x y z value-type species', &
          '] at line ', j, ' !!!'
        WRITE(*,'(1A)') &
          '  absolute cell position: x,y,z are floating point values'
        WRITE(*,'(3A)') '  given:"', line, '"'
        STOP
      END
