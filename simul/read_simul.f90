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

!> @brief read simulation parameters and allocate variables
!> @param[in] filename SIMUL parameter file name
!> @param[in] ismpl local sample index
!> @details
!> __ATTENTION__: IMPLICIT NONE was commented out before, if problems
!> in this subroutine arise, the reason could be that it needs to be
!> commented out again.
!>
!> `# subsample`: Set saa_sample, which has no obvious use. \n
!> `# parameter group`: Count `mpara`, set `max_gpara`. \n
!> `# bc group`: Count `mpara`, set `max_gpara`. \n
!>
!>  __Example__:
!>
!>        # parameter group, records=2
!>        4 parfile=visim_looms.par
!>          2 1
!>          3 1
!>          4 1
!>          5 1
!>        1 parfile=visim_looms2.par
!>          1 3
!>
!> Would lead to:
!>
!>      -> mpara = 4 + 1 = 5
!>      -> max_gpara = 2
!>
!> __Second loop__ `# parameter group`: Set `fnpara`, `gpara`, `seed_para`. \n
!> __Second loop__ `# bc group`: Set `fnpara`, `gpara`, `seed_para`. \n
!> The example would lead to
!>
!>       -> seed_para = [[1 1 1 1 3]
!>                       [2 3 4 5 1]]
!>       -> gpara=[1 1 1 1 2]
!>       -> fnpara=[visim_looms.par visim_looms2.par]
!>
!> __TODO__ under `# bc group`:
!>
!>            IF (found(79,key_char//' bc group',line,.FALSE.)) THEN
!>            aw-later        if (transient) then
!>            aw-later           write(*,'(1A)')
!>            aw-later     &   'error: bc optimization not allowed for time dependend models!'
!>            aw-later           write(*,'(2A)')
!>            aw-later     &      '  -> please, update this as tp optimization with the ',
!>            aw-later     &      'simulation begin as the start time.'
!>            aw-later           stop
!>            aw-later        endif
      SUBROUTINE read_simul(filename,ismpl)
        use arrays
        use simul_arrays
        use mod_genrl
        use mod_genrlc
        use mod_flow
        use mod_temp
        use mod_simul
        use mod_time
        use mod_data
        use mod_linfos
        IMPLICIT NONE
        character (len=80) ::filename
        character (len=5000) ::line
        character (len=4) :: ctmp
        character (len=4), dimension (max(nprop,nbc)) :: stmp
        DOUBLE PRECISION, ALLOCATABLE :: datmp(:,:)
        integer :: locstr, lblank, itmp, i2tmp, i3tmp, i4tmp, i5tmp, impara, &
          imax_gpara
        integer :: i, j, k, l,ident, ismpl
        DOUBLE PRECISION dtmp
        character (len=4), dimension (0:2) :: sll
        DATA sll/'    ', ' lin', ' log'/
        LOGICAL found, no_ext_link, no_ext_link_int
        EXTERNAL locstr, found, no_ext_link, no_ext_link_int, lblank
        integer :: param_enabled
        EXTERNAL param_enabled


!     open file
        OPEN(79,file=filename,status='old')

        WRITE(*,*)
        WRITE(*,*) '  reading simul parameter'
        WRITE(*,*)

!     init HDF5 support, when available
        CALL open_hdf5(' ')

!     sanity check
        IF (runmode>=2 .AND. nsmpl<=sm_max) THEN
          WRITE(*,'(1A)') 'error: for ENKF #samples needs &
            &to be greater than #simulate !'
          STOP
        END IF

        ! # subsample
        saa_sample = 1
        IF (found(79,key_char//' subsample',line,.FALSE.)) THEN
          WRITE(*,*) ' [R] : use sample average approximation'
          READ(79,*) saa_sample
        ELSE
          WRITE(*,*) ' <D> : no sample average approximation'
        END IF

!  --------------------

!     count 'mpara'
        mpara = 0
        impara = 0
        max_gpara = 0
        imax_gpara = 0

        ! # parameter group
        IF (found(79,key_char//' parameter group',line,.FALSE.)) THEN
          CALL get_arg('records',line,i,j)
          IF (i>0 .AND. j>=i) THEN
            READ(line(i:j),*) i4tmp
          ELSE
            READ(79,*) i4tmp
          END IF
          max_gpara = max_gpara + i4tmp
          DO j = 1, i4tmp
            READ(79,'(A)') line
            READ(line,*) itmp
            DO i = 1, itmp
              READ(79,'(A)') line
              mpara = mpara + 1
            END DO
          END DO
        END IF

        ! # bc group
        IF (found(79,key_char//' bc group',line,.FALSE.)) THEN
          CALL get_arg('records',line,i,j)
          IF (i>0 .AND. j>=i) THEN
            READ(line(i:j),*) i4tmp
          ELSE
            READ(79,*) i4tmp
          END IF
          max_gpara = max_gpara + i4tmp
          DO j = 1, i4tmp
            READ(79,'(A)') line
            READ(line,*) itmp
            DO i = 1, itmp
              READ(79,'(A)') line
              IF (locstr(line,'split')>=1) THEN
                DO k = 1, nbc_data
                  i2tmp = 0
                  IF (locstr(line,'head')>=1) i2tmp = pv_head
                  IF (locstr(line,'temp')>=1) i2tmp = pv_temp
                  IF (locstr(line,'conc')>=1) i2tmp = pv_conc
                  IF (locstr(line,'pres')>=1) i2tmp = pv_pres
                  IF (ibc_data(k,cbc_bcu)==itmp .AND. ibc_data(k,cbc_pv)==i2tmp) mpara = mpara + 1
                END DO
              ELSE
                mpara = mpara + 1
              END IF
            END DO
          END DO
        END IF

        ALLOCATE(fnpara(max(max_gpara,1)))
        ALLOCATE(gpara(max(mpara,1)))
        ALLOCATE(seed_para(2,max(mpara,1)))

        ! # parameter group
        DO j = 1, nunits
          DO i = 1, nprop
            opti_props(i,j) = 0
          END DO
        END DO
        IF (found(79,key_char//' parameter group',line,.FALSE.)) THEN
          CALL get_arg('records',line,i,j)
          IF (i>0 .AND. j>=i) THEN
            READ(line(i:j),*) i4tmp
          ELSE
             write(unit = *, fmt = *) '[E1] read_simul.f90: No records input for # parameter group'
          END IF
          WRITE(*,'(A,I3)') '  [R] : parameter groups, records=', &
            i4tmp
!       reading
          DO k = 1, i4tmp
            imax_gpara = imax_gpara + 1
            READ(79,'(A)') line
            READ(line,*) i5tmp
!          supported simulation functions
            CALL get_arg('parfile',line,i,j)
            READ(line(i:j),*) fnpara(imax_gpara)
            DO i = 1, i5tmp
              impara = impara + 1
              gpara(impara) = imax_gpara
              READ(79,'(1A)') line
              READ(line,*) seed_para(2,impara), seed_para(1,impara)
              IF (seed_para(1,impara)<firstidx .OR. &
                  seed_para(1,impara)>lastidx) THEN
                WRITE(*,'(1A,1I5,1A,1I5,1A)') &
                  'error: parameter index out of range, at group ', k, &
                  ', line ', i, '!'
                STOP
              END IF
!             default is linear (1)
              j = 1
              IF (locstr(line,'lin')>=1) j = 1
              IF (locstr(line,'log')>=1) j = 2
              opti_props(seed_para(1,impara),seed_para(2,impara)) = j
            END DO
          END DO
        ELSE
          WRITE(*,'(A,I3)') '  <D> : no parameter groups !'
        END IF

        ! # standard deviation
        IF (found(79,key_char//' standard deviation',line,.FALSE.)) THEN
          WRITE(*,*) ' [R] : standard deviation'
          ALLOCATE(datmp(nprop_load,maxunits))
          IF (no_ext_link(nprop_load,maxunits,1,datmp,'deviation', &
            line)) CALL read_array(79,nprop_load,maxunits,datmp, &
            key_char//' standard deviation',ismpl)
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
              'error: bug (1) in "read_simul.f", ask AW!'
            STOP
          END IF

          DO i = 1, maxunits
!          beCause of logarithimC sCale, suppress zeros
            DO j = firstidx, lastidx
              d_propunit(i,j) = max(1.0D-99,d_propunit(i,j))
            END DO
            IF (linfos(2)>=2) WRITE(*,'('//c_npropunit//'e12.4,1I8)') &
              (d_propunit(i,j),j=firstidx,lastidx), i
          END DO
          DEALLOCATE(datmp)
        ELSE
          WRITE(*,*) ' <D> : standard deviation !'
          DO i = 1, maxunits
            DO j = firstidx, lastidx
              d_propunit(i,j) = 0.5D0
            END DO
            IF (linfos(2)>=2) WRITE(*,'('//c_npropunit//'e12.4,1I8)') &
              (d_propunit(i,j),j=firstidx,lastidx), i
          END DO
        END IF

!  --------------------
!     init
        DO j = 1, bc_maxunits
          DO i = 1, nbc
            opti_bc(i,j) = 0
          END DO
        END DO
!     set
        ! # bc group
        IF (found(79,key_char//' bc group',line,.FALSE.)) THEN
          CALL get_arg('records',line,i,j)
          IF (i>0 .AND. j>=i) THEN
            READ(line(i:j),*) i4tmp
          ELSE
            READ(79,*) i4tmp
          END IF
          WRITE(*,'(A,I3)') '  [R] : bc groups, records=', i4tmp
!       reading
          DO k = 1, i4tmp
            imax_gpara = imax_gpara + 1
            READ(79,'(A)') line
            READ(line,*) l
!          supported simulation functions
            CALL get_arg('parfile',line,i,j)
            READ(line(i:j),*) fnpara(imax_gpara)
            DO i = 1, l
!             avaiting: [bc-unit  (lin|log)  (head|temp|conc)]
!                    == [unit-id   lin/log    pv-type        ]
              READ(79,'(1A)') line
              READ(line,*) itmp
              IF (itmp>bc_maxunits) THEN
                WRITE(*,'(A)') &
                  'error: "bc group" unit number out of range !!!'
                STOP
              END IF
              i2tmp = 1
              IF (locstr(line,'lin')>=1) i2tmp = 1
              IF (locstr(line,'log')>=1) i2tmp = 2
              i3tmp = 0
              IF (locstr(line,'head')>=1) i3tmp = pv_head
              IF (locstr(line,'temp')>=1) i3tmp = pv_temp
              IF (locstr(line,'conc')>=1) i3tmp = pv_conc
              IF (locstr(line,'pres')>=1) i3tmp = pv_pres
              IF (i3tmp==0) THEN
                WRITE(*,'(1A,1I5,1A,1I5,1A)') 'error: "bc &
                  &group" phys. value type not allowed, at &
                  &group ', k, ', line ', i, '!'
                STOP
              END IF
              opti_bc(i3tmp,itmp) = i2tmp

              IF (locstr(line,'split')>=1) THEN
                WRITE(*,'(1A,1I7)') '      splitting bc unit:', itmp
                DO j = 1, nbc_data
                  IF (ibc_data(j,cbc_bcu)==itmp .AND. &
                      ibc_data(j,cbc_pv)==i3tmp) THEN
                    impara = impara + 1
                    gpara(impara) = imax_gpara
                    seed_para(1,impara) = -2
                    seed_para(2,impara) = j
!                      copy bc-unit value into each bc-cell
                    dbc_data(j,1,ismpl) = propunit(itmp, &
                      bc_firstidx-1+i2tmp,ismpl)
                  END IF
                END DO
              ELSE
                impara = impara + 1
                gpara(impara) = imax_gpara
                seed_para(1,impara) = bc_firstidx - 1 + i3tmp
                seed_para(2,impara) = itmp
              END IF
            END DO
          END DO
        END IF

        ! # bc standard deviation
        IF (found(79,key_char//' bc standard deviation',line,.FALSE.)) THEN
          CALL get_arg('records',line,i,j)
          IF (i<1 .OR. j<i) THEN
            READ(79,*) itmp
          ELSE
            READ(line(i:j),*) itmp
          END IF

          DO i = 1, itmp
            READ(79,*) k, dtmp, ctmp
            IF ((k>bc_maxunits) .OR. (k<1)) THEN
              WRITE(*,'(2A,1I7,1A,1I7,1A)') 'error: bc standard &
                &deviation, unit number out of range or ', &
                'not used, (', k, ') at line ', i, '!'
              STOP
            END IF
            IF ((ctmp/='head') .AND. (ctmp/='temp') .AND. &
                (ctmp/='pres') .AND. (ctmp/='conc')) THEN
              WRITE(*,'(1A,1A1,1A,1I3,1A)') 'error: bc standard &
                &deviation, unit type not allowed, "', ctmp, &
                '" at line ', i, '!'
              STOP
            END IF

            IF (ctmp=='head') d_propunit(k,idx_hbc) = max(1.D-10,dtmp)
            IF (ctmp=='temp') d_propunit(k,idx_tbc) = max(1.D-10,dtmp)
            IF (ctmp=='conc') d_propunit(k,idx_cbc) = max(1.D-10,dtmp)
            IF (ctmp=='pres') d_propunit(k,idx_hbc) = max(1.D-10, dtmp) ! idx_hbc used for pressure
          END DO
          WRITE(*,'(A,I3)') '  [R] : bc standard deviation, records=' &
            , itmp
        ELSE
          WRITE(*,*) ' <D> : bc standard deviation !'
          DO i = 1, bc_maxunits
            DO j = bc_firstidx, bc_lastidx
              d_propunit(i,j) = 1.0D0
            END DO
          END DO
        END IF

!  --------------------
        IF (impara/=mpara) THEN
          WRITE(*,'(A)') 'error: "mpara" differs in "read_simul" !!!'
          STOP 1
        END IF

!  --------------------
!     print activity table
        IF (linfos(2)>=1) THEN
!        parameter units
          IF (maxunits>=1) THEN
            WRITE(*,*) ' '
            WRITE(*,'(6X,A)') &
              'activity matrix (param.-unit - group):'
            WRITE(*,'(8X,A2,A7,1X,'//c_npropunit//'(A4,1X),A1)') '| ', '   unit', &
              (properties(i),i=firstidx,lastidx), '|'
!           unit-loop
            DO j = 1, maxunits
              DO k = 1, lastidx - firstidx + 1
                stmp(k) = ' '
              END DO
!              activity list
              DO k = 1, mpara
!                 correct unit number?
                IF (seed_para(2,k)==j .AND. seed_para(1,k)>=firstidx &
                  .AND. seed_para(1,k)<=lastidx) &
                  WRITE(stmp(seed_para(1,k)-firstidx+1),'(I4)') gpara &
                  (k)
              END DO
              WRITE(*,'(8X,A2,I7,1X,13(A4,1X),A1)') '| ', j, &
                (stmp(i),i=1,nprop_load), '|'
            END DO
            WRITE(*,*) ' '
            WRITE(*,'(6X,A)') &
              'activity matrix (param.-unit - lin/log):'
            WRITE(*,'(8X,A2,A7,1X,'//c_npropunit//'(A4,1X),A1)') '| ', '   unit', &
              (properties(i),i=firstidx,lastidx), '|'
!           unit-loop
            DO j = 1, maxunits
              DO k = 1, lastidx - firstidx + 1
                stmp(k) = ' '
              END DO
!              activity list
              DO k = 1, mpara
!                 correct unit number?
                IF (seed_para(2,k)==j .AND. seed_para(1,k)>=firstidx &
                  .AND. seed_para(1,k)<=lastidx) &
                  WRITE(stmp(seed_para(1,k)-firstidx+1),'(A4)') sll( &
                  param_enabled(k,ismpl))
              END DO
              WRITE(*,'(8X,A2,I7,1X,13(A4,1X),A1)') '| ', j, &
                (stmp(i),i=1,nprop_load), '|'
            END DO
          END IF
!  --------------------
!        bc units
          IF (bc_maxunits>=1) THEN
            WRITE(*,*) ' '
            WRITE(*,'(6X,A)') 'activity matrix (bc-unit - group):'
            WRITE(*,'(8X,1A2,1A6,1X,'//c_nbcunit//'(A4,1X),1A1)') '| ', 'BCunit', &
              'Flow', 'Temp', 'Conc', '|'
            DO j = 1, bc_maxunits
              DO k = 1, bc_lastidx - bc_firstidx + 1
                stmp(k) = ' '
              END DO
!              activity list
              DO k = 1, mpara
!                 correct unit number?
                IF (seed_para(2,k)==j .AND. seed_para(1,k)>= &
                  bc_firstidx .AND. seed_para(1,k)<=bc_lastidx) &
                  WRITE(stmp(seed_para(1,k)-bc_firstidx+1),'(I4)') &
                  gpara(k)
                IF (seed_para(1,k)==-2) THEN
                  IF (ibc_data(seed_para(2,k),cbc_bcu)==j) WRITE(stmp &
                    (ibc_data(seed_para(2,k),cbc_pv)),'(I4)') gpara(k)
                END IF
              END DO
              WRITE(*,'(8X,1A2,1I6,1X,'//c_nbcunit//'(A4,1X),1A1)') '| ', j, &
                (stmp(i),i=1,nbc), '|'
            END DO
            WRITE(*,*) ' '
            WRITE(*,'(6X,A)') 'activity matrix (bc-unit - lin/log):'
            WRITE(*,'(8X,1A2,1A6,1X,'//c_nbcunit//'(A4,1X),1A1)') '| ', 'BCunit', &
              'Flow', 'Temp', 'Conc', '|'
            DO j = 1, bc_maxunits
              DO k = 1, bc_lastidx - bc_firstidx + 1
                stmp(k) = ' '
              END DO
!              activity list
              DO k = 1, mpara
!                 correct unit number?
                IF (seed_para(2,k)==j .AND. seed_para(1,k)>= &
                  bc_firstidx .AND. seed_para(1,k)<=bc_lastidx) &
                  WRITE(stmp(seed_para(1,k)-bc_firstidx+1),'(A4)') &
                  sll(param_enabled(k,ismpl))
                IF (seed_para(1,k)==-2) THEN
                  IF (ibc_data(seed_para(2,k),cbc_bcu)==j) WRITE(stmp &
                    (ibc_data(seed_para(2,k),cbc_pv)),'(A4)') sll( &
                    param_enabled(k,ismpl))
                END IF
              END DO
              WRITE(*,'(8X,1A2,1I6,1X,'//c_nbcunit//'(A4,1X),1A1)') '| ', j, &
                (stmp(i),i=1,nbc), '|'
            END DO
          END IF
        END IF
!  --------------------

        WRITE(*,*)
        CALL alloc_simul(ismpl)

!     finish HDF5 support, when available
        CALL close_hdf5()
        CLOSE(79)

          call enkf_make_output_dirs(.true.)

        RETURN
      END
