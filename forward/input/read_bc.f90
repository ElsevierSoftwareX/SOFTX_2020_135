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

!>    @brief read the boundary condition specifications
!>    @param[in] filen number of the opened file
!>    @param[in] line current character line
!>    @param[in] i_pv physical value / state variable - index
!>    @param[in] i_bt boundary type -index
!>    @param[in,out] posi number (position of the last one) of boundary points - before and after reading
!>    @param[out] ilost number of lost boundary points (reading ignored)
!>    @param[in] ismpl local sample index
!>    @details
!> read model boundary points\n
!> \n
!> Note: To be able to use input file parsing with hdf5, the
!> hdf5-input-files have to be generated using the script:
!> `convert_to_hdf5.py`. This script can be found in the repository
!> `SHEMAT-Suite_Scripts` under
!> `python/preprocessing/convert_to_hdf5.py`.
      SUBROUTINE read_bc(filen,line,i_pv,i_bt,posi,ilost,ismpl)
        use arrays
        use mod_genrl
        use mod_genrlc
        use mod_flow
        use mod_temp
        use mod_conc
        use mod_linfos
#ifndef noHDF
        use mod_input_file_parser_hdf5
#endif
        IMPLICIT NONE
        integer :: i, j, k
        integer :: ismpl
!
!       number of the opened file
        INTEGER filen, posi, ilost, i_errors
!       i_bcu    : bc-unit
!       i_bctp   : bc-time period
!       i_pv     : physical value index
!       i_bt     : boundary type (dirichlet, neumann)
!       i_si     : sub index (conc)
!       i_dir    : direction
!       i_records: number of entries for this section
        INTEGER i_bcu, i_bctp, i_pv, i_bt, i_si, i_dir, i_records
!       d_bcmy: bc-my
        DOUBLE PRECISION d_bcmy

        character (len=80) :: line

        LOGICAL, ALLOCATABLE :: tmpbl(:,:,:,:)
        INTEGER, ALLOCATABLE :: tmpind(:,:)
        DOUBLE PRECISION, ALLOCATABLE :: tmpval(:,:)

        INTEGER lblank, read_direction, i_b, i_e, j_b, j_e, k_b, k_e, ll
        INTEGER level_bcindex

        LOGICAL read_simple, read_bctp, read_species, read_bcval, &
          l_errign
        LOGICAL found, no_ext_link, no_ext_link_int
        EXTERNAL found, no_ext_link, no_ext_link_int, lblank, &
          read_direction

        character (len=5), dimension (0:6) :: c_dir
        DATA c_dir/'none', 'left', 'right', 'front', 'back', 'base', &
          'top'/
        character (len=8) :: c_name
        logical :: found_marker
        character (len=80) :: full_bc_name

!       if switch used, read one column more (bcindex)
        i_bcu = 0
!       0: do not read bc-units, but read pv-values
!       1: read bc-units
!       2: do not read bc-units and pv-values
        level_bcindex = 0
#ifndef noHDF
        if (h5parse_use_hdf5_datafile) then
            full_bc_name = trim(line)
            if (h5parse_check_attr_exist("bcindex","bc/i"//full_bc_name)) then
                call h5parse_read_integer_attribute("bcindex",level_bcindex,"bc/i"//full_bc_name)
            end if
        else
#endif
        CALL get_arg('bcindex',line,i,j)
        IF (i>=1 .AND. j>=i) THEN
!          prove "r" instead of "read"
          IF (line(i:i)=='r') THEN
            level_bcindex = 1
          ELSE
            level_bcindex = 2
            READ(line(i:j),*) i_bcu
!            write(*,'(1A,1I5)') '      bc unit=',i_bcu
          END IF
        END IF
#ifndef noHDF
        end if
#endif

        d_bcmy = 1.0D+18
#ifndef noHDF
        if (h5parse_use_hdf5_datafile) then
            if (h5parse_check_attr_exist("bcmy","bc/"//full_bc_name)) then
                call h5parse_read_double_attribute("bcmy",d_bcmy,"bc/"//full_bc_name)
            end if
        else
#endif
        CALL get_arg('bcmy',line,i,j)
        IF (i>=1 .AND. j>=i) THEN
          READ(line(i:j),*) d_bcmy
!          write(*,'(1A,1e16.8)') '      bcmy=',d_bcmy
        END IF
#ifndef noHDF
        end if
#endif

        i_si = 0
        read_species = .FALSE.
        IF (i_pv==pv_conc) read_species = .TRUE.
#ifndef noHDF
        if (h5parse_use_hdf5_datafile) then
            if (h5parse_check_attr_exist("species","bc/i"//full_bc_name)) then
                call h5parse_read_integer_attribute("species",i_si,"bc/i"//full_bc_name)
                read_species = .FALSE.
            end if
        else
#endif
        CALL get_arg('species',line,i,j)
        IF (i>=1 .AND. j>=i) THEN
          READ(line(i:j),*) i_si
          read_species = .FALSE.
!          write(*,'(1A,1I3)') '      speCies=',i_si
        END IF
#ifndef noHDF
        end if
#endif

        i_bctp = 0
        read_bctp = .TRUE.
#ifndef noHDF
        if (h5parse_use_hdf5_datafile) then
            if (h5parse_check_attr_exist("bctp","bc/i"//full_bc_name)) then
                call h5parse_read_integer_attribute("bctp",i_bctp,"bc/i"//full_bc_name)
                read_bctp = .False.
            end if
        else
#endif
        CALL get_arg('bctp',line,i,j)
        IF (i>=1 .AND. j>=i) THEN
          READ(line(i:j),*) i_bctp
          read_bctp = .FALSE.
!          write(*,'(1A,1I5)') '      bCtp=',i_bCtp
        END IF
#ifndef noHDF
        end if
#endif

#ifndef noHDF
        if (h5parse_use_hdf5_datafile) then
            read_bcval = .not. h5parse_check_attr_exist("value","bc/i"//full_bc_name)
        else
#endif
        read_bcval = .TRUE.
        CALL get_arg('value',line,i,j)
        IF (i>=1 .AND. j>=i) THEN
!         prove "i" instead of "init"
          IF (line(i:i)=='i') THEN
            read_bcval = .FALSE.
          ELSE
            WRITE(*,'(3A)') 'warning: option "value=', line(i:j), &
              '" ignored.'
          END IF
        END IF
#ifndef noHDF
        end if
#endif

#ifndef noHDF
        if (h5parse_use_hdf5_datafile) then
            l_errign = h5parse_check_attr_exist("error","bc/i"//full_bc_name)
        else
#endif
        l_errign = .FALSE.
        CALL get_arg('error',line,i,j)
        IF (i>=1 .AND. j>=i) THEN
!         prove "i" instead of "ignore"
          IF (line(i:i)=='i') THEN
            l_errign = .TRUE.
          ELSE
            WRITE(*,'(3A)') 'warning: option "error=', line(i:j), &
              '" ignored.'
          END IF
        END IF
#ifndef noHDF
        end if
#endif

        i_dir = -1
#ifndef noHDF
        if (h5parse_use_hdf5_datafile) then
            l_errign = h5parse_check_attr_exist("error","bc/i"//full_bc_name)
        else
#endif
        read_simple = .FALSE.
        CALL get_arg('simple',line,i,j)
        IF (i>=1 .AND. j>=i) THEN
          i_dir = read_direction(line(i:j))
          read_simple = .TRUE.
!         enable "error=ignore" as default in <simple> case
          l_errign = .TRUE.
        END IF
#ifndef noHDF
        end if
#endif

#ifndef noHDF
        if (h5parse_use_hdf5_datafile) then
            if (h5parse_check_attr_exist("direction","bc/i"//full_bc_name)) then
                call h5parse_read_integer_attribute("direction",i_dir,"bc/i"//full_bc_name)
            end if
        else
#endif
        CALL get_arg('direction',line,i,j)
        IF (i>=1 .AND. j>=i) THEN
          k = read_direction(line(i:j))
          IF (i_dir==-1) THEN
            i_dir = k
!debug            write(*,'(2A)') '      direCtion=',C_dir(i_dir)
          ELSE IF (i_dir/=k) THEN
            WRITE(*,'(1A)') 'error: "direction"-"simple" mismatch !'
            STOP
          END IF
        END IF
#ifndef noHDF
        end if
#endif

!       zero, should be forbidden !!!
        i_dir = max(0,i_dir)
#ifndef noHDF
        if (h5parse_use_hdf5_datafile) then
            i_records = h5parse_read_dimension_size_for_dataset("bc/"//full_bc_name)
        else
#endif
        i_b = 1
        i_e = i0
        j_b = 1
        j_e = j0
        k_b = 1
        k_e = k0
        IF (read_simple) THEN
          IF (i_dir==0) THEN
            WRITE(*,*) 'error: simple=none not allowed !!!'
            STOP
          ELSE IF (i_dir==1) THEN
            i_records = j0*k0
            i_b = 1
            i_e = 1
          ELSE IF (i_dir==2) THEN
            i_records = j0*k0
            i_b = i0
            i_e = i0
          ELSE IF (i_dir==3) THEN
            i_records = i0*k0
            j_b = 1
            j_e = 1
          ELSE IF (i_dir==4) THEN
            i_records = i0*k0
            j_b = j0
            j_e = j0
          ELSE IF (i_dir==5) THEN
            i_records = i0*j0
            k_b = 1
            k_e = 1
          ELSE IF (i_dir==6) THEN
            i_records = i0*j0
            k_b = k0
            k_e = k0
          END IF
        ELSE
          CALL get_arg('records',line,i,j)
          IF (i<1 .OR. j<i) THEN
            READ(filen,*) i_records
          ELSE
            READ(line(i:j),*) i_records
          END IF
        END IF
!debug       write(*,'(1A,1I5)') '      reCords=',i_reCords
#ifndef noHDF
        end if
#endif

        ALLOCATE(tmpind(i_records,nibc))
        ALLOCATE(tmpval(i_records,ndbc))

found_marker = .false.
#ifndef noHDF
        if (h5parse_use_hdf5_datafile) then
            call h5parse_read_2d_double_dataset("bc/"//full_bc_name,tmpval)
            call h5parse_read_2d_integer_dataset("bc/i"//full_bc_name,tmpind)
            found_marker = .true.
        else
#endif

        c_name = pv_name(i_pv) // '_' // bc_name(i_bt)
!       sainty check
        IF (i_pv==pv_conc .AND. i_si<1 .AND. .NOT. read_species) THEN
          WRITE(*,'(1A)') 'error: "species" index invalid!'
          STOP
        END IF

        IF (no_ext_link_int(i_records,nibc,1,tmpind,'i'//c_name,line) &
            .AND. no_ext_link(i_records,ndbc,1,tmpval,c_name,line)) &
            THEN
            found_marker = .true.
          IF (read_simple) THEN
!             e.g. [temp=10.d0]
!             need : bcindex = ?
!                    bctp    = ?
!                    species = ?
!                    simple  = ? (direction)
!                    bcmy    = ?
            IF (read_species) THEN
              WRITE(*,'(1A)') 'error: "species" index missing!!!'
              STOP
            END IF
            IF (level_bcindex==0 .AND. read_bcval) THEN
              READ(filen,*,err=1001,end=1001) (tmpval(ll,1),ll=1, &
                i_records)
            ELSE IF (level_bcindex==2) THEN
              DO ll = 1, i_records
                tmpval(ll,1) = 0.D0
              END DO
            ELSE IF (level_bcindex==1) THEN
              WRITE(*,'(1A)') 'error: when "simple=?", then "bcindex=read" not supported!!!'
              STOP
            END IF

            ll = 0
            DO k = k_b, k_e
              DO j = j_b, j_e
                DO i = i_b, i_e
                  ll = ll + 1
                  tmpind(ll,cbc_i) = i
                  tmpind(ll,cbc_j) = j
                  tmpind(ll,cbc_k) = k
                  tmpind(ll,cbc_bcu) = i_bcu
                  tmpind(ll,cbc_bctp) = i_bctp
                  tmpind(ll,cbc_pv) = i_pv
                  tmpind(ll,cbc_bt) = i_bt
                  tmpind(ll,cbc_si) = i_si
                  tmpind(ll,cbc_dir) = i_dir
                  tmpval(ll,2) = d_bcmy
                END DO
              END DO
            END DO
          ELSE
!           init, preset values
            DO ll = 1, i_records
              tmpind(ll,cbc_bcu) = i_bcu
              tmpind(ll,cbc_bctp) = i_bctp
              tmpind(ll,cbc_pv) = i_pv
              tmpind(ll,cbc_bt) = i_bt
              tmpind(ll,cbc_si) = i_si
              tmpind(ll,cbc_dir) = i_dir
              tmpval(ll,1) = 0.0D0
              tmpval(ll,2) = d_bcmy
              READ(filen,'(1A)',err=1000,end=1000) line
!             read, overwrite values
              IF (level_bcindex==1 .AND. read_bctp .AND. read_species) &
                  THEN
!               e.g. [i=1,j=10,k=3, bc-unit=4, bctp-id=0, species=1]
                READ(line,*,err=1002,end=1002) (tmpind(ll,j),j=1,3), &
                  tmpind(ll,cbc_bcu), tmpind(ll,cbc_bctp), &
                  tmpind(ll,cbc_si)
              ELSE IF (level_bcindex==1 .AND. .NOT. read_bctp .AND. &
                  read_species) THEN
!               e.g. [i=1,j=10,k=3, bc-unit=4, species=1]
                READ(line,*,err=1003,end=1003) (tmpind(ll,j),j=1,3), &
                  tmpind(ll,cbc_bcu), tmpind(ll,cbc_si)
              ELSE IF (level_bcindex==1 .AND. read_bctp .AND. &
                  .NOT. read_species) THEN
!               e.g. [i=1,j=10,k=3, bc-unit=4, bctp-id=0]
                READ(line,*,err=1004,end=1004) (tmpind(ll,j),j=1,3), &
                  tmpind(ll,cbc_bcu), tmpind(ll,cbc_bctp)
              ELSE IF (level_bcindex==1 .AND. .NOT. read_bctp .AND. &
                  .NOT. read_species) THEN
!               e.g. [i=1,j=10,k=3, bc-unit=4]
                READ(line,*,err=1005,end=1005) (tmpind(ll,j),j=1,3), &
                  tmpind(ll,cbc_bcu)
              ELSE IF (level_bcindex==0 .AND. read_bctp .AND. &
                  read_species .AND. read_bcval) THEN
!               e.g. [i=1,j=10,k=3, temp=10.d0, bctp-id=0, species=1]
                READ(line,*,err=1006,end=1006) (tmpind(ll,j),j=1,3), &
                  tmpval(ll,1), tmpind(ll,cbc_bctp), tmpind(ll,cbc_si)
              ELSE IF (level_bcindex==0 .AND. .NOT. read_bctp .AND. &
                  read_species .AND. read_bcval) THEN
!               e.g. [i=1,j=10,k=3, temp=10.d0, species=1]
                READ(line,*,err=1007,end=1007) (tmpind(ll,j),j=1,3), &
                  tmpval(ll,1), tmpind(ll,cbc_si)
              ELSE IF (level_bcindex==0 .AND. read_bctp .AND. &
                  .NOT. read_species .AND. read_bcval) THEN
!               e.g. [i=1,j=10,k=3, temp=10.d0, bctp-id=0]
                READ(line,*,err=1008,end=1008) (tmpind(ll,j),j=1,3), &
                  tmpval(ll,1), tmpind(ll,cbc_bctp)
              ELSE IF (level_bcindex==0 .AND. .NOT. read_bctp .AND. &
                  .NOT. read_species .AND. read_bcval) THEN
!               e.g. [i=1,j=10,k=3, temp=10.d0]
                READ(line,*,err=1009,end=1009) (tmpind(ll,j),j=1,3), &
                  tmpval(ll,1)
              ELSE IF (level_bcindex==2 .AND. read_bctp .AND. &
                  read_species) THEN
!               e.g. [i=1,j=10,k=3, bctp-id=0, species=1]
                READ(line,*,err=1010,end=1010) (tmpind(ll,j),j=1,3), &
                  tmpind(ll,cbc_bctp), tmpind(ll,cbc_si)
              ELSE IF (level_bcindex==2 .AND. .NOT. read_bctp .AND. &
                  read_species) THEN
!               e.g. [i=1,j=10,k=3, species=1]
                READ(line,*,err=1011,end=1011) (tmpind(ll,j),j=1,3), &
                  tmpind(ll,cbc_si)
              ELSE IF (level_bcindex==2 .AND. read_bctp .AND. &
                  .NOT. read_species) THEN
!               e.g. [i=1,j=10,k=3, bctp-id=0]
                READ(line,*,err=1012,end=1012) (tmpind(ll,j),j=1,3), &
                  tmpind(ll,cbc_bctp)
              ELSE IF (level_bcindex==2 .AND. .NOT. read_bctp .AND. &
                  .NOT. read_species) THEN
!               e.g. [i=1,j=10,k=3]
                READ(line,*,err=1013,end=1013) (tmpind(ll,j),j=1,3)
              ELSE IF (level_bcindex==0 .AND. read_bctp .AND. &
                  read_species .AND. .NOT. read_bcval) THEN
!               e.g. [i=1,j=10,k=3, bctp-id=0, species=1]
                READ(line,*,err=1014,end=1006) (tmpind(ll,j),j=1,3), &
                  tmpind(ll,cbc_bctp), tmpind(ll,cbc_si)
              ELSE IF (level_bcindex==0 .AND. .NOT. read_bctp .AND. &
                  read_species .AND. .NOT. read_bcval) THEN
!               e.g. [i=1,j=10,k=3, species=1]
                READ(line,*,err=1015,end=1007) (tmpind(ll,j),j=1,3), &
                  tmpind(ll,cbc_si)
              ELSE IF (level_bcindex==0 .AND. read_bctp .AND. &
                  .NOT. read_species .AND. .NOT. read_bcval) THEN
!               e.g. [i=1,j=10,k=3, bctp-id=0]
                READ(line,*,err=1016,end=1008) (tmpind(ll,j),j=1,3), &
                  tmpind(ll,cbc_bctp)
              ELSE IF (level_bcindex==0 .AND. .NOT. read_bctp .AND. &
                  .NOT. read_species .AND. .NOT. read_bcval) THEN
!               e.g. [i=1,j=10,k=3]
                READ(line,*,err=1017,end=1009) (tmpind(ll,j),j=1,3)
              END IF
 end do
          end if
        end if
#ifndef noHDF
        end if
#endif

        if (found_marker) then
          found_marker = .false.
          if (read_simple) then
            do ll = 1, i_records
            ! copy value, when "value=init" was speCified
                IF ( .NOT. read_bcval) THEN
                    i = tmpind(ll,cbc_i)
                    j = tmpind(ll,cbc_j)
                    k = tmpind(ll,cbc_k)
                    IF (i_pv==pv_head) tmpval(ll,1) = head(i,j,k, ismpl)
                    IF (i_pv==pv_pres) tmpval(ll,1) = pres(i,j,k, ismpl)
                    IF (i_pv==pv_temp) tmpval(ll,1) = temp(i,j,k, ismpl)
                    IF (i_pv==pv_conc) tmpval(ll,1) = conc(i,j,k,i_si, ismpl)
                END IF
            end do
          else
            do ll = 1, i_records

!
!     sanity checks
!             i,j,k index check
              IF (tmpind(ll,cbc_i)<1 .OR. tmpind(ll,cbc_i)>i0 .OR. &
                  tmpind(ll,cbc_j)<1 .OR. tmpind(ll,cbc_j)>j0 .OR. &
                  tmpind(ll,cbc_k)<1 .OR. tmpind(ll,cbc_k)>k0) THEN
                WRITE(*,'(1A,1I6,1A)') &
                  'error: index out of range, at line ', ll, '!'
                WRITE(*,'(3A)') '  given:"', line, '"'
                STOP
              END IF
!             pv index check
              IF (tmpind(ll,cbc_pv)<1 .OR. tmpind(ll,cbc_pv)>npv) THEN
                WRITE(*,'(1A,1I6,1A)') 'error: physical value index out of range [1..3], at line ', ll, '!'
                WRITE(*,'(3A)') '  given:"', line, '"'
                STOP
              END IF
!             species index check
              IF (tmpind(ll,cbc_si)<0 .OR. tmpind(ll,cbc_si)>ntrans) THEN
                WRITE(*,'(1A,1I6,1A)') &
                  'error: species index out of range, at line ', ll,'!'
                WRITE(*,'(3A)') '  given:"', line, '"'
                STOP
              END IF
!             direction index check
              IF (tmpind(ll,cbc_dir)<0 .OR. tmpind(ll,cbc_dir)>6) THEN
                WRITE(*,'(1A,1I6,1A)') 'error: direction index out of range [0..6], at line ', ll, '!'
                WRITE(*,'(3A)') '  given:"', line, '"'
                STOP
              END IF
!             wells - direction check
              IF (tmpind(ll,cbc_dir)/=0 .AND. tmpind(ll,cbc_bt)==bt_neuw) THEN
                WRITE(*,'(1A,1I6,1A)') 'error: well function needs direction index 0, at line ', ll, '!'
                WRITE(*,'(3A)') '  given:"', line, '"'
                STOP
              END IF
!             bc-type check
              IF (tmpind(ll,cbc_bt)/=bt_diri .AND. &
                  tmpind(ll,cbc_bt)/=bt_neum .AND. &
                  tmpind(ll,cbc_bt)/=bt_neuw) THEN
                WRITE(*,'(1A,1I6,1A)') &
                  'error: BC type out of range [1..3], at line ', ll, '!'
                WRITE(*,'(3A)') '  given:"', line, '"'
                STOP
              END IF
!             bc time period (BCTP) check
!               For deeper BCTP checks, a proof for 'tmpind(ll,cbc_bctp)' may be a good idea,
!               but the time period table is needed and readed later -> check then !
!
!             copy value, when "value=init" was speCified
              IF ( .NOT. read_bcval) THEN
                IF (tmpind(ll,cbc_pv)==pv_head) tmpval(ll,1) &
                  = head(tmpind(ll,cbc_i),tmpind(ll,cbc_j), &
                  tmpind(ll,cbc_k),ismpl)
                IF (tmpind(ll,cbc_pv)==pv_pres) tmpval(ll,1) &
                  = pres(tmpind(ll,cbc_i),tmpind(ll,cbc_j), &
                  tmpind(ll,cbc_k),ismpl)
                IF (tmpind(ll,cbc_pv)==pv_temp) tmpval(ll,1) &
                  = temp(tmpind(ll,cbc_i),tmpind(ll,cbc_j), &
                  tmpind(ll,cbc_k),ismpl)
                IF (tmpind(ll,cbc_pv)==pv_conc) tmpval(ll,1) &
                  = conc(tmpind(ll,cbc_i),tmpind(ll,cbc_j), &
                  tmpind(ll,cbc_k),tmpind(ll,cbc_si),ismpl)
              END IF
            END DO
          END IF
        END IF

!     sanity check
        ALLOCATE(tmpbl(I0,J0,K0,max(1,ntrans)))
        CALL set_lval(I0*J0*K0*max(1,ntrans), .FALSE., tmpbl)
!       mark existing boundary points
        DO i = 1, posi
          IF (ibc_data(i,cbc_pv)==i_pv) THEN
            tmpbl(ibc_data(i,cbc_i),ibc_data(i,cbc_j),ibc_data(i,cbc_k),max(1,ibc_data(i,cbc_si))) = .TRUE.
          END IF
        END DO
!       check for double definition of boundary points
        IF (l_errign) THEN
!         generate 'ignore' warnings
          i_errors = 0
          DO j = 1, i_records
            IF (tmpbl(tmpind(j,cbc_i),tmpind(j,cbc_j),tmpind(j,cbc_k),max(1,tmpind(j,cbc_si)))) THEN
              IF (read_simple) THEN
                i_errors = i_errors +1
              ELSE
                WRITE(*,'(1A,3I6,1A)') '        BC ignored at [', &
                  tmpind(j,cbc_i), tmpind(j,cbc_j), tmpind(j,cbc_k),']'
              END IF
              tmpind(j,cbc_pv) = -1
              ilost = ilost + 1
            END IF
          END DO
          IF (read_simple .AND. i_errors>0) THEN
            WRITE(*,'(1A,1I6,1A)') '        ',i_errors," BC's ignored (SIMPLE mode)"
          END IF
        ELSE
!         generate errors
          i_errors = 0
          DO j = 1, i_records
            IF (tmpbl(tmpind(j,cbc_i),tmpind(j,cbc_j),tmpind(j,cbc_k),max(1,tmpind(j,cbc_si)))) THEN
              WRITE(*,'(1A,3(1X,1I4),1A,1I1,1A,1I3,1A)') &
                '  double defined type at [', tmpind(j,cbc_i), tmpind(j,cbc_j), tmpind(j,cbc_k), &
                ' ],type=', i_pv,',species=', tmpind(j,cbc_si),', (different BC type ?) !'
              i_errors = i_errors +1
            END IF
          END DO
          IF (i_errors.gt.0) THEN
            WRITE(*,'(1A)') 'error: to many invalid declarations above!!!'
            STOP
          END IF
        END IF
        DEALLOCATE(tmpbl)
!
        DO i = 1, i_records
          IF (tmpind(i,cbc_pv)>=0) THEN
            posi = posi + 1
            DO j = 1, nibc
              ibc_data(posi,j) = tmpind(i,j)
            END DO
            DO j = 1, ndbc
              dbc_data(posi,j,ismpl) = tmpval(i,j)
            END DO
!           convert [MPa] into [Pa]
            IF (ibc_data(posi,cbc_pv)==pv_pres .AND. read_bcval) &
              dbc_data(posi,1,ismpl) = dbc_data(posi,1,ismpl)*pa_conv
!
            bc_maxunits = max(bc_maxunits,tmpind(i,cbc_bcu))
            IF ((tmpind(i,cbc_bcu)>nunits) .OR. (tmpind(i, &
                cbc_bcu)<0)) THEN
              WRITE(*,'(1A,1I7,1A)') &
                'error: bc-unit number out of range, at line ', i, '!'
              STOP
            END IF
          END IF
        END DO
!
        WRITE(*,'(4A)',advance='NO') '  [R] : ', pv_name(i_pv), ' ', &
          bc_name(i_bt)
        IF (read_simple) THEN
          WRITE(*,'(3A,1I6,1A)',advance='NO') ', simple=', &
            c_dir(i_dir), ', (size=', i_records, ')'
        ELSE
          WRITE(*,'(3A,1I6)',advance='NO') ', direction=', &
            c_dir(i_dir), ', records=', i_records
        END IF
        WRITE(*,'(1A,1e12.4)',advance='NO') ', bcmy=', d_bcmy
        IF (level_bcindex==1) THEN
          WRITE(*,'(1A)',advance='NO') ', bcindex=read'
        END IF
        IF (level_bcindex==2) THEN
          WRITE(*,'(1A,1I4)',advance='NO') ', bcindex=', i_bcu
        END IF
        IF (i_bctp/=0) WRITE(*,'(1A,1I4)',advance='NO') ', bctp=', &
          i_bctp
        IF (i_si/=0) WRITE(*,'(1A,1I4)',advance='NO') ', species=', &
          i_si
        IF (l_errign) WRITE(*,'(1A)',advance='NO') ', error=ignore'
        WRITE(*,'(1A,1I6,1A)') ', (offset=', posi - i_records + 1, &
          ')'
!

        DEALLOCATE(tmpval)
        DEALLOCATE(tmpind)
!
        RETURN

!       ERROR handler
1000    WRITE(*,'(1A)') 'error: expect a BC data line'
        STOP
1001    WRITE(*,'(1A)') 'error: to few values - "simple=?" specified'
        WRITE(*,'(1A)') '   expected: I0*J0*K0 x [bc-value]'
        STOP
1002    WRITE(*,'(1A)') &
          'error: to few values - "bcindex=read" specified'
        WRITE(*,'(1A)') &
          '   expected: records x [i,j,k,bc-unit,tpbc-id,species]'
        WRITE(*,'(3A)') '  given:"', line, '"'
        STOP
1003    WRITE(*,'(1A)') &
          'error: to few values - "bcindex=read,bctp=?" specified'
        WRITE(*,'(1A)') &
          '   expected: records x [i,j,k,bc-unit,species]'
        WRITE(*,'(3A)') '  given:"', line, '"'
        STOP
1004    WRITE(*,'(1A)') &
          'error: to few values - "bcindex=read,species=?" specified'
        WRITE(*,'(1A)') &
          '   expected: records x [i,j,k,bc-unit,tpbc-id]'
        WRITE(*,'(3A)') '  given:"', line, '"'
        STOP
1005    WRITE(*,'(1A)') 'error: to few values - "bcindex=read,bctp=?,species=?" specified'
        WRITE(*,'(1A)') '   expected: records x [i,j,k,bc-unit]'
        WRITE(*,'(3A)') '  given:"', line, '"'
        STOP
1006    WRITE(*,'(1A)') 'error: to few values - no option specified'
        WRITE(*,'(1A)') &
          '   expected: records x [i,j,k,bc-value,tpbc-id,species]'
        WRITE(*,'(3A)') '  given:"', line, '"'
        STOP
1007    WRITE(*,'(1A)') 'error: to few values - "bctp=?" specified'
        WRITE(*,'(1A)') &
          '   expected: records x [i,j,k,bc-value,species]'
        WRITE(*,'(3A)') '  given:"', line, '"'
        STOP
1008    WRITE(*,'(1A)') &
          'error: to few values - "species=?" specified'
        WRITE(*,'(1A)') &
          '   expected: records x [i,j,k,bc-value,tpbc-id]'
        WRITE(*,'(3A)') '  given:"', line, '"'
        STOP
1009    WRITE(*,'(1A)') &
          'error: to few values - "bctp=?,species=?" specified'
        WRITE(*,'(1A)') '   expected: records x [i,j,k,bc-value]'
        WRITE(*,'(3A)') '  given:"', line, '"'
        STOP
1010    WRITE(*,'(1A)') &
          'error: to few values - "bcindex=?" specified'
        WRITE(*,'(1A)') &
          '   expected: records x [i,j,k,tpbc-id,species]'
        WRITE(*,'(3A)') '  given:"', line, '"'
        STOP
1011    WRITE(*,'(1A)') &
          'error: to few values - "bcindex=?,bctp=?" specified'
        WRITE(*,'(1A)') '   expected: records x [i,j,k,species]'
        WRITE(*,'(3A)') '  given:"', line, '"'
        STOP
1012    WRITE(*,'(1A)') &
          'error: to few values - "bcindex=?,species=?" specified'
        WRITE(*,'(1A)') '   expected: records x [i,j,k,tpbc-id]'
        WRITE(*,'(3A)') '  given:"', line, '"'
        STOP
1013    WRITE(*,'(1A)') 'error: to few values - "bcindex=?,bctp=?,species=?" specified'
        WRITE(*,'(1A)') '   expected: records x [i,j,k]'
        WRITE(*,'(3A)') '  given:"', line, '"'
        STOP
1014    WRITE(*,'(1A)') &
          'error: to few values - "value=init" specified'
        WRITE(*,'(1A)') &
          '   expected: records x [i,j,k,tpbc-id,species]'
        WRITE(*,'(3A)') '  given:"', line, '"'
        STOP
1015    WRITE(*,'(1A)') &
          'error: to few values - "value=init,bctp=?" specified'
        WRITE(*,'(1A)') '   expected: records x [i,j,k,species]'
        WRITE(*,'(3A)') '  given:"', line, '"'
        STOP
1016    WRITE(*,'(1A)') &
          'error: to few values - "value=init,species=?" specified'
        WRITE(*,'(1A)') '   expected: records x [i,j,k,tpbc-id]'
        WRITE(*,'(3A)') '  given:"', line, '"'
        STOP
1017    WRITE(*,'(1A)') 'error: to few values - "value=init,bctp=?,species=?" specified'
        WRITE(*,'(1A)') '   expected: records x [i,j,k]'
        WRITE(*,'(3A)') '  given:"', line, '"'
        STOP
      END

!>    @brief sort the boundary points, physical value order
!>    @param[in] ismpl local sample index
!>    @details
!>sort the boundary points, physical value order\n
      SUBROUTINE sort_bc(ismpl)
        use arrays
        use mod_genrl
        use mod_genrlc
        use mod_flow
        use mod_temp
        use mod_linfos
        IMPLICIT NONE
        integer :: i, j
        integer :: ismpl
        INTEGER, ALLOCATABLE :: tmpind(:,:)
        DOUBLE PRECISION, ALLOCATABLE :: tmpval(:,:)


        ALLOCATE(tmpind(nbc_data,nibc))
        ALLOCATE(tmpval(nbc_data,ndbc))

        first_flow = nbc_data + 1
        last_flow = 0
        IF (head_active .OR. pres_active) THEN
          DO i = 1, nbc_data
            IF (ibc_data(i,cbc_pv)==pv_head.OR.ibc_data(i,cbc_pv)==pv_pres) THEN
              last_flow = last_flow + 1
              first_flow = min(first_flow,last_flow)
              DO j = 1, nibc
                tmpind(last_flow,j) = ibc_data(i,j)
              END DO
              DO j = 1, ndbc
                tmpval(last_flow,j) = dbc_data(i,j,ismpl)
              END DO
            END IF
          END DO
        END IF

        first_temp = nbc_data + 1
        last_temp = last_flow
        IF (temp_active) THEN
          DO i = 1, nbc_data
            IF (ibc_data(i,cbc_pv)==pv_temp) THEN
              last_temp = last_temp + 1
              first_temp = min(first_temp,last_temp)
              DO j = 1, nibc
                tmpind(last_temp,j) = ibc_data(i,j)
              END DO
              DO j = 1, ndbc
                tmpval(last_temp,j) = dbc_data(i,j,ismpl)
              END DO
            END IF
          END DO
        END IF

        first_conc = nbc_data + 1
        last_conc = last_temp
        IF (trans_active) THEN
          DO i = 1, nbc_data
            IF (ibc_data(i,cbc_pv)==pv_conc) THEN
              last_conc = last_conc + 1
              first_conc = min(first_conc,last_conc)
              DO j = 1, nibc
                tmpind(last_conc,j) = ibc_data(i,j)
              END DO
              DO j = 1, ndbc
                tmpval(last_conc,j) = dbc_data(i,j,ismpl)
              END DO
            END IF
          END DO
        END IF

!       copy back sorted values
        DO j = 1, nibc
          DO i = 1, nbc_data
            ibc_data(i,j) = tmpind(i,j)
          END DO
        END DO
        DO j = 1, ndbc
          DO i = 1, nbc_data
            dbc_data(i,j,ismpl) = tmpval(i,j)
          END DO
        END DO

        DEALLOCATE(tmpval)
        DEALLOCATE(tmpind)

        RETURN
      END
