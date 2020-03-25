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

!>    @brief read model parameter
!>    @param[in] filename model file name
!>    @param[in] ismpl local sample index
!>    @details
!> Note: To be able to use input file parsing with hdf5, the
!> hdf5-input-files have to be generated using the script:
!> `convert_to_hdf5.py`. This script can be found in the repository
!> `SHEMAT-Suite_Scripts` under
!> `python/preprocessing/convert_to_hdf5.py`.
      SUBROUTINE read_model(filename,ismpl)
        use arrays
        use mod_genrl
        use mod_genrlc
        use mod_flow
        use mod_temp
        use mod_conc
        use mod_time
        use mod_data
        use mod_blocking_size
        use mod_OMP_TOOLS
        use mod_linfos
#ifndef noHDF
        use mod_input_file_parser_hdf5
#endif
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k, l

        INCLUDE 'OMP_TOOLS.inc'
!
        character (len=80) :: filename
        character (len=80) :: line
        character (len=320) :: longline
        character (len=4) :: ctmp
        character (len=32) :: strng
        character (len=1) :: sbc
!
        DOUBLE PRECISION dtmp
        DOUBLE PRECISION, ALLOCATABLE :: datmp(:,:)
        LOGICAL, ALLOCATABLE :: ltmp(:,:)
!       is head needed - head based computation
        LOGICAL head_needed
!
        INTEGER tmplen, minunits, posi, ilost
!
        INTEGER ijk, omp_inner, omp_outer
        INTEGER i1, i2, nbc_sections
        INTEGER tracer
        INTEGER sm_max
!
        INTEGER locstr, lblank, read_direction, get_ioptval
        EXTERNAL locstr, lblank, read_direction, get_ioptval
        LOGICAL found, no_ext_link, no_ext_link_int, test_null, &
          test_option
        EXTERNAL found, no_ext_link, no_ext_link_int, test_null, &
          test_option
        logical :: found_marker
        character(len=80) :: full_bc_name

        found_marker = .false.

!
        CALL read_check(filename)
!
        WRITE(*,*) ' '
        WRITE(*,*) '  reading model input parameter:'
        WRITE(*,*) '    from file "', filename(:lblank(filename)),'"'
        WRITE(*,*) ' '

! ------------------
!       default sample index
        ismpl = 1
        nsmpl = 1
!       generic init staff equal for all models
        CALL model_init(ismpl)
!       setup string constants for the number of property-units and bc-units
        WRITE(c_npropunit,'(I2)') lastidx -firstidx +1
        WRITE(c_nbcunit,'(I2)') bc_lastidx -bc_firstidx +1
!       setup string constant for the number of variables
        WRITE(c_npv,'(I2)') npv
! ------------------
!     read file
        OPEN(79,file=filename,status='old')

!     init HDF5 support, when available
        CALL open_hdf5(' ')

        title = 'NO TITLE'
        IF (found(79,key_char//' title',line,.FALSE.)) THEN
          READ(79,'(1A)',err=200,end=200) title
          WRITE(*,*) ' [R] : title'
        ELSE
          WRITE(*,*) ' <D> : no title '
        END IF

        runmode = 0
        IF (found(79,key_char//' runmode',line,.FALSE.)) THEN
          READ(79,*,err=202,end=202) runmode
          WRITE(*,*) ' [R] : runmode'
          IF (runmode==0) WRITE(*,*) '    >>> forward modeling'
          IF (runmode==1) WRITE(*,*) &
            '    >>> forward modeling and data fit'
          IF (runmode==2) WRITE(*,*) &
            '    >>> inverse modeling (extra steady state)'
          IF (runmode==3) WRITE(*,*) '    >>> inverse modeling '
        ELSE
          WRITE(*,*) ' <D> : runmode=0, forward modeling !'
        END IF

#ifdef DEBUG
        n_debugout = 0
        CALL read_debugout(ismpl)
#endif

        write_smonitor = .FALSE.
        transient = .FALSE.

        ! Default time unit [s]
        tunit = tunit_const

        linfos(1) = 2
        linfos(2) = 1
        linfos(3:4) = 0
        IF (found(79,key_char//' linfo',line,.FALSE.)) THEN
          READ(79,*) linfos
          WRITE(*,*) ' [R] : linfo'
        END IF

#ifdef BENCH
        IF (found(79,key_char//' ilu block size',line,.FALSE.)) THEN
          READ(79,*) block_i,block_j,block_k
          linfos(4) = -1000
          WRITE(*,*) ' [R] : ILU block size (benchmark mode)'
        END IF
#endif

!       (pseudo) realisation number
        sm_max = 1
        IF (found(79,key_char//' samples',line,.FALSE.)) THEN
          READ(79,*,err=201,end=201) nsmpl
          WRITE(*,'(1A,1I6)') '  [R] : samples =',nsmpl
        END IF
        sm_max = nsmpl

        IF (found(79,key_char//' PROPS',line,.FALSE.)) THEN
          CALL get_arg('PROPS',line,i,j)
          IF (i<1 .OR. j<i) THEN
            READ(79,*,err=206,end=206) def_props
          ELSE
            READ(line(i:j),*) def_props
          END IF
        ELSE
          WRITE(*,'(2A)') &
            'error: can not find section "'//key_char//' PROPS=<...>",', &
            ' must be defined!'
          def_props = '<name>'
        END IF
        CALL props_check(ismpl)

        IF (found(79,key_char//' USER',line,.FALSE.)) THEN
          CALL get_arg('USER',line,i,j)
          IF (i<1 .OR. j<i) THEN
            READ(79,*,err=207,end=207) def_user
          ELSE
            READ(line(i:j),*) def_user
          END IF
        ELSE
          WRITE(*,'(2A)') &
            'error: can not find section "'//key_char//' USER=<...>",', &
            ' must be defined!'
          def_user = '<name>'
        END IF
        CALL user_check(ismpl)

!     test for command line option given OpenMP parallelisation
        omp_outer = 0
        omp_inner = 0
        IF (test_option('-tsample')) THEN
          omp_outer = get_ioptval('-tsample')
        END IF
        IF (test_option('-tsolve')) THEN
          omp_inner = get_ioptval('-tsolve')
        END IF
        IF (omp_outer>0 .AND. omp_inner>0) WRITE(*,'(2(1A,1I3),1A)') &
          '  [R] : command line OpenMP thread configuration (', &
          max(omp_outer,1), 'x', max(omp_inner,1), ' threads)'
        IF (omp_outer==0 .AND. omp_inner>0) WRITE(*,'(2A,1I3,1A)') &
          '  [R] : command line OpenMP thread configuration (', &
          '  ?x', max(omp_inner,1), ' threads)'
        IF (omp_outer>0 .AND. omp_inner==0) WRITE(*,'(1A,1I3,1A)') &
          '  [R] : command line OpenMP thread configuration (', &
          max(omp_outer,1), 'x  ? threads)'

!       test for default environment given OpenMP parallelisation
        IF (omp_inner>0) THEN
!$OMP     parallel default(none) shared(Tlevel_0, omp_inner)&
!$OMP       num_threads(omp_inner)
!$OMP       master
!             test thread configuartion for thread-level 0
!             (max number of sample threads)
              tlevel_0 = omp_get_num_of_threads()
!$OMP       end master
!$OMP     end parallel
        ELSE
!$OMP     parallel default(none) shared(Tlevel_0)
!$OMP       master
!             get the number threads for thread-level 0
!             (max number of sample threads)
              tlevel_0 = omp_get_num_of_threads()
!$OMP       end master
!$OMP     end parallel
        END IF
!$OMP   parallel default(none) num_threads(Tlevel_0)&
!$OMP     shared(Tlevel_0,Tlevel_1,omp_inner)
        IF (omp_get_his_thread_num()==0) THEN
          tlevel_1 = tlevel_0
          tlevel_0 = 1
        END IF
!$OMP   end parallel

!       use "sm_max" instead of "nsmpl", because of the ENKF case nsmpl = sm_max+1
        IF (omp_outer>1) WRITE(*,'(2A)') &
          '  [I] : OpenMP parallelisation limited, ', &
          'target build not support nesting !'
#ifdef PROPS_IAPWS
        WRITE(*,'(2A)') '  [I] : OpenMP parallelisation limited, ', &
          'IAPWS target specific behaviour.'
#endif
        IF (tlevel_0>1 .OR. tlevel_1>1) WRITE(*,'(1A,1I3,1A,1I3,1A)') &
          '  [I] : OpenMP parallelisation enabled (', tlevel_0, 'x', &
          tlevel_1, ' threads)'

        IF (test_option('-scalemp')) THEN
          CALL init_scalemp_binding()
        ELSE IF (test_option('-libnuma')) THEN
          CALL init_scalemp_binding()
        ELSE IF (test_option('-tbind')) THEN
          CALL get_coptval('-tbind',line)
          CALL load_binding(line)
        ELSE
          CALL load_binding('default')
        END IF

#ifndef noHDF
        if (h5parse_use_hdf5_datafile) then
          i0 = h5parse_read_dimension_size_for_dataset("grid/delx")
          j0 = h5parse_read_dimension_size_for_dataset("grid/dely")
          k0 = h5parse_read_dimension_size_for_dataset("grid/delz")
          found_marker = .true.
        else
#endif
        IF (found(79,key_char//' grid',line,.FALSE.)) THEN
          READ(79,*,err=203,end=203) i0, j0, k0
          found_marker = .true.
        ELSE
          WRITE(*,'(1A)') &
            'error: no grid dimensions i0, j0, k0 defined !'
          STOP
        END IF
#ifndef noHDF
        endif
#endif
        if (found_marker) then
          found_marker = .false.
          write(*,'(A,3(I4,A))') '  [R] : [i0, j0, k0] = [', i0, ',', &
            j0, ',', k0, ']'
        end if


        maxiter_nl = 100
        nladapt = 0
        nlconverge = 0
        IF (found(79,key_char//' nlsolve',line,.FALSE.)) THEN
          READ(79,*,err=204,end=204) maxiter_nl, nladapt
          READ(79,*,err=1432,end=1432) nlconverge
          GOTO 1433
1432      nlconverge=0
1433      WRITE(*,*) ' [R] : nonlinear solver parameter'
          if (.NOT. nlconverge .eq. 0) write(*,*) ' [W] : Nonlinear Convergence test disabled!'
          IF (nladapt==0) THEN
            WRITE(*,*) ' [I] : fixed relxation factor assumed'
          ELSE IF (nladapt==1) THEN
            WRITE(*,*) ' [I] : adaptive relaxation type ', nladapt
            WRITE(*,*) &
              ' [E] : adaptive relaxation disabled manually!', ' &
              &-> please select the fixed mode (set 0) or &
              &remove this STOP!'
!         Please ask Volker Rath or Andreas Wolf !!!
            STOP
          ELSE
            WRITE(*,*) 'error: relaxation type ', nladapt, &
              ' not defined'
            STOP
          END IF
        END IF

!     --- read switches ---
        IF (found(79,key_char//' active',line,.FALSE.)) THEN
          head_active = .FALSE.
          pres_active = .FALSE.
          temp_active = .FALSE.
          trac_active = .FALSE.
          chem_active = .FALSE.
          trans_active = .FALSE.
!        test for each allowed type of calculation
          i = locstr(line,'head')
          IF (i>=1) THEN
            head_active = .TRUE.
            pres_active = .TRUE.
          END IF
          i = locstr(line,'pres')
          IF (i>=1) THEN
            pres_active = .TRUE.
#ifdef head_base
            head_active = .TRUE.
#endif
          END IF
          i = locstr(line,'temp')
          IF (i>=1) THEN
            temp_active = .TRUE.
          END IF
          i = locstr(line,'trac')
          IF (i>=1) THEN
            trac_active = .TRUE.
          END IF
          i = locstr(line,'chem')
          IF (i>=1) THEN
            chem_active = .TRUE.
          END IF
          IF (trac_active .OR. chem_active) trans_active = .TRUE.
        ELSE
          head_active = .TRUE.
          pres_active = .TRUE.
          temp_active = .TRUE.
          trac_active = .FALSE.
          chem_active = .FALSE.
          trans_active = .FALSE.
        END IF

!     default file output compression: plain, compress_out=1
!     default file output compression: bzip2 -> "*.bz2", compress_out=2
        compress_out = 1
        IF (found(79,key_char//' file output',line,.FALSE.)) THEN
          hdf_out = .FALSE.
          tec_out = .FALSE.
          vtk_out = .FALSE.
          txt_out = .FALSE.
          ctmp = '    '

!        test for each allowed type of output files
          i = locstr(line,'h5')
          IF (i>=1) THEN
            hdf_out = .TRUE.
            ctmp(1:1) = 'X'
          END IF
!-
          i = locstr(line,'plt')
          IF (i>=1) THEN
            tec_out = .TRUE.
            ctmp(2:2) = 'X'
          END IF
!-
          i = locstr(line,'vtk')
          IF (i>=1) THEN
            vtk_out = .TRUE.
            ctmp(3:3) = 'X'
          END IF
!-
          i = locstr(line,'txt')
          IF (i>=1) THEN
            txt_out = .TRUE.
            ctmp(4:4) = 'X'
          END IF
!- may be obsolete later, use "h5" instead
          i = locstr(line,'hdf')
          IF (i>=1) THEN
            hdf_out = .TRUE.
            ctmp(1:1) = 'X'
          END IF
!- may be obsolete later, use "plt" instead
          i = locstr(line,'tec')
          IF (i>=1) THEN
            tec_out = .TRUE.
            ctmp(2:2) = 'X'
          END IF
!-
          WRITE(*,'(1A,4(1A,1A1),1A)') '  [R] : output suffix: ', &
            '.h5:[', ctmp(1:1), '], .plt:[', ctmp(2:2), '], .vtk:[', &
            ctmp(3:3), '], .txt:[', ctmp(4:4), ']'

          DO j = 1, ncompress
!           looking for suffix names -> enable compression tool
            i = locstr(line,compress_suffix(j))
            IF (i>=1) THEN
              compress_out = j
            END IF
          END DO
          IF (compress_out>1) WRITE(*,'(3A)') &
            '  [R] : output file compression enabled [*.', &
            compress_suffix(compress_out), ']'
        ELSE
          hdf_out = .TRUE.
          tec_out = .FALSE.
          vtk_out = .FALSE.
          txt_out = .FALSE.
        END IF

        DO i = 1, nprop
          out_prop(i) = .TRUE.
        END DO
        DO i = 1, npv
          out_pv(i) = .TRUE.
        END DO
        DO i = 1, nout_ijk
          out_ijk(i) = .TRUE.
        END DO
!       disable specific hdf5-outputs
        IF (found(79,key_char//' disable output',line,.FALSE.)) THEN
          READ(79,'(1A)',end=215) longline
          CALL read_oprop(longline)
          CALL read_opv(longline)
          CALL read_oijk(longline)
        END IF

        write_param = .TRUE.
        IF (found(79,key_char//' disable small output',line,.FALSE.)) THEN
!         disable the additional parameter output file
          write_param = .FALSE.
        END IF

        write_disable = .FALSE.
        IF (found(79,key_char//' set write disable',line,.FALSE.)) THEN
          write_disable = .TRUE.
            write(*,'(1A)') &
              '  [R] : Set write disable, write_disable = .TRUE.'
        END IF

        write_iter_disable = .FALSE.
        IF (found(79,key_char//' set write iter disable',line,.FALSE.)) THEN
          write_iter_disable = .TRUE.
            write(*,'(1A)') &
              '  [R] : Set write iter disable, write_iter_disable = .TRUE.'
        END IF

        read_external_input = .TRUE.
        IF (found(79,key_char//' read external input',line,.FALSE.)) THEN
          read(unit = 79, fmt = '(l1)') read_external_input
          write(unit = *, fmt = *) '  [R] : read_external_input = ', read_external_input
        END IF

!       count the number of bc-entries (sections and entries)
        nbc_sections = 0
        nbc_data = 0
#ifndef noHDF
        if (h5parse_use_hdf5_datafile) then
          do i = 1, size(pv_name)
          do k = 1, size(bc_name)
              j=1
              write(full_bc_name, '(A,"_",A,"_",I0)') pv_name(i), bc_name(k), j
                do while (h5parse_check_dataset_exist("bc/"//full_bc_name))
                    nbc_sections = nbc_sections + 1
                    nbc_data = nbc_data + h5parse_read_dimension_size_for_dataset("bc/"//full_bc_name)
                    j = j+1
                    write(full_bc_name, '(A,"_",A,"_",I0)') pv_name(i), bc_name(k), j
                end do
          end do
          end do
        else
#endif
        REWIND 79
10      READ(79,'(1A)',end=11) line
        i = 0
        IF (head_active) i = i + locstr(line,key_char//' head bcd')
        IF (head_active) i = i + locstr(line,key_char//' head bcn')
        IF (head_active) i = i + locstr(line,key_char//' head bcw')
        IF (pres_active) i = i + locstr(line,key_char//' pres bcd')
        IF (pres_active) i = i + locstr(line,key_char//' pres bcn')
        IF (pres_active) i = i + locstr(line,key_char//' pres bcw')
        IF (temp_active) i = i + locstr(line,key_char//' temp bcd')
        IF (temp_active) i = i + locstr(line,key_char//' temp bcn')
        IF (trans_active) i = i + locstr(line,key_char//' conc bcd')
        IF (trans_active) i = i + locstr(line,key_char//' conc bcn')
!       bc line?
        IF (i==1) THEN
!         increase number of sections (reading later)
          nbc_sections = nbc_sections + 1
          CALL get_arg('simple',line,i,j)
          IF (i>=1 .AND. j>=i) THEN
            k = read_direction(line(i:j))
            tmplen = 0
            IF (k==1 .OR. k==2) tmplen = j0*k0
            IF (k==3 .OR. k==4) tmplen = i0*k0
            IF (k==5 .OR. k==6) tmplen = i0*j0
          ELSE
            CALL get_arg('records',line,i,j)
            IF (i<1 .OR. j<i) THEN
              READ(79,*,err=205,end=205) tmplen
            ELSE
              READ(line(i:j),*) tmplen
            END IF
          END IF
!         increase number of entries
          nbc_data = nbc_data + tmplen
!debug          if (linfos(1).ge.2) write(*,'(1A,1I6)')
!debug     &      '  [I] : BC block found, but reading later, size=',tmplen
        END IF
!     read next line, up to the end of file
        GO TO 10
!     restart file, because of the end of file here
#ifndef noHDF
        end if
#endif
11      REWIND 79
!     --- end switches ---

        IF (head_active .OR. pres_active) THEN
!        control environment for head/pres based computation
          WRITE(*,*) ' '
          WRITE(*,*) '  reading flow parameters'
          WRITE(*,*) ' '

          aparf = 1.0D0
          IF (found(79,key_char//' lsolvef',line,.FALSE.)) THEN
            READ(79,'(A)') line
            CALL read_solvpar(line,errf,controlf,lmaxitf,ismpl)
            WRITE(*,'(1A,1I4,1A)') &
              '  [R] : flow linear solver (control=', controlf, ')'
          ELSE
            IF (found(79,key_char//' error lsolvef',line,.TRUE.)) THEN
              READ(79,*,err=210,end=210) errf
              WRITE(*,*) ' [R] : error for flow linear solver'
            END IF
            IF (found(79,key_char//' maxiter lsolvef',line,.TRUE.)) THEN
              READ(79,*,err=211,end=211) lmaxitf
              WRITE(*,*) ' [R] : max-iter for flow linear solver'
            END IF
            IF (found(79,key_char//' name lsolvef',line,.TRUE.)) THEN
              READ(79,'(A)') line
              CALL read_solver(line,i,ismpl)
              IF (i==-1) READ(line,*,err=212,end=212) i
              WRITE(*,*) ' [R] : solver name for flow linear solver'
            END IF
            IF (found(79,key_char//' criteria lsolvef',line,.TRUE.)) THEN
              READ(79,'(A)') line
              CALL read_criteria(line,j,ismpl)
              IF (j==-1) READ(line,*,err=213,end=213) j
              WRITE(*,*) ' [R] : criteria for flow linear solver'
            END IF
            IF (found(79,key_char//' precondition lsolvef',line,.TRUE.)) THEN
              READ(79,'(A)') line
              CALL read_preco(line,k,ismpl)
              IF (k==-1) READ(line,*,err=214,end=214) k
              WRITE(*,*) ' [R] : precondition for flow linear solver'
            END IF
            CALL encntrl3(controlf,i,j,k)
          END IF

          IF (found(79,key_char//' nliterf',line,.TRUE.)) THEN
            dtmp = -1.0D0
            READ(79,*,err=101,end=101) nltolf, dtmp
101         WRITE(*,*) ' [R] : flow nonlinear iteration tolerance'
            IF (dtmp<0.0D0) THEN
              WRITE(*,'(1A)') 'error: not reading relaxation factor!'
              STOP
            END IF
            IF (nladapt==1) THEN
              nlmaxf = dtmp
              WRITE(*,*) ' [I] : max flow change/iteration'
            ELSE
              nlrelaxf = dtmp
              WRITE(*,*) ' [I] : flow nonlinear relaxation factor'
            END IF
          END IF
        ELSE
          nltolf = 1.D-30
          nltols = 1.D-30
          WRITE(*,*) ' '
          WRITE(*,*) ' <D> : no flow [disabled]'
          WRITE(*,*) ' '
        END IF

        IF (temp_active) THEN
          WRITE(*,*) ' '
          WRITE(*,*) '  reading heat transport parameters'
          WRITE(*,*) ' '
          apart = 1.0D0
          IF (found(79,key_char//' lsolvet',line,.FALSE.)) THEN
            READ(79,'(A)') line
            CALL read_solvpar(line,errt,controlt,lmaxitt,ismpl)
            WRITE(*,'(1A,1I4,1A)') &
              '  [R] : temperature linear solver (control=', controlt,')'
          ELSE
            IF (found(79,key_char//' error lsolvet',line,.TRUE.)) THEN
              READ(79,*,err=220,end=220) errt
              WRITE(*,*) ' [R] : error for temperature linear solver'
            END IF
            IF (found(79,key_char//' maxiter lsolvet',line,.TRUE.)) THEN
              READ(79,*,err=221,end=221) lmaxitt
              WRITE(*,*) &
                ' [R] : max-iter for temperature linear solver'
            END IF
            IF (found(79,key_char//' name lsolvet',line,.TRUE.)) THEN
              READ(79,'(A)') line
              CALL read_solver(line,i,ismpl)
              IF (i==-1) READ(line,*,err=222,end=222) i
              WRITE(*,*) &
                ' [R] : solver name for temperature linear solver'
            END IF
            IF (found(79,key_char//' criteria lsolvet',line,.TRUE.)) THEN
              READ(79,'(A)') line
              CALL read_criteria(line,j,ismpl)
              IF (j==-1) READ(line,*,err=223,end=223) j
              WRITE(*,*) &
                ' [R] : criteria for temperature linear solver'
            END IF
            IF (found(79,key_char//' precondition lsolvet',line,.TRUE.)) THEN
              READ(79,'(A)') line
              CALL read_preco(line,k,ismpl)
              IF (k==-1) READ(line,*,err=224,end=224) k
              WRITE(*,*) &
                ' [R] : precondition for temperature linear solver'
            END IF
            CALL encntrl3(controlt,i,j,k)
          END IF

          IF (found(79,key_char//' nlitert',line,.TRUE.)) THEN
            dtmp = -1.0D0
            READ(79,*,err=102,end=102) nltolt, dtmp
102         WRITE(*,*) &
              ' [R] : temperature nonlinear iteration tolerance'
            IF (dtmp<0.0D0) THEN
              WRITE(*,'(1A)') 'error: not reading relaxation factor!'
              STOP
            END IF
            IF (nladapt==1) THEN
              nlmaxt = dtmp
              WRITE(*,*) ' [I] : max temperature change/iteration'
            ELSE
              nlrelaxt = dtmp
              WRITE(*,*) &
                ' [I] : temperature nonlinear relaxation factor'
            END IF
          END IF
        ELSE
          nltolt = 1.D-30
          WRITE(*,*) ' '
          WRITE(*,*) ' <D> : no temperature [disabled]'
          WRITE(*,*) ' '
        END IF

        ntrac = 0
        nchem = 0
        ntrans = ntrac + nchem
        IF (trans_active) THEN
          WRITE(*,*) ' '
          WRITE(*,*) '  reading chemical transport parameters'
          WRITE(*,*) ' '

          IF (found(79,key_char//' ntrans',line,.FALSE.)) THEN
            READ(79,*,err=250,end=250) ntrac, nchem
            ntrans = ntrac + nchem
            WRITE(*,'(1A,2I4)') &
              '  [R] : tracers, reactive components=', ntrac, nchem
          ELSE
            WRITE(*,*) ' <D> : no tracers or reactive components'
          END IF
          IF (ntrac<=0) trac_active = .FALSE.
          IF (nchem<=0) chem_active = .FALSE.
          trans_active = .FALSE.
          IF (trac_active .OR. chem_active) trans_active = .TRUE.
        END IF
!       [head/pres,temp, conc...]
        conv_hmax = 3 + ntrans

        IF (trans_active) THEN
          aparc = 1.0D0
          IF (found(79,key_char//' lsolvec',line,.FALSE.)) THEN
            READ(79,'(A)') line
            CALL read_solvpar(line,errc,controlc,lmaxitc,ismpl)
            WRITE(*,'(1A,1I4,1A)') &
              '  [R] : transport linear solver (control=', controlc, &
              ')'
          ELSE
            IF (found(79,key_char//' error lsolvec',line,.TRUE.)) THEN
              READ(79,*,err=230,end=230) errc
              WRITE(*,*) ' [R] : error for transport linear solver'
            END IF
            IF (found(79,key_char//' maxiter lsolvec',line,.TRUE.)) THEN
              READ(79,*,err=231,end=231) lmaxitc
              WRITE(*,*) &
                ' [R] : max-iter for transport linear solver'
            END IF
            IF (found(79,key_char//' name lsolvec',line,.TRUE.)) THEN
              READ(79,'(A)') line
              CALL read_solver(line,i,ismpl)
              IF (i==-1) READ(line,*,err=232,end=232) i
              WRITE(*,*) &
                ' [R] : solver name for transport linear solver'
            END IF
            IF (found(79,key_char//' criteria lsolvec',line,.TRUE.)) THEN
              READ(79,'(A)') line
              CALL read_criteria(line,j,ismpl)
              IF (j==-1) READ(line,*,err=233,end=233) j
              WRITE(*,*) &
                ' [R] : criteria for transport linear solver'
            END IF
            IF (found(79,key_char//' precondition lsolvec',line,.TRUE.)) THEN
              READ(79,'(A)') line
              CALL read_preco(line,k,ismpl)
              IF (k==-1) READ(line,*,err=234,end=234) k
              WRITE(*,*) &
                ' [R] : precondition for transport linear solver'
            END IF
            CALL encntrl3(controlc,i,j,k)
          END IF

          IF (found(79,key_char//' nliterc',line,.TRUE.)) THEN
            dtmp = -1.0D0
            READ(79,*,err=103,end=103) nltolc, dtmp
103         WRITE(*,*) &
              ' [R] : transport nonlinear iteration tolerance'
            IF (dtmp<0.0D0) THEN
              WRITE(*,'(1A)') 'error: not reading relaxation factor!'
              STOP
            END IF
            IF (nladapt==1) THEN
              nlmaxc = dtmp
              WRITE(*,*) ' [I] : max transport change/iteration'
            ELSE
              nlrelaxc = dtmp
              WRITE(*,*) &
                ' [I] : transport nonlinear relaxation factor'
            END IF
          END IF
        ELSE
          nltolc = 1.D-30
          WRITE(*,*) ' '
          WRITE(*,*) ' <D> : no transport [disabled]'
          WRITE(*,*) ' '
        END IF


        WRITE(*,*) ' '
        grav = 9.81D0
        IF (found(79,key_char//' grav',line,.FALSE.)) THEN
          READ(79,*,err=251,end=251) grav
          WRITE(*,'(1a,1e12.4,1a)') '  [R] : grav', grav, &
            ' (>1.0d-30)'
          grav = max(grav,1.D-30)
        ELSE
          WRITE(*,'(1a,1e12.4)') '  <D> : grav = ', grav
        END IF

        hpf = 0.0D0
        IF (found(79,key_char//' hpf',line,.FALSE.)) THEN
          READ(79,*,err=252,end=252) hpf
          WRITE(*,*) ' [R] : hpf, fluid heat production'
        ELSE
          WRITE(*,*) ' <D> : no fluid heat production'
        END IF

!#ifdef head_base
        IF (found(79,key_char//' rref',line,.FALSE.)) THEN
          READ(79,*) rref
          WRITE(*,'(1a,1e12.4)') '  [R] : rref, reference density = ' &
            , rref
        ELSE
          rref = 998.D0
          WRITE(*,'(1a,1e12.4)') '  <D> : rref, reference density = ' &
            , rref
        END IF
!#endif

!      if (found(79,'? tref',line,.false.)) then
!        read(79,*) tref
!        write(*,'(1a,1e12.4)')
!     &     '  [R] : tref, referenCe temperature = ',tref
!      else
        tref = 20.D0
        WRITE(*,'(a,e12.4)') '  <D> : tref, reference temperature = ' &
          , tref
!      endif

        rhom = 2500D0
        cma1 = 1.D0
        cma2 = 0.D0
        cma3 = 0.D0
        IF (found(79,key_char//' rhocm',line,.FALSE.)) THEN
          READ(79,*,err=253,end=253) rhom, cma1, cma2, cma3
          WRITE(*,*) ' [R] : rhom, rock heat capacity model'
        ELSE
          WRITE(*,'(1A)') '  <D> : fixed rhocm, defined by unit'
        END IF
!     normalise
        IF (test_null(cma1)) THEN
          WRITE(*,'(1A)') 'error: "cma1" equals to zero !!!'
          STOP
        ELSE
          cma2 = cma2/cma1
          cma3 = cma3/cma1
          cma1 = 1.D0
        END IF

!       before memory allocating, counting bc-tp entries ('ngsmax') needed
        ngsmax = 1
        nbctp = 0
        IF (found(79,key_char//' bc time periods',line,.FALSE.)) THEN
          CALL get_arg('records',line,i,j)
          IF (i<1 .OR. j<i) THEN
            READ(79,*,err=254,end=254) l
          ELSE
            READ(line(i:j),*) l
          END IF
          nbctp = l
          DO i = 1, l
!           [bctp-unit, number of entries]
            READ(79,*,err=255,end=255) k, j
            nbctp = max(nbctp,k)
            ngsmax = max(ngsmax,j)
            DO k = 1, j
!              skip 'j' lines
              READ(79,'(1A)') line
            END DO
          END DO
        END IF

!       before memory allocating, counting bc-units entries ('nunits') nee
        nunits = 0
        IF (found(79,key_char//' bcunits',line,.FALSE.)) THEN
          CALL get_arg('records',line,i,j)
          IF (i<1 .OR. j<i) THEN
            READ(79,*,err=256,end=256) tmplen
          ELSE
            READ(line(i:j),*) tmplen
          END IF
          DO i = 1, tmplen
            READ(79,*,err=257,end=257) k, dtmp, ctmp
            nunits = max(k,nunits)
          END DO
        END IF

!----------------------------------------------------------------
!       zones
        WRITE(*,*) ' '
        WRITE(*,*) '  reading zone parameter '
        WRITE(*,*) ' '

!        allocation here, because of the need for "nunits"
          ALLOCATE(uindex(i0,j0,k0))
          memory = memory + i0*j0*k0
#ifndef noHDF
        if (h5parse_use_hdf5_datafile) then
            found_marker = .true.
            call h5parse_read_3d_integer_dataset("uindex", uindex)
        else
#endif
        IF (found(79,key_char//' uindex',line,.TRUE.)) THEN
          IF (no_ext_link_int(i0,j0,k0,uindex,'uindex',line)) &
            READ(79,*,err=110,end=110) (((uindex(i,j,k),i=1, &
            i0),j=1,j0),k=1,k0)
          found_marker = .true.
        end if
#ifndef noHDF
        end if
#endif
        if (found_marker) then
          found_marker = .false.
          WRITE(*,*) ' [R] : "uindex" - unit index number, unit-cell assignment'

          maxunits = uindex(1,1,1)
          minunits = uindex(1,1,1)
          DO k = 1, k0
            DO j = 1, j0
              DO i = 1, i0
                maxunits = max(uindex(i,j,k),maxunits)
                minunits = min(uindex(i,j,k),minunits)
              END DO
            END DO
          END DO

          IF ((maxunits>i0*j0*k0) .OR. (minunits<1)) THEN
            WRITE(*,'(A,I7,A,I7,A)') &
              'error: unit index out of range [', minunits, ',', &
              maxunits, '] !'
            STOP
          END IF

!        also used as maximum for BC units ("bc_maxunits")
          nunits = max(maxunits,nunits)

!        check thickness for each unit layer
          CALL check_units(ismpl)
        END IF

        nbh_logs = 0
        IF (found(79,key_char//' borehole log',line,.FALSE.)) THEN
          CALL get_arg('records',line,i,j)
          IF (i<1 .OR. j<i) THEN
            READ(79,*,err=254,end=254) nbh_logs
          ELSE
            READ(line(i:j),*) nbh_logs
          END IF
        END IF

! ------------------
!     initialisation for linear system solver
!     init array for prozessor grid, see more in 'solve/omp_preconditioniers.f'
        CALL par_init2(I0,J0,K0)

! ------------------
!     memory managment
        WRITE(*,*) ' '
        CALL alloc_arrays(ismpl)

! ------------------
!     reading later
        ndata = 0
        WRITE(*,*) ' [R] : "units" - unit (rock) properties'
        ALLOCATE(datmp(nprop_load,maxunits))
!       init property vaues to the default, needed when not enough values readed
        DO j = 1, maxunits
          DO i = 1, nprop_load
            datmp(i,j) = prop_default(i)
          END DO
        END DO
 
#ifndef noHDF
       if (h5parse_use_hdf5_datafile) then
          j = h5parse_read_dimension_size_for_dataset("units")
          call h5parse_read_2d_double_dataset("units",datmp(1:j,:))
          found_marker = .true.
       else
#endif
       if (found(79,key_char//' units',line,.TRUE.)) then
          IF (no_ext_link(nprop_load,maxunits,1,datmp,'units',line)) &
              THEN
            CALL read_array(79,nprop_load,maxunits,datmp,key_char//' units', &
              ismpl)
          END IF
          found_marker = .true.
       end if
#ifndef noHDF
       end if
#endif
       IF (found_marker) THEN
          found_marker = .false.
          IF (nprop_load==lastidx-firstidx+1) THEN
!          load all, dense entries (no specific index)
            DO i = 1, maxunits
              DO j = 1, nprop_load
                propunit(i,firstidx-1+j,ismpl) = datmp(j,i)
              END DO
            END DO
          ELSE
!          needs to handle manual reading with specific index
            WRITE(*,'(1A)') &
              'error: bug (1) in "read_model.f", ask AW!'
            STOP
          END IF

          IF (linfos(1)>=2 .AND. maxunits<=64) WRITE(*, &
            '(1A18,'//c_npropunit//'("    ",1A4,"    "),1A6)') '  unit properties:', &
            (properties(i),i=firstidx,lastidx), ' unit#'
          DO i = 1, maxunits
!         because of logarithimc scale, suppress zeros
            propunit(i,idx_por,ismpl) = max(prop_min(idx_por), &
              propunit(i,idx_por,ismpl))
            propunit(i,idx_an_kx,ismpl) = max(prop_min(idx_an_kx), &
              propunit(i,idx_an_kx,ismpl))
            propunit(i,idx_an_ky,ismpl) = max(prop_min(idx_an_ky), &
              propunit(i,idx_an_ky,ismpl))
            propunit(i,idx_kz,ismpl) = max(prop_min(idx_kz), &
              propunit(i,idx_kz,ismpl))
!?          propunit(i,idx_comp,ismpl) = max(prop_min(idx_comp), &
!?            propunit(i,idx_comp,ismpl))
            propunit(i,idx_an_lx,ismpl) = max(prop_min(idx_an_lx), &
              propunit(i,idx_an_lx,ismpl))
            propunit(i,idx_an_ly,ismpl) = max(prop_min(idx_an_ly), &
              propunit(i,idx_an_ly,ismpl))
            propunit(i,idx_lz,ismpl) = max(prop_min(idx_lz), &
              propunit(i,idx_lz,ismpl))
            propunit(i,idx_q,ismpl) = max(prop_min(idx_q), &
              propunit(i,idx_q,ismpl))
            propunit(i,idx_rc,ismpl) = max(prop_min(idx_rc), &
              propunit(i,idx_rc,ismpl))
            propunit(i,idx_df,ismpl) = max(prop_min(idx_df), &
              propunit(i,idx_df,ismpl))
            propunit(i,idx_ec,ismpl) = max(prop_min(idx_ec), &
              propunit(i,idx_ec,ismpl))
!?          propunit(i,idx_lc,ismpl) = max(prop_min(idx_lc), &
!?            propunit(i,idx_lc,ismpl))
            propunit(i,idx_s_nr,ismpl) = max(prop_min(idx_s_nr), &
              propunit(i,idx_s_nr,ismpl))
            propunit(i,idx_s_wr,ismpl) = max(prop_min(idx_s_wr), &
              propunit(i,idx_s_wr,ismpl))
            IF (linfos(1)>=2) WRITE(*,'(18X,'//c_npropunit//'e12.4,i6)') &
              (propunit(i,j,ismpl),j=firstidx,lastidx), i
          END DO

        END IF
        DEALLOCATE(datmp)

! ------------------

        IF (nbh_logs>0) THEN
          IF (found(79,key_char//' borehole log',line,.FALSE.)) THEN
            CALL get_arg('records',line,i,j)
            IF (i<1 .OR. j<i) THEN
              READ(79,*,err=254,end=254) i
            ELSE
              READ(line(i:j),*) i
            END IF
            IF (nbh_logs/=i) THEN
              WRITE(*,*) 'error: more than one "borehole logs" section defined'
              STOP
            END IF
            DO i = 1, nbh_logs
!             [x-position, y-position, bh-name]
              READ(79,*,err=255,end=255) ibh_pos(1,i),ibh_pos(2,i),cbh_name(i)
            END DO
          END IF
        END IF

! ------------------

        WRITE(*,*) ' '
        WRITE(*,*) '  reading/preprocessing arrays'
        WRITE(*,*) ' '

        ijk = i0*j0*k0

#ifndef noHDF
        if (h5parse_use_hdf5_datafile) then
          call h5parse_read_1d_double_dataset("grid/delx",delx,i0)
          found_marker = .true.
        else
#endif
        IF (found(79,key_char//' delx',line,.TRUE.)) THEN
          IF (no_ext_link(i0,1,1,delx,'delx',line)) READ(79,*, &
            err=260,end=260) (delx(i),i=1,i0)
          found_marker = .true.
        END IF
#ifndef noHDF
        end if
#endif
        if (found_marker) then
            found_marker = .false.
            WRITE(*,*) ' [R] : delx'
        end if

#ifndef noHDF
        if (h5parse_use_hdf5_datafile) then
          call h5parse_read_1d_double_dataset("grid/dely",dely,j0)
          found_marker = .true.
        else
#endif
        IF (found(79,key_char//' dely',line,.TRUE.)) THEN
          IF (no_ext_link(1,j0,1,dely,'dely',line)) READ(79,*, &
            err=261,end=261) (dely(j),j=1,j0)
          found_marker = .true.
        END IF
#ifndef noHDF
        end if
#endif
        if (found_marker) then
            found_marker = .false.
            WRITE(*,*) ' [R] : dely'
        end if
#ifndef noHDF
        if (h5parse_use_hdf5_datafile) then
          call h5parse_read_1d_double_dataset("grid/delz",delz,k0)
          found_marker = .true.
        else
#endif
        IF (found(79,key_char//' delz',line,.TRUE.)) THEN
          IF (no_ext_link(1,1,k0,delz,'delz',line)) READ(79,*, &
            err=262,end=262) (delz(k),k=1,k0)
          found_marker = .true.
        END IF
#ifndef noHDF
        end if
#endif
        if (found_marker) then
            found_marker = .false.
            WRITE(*,*) ' [R] : delz'
        end if

!     compute absolute xyz-positions
        delxa(1) = 0.5D0*delx(1)
        DO i = 2, i0
          delxa(i) = delxa(i-1) + 0.5D0*(delx(i-1)+delx(i))
        END DO
        delya(1) = 0.5D0*dely(1)
        DO j = 2, j0
          delya(j) = delya(j-1) + 0.5D0*(dely(j-1)+dely(j))
        END DO
        delza(1) = 0.5D0*delz(1)
        DO k = 2, k0
          delza(k) = delza(k-1) + 0.5D0*(delz(k-1)+delz(k))
        END DO

!     initial values
        WRITE(*,*) ' '
        WRITE(*,*) '  reading initial values '
        WRITE(*,*) ' '

        head_needed = .false.
#ifdef head_base
        ! Head input forced for head computation
        head_needed = .true.
#endif
        ! Is pres2head or head2pres needed for initial data
        is_init_flow_trafo_needed = .true.

#ifndef noHDF
        if (h5parse_use_hdf5_datafile) then
          if (h5parse_check_dataset_exist("init/head") .or. head_needed) then
            call h5parse_read_3d_double_dataset("init/head", head(:,:,:,ismpl))
            found_marker = .true.
          end if
        else
#endif
        IF (found(79,key_char//' head init',line,head_needed)) THEN
          found_marker = .true.
          IF (no_ext_link(i0,j0,k0,head(1,1,1, &
              ismpl),'head',line)) READ(79,*,err=265,end=265) ((( &
              head(i,j,k,ismpl),i=1,i0),j=1,j0),k=1,k0)
        END IF
#ifndef noHDF
        end if
#endif
        if (found_marker) then
          found_marker = .false.

          ! Head input for pres_base: No trafo needed
          if (.not.head_needed) is_init_flow_trafo_needed = .false.
          write(*,*) ' [R] : head, in [m]'
        end if

#ifndef noHDF
        if (h5parse_use_hdf5_datafile) then
            call h5parse_read_3d_double_dataset("init/temp", temp(:,:,:,ismpl))
        else
#endif
        IF (found(79,key_char//' temp init',line,.TRUE.)) THEN
          IF (no_ext_link(i0,j0,k0,temp(1,1,1, &
              ismpl),'temp',line)) READ(79,*,err=266,end=266) ((( &
              temp(i,j,k,ismpl),i=1,i0),j=1,j0),k=1,k0)
          END IF
#ifndef noHDF
        end if
#endif
        WRITE(*,*) ' [R] : temp, in [degree Celsius]'

#ifndef noHDF
        if (h5parse_use_hdf5_datafile) then
          if (h5parse_check_dataset_exist("init/pres") .or. .not. head_needed) then
            call h5parse_read_3d_double_dataset("init/pres", pres(:,:,:,ismpl))
            found_marker = .true.
          end if
        else
#endif
        IF (found(79,key_char//' pres init',line,.NOT.head_needed)) THEN
          found_marker = .true.
          IF (no_ext_link(i0,j0,k0,pres(1,1,1, &
            ismpl),'pres',line)) READ(79,*,err=267,end=267) ((( &
            pres(i,j,k,ismpl),i=1,i0),j=1,j0),k=1,k0)
        END IF
#ifndef noHDF
        end if
#endif
        if (found_marker) then
          found_marker = .false.

          ! Pres input for head_base: no trafo needed
          if (head_needed) is_init_flow_trafo_needed = .false.
        end if

!       convert [MPa] into [Pa]
        CALL dscal(i0*j0*k0,pa_conv,pres(1,1,1,ismpl),1)
        WRITE(*,*) ' [R] : pres, in [MPa]'

        IF (trac_active) THEN
          DO tracer = 1, ntrac
            WRITE(strng,'(1A,1I4.4,1A)') key_char//' tracer', tracer, &
              ' init'
            CALL chln(strng,i1,i2)
            IF (found(79,strng(i1:i2),line,.TRUE.)) THEN
              IF (no_ext_link(i0,j0,k0,conc(1,1,1,tracer, &
                ismpl),strng(i1+2:i2-5),line)) READ(79,*,err=268, &
                end=268) (((conc(i,j,k,tracer,ismpl),i=1, &
                i0),j=1,j0),k=1,k0)
              WRITE(*,'(1A,1I4.4,1A)') '  [R] : tracer', tracer, &
                ', in [Mol/l]'
            END IF
          END DO
        END IF

!----------------------------------------------------------------

!       boundary conditions if not on computional domain boundaries
        WRITE(*,*) ' '
        WRITE(*,*) '  reading boundary conditions '
        WRITE(*,*) ' '

!       initialize current number of bc-units
        bc_maxunits = 0
        DO i = 1, nunits
!         init undefined max values
          propunit(i,idx_hbc,ismpl) = const_dble(3)
          propunit(i,idx_tbc,ismpl) = const_dble(3)
          propunit(i,idx_cbc,ismpl) = const_dble(3)
          propunit(i,idx_snbc,ismpl) = const_dble(3)
          propunit(i,idx_ebc,ismpl) = const_dble(3)
        END DO

!       read the entries of irregular boundary conditions
            REWIND 79
            posi = 0
            ilost = 0
#ifndef noHDF
            if (h5parse_use_hdf5_datafile) then
                do ijk = 1, 2
                do i = 1, size(pv_name)
                do k = 1, size(bc_name)
                    j=1
                    write(full_bc_name, '(A,"_",A,"_",I0)') pv_name(i), bc_name(k), j
                    do while (h5parse_check_dataset_exist("bc/"//full_bc_name))
                        if ((h5parse_check_attr_exist("simple", "bc/i"//full_bc_name) .and. ijk == 2) .or. &
                          & (.not. h5parse_check_attr_exist("simple", "bc/i"//full_bc_name) .and. ijk == 1)) then
                            CALL read_bc(79,full_bc_name,i,k,posi,ilost,ismpl)
                        end if
                        j = j+1
                        write(full_bc_name, '(A,"_",A,"_",I0)') pv_name(i), bc_name(k), j
                    end do
                end do
                end do
                end do
            else
#endif
            DO k = 1, nbc_sections
20            READ(79,'(1A)',err=21,end=21) line
              i = 0
              j = 0
              IF (locstr(line,key_char//' head bcd')==1 .AND. head_active) THEN
                i = pv_head
                j = bt_diri
              ELSE IF (locstr(line,key_char//' head bcn')==1 .AND. head_active) THEN
                i = pv_head
                j = bt_neum
              ELSE IF (locstr(line,key_char//' head bcw')==1 .AND. head_active) THEN
                i = pv_head
                j = bt_neuw
              ELSE IF (locstr(line,key_char//' pres bcd')==1 .AND. pres_active) THEN
                i = pv_pres
                j = bt_diri
              ELSE IF (locstr(line,key_char//' pres bcn')==1 .AND. pres_active) THEN
                i = pv_pres
                j = bt_neum
              ELSE IF (locstr(line,key_char//' pres bcw')==1 .AND. pres_active) THEN
                i = pv_pres
                j = bt_neuw
              ELSE IF (locstr(line,key_char//' temp bcd')==1 .AND. temp_active) THEN
                i = pv_temp
                j = bt_diri
              ELSE IF (locstr(line,key_char//' temp bcn')==1 .AND. temp_active) THEN
                i = pv_temp
                j = bt_neum
              ELSE IF (locstr(line,key_char//' conc bcd')==1 .AND. trans_active) THEN
                i = pv_conc
                j = bt_diri
              ELSE IF (locstr(line,key_char//' conc bcn')==1 .AND. trans_active) THEN
                i = pv_conc
                j = bt_neum
              END IF
!             bc line?
              IF (i>0 .AND. j>0) THEN
                CALL read_bc(79,line,i,j,posi,ilost,ismpl)
              ELSE
!               try next line
                GO TO 20
              END IF
            END DO
#ifndef noHDF
            end if
#endif
!           sanity check
21          IF (posi+ilost/=nbc_data) THEN
              WRITE(*,'(1A)') &
                'error: lost some BC entries, in "read_model"!'
              STOP
            END IF
!           update 'nbc_data'
            nbc_data = nbc_data - ilost

!     sort "*bc_data" and setup the first/last-index for all physical va
            CALL sort_bc(ismpl)

!     print out dependency betwee bc and bctp
            IF (linfos(1)>=1) THEN
              CALL show_bcdep(ismpl)
              WRITE(*,*) ' '
            END IF
!----------------------------------------------------------------

            IF (bc_maxunits>0) THEN
              IF (found(79,key_char//' bcunits',line,.TRUE.)) THEN

                CALL get_arg('records',line,i,j)
                IF (i<1 .OR. j<i) THEN
                  READ(79,*,err=256,end=256) tmplen
                ELSE
                  READ(line(i:j),*) tmplen
                END IF

!         mark the needed/used bc-units
                ALLOCATE(ltmp(bc_maxunits,2))
                DO i = 1, bc_maxunits
!            needed
                  ltmp(i,1) = .FALSE.
!            used
                  ltmp(i,2) = .FALSE.
                END DO
                DO i = 1, nbc_data
                  k = ibc_data(i,cbc_bcu)
                  IF (k>0) ltmp(k,1) = .TRUE.
                END DO

                IF (linfos(1)>=2) WRITE(*,'(A/2A12,A6)') &
                  '  bc unit properties:', '      unit# ', &
                  '      value ', ' type '
                DO i = 1, tmplen
                  READ(79,*,err=257,end=257) k, dtmp, ctmp
                  IF ((k>bc_maxunits) .OR. (k<1)) THEN
                    WRITE(*,'(2A,1I3,1A,1I3,1A)') 'error: &
                      &in section "bcunits", unit number &
                      &out of', ' range, (', k, ') at line ', i, &
                      ' !!!'
                    STOP
                  END IF
                  IF ( .NOT. ltmp(k,1)) THEN
                    WRITE(*,'(2A,1I3,1A,1I3,1A)') &
                      'warning: in section "bcunits", unit number', &
                      ' not used, (', k, ') at line ', i, '.'
                  END IF
                  ltmp(k,2) = .TRUE.
                  IF ((ctmp/='head') .AND. (ctmp/='pres') .AND. &
                      (ctmp/='temp') .AND. (ctmp/='conc')) THEN
                    WRITE(*,'(1A,1A4,1A,1I3,1A)') 'error: &
                      &in section "bcunits", bc type not &
                      &allowed, "', ctmp, '" at line ', i, ' !!!'
                    STOP
                  END IF

!            sanity checks
                  IF (ctmp=='head' .AND. propunit(k,idx_hbc,ismpl)/= &
                      const_dble(3)) THEN
                    WRITE(*,'(1A,1I7,1A)') 'error: head (/ pres) bc unit ', k, &
                      ' multi defined !!!'
                    STOP
                  END IF
                  IF (ctmp=='pres' .AND. propunit(k,idx_hbc,ismpl)/= &
                      const_dble(3)) THEN
                    WRITE(*,'(1A,1I7,1A)') 'error: pres (/ head) bc unit ', k, &
                      ' multi defined !!!'
                    STOP
                  END IF
                  IF (ctmp=='temp' .AND. propunit(k,idx_tbc,ismpl)/= &
                      const_dble(3)) THEN
                    WRITE(*,'(1A,1I7,1A)') 'error: temp bc unit ', k, &
                      ' multi defined !!!'
                    STOP
                  END IF
                  IF (ctmp=='conc' .AND. propunit(k,idx_cbc,ismpl)/= &
                      const_dble(3)) THEN
                    WRITE(*,'(1A,1I7,1A)') 'error: conc bc unit ', k, &
                      ' multi defined !!!'
                    STOP
                  END IF
!            update values
                  IF (ctmp=='head') propunit(k,idx_hbc,ismpl) = max(0.D0,dtmp)
!                 convert [MPa] into [Pa]
                  IF (ctmp=='pres') propunit(k,idx_hbc,ismpl) = max(0.D0,dtmp*pa_conv)
                  IF (ctmp=='temp') propunit(k,idx_tbc,ismpl) = dtmp
                  IF (ctmp=='conc') propunit(k,idx_cbc,ismpl) = max(0.D0,dtmp)
                  IF (linfos(1)>=2) WRITE(*,'(I11,1X,1e12.4,1X,A4)') &
                    k, dtmp, ctmp
                END DO
                WRITE(*,'(A,I7)') &
                  '  [R] : bc unit properties, records=', tmplen
                DO k = 1, bc_maxunits
                  IF (ltmp(k,1) .AND. .NOT. ltmp(k,2)) THEN
                    WRITE(*,'(2A,1I3,1A,1I3,1A)') &
                      'error: in section "bcunits", unit number', &
                      ' not defined, (', k, ') !!!'
                    STOP
                  END IF
                END DO
                DEALLOCATE(ltmp)
              END IF
            END IF
            DO i = 1, nunits
!       clear unread values
              IF (propunit(i,idx_hbc,ismpl)==const_dble(3)) propunit(i,idx_hbc,ismpl) = 0.D0
              IF (propunit(i,idx_tbc,ismpl)==const_dble(3)) propunit(i,idx_tbc,ismpl) = 0.D0
              IF (propunit(i,idx_cbc,ismpl)==const_dble(3)) propunit(i,idx_cbc,ismpl) = 0.D0
              IF (propunit(i,idx_snbc,ismpl)==const_dble(3)) propunit(i,idx_snbc,ismpl) = 0.D0
              IF (propunit(i,idx_ebc,ismpl)==const_dble(3)) propunit(i,idx_ebc,ismpl) = 0.D0
            END DO

            IF (trans_active) THEN
              IF (found(79,key_char//' transpar',line,.FALSE.)) THEN
                DO i = 1, ntrans
                  READ(79,*,err=270,end=270) diff_c(i), mmas_c(i)
                END DO
                WRITE(*,'(1A,1I6,1A)') '  [R] : ', ntrans, &
                  ' diffusion constants and mol masses'
              ELSE
                dtmp = 1.D-9
                DO i = 1, ntrans
                  diff_c(i) = dtmp
                  mmas_c(i) = mmas_nacl
                END DO
                WRITE(*,*) ' <D> : uniform diffusion constant = ', &
                  dtmp
              END IF
            END IF


          IF (found(79,key_char//' prop limit',line,.FALSE.)) THEN
            WRITE(*,*) ' [R] : property limitations'
            CALL get_arg('records',line,i,j)
            IF (i<1 .OR. j<i) THEN
              READ(79,*,err=280,end=280) tmplen
            ELSE
              READ(line(i:j),*) tmplen
            END IF
            DO i = 1, tmplen
              READ(79,'(a)',err=281,end=281) line
              CALL read_limit(line)
            END DO
          ELSE
            WRITE(*,*) ' <D> : property limitations'
          END IF

          IF (found(79,key_char//' velocity',line,.FALSE.)) THEN
             ALLOCATE(vdefault(3,nsmpl))
             
             WRITE(*,*) ' [R] : default Darcy velocity'
             READ(79,*,err=282,end=282) vdefault(1,1), vdefault(2,1), vdefault(3,1)
             WRITE(*,*) '  [R] : [x:vdefault(1,1), y:vdefault(2,1), z:vdefault(3,1)] = [', &
                  vdefault(1,1), ',', vdefault(2,1), ',', vdefault(3,1), ']'
             vdefaultswitch = .true.
             
             DO i = 2, nsmpl
                vdefault(1,i) = vdefault(1,1)
                vdefault(2,i) = vdefault(2,1)
                vdefault(3,i) = vdefault(3,1)
             END DO
          ELSE
             WRITE(*,'(1A)') &
                  '  [R] : Zero default Darcy velocity -  no # velocity!'
             vdefaultswitch = .false.
          END IF
          
!      if (found(79,key_char//' pumping',line,.false.)) then
!          Call get_arg('reCords',line,i,j)
!          if (i.lt.1.or.j.lt.i) then
!             read(79,*) templen
!          else
!             read(line(i:j),*) templen
!          endif
!         npump=templen
!       do i = 1, npump
!               read(79,*) (tocell(j,i),j=1,3), (fromcell(j,i),j=1,3), howmuch(i)
!         end do
!          write(*,'(1A,1I6,1A)')
!     &        '  [R] : ',npump,' pumping/reinjeCtion setups'
!      else
!          write(*,*) ' <D> : no pumping/reinjeCtion setups '
!      endif

!     clean node info
            DO k = 1, k0
              DO j = 1, j0
                DO i = 1, i0
                  node_info(i,j,k) = '     '
!        node_info(i,j,k)(pv_head:pv_head) = ' '
!        node_info(i,j,k)(pv_temp:pv_temp) = ' '
!        node_info(i,j,k)(pv_conc:pv_conc) = ' '
!        node_info(i,j,k)(pv_pres:pv_pres) = ' '
                END DO
              END DO
            END DO
!     init node info
            DO l = 1, nbc_data
              i = ibc_data(l,cbc_i)
              j = ibc_data(l,cbc_j)
              k = ibc_data(l,cbc_k)
              sbc = ' '
              IF (ibc_data(l,cbc_bt)==bt_diri) sbc = 'd'
              IF (ibc_data(l,cbc_bt)==bt_neum) sbc = 'n'
              IF (ibc_data(l,cbc_bt)==bt_neuw) sbc = 'w'
              node_info(i,j,k) (ibc_data(l,cbc_pv):ibc_data(l,cbc_pv)) &
                = sbc
            END DO

! ------------------

!     finish HDF5 support, when available
            CALL close_hdf5()

!     close project config file
            CLOSE(79)

            RETURN

!     error handler
110         WRITE(*,'(1A,3I7,1A)') &
              'error: reading section "uindex", at index [', i, j, k, &
              ']!'
            STOP
200         WRITE(*,'(1A)') 'error: can not read "title"!'
            STOP
201         WRITE(*,'(1A)') &
              'error: can not read number of "samples"!'
            STOP
202         WRITE(*,'(1A)') 'error: can not read "runmode"!'
            STOP
203         WRITE(*,'(1A)') &
              'error: can not read dimensions in section "grid"!'
            STOP
204         WRITE(*,'(2A)') 'error: in section "nlsolve" expecting:', &
              ' <max iter> <adapting switch = {0,1}>!'
            STOP
205         WRITE(*,'(2A)') 'error: no size specified for some ', &
              '"boundary condition definition"!'
            STOP
206         WRITE(*,'(1A)') &
              'error: can not find section "'//key_char//' PROPS=<...>"!'
            STOP
207         WRITE(*,'(1A)') &
              'error: can not find section "'//key_char//' USER=<...>"!'
            STOP
210         WRITE(*,'(1A)') &
              'error: no definition for "error lsolvef"!'
            STOP
211         WRITE(*,'(1A)') &
              'error: no definition for "maxiter lsolvef"!'
            STOP
212         WRITE(*,'(1A)') &
              'error: no definition for "name lsolvef"!'
            STOP
213         WRITE(*,'(1A)') &
              'error: no definition for "criteria lsolvef"!'
            STOP
214         WRITE(*,'(1A)') &
              'error: no definition for "precondition lsolvef"!'
            STOP
215         WRITE(*,'(1A)') &
              'error: file end after section "disable output"!'
            STOP
220         WRITE(*,'(1A)') &
              'error: no definition for "error lsolvet"!'
            STOP
221         WRITE(*,'(1A)') &
              'error: no definition for "maxiter lsolvet"!'
            STOP
222         WRITE(*,'(1A)') &
              'error: no definition for "name lsolvet"!'
            STOP
223         WRITE(*,'(1A)') &
              'error: no definition for "criteria lsolvet"!'
            STOP
224         WRITE(*,'(1A)') &
              'error: no definition for "precondition lsolvet"!'
            STOP
250         WRITE(*,'(1A)') 'error: in section "ntrans"', &
              ' expecting: <#tracer> <#chemicals>!'
            STOP
230         WRITE(*,'(1A)') &
              'error: no definition for "error lsolvec"!'
            STOP
231         WRITE(*,'(1A)') &
              'error: no definition for "maxiter lsolvec"!'
            STOP
232         WRITE(*,'(1A)') &
              'error: no definition for "name lsolvec"!'
            STOP
233         WRITE(*,'(1A)') &
              'error: no definition for "criteria lsolvec"!'
            STOP
234         WRITE(*,'(1A)') &
              'error: no definition for "precondition lsolvec"!'
            STOP
251         WRITE(*,'(1A)') &
              'error: can not read gravity in section "grav"!'
            STOP
252         WRITE(*,'(1A)') &
              'error: can not read the value in section "hpf"!'
            STOP
253         WRITE(*,'(1A)') 'error: in section "rhocm"', &
              'expecting: "rhom cma1 cma2 cma3"!'
            STOP
254         WRITE(*,'(1A)') &
              'error: no size specified for "bc time periods"!'
            STOP
255         WRITE(*,'(2A)') 'error: in data section of "bc time periods",',&
               ' expecting:',' <unit-id> <number of entries>!'
            STOP
256         WRITE(*,'(1A)') 'error: no size specified for "bcunits"!'
            STOP
257         WRITE(*,'(1A)') 'error: in data section of "bcunits"', &
              'please read manual!'
            STOP
260         WRITE(*,'(1A)') 'error: to few values in section "delx"!'
            STOP
261         WRITE(*,'(1A)') 'error: to few values in section "dely"!'
            STOP
262         WRITE(*,'(1A)') 'error: to few values in section "delz"!'
            STOP
265         WRITE(*,'(1A)') &
              'error: to few values in section "head init"!'
            STOP
266         WRITE(*,'(1A)') &
              'error: to few values in section "temp init"!'
            STOP
267         WRITE(*,'(1A)') &
              'error: to few values in section "pres init"!'
            STOP
268         WRITE(*,'(1A)') &
              'error: to few values in section "tracer**** init"!'
            STOP
270         WRITE(*,'(1A)') &
              'error: to few values in section "transpar"!'
            STOP
280     WRITE(*,'(1A)') &
          'error: can not read number of limits in section "prop limit"!'
        STOP
281     WRITE(*,'(1A)') &
          'error: to few lines in section "prop limit"!'
        STOP
282     WRITE(*,'(1A)') &
          'error: wrong specification of default velocity!'
        STOP
      END

!>    @brief print out the dependency between bc and tpbc
!>    @param[in] ismpl local sample index
      SUBROUTINE show_bcdep(ismpl)
        use arrays
        use mod_genrl
        use mod_genrlc
        use mod_linfos
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k
        character (len=800) :: line
        INTEGER i_bcu, i_bctp

!     print bc and bctp dependency
        i_bcu = -1
        i_bctp = -1
        DO i = 1, nbc_data
!          max bc-unit index
          i_bcu = max(i_bcu,ibc_data(i,cbc_bcu))
!          max bctp index
          i_bctp = max(i_bctp,ibc_data(i,cbc_bctp))
        END DO
        IF (i_bcu>=0 .AND. i_bctp>=1) THEN
          WRITE(*,*) ' '
          WRITE(*,'(6X,A)') 'dependency matrix [time &
            &depended boundary condition]:'
          WRITE(line,'(100(I8,1X))') (i,i=1,i_bctp)
          WRITE(*,'(8X,1A14,1A,1A1)') '|BCunit| bctp:', &
            line(1:i_bctp*11), '|'
          DO i = 0, i_bcu
            line = ' '
            DO j = 1, nbc_data
              k = (ibc_data(j,cbc_bctp)-1)*9 + 1
              IF (ibc_data(j,cbc_bcu)==i .AND. &
                ibc_data(j,cbc_bctp)>=1 .AND. &
                ibc_data(j,cbc_pv)==pv_head .AND. &
                ibc_data(j,cbc_bt)==bt_diri) line(k+1:k+2) = 'Hd'
              IF (ibc_data(j,cbc_bcu)==i .AND. &
                ibc_data(j,cbc_bctp)>=1 .AND. &
                ibc_data(j,cbc_pv)==pv_head .AND. &
                ibc_data(j,cbc_bt)==bt_neum) line(k+1:k+2) = 'Hn'
              IF (ibc_data(j,cbc_bcu)==i .AND. &
                ibc_data(j,cbc_bctp)>=1 .AND. &
                ibc_data(j,cbc_pv)==pv_head .AND. &
                ibc_data(j,cbc_bt)==bt_neuw) line(k+1:k+2) = 'Hw'
              IF (ibc_data(j,cbc_bcu)==i .AND. &
                ibc_data(j,cbc_bctp)>=1 .AND. &
                ibc_data(j,cbc_pv)==pv_pres .AND. &
                ibc_data(j,cbc_bt)==bt_diri) line(k+1:k+2) = 'Pd'
              IF (ibc_data(j,cbc_bcu)==i .AND. &
                ibc_data(j,cbc_bctp)>=1 .AND. &
                ibc_data(j,cbc_pv)==pv_pres .AND. &
                ibc_data(j,cbc_bt)==bt_neum) line(k+1:k+2) = 'Pn'
              IF (ibc_data(j,cbc_bcu)==i .AND. &
                ibc_data(j,cbc_bctp)>=1 .AND. &
                ibc_data(j,cbc_pv)==pv_pres .AND. &
                ibc_data(j,cbc_bt)==bt_neuw) line(k+1:k+2) = 'Pw'
              IF (ibc_data(j,cbc_bcu)==i .AND. &
                ibc_data(j,cbc_bctp)>=1 .AND. &
                ibc_data(j,cbc_pv)==pv_temp .AND. &
                ibc_data(j,cbc_bt)==bt_diri) line(k+3:k+4) = 'Td'
              IF (ibc_data(j,cbc_bcu)==i .AND. &
                ibc_data(j,cbc_bctp)>=1 .AND. &
                ibc_data(j,cbc_pv)==pv_temp .AND. &
                ibc_data(j,cbc_bt)==bt_neum) line(k+3:k+4) = 'Tn'
              IF (ibc_data(j,cbc_bcu)==i .AND. &
                ibc_data(j,cbc_bctp)>=1 .AND. &
                ibc_data(j,cbc_pv)==pv_conc .AND. &
                ibc_data(j,cbc_bt)==bt_diri) line(k+5:k+6) = 'Cd'
              IF (ibc_data(j,cbc_bcu)==i .AND. &
                ibc_data(j,cbc_bctp)>=1 .AND. &
                ibc_data(j,cbc_pv)==pv_conc .AND. &
                ibc_data(j,cbc_bt)==bt_neum) line(k+5:k+6) = 'Cn'
            END DO
            WRITE(*,'(8X,1A2,1I4,1A2,6X,1A,1A1)') '| ', i, ' |', &
              line(1:i_bctp*11), '|'
          END DO
        END IF
!
        RETURN
      END

!>    @brief read all solver parameter
!>    @param[in] line current character line
!>    @param[out] errft break value (error)
!>    @param[out] controlft solver code
!>    @param[out] lmaxitft max iteration number for lin. system solver
!>    @param[in] ismpl local sample index
      SUBROUTINE read_solvpar(line,errft,controlft,lmaxitft,ismpl)
        IMPLICIT NONE
        character (len=80) :: line
        DOUBLE PRECISION errft
        INTEGER controlft, lmaxitft, ismpl
        INTEGER i, j, c1, c2, c3, clast, beginlast
        EXTERNAL clast, beginlast

        CALL read_solver(line,c1,ismpl)
        CALL read_criteria(line,c2,ismpl)
!     default: auto
        IF (c2==-1) c2 = 4
        CALL read_preco(line,c3,ismpl)
!     default: ILU
        IF (c3==-1) c3 = 0
!
        IF (c1==-1) THEN
          READ(line,*) errft, controlft, lmaxitft
        ELSE
          CALL encntrl3(controlft,c1,c2,c3)
          READ(line,*) errft
          j = clast(line)
          i = beginlast(line(1:j))
          READ(line(i:j),*) lmaxitft
        END IF
!
        RETURN
      END

!>    @brief read the solver type
!>    @param[in] line current character line
!>    @param[out] c1 linear system solver code
!>    @param[in] ismpl local sample index
      SUBROUTINE read_solver(line,c1,ismpl)
        IMPLICIT NONE
        character (len=80) :: line
        INTEGER i, c1, ismpl, locstr
        EXTERNAL locstr

        c1 = -1
!     check CG before BiCGStab, because CG is a substring
        i = locstr(line,'cg')
        IF (i>=1) c1 = 2
        i = locstr(line,'bicg')
        IF (i>=1) c1 = 0
        i = locstr(line,'nag')
        IF (i>=1) c1 = 1
        i = locstr(line,'plu')
        IF (i>=1) c1 = 3
!
        RETURN
      END

!>    @brief read the solver break criteria
!>    @param[in] line current character line
!>    @param[out] c2 break criteria code
!>    @param[in] ismpl local sample index
      SUBROUTINE read_criteria(line,c2,ismpl)
        IMPLICIT NONE
        character (len=80) :: line
        INTEGER i, c2, ismpl, locstr
        EXTERNAL locstr

        c2 = -1
        i = locstr(line,'rel')
        IF (i>=1) c2 = 0
        i = locstr(line,'abs')
        IF ((i>=1) .AND. (c2/=0)) c2 = 1
!     Rel + Abs
        IF ((i>=1) .AND. (c2==0)) c2 = 3
        i = locstr(line,'max')
        IF (i>=1) c2 = 2
        i = locstr(line,'auto')
        IF (i>=1) c2 = 4
!
        RETURN
      END

!>    @brief read the preconditioner type
!>    @param[in] line current character line
!>    @param[out] c3 preconditioner code
!>    @param[in] ismpl local sample index
      SUBROUTINE read_preco(line,c3,ismpl)
        IMPLICIT NONE
        character (len=80) :: line
        INTEGER i, c3, ismpl, locstr
        EXTERNAL locstr

        c3 = -1
        i = locstr(line,'ilu')
        IF (i>=1) c3 = 0
        i = locstr(line,'ssor')
        IF (i>=1) c3 = 1
        i = locstr(line,'diag')
        IF (i>=1) c3 = 2
        i = locstr(line,'none')
        IF (i>=1) c3 = 3
!
        RETURN
      END

!>    @brief read the direction-type
!>    @param[in] richtung direction string
!>    @return direction code
      INTEGER FUNCTION read_direction(richtung)
        IMPLICIT NONE
        character (len=*) :: richtung
        INTEGER locstr
        EXTERNAL locstr

        read_direction = 0
        IF (locstr(richtung,'none')==1) read_direction = 0
        IF (locstr(richtung,'0')==1) read_direction = 0
        IF (locstr(richtung,'left')==1) read_direction = 1
        IF (locstr(richtung,'1')==1) read_direction = 1
        IF (locstr(richtung,'right')==1) read_direction = 2
        IF (locstr(richtung,'2')==1) read_direction = 2
        IF (locstr(richtung,'front')==1) read_direction = 3
        IF (locstr(richtung,'3')==1) read_direction = 3
        IF (locstr(richtung,'back')==1) read_direction = 4
        IF (locstr(richtung,'4')==1) read_direction = 4
        IF (locstr(richtung,'base')==1) read_direction = 5
        IF (locstr(richtung,'5')==1) read_direction = 5
        IF (locstr(richtung,'top')==1) read_direction = 6
        IF (locstr(richtung,'6')==1) read_direction = 6
!
        RETURN
      END

#ifdef DEBUG
!>    @brief debug: read additional output times
!>    @param[in] ismpl local sample index
      SUBROUTINE read_debugout(ismpl)
        use arrays
        use mod_genrl
        IMPLICIT NONE
        character (len=80) :: line
        INTEGER i, j, ismpl
        LOGICAL found
        EXTERNAL found

        IF (found(79,key_char//' debug output times',line,.FALSE.)) THEN
          CALL get_arg('records',line,i,j)
          IF (i<1 .OR. j<i) THEN
            READ(79,*) n_debugout
          ELSE
            READ(line(i:j),*) n_debugout
          END IF
          ALLOCATE(debugout(2,n_debugout))
          DO i = 1, n_debugout
            READ(79,*) (debugout(j,i),j=1,2)
          END DO
          WRITE(*,'(1A,1I4)') '  [R] : debug output, records=', &
            n_debugout
        ELSE
          WRITE(*,'(1A,1I4)') &
            '  [I] : no additional debug output'
!        dummy allocation
          ALLOCATE(debugout(1,1))
        END IF
!
        RETURN
      END
#endif

!>    @brief read the property output switch
!>    @param[in] line current character line
      SUBROUTINE read_oprop(line)
        use arrays
        IMPLICIT NONE
        character (len=*) :: line
        INTEGER i, l, locstr, ibegin, iend
        EXTERNAL locstr

        DO i = firstidx, bc_lastidx
          CALL chln(properties(i),ibegin,iend)
          l = locstr(line,properties(i)(ibegin:iend))
          IF (l>=1) out_prop(i) = .FALSE.
        END DO
!
        WRITE(*,*) ' [I] : HDF5 output configuration'
        WRITE(*,'(1A,'//c_npropunit//'A5)') '     ', (properties(i),i=firstidx,lastidx)
        WRITE(*,'(1A,'//c_npropunit//'L5)') '     ', (out_prop(i),i=firstidx,lastidx)
!
        WRITE(*,'(1A,'//c_nbcunit//'A5)') '     ', (properties(i),i=bc_firstidx,bc_lastidx)
        WRITE(*,'(1A,'//c_nbcunit//'L5)') '     ', (out_prop(i),i=bc_firstidx,bc_lastidx)
!
        RETURN
      END

!>    @brief read the physical-value output switch
!>    @param[in] line current character line
      SUBROUTINE read_opv(line)
        use arrays
        IMPLICIT NONE
        character (len=*) :: line
        INTEGER i, l, locstr, ibegin, iend
        EXTERNAL locstr

        DO i = 1, npv
          CALL chln(pv_name(i),ibegin,iend)
          l = locstr(line,pv_name(i)(ibegin:iend))
          IF (l>=1) out_pv(i) = .FALSE.
        END DO
!
        WRITE(*,'(1A,'//c_npv//'A5)') '     ', (pv_name(i),i=1,npv)
        WRITE(*,'(1A,'//c_npv//'L5)') '     ', (out_pv(i),i=1,npv)
!
        RETURN
      END

!>    @brief read the xyz-position output switch
!>    @param[in] line current character line
      SUBROUTINE read_oijk(line)
        use arrays
        IMPLICIT NONE
        character (len=*) :: line
        INTEGER i, l, locstr, ibegin, iend
        EXTERNAL locstr
        character (len=6) :: ijk_name(nout_ijk)
        DATA ijk_name/'    px', '    py', '    pz', '    vx', &
          '    vy', '    vz', '  rhof', '  visf', 'uindex'/

        DO i = 1, nout_ijk
          CALL chln(ijk_name(i),ibegin,iend)
          l = locstr(line,ijk_name(i)(ibegin:iend))
          IF (l>=1) out_ijk(i) = .FALSE.
        END DO
!
        WRITE(*,'(1A,9A7)') '     ', (ijk_name(i),i=1,nout_ijk)
        WRITE(*,'(1A,9L7)') '     ', (out_ijk(i),i=1,nout_ijk)
!
        RETURN
      END

!>    @brief read property limits
!>    @param[in] line current character line
      SUBROUTINE read_limit(line)
        USE ARRAYS
        IMPLICIT NONE
        character (len=*) :: line
        character (len=8) :: skey
        INTEGER i, j, locstr, lblank
        EXTERNAL locstr, lblank
        INTRINSIC adjustl, trim

!       search the parameter key-name
        j = 0
        DO i = 1, nprop_load
          skey = adjustl(properties(i))
          IF (locstr(line,skey(1:lblank(skey)))>=1) j = i
        END DO
        IF (j>=1) THEN
!         MIN or MAX limitation
          IF (locstr(line,'min')>=1) THEN
            READ(line,*) prop_min(j)
          ELSE IF (locstr(line,'max')>=1) THEN
            READ(line,*) prop_max(j)
          ELSE
            WRITE(*,'(4A)') 'error: wrong limit specification (min|max)', &
              ' in section "prop limit" line: "',trim(line),'"!'
            STOP
          ENDIF
        ELSE
          WRITE(*,'(4A)') 'error: wrong parameter keyword in section', &
            ' "prop limit" line: "',trim(line),'"!'
          STOP
        ENDIF
!
        RETURN
      END
