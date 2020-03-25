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

!>    @brief read time parameter
!>    @param[in] filename model file name
!>    @param[in] ismpl local sample index
!>    @details
!> Note: To be able to use input file parsing with hdf5, the
!> hdf5-input-files have to be generated using the script:
!> `convert_to_hdf5.py`. This script can be found in the repository
!> `SHEMAT-Suite_Scripts` under
!> `python/preprocessing/convert_to_hdf5.py`.
      SUBROUTINE read_time(filename,ismpl)
        use arrays
        use mod_genrl
        use mod_genrlc
        use mod_flow
        use mod_temp
        use mod_conc
        use mod_time
        use mod_data
        use mod_linfos
#ifndef noHDF
        use mod_input_file_parser_hdf5
#endif
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k, l
!
        character (len=80) :: filename
        character (len=80) :: line
!
        INTEGER lblank, id, num, step_type, locstr
        LOGICAL found, no_ext_link, no_ext_link_int
        EXTERNAL found, no_ext_link, no_ext_link_int, lblank, locstr

        logical :: found_marker = .false.


        WRITE(*,*)
        WRITE(*,*) '  reading time input parameter:'
        WRITE(*,*) '    from file "', filename(:lblank(filename)), &
          '"'
        WRITE(*,*)

!     read file
        OPEN(79,file=filename,status='old')

!     init HDF5 support, when available
        CALL open_hdf5(' ')

! ------------------
        transient = .FALSE.
        thetaf = 1.D0
        thetat = 1.D0
!     start of the simulation time
        simtime(ismpl) = 0.D0
#ifndef noHDF
        if (h5parse_use_hdf5_datafile) then
            if (h5parse_check_attr_exist("thetaf","time")) then
                call h5parse_read_double_attribute("thetaf", thetaf, "time")
                call h5parse_read_double_attribute("thetat", thetat, "time")
                call h5parse_read_double_attribute("thetac", thetac, "time")
                call h5parse_read_double_attribute("tstart", simtime(ismpl), "time")
                found_marker = .true.
            end if
        else
#endif
        IF (found(79,key_char//' timestep control',line,.FALSE.)) THEN
          READ(79,*) i
          IF (i/=0) THEN
            READ(79,*,err=1001) thetaf, thetat, thetac, &
              simtime(ismpl)
              found_marker = .true.
          ELSE
            WRITE(*,'(A)') '  [R] : Steady state !'
          END IF
        ELSE
          WRITE(*,'(A)') '  <D> : Steady state !'
        END IF
#ifndef noHDF
        endif
#endif
        if (found_marker) then
            found_marker = .false.
            transient = .TRUE.
            WRITE(*,'(A)') '  [R] : timestep control'
            IF (thetaf<0.5D0) THEN
              WRITE(*,'(A)') '  theta flow < .5, set to .5!'
              thetaf = 0.5D0
            END IF
            IF (thetat<0.5D0) THEN
              WRITE(*,'(A)') '  theta temp < .5, set to .5!'
              thetat = 0.5D0
            END IF
            IF (thetac<0.5D0) THEN
              WRITE(*,'(A)') '  theta conc < .5, set to .5!'
              thetat = 0.5D0
            END IF
        END IF

        ! Default time unit [s]
        tunit = tunit_const
#ifndef noHDF
        if (h5parse_use_hdf5_datafile) then
            if (h5parse_check_attr_exist("tunit","time")) then
                call h5parse_read_double_attribute("tunit", tunit, "time")
                found_marker = .true.
            end if
        else
#endif
        IF (found(79,key_char//' tunit',line,.FALSE.)) THEN
          READ(79,*) tunit
          found_marker = .true.
        END IF
#ifndef noHDF
        end if
#endif
        if (found_marker) then
            found_marker = .false.
            WRITE(*,'(A,1e24.8)') '  [R] : tunit =', tunit
        else
            WRITE(*,'(A,1e24.8,A)') '  <D> : tunit =', tunit, ' !'
        end if

!        simulation start time (overwrite the reading in section "# time")
#ifndef noHDF
        if (h5parse_use_hdf5_datafile) then
            if (h5parse_check_attr_exist("tstart","time")) then
                found_marker = .true.
                ! tstart was already read before
            end if
        else
#endif
        IF (found(79,key_char//' tstart',line,.FALSE.)) THEN
          READ(79,*) simtime(ismpl)
          found_marker = .true.
        END IF
#ifndef noHDF
        end if
#endif
        if (found_marker) then
            found_marker = .false.
            WRITE(*,'(A,1e24.8)') '  [R] : tstart =', simtime(ismpl)
        else
            WRITE(*,'(A,1e24.8)') '  <D> : tstart=', simtime(ismpl)
        end if

!        offset of the simulation time index
        itimestep_0 = 0
#ifndef noHDF
        if (h5parse_use_hdf5_datafile) then
            if (h5parse_check_attr_exist("titer","time")) then
                call h5parse_read_integer_attribute("titer", itimestep_0, "time")
                found_marker = .true.
            end if
        else
#endif
        IF (found(79,key_char//' titer',line,.FALSE.)) THEN
          READ(79,*) itimestep_0
          found_marker = .true.
        END IF
#ifndef noHDF
        end if
#endif
        if (found_marker) then
            found_marker = .false.
            WRITE(*,'(A,1I7)') '  [R] : titer =', itimestep_0
        else
            WRITE(*,'(A,1I7)') '  <D> : titer =', itimestep_0
        end if

!        convert into sec.
        simtime(ismpl) = simtime(ismpl)*tunit
        simtime_0 = simtime(ismpl)

!       check for variable time step size
#ifndef noHDF
        if (h5parse_use_hdf5_datafile) then
            if (h5parse_check_attr_exist("delt_start","time")) then
                found_marker = .true.
                call h5parse_read_double_attribute("delt_start", delt_start, "time")
                call h5parse_read_double_attribute("delt_min", delt_min, "time")
                call h5parse_read_double_attribute("delt_max", delt_max, "time")
                call h5parse_read_double_attribute("max_simtime", max_simtime, "time")
                call h5parse_read_integer_attribute("delt_double", delt_double, "time")
            end if
        else
#endif
        IF (found(79,key_char//' variable step size',line,.FALSE.)) THEN
          READ(79,*) i
          IF (i/=0) THEN
            READ(79,*,err=1002) delt_start, delt_min, delt_max, max_simtime, delt_double
            found_marker = .true.
          END IF
        END IF
#ifndef noHDF
        end if
#endif
        ! Set flag for variable step size
        delt_vary = found_marker

        if (found_marker) then
            found_marker = .false.
            WRITE(*,*)
            WRITE(*,'(1A)') '  [R] : Variable time stepping: '
            WRITE(*,'(1A,1e12.4)') '    Starting step size: ', delt_start
            WRITE(*,'(1A,1e12.4)') '    Minimum  step size: ', delt_min
            WRITE(*,'(1A,1e12.4)') '    Maximum  step size: ', delt_max
            WRITE(*,'(1A,1e12.4)') '    Maximum  simulation time: ', max_simtime
            WRITE(*,'(1A,1I7)') '    Doubling time: ', delt_double
            WRITE(*,*)
        end if

        ! Default number of time periods
        nperiod = 0

        ! Default total number of time steps
        ntimestep = 0

        ! Define time periods from input
        IF (.not. delt_vary) THEN

          max_simtime = simtime(ismpl)
#ifndef noHDF
          if (h5parse_use_hdf5_datafile) then
            if (h5parse_check_dataset_exist("time/periods")) then
                found_marker = .true.
                nperiod = h5parse_read_dimension_size_for_dataset("time/periods")
                allocate(iperiod(nperiod,2))
                call h5parse_read_2d_double_dataset("time/periods", dperiod(:nperiod,:))
                call h5parse_read_2d_integer_dataset("time/iperiods", iperiod(:nperiod,:2))
            end if
          else
#endif
          IF (found(79,key_char//' time periods',line,.FALSE.)) THEN
            found_marker = .true.

            ! Read number of time periods
            CALL get_arg('records',line,i,j)
            IF (i<1 .OR. j<i) THEN
              READ(79,*) nperiod
            ELSE
              READ(line(i:j),*) nperiod
            END IF

            ! Test that number of periods is a positive integer
            if (nperiod < 1) then
              write(unit = *, fmt = *) "[E]: Number of time periods in ", &
                  "'# time periods' must be positive integer.\n", &
                  "Input was nperiod=", nperiod
              stop
            end if

            allocate(iperiod(nperiod,2))
            allocate(dperiod(nperiod,2))

            ! Read time period information
            DO k = 1, nperiod
! ---------------
              READ(79,'(1A)') line
!               step_type : 0=not set, 1=linear, 2=logarithmic
              step_type = 0
              i = locstr(line,'lin')
              IF (i>=1) step_type = 1
              i = locstr(line,'log')
              IF (i>=1) step_type = 2
!              JK: Logarithm
              i = locstr(line,'jlo')
              IF (i>=1) step_type = 3
!               increase
              i = locstr(line,'inc')
              IF (i>=1) step_type = abs(step_type)
!               decrease
              i = locstr(line,'dec')
              IF (i>=1) step_type = -abs(step_type)

              ! For time period k:
              ! - Read 'dperiod': [ start of time period, end of time period ]
              ! - Read 'iperiod': [ # time steps, time step type]
              IF (step_type==0) THEN
!                  need a third value as the time step type
                READ(line,*,err=2001) dperiod(k,1), dperiod(k,2), &
                  (iperiod(k,j),j=1,2)
              ELSE
                READ(line,*,err=2002) dperiod(k,1), dperiod(k,2), &
                  (iperiod(k,j),j=1,1)
                iperiod(k,2) = step_type
              END IF
            end do
          end if
#ifndef noHDF
          end if
#endif

          if (found_marker) then
            WRITE(*,'(A,I4)') '  [R] : periods, records=', nperiod
            found_marker = .false.
            DO k = 1, nperiod
              ntimestep = ntimestep + iperiod(k,1)
              dperiod(k,1) = dperiod(k,1)*tunit
              dperiod(k,2) = dperiod(k,2)*tunit
              max_simtime = max(max_simtime,dperiod(k,2))

              IF (iperiod(k,2)==1) THEN
                WRITE(*,'(1A,1e12.4,1A,1e12.4,1A,1I8,1A)') &
                  '    time :', dperiod(k,1)/tunit, ' - ', &
                  dperiod(k,2)/tunit, ', '//key_char//' ', iperiod(k,1), ', linear'
              ELSE IF (iperiod(k,2)==2) THEN
                WRITE(*,'(1A,1e12.4,1A,1e12.4,1A,1I8,1A)') &
                  '    time :', dperiod(k,1)/tunit, ' - ', &
                  dperiod(k,2)/tunit, ', '//key_char//' ', iperiod(k,1), &
                  ', logarithmic ascending'
              ELSE IF (iperiod(k,2)==3) THEN
                WRITE(*,'(1A,1e12.4,1A,1e12.4,1A,1I8,1A)') &
                  '    time :', dperiod(k,1)/tunit, ' - ', &
                  dperiod(k,2)/tunit, ', '//key_char//' ', iperiod(k,1), &
                  ', varied logarithmic ascending'
              ELSE IF (iperiod(k,2)==-2) THEN
                WRITE(*,'(1A,1e12.4,1A,1e12.4,1A,1I8,1A)') &
                  '    time :', dperiod(k,1)/tunit, ' - ', &
                  dperiod(k,2)/tunit, ', '//key_char//' ', iperiod(k,1), &
                  ', logarithmic descending'
              ELSE
                WRITE(*,'(1A,1I2,1A,1I2,1A)') &
                  'error: wrong time step type ', iperiod(k,2), ' at ', &
                  k, ' !'
                STOP
              END IF

              ! Sanity checks
              ! -------------

              ! start time < end time
              if (dperiod(k,2) <= dperiod(k,1)) then
                write(*,'(1A,1I3,1A)') &
                  'error: time period START smaller then the END at ', &
                  k, ' !'
                write(*,'(1A,1e12.4,1A,1e12.4,1A)') 't_start = ', &
                    dperiod(k,1)/tunit, ' > t_end = ', dperiod(k,2)/tunit, &
                    ' in "read_time.f" !'
                stop
              end if
              ! no negative start / end times
              if (dperiod(k,1) < 0.0d0 .or. dperiod(k,2) < 0.0d0) then
                write(*,'(1A,1I4,1A)') &
                    'error: no negative times allowed (period', k, &
                    ') in "read_time.f" !'
                write(*,'(1A,1e12.4,1A,1e12.4,1A)') 't_start = ', &
                    dperiod(k,1)/tunit, ' < 0, or  t_end', dperiod(k,2)/tunit, &
                    ' < 0 in "read_time.f" !'
                stop
              end if
              ! end time = next start time
              if (k>=2) then
                if (1.0d0-dperiod(k-1,2)/dperiod(k,1) > 1.0d-14) then
                  write(*,'(1A,1I3,1A,1I3,1A)') &
                    'error: time gap between period ', k - 1, ' and ', &
                    k, ' !'
                  stop
                end if
                if (1.0d0-dperiod(k,1)/dperiod(k-1,2) > 1.0d-14) then
                  write(*,'(1A,1I3,1A,1I3,1A)') &
                    'error: time overlapping between period ', k - 1, &
                    ' and ', k, ' !'
                  stop
                end if
                ! first start time <= beginning of simulation
              else if (dperiod(1,1)>simtime(ismpl)) then
                write(*,'(2A)') 'error: first time period begins', &
                  ' later than the simulation start time !'
                stop
              end if
              ! -------------

            END DO

            ! Check the number of time periods
            if (ntimestep < 1) then
              write(unit = *, fmt = *) "Error: number of time steps (ntimestep = ", &
                  ntimestep, ") should be larger than zero."
              stop
            end if

            ! Compute time period table
            CALL calc_deltatime(ismpl)

            deallocate(iperiod)
            deallocate(dperiod)

          END IF

        END IF

!     read bc time periods
        bctp = .FALSE.
        IF (found(79,key_char//' bc time periods',line,.FALSE.)) THEN
          CALL get_arg('records',line,i,j)
          IF (i<1 .OR. j<i) THEN
            READ(79,*) l
          ELSE
            READ(line(i:j),*) l
          END IF
          bctp = .TRUE.
!         clean counter table
          DO i = 1, nbctp
            ibcperiod(i) = 0
          END DO
!           "bcperiod": (period-index,value-type,TP-ID,sample)
!           - value-type: time, Alpha, Beta
!           "ibcperiod" - number of periods: (TP-ID)
!           "lbcperiod" - on/off switch: (period-index,TP-ID)
          DO i = 1, l
            READ(79,*,err=3000,end=3000) id, num
!           sanity check
            IF (num>ngsmax) THEN
              WRITE(*,'(1A,1I3,1A,1I5,1A)') 'error: number of periods in "bc time periods" to big (max. ', &
                ngsmax, ', ID=', id, ') !'
              STOP
            END IF
            IF (id>nbctp) THEN
              WRITE(*,'(1A)') 'error: something goes wrong with BC-TP indexing (software bug)!'
              STOP
            END IF
!
            ibcperiod(id) = num
            DO k = 1, num
              READ(79,'(1A)',err=3001,end=3001) line
              READ(line,*,err=3002,end=3002) (bcperiod(k,j,id,ismpl), &
                j=1,3)
              lbcperiod(k,id) = .TRUE.
              j = locstr(line,'off')
              IF (j>=1) lbcperiod(k,id) = .FALSE.
            END DO
            DO k = 1, num
              bcperiod(k,1,id,ismpl) = tunit*bcperiod(k,1,id,ismpl)
              IF (k>=2) THEN
                IF (bcperiod(k,1,id,ismpl)<=bcperiod(k-1,1,id,ismpl)) &
                    THEN
                  WRITE(*,'(1A,1I2,1A,1I4,2A)') 'error: BC time period record=', i, &
                    ' and time changing=', k, ' has a smaller start time than before', &
                    ' (sorting needed)!'
                  STOP
                END IF
              END IF
            END DO
          END DO
          WRITE(*,'(1A,1I4,1A,1I4,1A)') '  [R] : BC time periods, records=', l,' (max index =',nbctp,')'
!         sanity check, full check for time period table - usage
          DO i = 1, nbc_data
            j = ibc_data(i,cbc_bctp)
            IF (j>nbctp) THEN
              WRITE(*,'(5A,1I4,1A)') 'error: "',pv_name(ibc_data(i,cbc_pv)), ' ', &
                bc_name(ibc_data(i,cbc_bt)),'" time period mismatch, use of a BCTP table index ', &
                j,' in the declaration, but not defined under "time periods"!'
              STOP
            END IF
            IF (j>0) THEN
              IF (ibcperiod(j)<1) THEN
                WRITE(*,'(5A,1I4,1A)') 'error: "',pv_name(ibc_data(i,cbc_pv)), ' ', &
                  bc_name(ibc_data(i,cbc_bt)),'" time period mismatch, use of a BCTP table index ', &
                  j,' in the declaration, but not defined under "time periods" or zero!'
                STOP
              END IF
            END IF
          END DO
        END IF

!     monitoring time steps
        nmon = 0
        monitor = .FALSE.
        out_orientation = 0
        IF (found(79,key_char//' monitor',line,.FALSE.)) THEN
          monitor = .TRUE.
          CALL get_arg('records',line,i,j)
          IF (i<1 .OR. j<i) THEN
            WRITE(*,'(1A)') &
              'error: in section "monitor", records=<number> needed!'
            STOP
          ELSE
            READ(line(i:j),*) nmon
          END IF
          i = locstr(line,'col')
          IF (i>=1) THEN
!           column wise orientation
            out_orientation = 0
          END IF
          i = locstr(line,'row')
          IF (i>=1) THEN
!           row wise orientation
            out_orientation = 1
          END IF
          i = locstr(line,'new')
          IF (i>=1) THEN
            CALL get_arg('new',line,i,j)
            IF (i>=1 .AND. j>=i) THEN
!             sainty check
              IF (line(i:i)=='p' .AND. out_orientation==1) THEN
                WRITE(*,'(1A)') 'error: "new=position" defined, &
                  &but "row" not supported for it!'
                STOP
              END IF
!             "time" defined, time named files
              IF (line(i:i)=='t') out_orientation = out_orientation + &
                2
!             "position" defined, position named files
              IF (line(i:i)=='p') out_orientation = out_orientation + &
                4
            ELSE
!             default: time files
              out_orientation = out_orientation + 2
            END IF
          END IF
          WRITE(*,'(A,I4,1A,1I2)') &
            '  [R] : monitoring points, records=', nmon, &
            ', data orientation=', out_orientation
          READ(79,*,err=4001) ((imon(i,j),j=1,4),i=1,nmon)
!        sanity checks
          DO i = 1, nmon
            IF (imon(i,1)>i0 .OR. imon(i,1)<1 .OR. imon(i,2)>j0 .OR. &
                imon(i,2)<1 .OR. imon(i,3)>k0 .OR. imon(i,3)<1) THEN
              WRITE(*,'(1A,1I5,1A)') 'error: monitor point ', i, &
                ' out of range !'
              STOP
            END IF
          END DO
        ELSE
          WRITE(*,'(1A)') '  <D> : no monitoring points !'
        END IF

!     output times
        noutt = 0
        IF (found(79,key_char//' output times',line,.FALSE.)) THEN
          CALL get_arg('records',line,i,j)
          IF (i<1 .OR. j<i) THEN
            READ(79,*) noutt
          ELSE
            READ(line(i:j),*) noutt
          END IF
          DEALLOCATE(outt)
          memory = memory - 1
          ALLOCATE(outt(noutt+1))
          memory = memory + noutt + 1
          IF (noutt > 0) READ(79,*) (outt(i),i=1,noutt)
          DO i = 1, noutt
            outt(i) = outt(i)*tunit
          END DO
          outt(noutt+1) = max_simtime*2.D0 +1.0d0
          WRITE(*,'(A,I4)') '  [R] : output times, records=', noutt
        ELSE
          outt(1) = max_simtime*2.D0 +1.0d0
          WRITE(*,'(A)') '  <D> : no output times !'
        END IF

!aw??      endif
! ------------------

!     finish HDF5 support, when available
        CALL close_hdf5()

!     close project config file
        CLOSE(79)

        WRITE(*,*)
        WRITE(*,*)

        RETURN

!     error handler
1001    WRITE(*,'(2A)') 'error: in section "timestep control",', &
          ' awaiting [thetaf thetat thetac simtime] !'
        STOP
1002    WRITE(*,'(2A)') 'error: in section "variable step size",', &
          ' awaiting [delt_start, delt_min, delt_max, max_simtime, delt_double] !'
        STOP
2001    WRITE(*,'(1A)') 'error: awaiting four values per line !'
        STOP
2002    WRITE(*,'(1A)') 'error: awaiting three values per line !'
        STOP
3000    WRITE(*,'(2A)') 'error: in section "bc time periods",', &
          ' awaiting [tp-ID tp-number] !'
        STOP
3001    WRITE(*,'(1A,1I3,1A)') &
          'error: awaiting BC time period information at line ', k, &
          '!'
        STOP
3002    WRITE(*,'(1A,1I3,1A)') 'error: awaiting three values for ', &
          num, ' lines !'
        STOP
4001    WRITE(*,'(1A,1I3,1A)') 'error: awaiting four values for ', &
          nmon, ' lines !'
        STOP
      END
