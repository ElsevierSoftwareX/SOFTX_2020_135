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

!>    @brief read data & parameter for restart
!>    @param[in] filename file name with restart information
!>    @param[out] r_time simulation time
!>    @param[out] r_itim time step iteration number
!>    @param[in] ismpl local sample index
      SUBROUTINE read_restartfw(filename,r_time,r_itim,ismpl)
        use arrays
        use mod_genrl
        use mod_flow
        use mod_temp
        use mod_conc
        use mod_time
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k
!
        character (len=80) :: filename
        character (len=80) :: line
        character (len=32) :: strng
        DOUBLE PRECISION r_time
        INTEGER lblank, r_itim, lout, tracer, i1, i2
        LOGICAL found, no_ext_link
        EXTERNAL found, no_ext_link, lblank

!
        WRITE(*,*)
        WRITE(*,*) &
          '  reading forward model input parameter for restart:'
        WRITE(*,*) '    from file "', filename(:lblank(filename)), &
          '"'
        WRITE(*,*)

        CALL omp_new_file_handler(lout,1)
!     read restart file
        OPEN(lout,file=filename,status='old')

        IF (found(lout,key_char//' simtime',line,.TRUE.)) THEN
          READ(lout,*) r_time
          WRITE(*,*) ' [R] : simtime - for restart'
          simtime_0 = r_time
        END IF
        IF (found(lout,key_char//' itimestep',line,.TRUE.)) THEN
          READ(lout,*) r_itim
          WRITE(*,*) ' [R] : itimestep - for restart'
          itimestep_0 = r_itim
        END IF

!     init HDF5 support, when available
        CALL open_hdf5(' ')

        IF (found(lout,key_char//' head init',line,.TRUE.)) THEN
          IF (no_ext_link(i0,j0,k0,head(1,1,1, &
            ismpl),'head',line)) READ(lout,*) (((head(i,j,k,ismpl), &
            i=1,i0),j=1,j0),k=1,k0)
          WRITE(*,*) ' [R] : head - for restart'
        END IF

        IF (found(lout,key_char//' temp init',line,.TRUE.)) THEN
          IF (no_ext_link(i0,j0,k0,temp(1,1,1, &
            ismpl),'temp',line)) READ(lout,*) (((temp(i,j,k,ismpl), &
            i=1,i0),j=1,j0),k=1,k0)
          WRITE(*,*) ' [R] : temp - for restart'
        END IF

        IF (found(lout,key_char//' pres init',line,.TRUE.)) THEN
          IF (no_ext_link(i0,j0,k0,pres(1,1,1, &
            ismpl),'pres',line)) READ(lout,*) (((pres(i,j,k,ismpl), &
            i=1,i0),j=1,j0),k=1,k0)
          WRITE(*,*) ' [R] : pres - for restart'
        END IF
!       convert [MPa] into [Pa]
        CALL dscal(i0*j0*k0,pa_conv,pres(1,1,1,ismpl),1)

        IF (trac_active) THEN
          DO tracer = 1, ntrac
            WRITE(strng,'(1A,1I4.4,1A)') key_char//' tracer', tracer, ' init'
            CALL chln(strng,i1,i2)
            IF (found(79,strng(i1:i2),line,.TRUE.)) THEN
              IF (no_ext_link(i0,j0,k0,conc(1,1,1,tracer, &
                ismpl),strng(i1+2:i2-5),line)) READ(79,*) (((conc(i,j &
                ,k,tracer,ismpl),i=1,i0),j=1,j0),k=1,k0)
              WRITE(*,'(1A,1I4.4,1A)') '  [R] : tracer', tracer, &
                ' - for restart'
            END IF
          END DO
        END IF

        is_init_flow_trafo_needed = .false.

!     finish HDF5 support, when available
        CALL close_hdf5()

!     close restart file
        CLOSE(lout)
        CALL omp_del_file_handler(lout)

        RETURN
      END

!>    @brief write data & parameter for restart
!>    @param[in] filename file name with restart information
!>    @param[in] r_time simulation time
!>    @param[in] r_itim time step iteration number
!>    @param[in] ismpl local sample index
      SUBROUTINE write_restartfw(filename,r_time,r_itim,ismpl)
        use arrays
        use mod_genrl
        use mod_flow
        use mod_temp
        use mod_conc
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k
!
        character (len=80) :: filename
        character (len=80) :: hfilename
        character (len=32) :: strng
        DOUBLE PRECISION r_time
        INTEGER lblank, r_itim, lout, tracer, i1, i2
        LOGICAL test_hdf5
        EXTERNAL lblank, test_hdf5

!
!$OMP master
        CALL omp_new_file_handler(lout,1)
!     write restart file
        OPEN(lout,file=filename,status='unknown')

        WRITE(lout,'(A)') key_char//' simtime'
        WRITE(lout,'(e24.16)') r_time

        WRITE(lout,'(A)') key_char//' itimestep'
        WRITE(lout,'(I8)') r_itim

!AW-n.parallel      write(lout,'(A)') '# rms'
!AW-n.parallel      write(lout,'(1e24.16)') difrmsf, difrmst, difrmse

        WRITE(*,'(3A)') '  [W] : control parameter to "', &
          filename(:lblank(filename)), '" - for restart'

        IF (test_hdf5()) THEN

!       init HDF5 support, when available
          CALL open_hdf5(' ')

          hfilename = filename(:lblank(filename)) // '.h5'
          CALL create_hdf5(hfilename)

          WRITE(lout,'(A)') key_char//' head init, HDF5 = ' // hfilename
          CALL write_hdf5(i0,j0,k0,head(1,1,1,ismpl),'head', &
            hfilename)
!       write(*,*) ' [W] : head - for restart'

          WRITE(lout,'(A)') key_char//' temp init, HDF5 = ' // hfilename
          CALL write_hdf5(i0,j0,k0,temp(1,1,1,ismpl),'temp', &
            hfilename)
!       write(*,*) ' [W] : temp - for restart'

          WRITE(lout,'(A)') key_char//' pres init, HDF5 = ' // hfilename
!       convert [Pa] into [MPa]
          DO k = 1, k0
            DO j = 1, j0
              DO i = 1, i0
                x(i,j,k,ismpl) = pres(i,j,k,ismpl)*pa_conv1
              END DO
            END DO
          END DO
          CALL write_hdf5(i0,j0,k0,x(1,1,1,ismpl),'pres',hfilename)
!       write(*,*) ' [W] : pres - for restart'

          IF (trac_active) THEN
            DO tracer = 1, ntrac
              WRITE(strng,'(1A,1I4.4)') 'tracer', tracer
              CALL chln(strng,i1,i2)
              CALL write_hdf5(i0,j0,k0,conc(1,1,1,tracer,ismpl), &
                strng(i1:i2),hfilename)
!            write(*,'(1A,1I4.4,1A)') '  [W] : tracer',tracer,' - for restart'
            END DO
          END IF

          WRITE(*,'(3A)') '  [W] : physical state to "', &
            filename(:lblank(filename)), '.h5" - for restart'

!       finish HDF5 support, when available
          CALL close_hdf5()

        ELSE

          WRITE(lout,'(A)') key_char//' head init'
          CALL write_dense3d(head(1,1,1,ismpl),i0,j0,k0,lout,ismpl)
!       write(*,*) ' [W] : head - for restart'

          WRITE(lout,'(A)') key_char//' temp init'
          CALL write_dense3d(temp(1,1,1,ismpl),i0,j0,k0,lout,ismpl)
!       write(*,*) ' [W] : temp - for restart'

          WRITE(lout,'(A)') key_char//' pres init'
!       convert [Pa] into [MPa]
          DO k = 1, k0
            DO j = 1, j0
              DO i = 1, i0
                x(i,j,k,ismpl) = pres(i,j,k,ismpl)*pa_conv1
              END DO
            END DO
          END DO
          CALL write_dense3d(x(1,1,1,ismpl),i0,j0,k0,lout,ismpl)
!       write(*,*) ' [W] : pres - for restart'

          IF (trac_active) THEN
            DO tracer = 1, ntrac
              WRITE(strng,'(1A,1I4.4)') 'tracer', tracer
              CALL chln(strng,i1,i2)
              WRITE(lout,'(3A)') key_char//' ', strng(i1:i2), ' init'
              CALL write_dense3d(conc(1,1,1,tracer,ismpl),i0,j0,k0, &
                lout,ismpl)
!           write(*,'(1A,1I4.4,1A)') '  [W] : tracer',tracer,' - for restart'
            END DO
          END IF

          WRITE(*,'(3A)') '  [W] : physical state to "', &
            filename(:lblank(filename)), '" - for restart'

        END IF

!     close restart file
        CLOSE(lout)
        CALL omp_del_file_handler(lout)
!$OMP end master

        RETURN
      END
