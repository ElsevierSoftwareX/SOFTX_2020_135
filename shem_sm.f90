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

!>    @brief main program stochstic simulation/inversion (based on SGSIM, VISIM, ENKF)
!>    @details
!>
!> **SHEMAT-Suite (Simulator for HEat and MAss Transport)** is a
!> numerical code for computing flow, heat and species transport
!> equations in porous media. The governing equations of the code are
!> the groundwater flow equation, the heat transport equation and the
!> species transport equation. \n\n
!>
!> SHEMAT-Suite includes parameter estimation and data assimilation
!> approaches, both stochastic (Monte Carlo, ensemble Kalman filter)
!> and deterministic (Bayesian inversion using automatic
!> differentiation for calculating derivatives).\n\n
!>
!> Note: To be able to use input file parsing with hdf5, the
!> hdf5-input-files have to be generated using the script:
!> `convert_to_hdf5.py`. This script can be found in the repository
!> `SHEMAT-Suite_Scripts` under
!> `python/preprocessing/convert_to_hdf5.py`.
      PROGRAM shem_sm
        use arrays
        use mod_genrl
        use mod_genrlc
        use mod_data
        use mod_time
        use mod_linfos
        use mod_OMP_TOOLS
        use mod_simul
#ifndef noHDF
        use mod_input_file_parser_hdf5
#endif
        IMPLICIT NONE
        integer :: ismpl
        integer :: i
        double precision :: tsglobal, tfglobal
        INCLUDE 'version.inc'
        INCLUDE 'OMP_TOOLS.inc'
!
        character (len=80) :: filename
!
        INTEGER lblank, iter_enkf
        EXTERNAL lblank
!     if "restart" option is used
        LOGICAL test_option
        EXTERNAL test_option
!     time of the last saving of the restart data
        DOUBLE PRECISION backuptime

        CALL sys_cputime(tsglobal)
        backuptime = tsglobal
        write_disable = .FALSE.
        write_iter_disable = .FALSE.
        nested_build = .TRUE.
        ismpl = 1
        def_binary = 'simul'
        
        WRITE(*,*) ' '
        WRITE(*,*) '======================================'
        WRITE(*,*) version
        WRITE(*,*) '======================================'
        WRITE(*,*) ' '
        OPEN(66,file='shemade.job')

! -----------------------------------------

!     read new model
10      READ(66,'(A)',end=99999) filename
        CALL sys_cputime(tslocal)


        IF (filename(1:1)=='!') GO TO 10
        IF (filename(1:1)=='%') GO TO 10
        IF (filename(1:1)=='#') GO TO 10

        WRITE(*,*) ' '
        WRITE(*,*) ' '
        WRITE(*,*) ' *** NEW MODEL '
        project = filename
        runmode = 0
        transient = .FALSE.
        nperiod = 0


#ifndef noHDF
        CALL h5parse_open_datafile(filename)
#endif
!     read forward model
        CALL read_model(filename,ismpl)
        CALL read_control(filename,ismpl)
!     read timestep parameter
        CALL read_time(filename,ismpl)
!     read simulation parameters
        CALL read_simul(filename_simul,ismpl)
!     read ENKF parameters
        CALL read_enkf(filename_enkf,ismpl)
!     read data
        IF (runmode>0) CALL read_data(filename_data,ismpl)
!     split units
        CALL read_split_sm(filename,ismpl)

        status_log = filename(1:lblank(filename)) // '_status.log'
        restart_name = filename(1:lblank(filename)) // '_restart'

        IF (test_option('restart')) THEN
!       read all necessary data from last save point
          CALL read_restartfw(restart_name,simtime(ismpl),itimestep_0, &
            ismpl)

          WRITE(*,'(3A)') '  [W] : "', status_log(1:lblank(status_log &
            )), '"'
          OPEN(76,file=status_log,status='unknown',position='append')
          WRITE(76,'(3A)') '%    restart Project: "', &
            filename(1:lblank(filename)), '"'
          WRITE(76,'(1A)') '%'
          IF (transient) THEN
            WRITE(76,'(2A)') '% <time step>', ' <cum time>'
          END IF
#ifdef head_base
          WRITE(76,'(4A)') '%    <iteration>',' <delta head>',' <delta temp>',' (<delta conc> ...)'
#endif
#ifdef pres_base
          WRITE(76,'(4A)') '%    <iteration>',' <delta pres>',' <delta temp>',' (<delta conc> ...)'
#endif
          CLOSE(76)

        ELSE
          WRITE(*,'(3A)') '  [W] : "', status_log(1:lblank(status_log &
            )), '"'
          OPEN(76,file=status_log)
          WRITE(76,'(2A)') '% Shemat-Suite version: ', version
          WRITE(76,'(2A)') '%          build: ', datum
          WRITE(76,'(2A)') '%   build command line: ', makecmd
          WRITE(76,'(3A)') '%        Project: "', &
            filename(1:lblank(filename)), '"'
          WRITE(76,'(1A)') '%'
          IF (transient) THEN
            WRITE(76,'(2A)') '% <time step>', ' <cum time>'
          END IF
#ifdef head_base
          WRITE(76,'(4A)') '%    <iteration>',' <delta head>',' <delta temp>',' (<delta conc> ...)'
#endif
#ifdef pres_base
          WRITE(76,'(4A)') '%    <iteration>',' <delta pres>',' <delta temp>',' (<delta conc> ...)'
#endif
          CLOSE(76)
        END IF

#ifndef noHDF
        CALL h5parse_close_datafile()
#endif

!
!     initialize

        CALL forward_init(ismpl)

! ---------
        WRITE(*,*) ' '
        DO iter_enkf = 1, maxiter_enkf
          WRITE(*,*) ' '
          WRITE(*,'(1A,1I2,1A)') ' ***   ENKF iteration ',iter_enkf,' ***'
          WRITE(*,*) ' '
!         compute realisations
          CALL enkf_iter(iter_enkf)
        END DO
! ---------

!     output - not needed for SIMUL/ENKF, because of the samples output
!AW        CALL forward_write(-1,ismpl)
!AW        IF (runmode>0) CALL write_data(-1,ismpl)

        CALL dealloc_simul(ismpl)
        CALL dealloc_arrays(ismpl)
        CALL props_end(ismpl)
        CALL dealloc_data(ismpl)

!     Another file to load?
        GO TO 10

! -----------------------------------------

!     finis terrae
99999   CONTINUE


        CLOSE(66)
        CALL sys_cputime(tfglobal)
        tslocal = tfglobal - tsglobal
        i = int(tslocal/60.D0)
        WRITE(*,'(1A,1I4,1A,1F5.2,1A)') ' total cpu time: ', i, ':', &
          tslocal - dble(i)*60.D0, ' min'
        WRITE(*,*) 'RUN O.K.'
      END
