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

!>    @brief main program for forward simulation (only)
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
program shem_fw
  use arrays
  use mod_genrl
  use mod_genrlc
  use mod_time
  use mod_OMP_TOOLS
#ifndef noHDF
  use mod_input_file_parser_hdf5
#endif
  implicit none

  ! Global sampling index
  integer :: ismpl

  ! Starting time of SHEMAT-Suite execution
  double precision :: tsglobal

  ! Finishing time of SHEMAT-Suite execution
  double precision :: tfglobal

  ! Runtime of SHEMAT-Suite execution
  double precision ::ttglobal

  ! Number of full minutes in runtime
  integer :: tmins

  ! Include version information from compilation
  INCLUDE 'version.inc'

  ! SHEMAT-Suite input filename read from shemade.job
  character (len=80) :: filename


  CALL sys_cputime(tsglobal)
  nested_build = .FALSE.
  ismpl = 1
  def_binary = 'forward'

  write(*,*) ' '
  write(*,*) '======================================'
  write(*,*) version
  write(*,*) '======================================'
  write(*,*) ' '

  open(unit = 66, file='shemade.job')

  ! -----------------------------------------

  ! read new input file name
10 read(unit = 66, fmt = '(A)', end = 99999) filename

  if (filename(1:1)=='!') go to 10
  if (filename(1:1)=='%') go to 10
  if (filename(1:1)=='#') go to 10

  write(unit = *, fmt = *) ' '
  write(unit = *, fmt = *) ' '
  write(unit = *, fmt = *) ' *** NEW MODEL '

  project = filename
  runmode = 0
  transient = .false.

  ! reading input
  call forward_preparation(filename, ismpl)

  !     initialize
  call forward_init(ismpl)

  ! ---------
  call forward_iter(simtime_0,max_simtime,-1,0,ismpl)
  ! ---------

  !     output
  call forward_write(-1,ismpl)
  if (runmode > 0) call write_data(-1,ismpl)

  call dealloc_arrays(ismpl)
  call props_end(ismpl)
  call dealloc_data(ismpl)

  !     Another file to load?
  go to 10

  ! -----------------------------------------

  !     finis terrae
99999 continue
  !
  close(unit = 66)
  call sys_cputime(tfglobal)
  ttglobal = tfglobal - ttglobal
  tmins = int(ttglobal/60.D0)
  write(unit = *, fmt = '(1A,1I4,1A,1F5.2,1A)') ' total cpu time: ', tmins, ':', &
      ttglobal - dble(tmins)*60.D0, ' min'
  write(unit = *, fmt = *) 'RUN O.K.'
  !
end program shem_fw
