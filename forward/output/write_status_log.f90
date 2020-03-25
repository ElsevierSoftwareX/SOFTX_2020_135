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

!>    @brief naming status.log file and writing header
!>    @param[in] filename input filename as read from shemade.job
!>    @param[inout] ismpl local sample index
!>    @details
!>    The status.log file is named and the header is written depending
!>    on whether the restart option was given in the command line.
subroutine write_status_log(filename, ismpl)

  use mod_genrlc
  use mod_time

  implicit none

  integer :: ismpl
  
  character (len=80), intent (in) :: filename
  
  INCLUDE 'version.inc'

  integer, external :: lblank
  logical, external :: test_option

  status_log = filename(1:lblank(filename)) // '_status.log'
  status_log_inv = filename(1:lblank(filename)) // &
      '_status-inv.log'

  IF (test_option('restart')) THEN
    WRITE(*,'(3A)') '  [W] : "', status_log(1:lblank(status_log)), '"'
    OPEN(76,file=status_log,status='unknown',position='append')
    WRITE(76,'(3A)') '%    restart Project: "',filename(1:lblank(filename)), '"'
    WRITE(76,'(1A)') '%'
    IF (transient) THEN
      WRITE(76,'(3A)') '% <time step>', ' <deltat>', ' <cum time>'
    END IF
#ifdef head_base
    WRITE(76,'(4A)') '%    <iteration>',' <delta head>',' <delta temp>',' (<delta conc> ...)'
#endif
#ifdef pres_base
    WRITE(76,'(4A)') '%    <iteration>',' <delta pres>',' <delta temp>',' (<delta conc> ...)'
#endif
    CLOSE(76)

  ELSE
    WRITE(*,'(3A)') '  [W] : "', status_log(1:lblank(status_log)), '"'
    OPEN(76,file=status_log)
    WRITE(76,'(2A)') '% Shemat-Suite version: ', version
    WRITE(76,'(2A)') '%          build: ', datum
    WRITE(76,'(2A)') '%   build command line: ', makecmd
    WRITE(76,'(3A)') '%        Project: "', filename(1:lblank(filename)), '"'
    WRITE(76,'(1A)') '%'
    IF (transient) THEN
      WRITE(76,'(3A)') '% <time step>', ' <deltat>', ' <cum time>'
    END IF
#ifdef head_base
    WRITE(76,'(4A)') '%    <iteration>',' <delta head>',' <delta temp>',' (<delta conc> ...)'
#endif
#ifdef pres_base
    WRITE(76,'(4A)') '%    <iteration>',' <delta pres>',' <delta temp>',' (<delta conc> ...)'
#endif
    CLOSE(76)
  END IF
  
end subroutine write_status_log
