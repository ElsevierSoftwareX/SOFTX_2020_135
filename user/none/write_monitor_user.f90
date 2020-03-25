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

!> @brief writes data in tecplot-format (example routine)
!> @param[in] otype orientation type:\n
!> - 1 new file\n
!> - 2 append\n
!> @param[in] ismpl local sample index
      subroutine write_monitor_user(otype,ismpl)

        ! use mod_genrl, only: out_orientation
        ! use mod_genrlc, only: project
        ! use mod_time, only: simtime, tunit

        implicit none

        ! otype:
        ! 1  new file
        ! 2  append
        integer, intent (in) :: otype

        ! Samples index
        integer :: ismpl

        ! logical :: ltimes
        ! logical :: lmonip
        ! character (len=256) :: filename
        ! integer, external :: lblank


        ! ! Orientation of monitoring output
        ! if (out_orientation==2 .or. out_orientation==3) then
        !   ltimes = .true.
        ! else
        !   ltimes = .false.
        ! end if

        ! if (out_orientation==4) then
        !   lmonip = .true.
        ! else
        !   lmonip = .false.
        ! end if

        ! ! Output Filename
        ! call chln(project,i1,i2)
        ! filename = project(i1:i2)//'_monitor_user.dat'

        ! if (ltimes) then
        !   write(filename,'(2A,1e14.8,1A)') &
        !       project(i1:i2),'_', &
        !       (simtime(ismpl))/tunit,'_monitor_user.dat'
        ! endif

        ! if (linfos(3).ge.2.and..not.lmonip) then
        !   write(*,'(3A)') '  [W] : User Monitor points to "', &
        !       filename(1:lblank(filename)),'"'
        ! endif


        return

      end subroutine write_monitor_user
