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

!> @brief Simulation time related manipulations
!> @param[in] i_simtime local index for different manipulations
!> @param[in] iter_enkf ENKF iteration counter
!> @details
!> Collection of simulation time related manipulations.
subroutine enkf_simtime(i_simtime, iter_enkf)

  use arrays, only:&
       project_sfx,&
       simtime
  
  use mod_time, only:&
       itimestep_0,&
       simtime_0,&
       ntimestep,&
       max_simtime,&
       tunit

  use mod_genrl, only:&
       nsmpl

  use mod_enkf, only:&
       itimestep_ismpl,&
       iter_out,&
       ttstart,&
       ttend,&
       nrens,&
       obst,&
       irobs

  implicit none

  integer, intent(in) :: i_simtime
  integer, intent(in) :: iter_enkf

  integer :: irens

  !---------------------------------------------------------------
  if (i_simtime == 1) then

     itimestep_ismpl = itimestep_0
     DO irens = 1, nsmpl
        WRITE(project_sfx(irens),'(1I5)') iter_enkf
        project_sfx(irens) = '_E' // adjustl(project_sfx(irens))
     END DO

     !---------------------------------------------------------------
  else if (i_simtime == 2) then

     ttstart = simtime_0
     DO irens = 1, nrens
        simtime(irens) = simtime_0
     END DO
     ! Outer iteration increade (-1) at every observation-interval
     iter_out = 0
     obst(0) = 0.D0

     !---------------------------------------------------------------
  else if  (i_simtime == 3) then

     itimestep_0 = max(itimestep_0 +1, &
          int( dble(ntimestep -itimestep_ismpl) *ttend /(max_simtime - simtime_0) ) )

     ! increase assimilation counter
     iter_out = iter_out +1

     !---------------------------------------------------------------
  else if (i_simtime == 4) then

     ttstart = simtime(1)

     !---------------------------------------------------------------
  else if (i_simtime == 5) then

     itimestep_0 = itimestep_ismpl

     !---------------------------------------------------------------
  else if (i_simtime == 6) then

     ttend = obst(irobs)*tunit

     !---------------------------------------------------------------
  else

     write(unit = *, fmt = *) 
     write(unit = *, fmt = *) "!!!!!!time!!!!!!!"
     write(unit = *, fmt = *) 
     write(unit = *, fmt = *) "WRONG i_simtime in enkf_iter"
     write(unit = *, fmt = *) 
     write(unit = *, fmt = *) "!!!!!!!!time!!!!!!"
     write(unit = *, fmt = *) 
     stop

  end if


end subroutine enkf_simtime
