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

!>    @brief dependency wrapper (parameter -> simulation) useful for INVERSE and SIMUL
!>    @param[in] dinput parameter vector for initial setup (dependency vector)
!>    @param[out] doutput simulated/computed data/result values
!>    @param[in] simtime_run start time of the simulation
!>    @param[in] simtime_end finish time of the simulation
!>    @param[in] iter_out inverse iteration, SM realisation
!>    @param[in] iseed 0: FW simulation, 1 .. <mpara>: AD seeding index
!>    @param[in] ismpl local sample index
      SUBROUTINE forward_compute(dinput, doutput, simtime_run, simtime_end, iter_out, iseed, ismpl)
        use arrays
        use mod_data
        use mod_genrl
        IMPLICIT NONE
        integer :: ismpl
        integer :: i
!       input vector
        DOUBLE PRECISION dinput(mpara)
!       output vector
        DOUBLE PRECISION doutput(ndata)
!       benchmark time meassurment stuff
#ifdef BENCH
        DOUBLE PRECISION trun, tend
#endif
!       parameter declaration
        INTEGER iter_out, iseed
        DOUBLE PRECISION simtime_run, simtime_end

!
#ifdef  BENCH
        CALL sys_cputime(trun)
#endif
!
!       dependency setup of the parameters
        DO i = 1, mpara
          CALL set_optip(i,dinput(i),ismpl)
        END DO
!       run forward iteration
        CALL forward_iter(simtime_run,simtime_end,iter_out,iseed,ismpl)
!       copy result values
        CALL dcopy(ndata,sdata(1,ismpl),1,doutput,1)
!
!
#ifdef  BENCH
        CALL sys_cputime(tend)
        WRITE(*,'(1A,I7,1A,F9.2,1A,2F20.2)') &
          '  [T] : component =', iseed, ', time =',tend - trun, ' sec, [run] [end]:', trun, tend
#endif
!
        RETURN
      END
