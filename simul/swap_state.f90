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

!> @brief parallelisation wrapper for "omp_swap_state"
!> @param[in] ismpl_a sample A to exchange with B
!> @param[in] ismpl_b sample B to exchange with A
!> @param[in] ismpl local sample index (here only for omp_binding)
      SUBROUTINE swap_state(ismpl_a, ismpl_b, ismpl)
        use mod_genrl
        use mod_OMP_TOOLS
        IMPLICIT NONE
        INCLUDE 'OMP_TOOLS.inc'
        integer :: ismpl
        integer :: ismpl_a, ismpl_b

#ifdef fOMP
!$OMP parallel num_threads(Tlevel_1)
!$      call omp_binding(ismpl)
#endif
        CALL omp_swap_state(ismpl_a, ismpl_b, ismpl)
#ifdef fOMP
!$OMP end parallel
#endif

        RETURN
      END

!> @brief exchange two (sample) states
!> @param[in] ismpl_a sample A to exchange with B
!> @param[in] ismpl_b sample B to exchange with A
!> @param[in] ismpl local sample index (dummy here)
      SUBROUTINE omp_swap_state(ismpl_a, ismpl_b, ismpl)
        use arrays
        USE simul_arrays
        use mod_genrl
        use mod_conc
        use mod_simul
        use mod_time
        use mod_data
        IMPLICIT NONE
!     realisation
        integer :: ismpl
        integer :: i, j, k, l
        integer :: ismpl_a, ismpl_b
        integer :: tpos, tanz
        integer :: iv
        DOUBLE PRECISION s_tmp, vdefault_tmp

        IF (ismpl_a==ismpl_b) RETURN
!
!$OMP master
        CALL dswap(nbc_data,dbc_data(1,1,ismpl_a),1,dbc_data(1,1,ismpl_b),1)
        CALL dswap(nunits*nprop,propunit(1,1,ismpl_a),1,propunit(1,1,ismpl_b),1)
        CALL dswap(ngsmax*3*nbctp,bcperiod(1,1,1,ismpl_a),1,bcperiod(1,1,1,ismpl_b),1)
!
        CALL dswap(mpara,main_input(1,ismpl_a),1,main_input(1,ismpl_b),1)
        CALL dswap(ndata,main_output(1,ismpl_a),1,main_output(1,ismpl_b),1)
!       swap simutaltion time
        s_tmp = simtime(ismpl_a)
        simtime(ismpl_a) = simtime(ismpl_b)
        simtime(ismpl_b) = s_tmp
!       swap default velocities
        if(vdefaultswitch) then
           do iv = 1, 3
              vdefault_tmp = vdefault(iv,ismpl_a)
              vdefault(iv,ismpl_a) = vdefault(iv,ismpl_b)
              vdefault(iv,ismpl_b) = vdefault_tmp
           end do
        end if
!$OMP end master
!
        CALL omp_part(i0*j0*k0,tpos,tanz)
        CALL ijk_m(tpos,i,j,k)
!
        CALL dswap(tanz,head(i,j,k,ismpl_a),1,head(i,j,k,ismpl_b),1)
        CALL dswap(tanz,temp(i,j,k,ismpl_a),1,temp(i,j,k,ismpl_b),1)
        DO l = 1, ntrans
          CALL dswap(tanz,conc(i,j,k,l,ismpl_a),1,conc(i,j,k,l,ismpl_b),1)
        END DO
        CALL dswap(tanz,pres(i,j,k,ismpl_a),1,pres(i,j,k,ismpl_b),1)
!
        CALL dswap(tanz,tsal(i,j,k,ismpl_a),1,tsal(i,j,k,ismpl_b),1)
!
        CALL dswap(tanz,w(i,j,k,ismpl_a),1,w(i,j,k,ismpl_b),1)
!
        RETURN
      END
