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

!>    @brief free the memory of all global main arrays (with dynamic size)
!>    @param[in] ismpl local sample index
      SUBROUTINE dealloc_arrays(ismpl)
        use arrays
        use mod_genrl
        use mod_genrlc
        use mod_data
        use mod_time
        use mod_conc
        use mod_linfos
        use mod_OMP_TOOLS
        use mod_blocking_size
        IMPLICIT NONE
        integer :: ismpl
        INCLUDE 'OMP_TOOLS.inc'
        INTEGER mfactor
        INTRINSIC max

!
        IF (linfos(1)>=2) WRITE(*,*) ' [I] : ... dealloc_arrays'
!
!       single system size (factor == 1)
        mfactor = 1

        DEALLOCATE(project_sfx)

        DEALLOCATE(propunit)
        memory = memory - nunits*nprop*nsmpl

!       temporary convergence list
        DEALLOCATE(conc_conv)
        memory = memory - ntrac*nsmpl

!       inversion
        DEALLOCATE(opti_props)
        memory = memory - nunits*nprop
        DEALLOCATE(opti_bc)
        memory = memory - nunits*nbc
        DEALLOCATE(a_propunit)
        DEALLOCATE(d_propunit)
        DEALLOCATE(e_propunit)
        memory = memory - 3*nunits*nprop

        DEALLOCATE(node_info)
        memory = memory - i0*j0*k0
!     additional global & private vectors for linear system solver
!       global buffer for boundary exchange (+ismpl)
        DEALLOCATE(lss_bound_block)
        memory = memory - (block_i*block_j+block_i*block_k+block_j*block_k)*bdim_i*bdim_j*bdim_k*2*nsmpl
        DEALLOCATE(lss_dnrm)
        memory = memory - i0*j0*k0*nsmpl
        DEALLOCATE(lss_tmp)
        memory = memory - i0*j0*k0*nsmpl
!       private copy for preconditioning (+Tlevel_1 +ismpl)
        DEALLOCATE(lss_lma)
        memory = memory - max_blocks*block_i*block_j*block_k*tlevel_1*nsmpl
        DEALLOCATE(lss_lmb)
        memory = memory - max_blocks*block_i*block_j*block_k*tlevel_1*nsmpl
        DEALLOCATE(lss_lmc)
        memory = memory - max_blocks*block_i*block_j*block_k*tlevel_1*nsmpl
        DEALLOCATE(lss_lmd)
        memory = memory - max_blocks*block_i*block_j*block_k*tlevel_1*nsmpl
        DEALLOCATE(lss_lme)
        memory = memory - max_blocks*block_i*block_j*block_k*tlevel_1*nsmpl
        DEALLOCATE(lss_lmf)
        memory = memory - max_blocks*block_i*block_j*block_k*tlevel_1*nsmpl
        DEALLOCATE(lss_lmg)
        memory = memory - max_blocks*block_i*block_j*block_k*tlevel_1*nsmpl
        DEALLOCATE(lss_lud)
        memory = memory - max_blocks*block_i*block_j*block_k*tlevel_1*nsmpl
        DEALLOCATE(lss_lx)
        memory = memory - max_blocks*block_i*block_j*block_k*tlevel_1*nsmpl
        DEALLOCATE(lss_lb)
        memory = memory - max_blocks*block_i*block_j*block_k*tlevel_1*nsmpl
!
        DEALLOCATE(lss_lloctmp)
        memory = memory - max_blocks*block_i*block_j*block_k*mfactor*max_loctmp*tlevel_1*nsmpl
!
        DEALLOCATE(lss_ldnrm)
        memory = memory - max_blocks*block_i*block_j*block_k*tlevel_1*nsmpl
        DEALLOCATE(lss_ud_block)
        memory = memory - (block_i*block_j+block_i*block_k+block_j*block_k)*max_blocks*tlevel_1*nsmpl
        DEALLOCATE(lss_lr0_hat)
        memory = memory - max_blocks*block_i*block_j*block_k*tlevel_1*nsmpl

!       variables for variable time stepping
        deallocate(flag_delt)
        memory = memory - nsmpl
        deallocate(delt_count)
        memory = memory - nsmpl
        deallocate(flag_1st_timestep)
        memory = memory - nsmpl
        deallocate(delt_old)
        memory = memory - nsmpl
!
!
!       coefficients for linear equations
        DEALLOCATE(a)
        memory = memory - i0*j0*k0*mfactor*nsmpl
        DEALLOCATE(b)
        memory = memory - i0*j0*k0*mfactor*nsmpl
        DEALLOCATE(c)
        memory = memory - i0*j0*k0*mfactor*nsmpl
        DEALLOCATE(d)
        memory = memory - i0*j0*k0*mfactor*nsmpl
        DEALLOCATE(e)
        memory = memory - i0*j0*k0*mfactor*nsmpl
        DEALLOCATE(f)
        memory = memory - i0*j0*k0*mfactor*nsmpl
        DEALLOCATE(g)
        memory = memory - i0*j0*k0*mfactor*nsmpl
        DEALLOCATE(x)
        memory = memory - i0*j0*k0*mfactor*nsmpl
!
        DEALLOCATE(w)
        memory = memory - i0*j0*k0*mfactor*nsmpl
        DEALLOCATE(r)
        memory = memory - i0*j0*k0*mfactor
        DEALLOCATE(bc_mask)
        memory = memory - i0*j0*k0*mfactor*nsmpl
!       only for ilu-precond. (shadow vectors)
        DEALLOCATE(ud)
        memory = memory - i0*j0*k0*mfactor*nsmpl
!
!
        DEALLOCATE(head)
        memory = memory - i0*j0*k0*nsmpl
        DEALLOCATE(temp)
        memory = memory - i0*j0*k0*nsmpl
        DEALLOCATE(conc)
        memory = memory - i0*j0*k0*max(ntrans,1)*nsmpl
        DEALLOCATE(pres)
        memory = memory - i0*j0*k0*nsmpl
        DEALLOCATE(tsal)
        memory = memory - i0*j0*k0*nsmpl
        DEALLOCATE(uindex)
        memory = memory - i0*j0*k0

        DEALLOCATE(headold)
        memory = memory - i0*j0*k0*ncgen*nsmpl
        DEALLOCATE(tempold)
        memory = memory - i0*j0*k0*ncgen*nsmpl
        DEALLOCATE(concold)
        memory = memory - i0*j0*k0*ncgen*nsmpl*max(ntrans,1)
        DEALLOCATE(presold)
        memory = memory - i0*j0*k0*ncgen*nsmpl

        DEALLOCATE(delx)
        memory = memory - i0
        DEALLOCATE(dely)
        memory = memory - j0
        DEALLOCATE(delz)
        memory = memory - k0
        DEALLOCATE(delxa)
        memory = memory - i0
        DEALLOCATE(delya)
        memory = memory - j0
        DEALLOCATE(delza)
        memory = memory - k0

        if (vdefaultswitch) then
           DEALLOCATE(vdefault)
        end if
        
!      dealloc. proz. grid array, see more in 'solve/omp_preconditioniers.f'
        CALL par_end2()

!     time periods
        DEALLOCATE(bcperiod)
        memory = memory - ngsmax*3*max(nbctp,1)*nsmpl
        DEALLOCATE(ibcperiod)
        memory = memory - max(nbctp,1)
        DEALLOCATE(lbcperiod)
        memory = memory - ngsmax*max(nbctp,1)
        DEALLOCATE(outt)
        memory = memory - max(noutt+1,1)
        DEALLOCATE(smon_idx)
        memory = memory - nsmpl
        deallocate(delta_time)
        memory = memory - max(ntimestep,1)
        DEALLOCATE(opti_tp)
        memory = memory - 3*ngsmax*mopti_tp*max(nbctp,1)

        DEALLOCATE(a_bcperiod)
        DEALLOCATE(d_bcperiod)
        DEALLOCATE(e_bcperiod)
        DEALLOCATE(simtime)
        memory = memory - nsmpl

        DEALLOCATE(tr_switch)
        memory = memory - nsmpl

        DEALLOCATE(fh_table)
        memory = memory - c_fhandler*tlevel_0

        DEALLOCATE(diff_c)
        memory = memory - max(ntrans,1)
        DEALLOCATE(mmas_c)
        memory = memory - max(ntrans,1)
        DEALLOCATE(beta_c)
        memory = memory - max(ntrans,1)

!     boundary structures
        DEALLOCATE(ibc_data)
        memory = memory - max(nbc_data,1)*nibc
        DEALLOCATE(dbc_data)
        memory = memory - max(nbc_data,1)*ndbc*nsmpl
        DEALLOCATE(dbc_dataold)
        memory = memory - max(nbc_data,1)

!       borehole logs
        DEALLOCATE(ibh_pos)
        memory = memory - 2*nbh_logs
        DEALLOCATE(cbh_name)
        memory = memory - nbh_logs*64

!     - convergency history buffer -
        DEALLOCATE(conv_history)
        memory = memory - conv_hlen*conv_hmax*nsmpl
        DEALLOCATE(conv_chlen)
        DEALLOCATE(conv_ipos)
        DEALLOCATE(lcon)
        memory = memory - 3*conv_hmax*nsmpl

#ifdef DEBUG
        DEALLOCATE(debugout)
        n_debugout = 0
#endif

!     OpenMP specific REDUCTION staff
        DEALLOCATE(omp_dglobal)
        DEALLOCATE(omp_iglobal)
        memory = memory - 12*tlevel_1*nsmpl


!       allocated in "forward_init"
        DEALLOCATE(main_input)
        memory = memory - max(mpara,1)*nsmpl
        DEALLOCATE(main_output)
        memory = memory - max(ndata,1)*nsmpl

        write(*,*) "Deallocated Memory, remaining ",memory*8/1024/1024, " MB"
        RETURN
      END SUBROUTINE dealloc_arrays
