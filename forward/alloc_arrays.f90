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

!>    @brief allocate all global main arrays (dynamic size)
!>    @param[in] ismpl local sample index
!>    @details
!> an additional ccNUMA initialisation is performed after the allocation,\n
!> see "numa_init"\n
      SUBROUTINE alloc_arrays(ismpl)
        use arrays
        use mod_genrl
        use mod_genrlc
        use mod_data
        use mod_time
        use mod_conc
        use mod_blocking_size
        use mod_linfos
        use mod_OMP_TOOLS
        IMPLICIT NONE
        integer :: ismpl
        integer :: i
        INCLUDE 'OMP_TOOLS.inc'
        INTEGER l2
!       thread stuff
        ! INTEGER tpos, tanz
        INTEGER mfactor
!       only for benchmarking
#ifdef BENCH
        DOUBLE PRECISION trun, tend
#endif

!
        IF (linfos(1)>=2) WRITE(*,*) ' [I] : ... alloc_arrays'
#ifdef BENCH
        CALL sys_cputime(trun)
#endif

!       single system size (factor == 1)
        mfactor = 1

        ALLOCATE(project_sfx(nsmpl))
        DO i = 1, nsmpl
          project_sfx(i) = ''
        END DO

        ALLOCATE(propunit(nunits,nprop,nsmpl))
        memory = memory + nunits*nprop*nsmpl
        CALL set_dval(nunits*nprop*nsmpl,0.D0,propunit)

!       temporary convergence list
        ALLOCATE(conc_conv(ntrac,nsmpl))
        memory = memory + ntrac*nsmpl

!       inversion
        ALLOCATE(opti_props(nprop,nunits))
        memory = memory + nunits*nprop
        ALLOCATE(opti_bc(nbc,nunits))
        memory = memory + nunits*nbc
        ALLOCATE(a_propunit(nunits,nprop))
        ALLOCATE(d_propunit(nunits,nprop))
        ALLOCATE(e_propunit(nunits,nprop))
        memory = memory + 3*nunits*nprop
        CALL set_dval(nunits*nprop,0.D0,a_propunit)
        CALL set_dval(nunits*nprop,0.D0,d_propunit)
        CALL set_dval(nunits*nprop,0.D0,e_propunit)

        ALLOCATE(node_info(i0,j0,k0))
        memory = memory + i0*j0*k0
!     additional global & private vectors for linear system solver
!       global buffer for boundary exchange (+ismpl)
        ALLOCATE(lss_bound_block(block_i*block_j+block_i*block_k+ block_j*block_k,bdim_i,bdim_j,bdim_k,2,nsmpl))
        memory = memory + (block_i*block_j+block_i*block_k+block_j*block_k)*bdim_i*bdim_j*bdim_k*2*nsmpl
        ALLOCATE(lss_dnrm(i0*j0*k0*mfactor,nsmpl))
        memory = memory + i0*j0*k0*mfactor*nsmpl
        ALLOCATE(lss_tmp(i0*j0*k0,nsmpl))
        memory = memory + i0*j0*k0*nsmpl
!       private copy for preconditioning (+Tlevel_1 +ismpl)
        ALLOCATE(lss_lma(max_blocks*block_i*block_j*block_k,tlevel_1, nsmpl))
        memory = memory + max_blocks*block_i*block_j*block_k*tlevel_1*nsmpl
        ALLOCATE(lss_lmb(max_blocks*block_i*block_j*block_k,tlevel_1, nsmpl))
        memory = memory + max_blocks*block_i*block_j*block_k*tlevel_1*nsmpl
        ALLOCATE(lss_lmc(max_blocks*block_i*block_j*block_k,tlevel_1, nsmpl))
        memory = memory + max_blocks*block_i*block_j*block_k*tlevel_1*nsmpl
        ALLOCATE(lss_lmd(max_blocks*block_i*block_j*block_k,tlevel_1, nsmpl))
        memory = memory + max_blocks*block_i*block_j*block_k*tlevel_1*nsmpl
        ALLOCATE(lss_lme(max_blocks*block_i*block_j*block_k,tlevel_1, nsmpl))
        memory = memory + max_blocks*block_i*block_j*block_k*tlevel_1*nsmpl
        ALLOCATE(lss_lmf(max_blocks*block_i*block_j*block_k,tlevel_1, nsmpl))
        memory = memory + max_blocks*block_i*block_j*block_k*tlevel_1*nsmpl
        ALLOCATE(lss_lmg(max_blocks*block_i*block_j*block_k,tlevel_1, nsmpl))
        memory = memory + max_blocks*block_i*block_j*block_k*tlevel_1*nsmpl
        ALLOCATE(lss_lud(max_blocks*block_i*block_j*block_k,tlevel_1, nsmpl))
        memory = memory + max_blocks*block_i*block_j*block_k*tlevel_1*nsmpl
        ALLOCATE(lss_lx(max_blocks*block_i*block_j*block_k,tlevel_1, nsmpl))
        memory = memory + max_blocks*block_i*block_j*block_k*tlevel_1*nsmpl
        ALLOCATE(lss_lb(max_blocks*block_i*block_j*block_k,tlevel_1, nsmpl))
        memory = memory + max_blocks*block_i*block_j*block_k*tlevel_1*nsmpl
!
        ALLOCATE(lss_lloctmp(max_blocks*block_i*block_j*block_k, max_loctmp,tlevel_1,nsmpl))
        memory = memory + max_blocks*block_i*block_j*block_k*max_loctmp*tlevel_1*nsmpl
!
        ALLOCATE(lss_ldnrm(max_blocks*block_i*block_j*block_k, tlevel_1,nsmpl))
        memory = memory + max_blocks*block_i*block_j*block_k*tlevel_1*nsmpl
        ALLOCATE(lss_ud_block(block_i*block_j+block_i*block_k+ block_j*block_k,max_blocks,tlevel_1,nsmpl))
        memory = memory + (block_i*block_j+block_i*block_k+block_j*block_k)*max_blocks*tlevel_1*nsmpl
        ALLOCATE(lss_lr0_hat(max_blocks*block_i*block_j*block_k, tlevel_1,nsmpl))
        memory = memory + max_blocks*block_i*block_j*block_k*tlevel_1*nsmpl
!     sanity check
        IF (max_blocks*block_i*block_j*block_k*tlevel_1<i0*j0*k0*mfactor) THEN
!        when fewer than I0*J0*K0 memory is allocated,
!          then we have lost some elements in "par_init2(I0,J0,K0)"
          WRITE(*,'(1A)') 'error: software bug, something goes wrong in "alloc_arrays"!'
          STOP
        END IF

!
!       variables for variable time stepping
        allocate(flag_delt(nsmpl))
        memory = memory + nsmpl
        ALLOCATE(delt_count(nsmpl))
        memory = memory + nsmpl
        allocate(flag_1st_timestep(nsmpl))
        memory = memory + nsmpl
        allocate(delt_old(nsmpl))
        memory = memory + nsmpl
!

!
!       coefficients for linear equations
        ALLOCATE(a(i0,j0,k0*mfactor,nsmpl))
        memory = memory + i0*j0*k0*mfactor*nsmpl
        ALLOCATE(b(i0,j0,k0*mfactor,nsmpl))
        memory = memory + i0*j0*k0*mfactor*nsmpl
        ALLOCATE(c(i0,j0,k0*mfactor,nsmpl))
        memory = memory + i0*j0*k0*mfactor*nsmpl
        ALLOCATE(d(i0,j0,k0*mfactor,nsmpl))
        memory = memory + i0*j0*k0*mfactor*nsmpl
        ALLOCATE(e(i0,j0,k0*mfactor,nsmpl))
        memory = memory + i0*j0*k0*mfactor*nsmpl
        ALLOCATE(f(i0,j0,k0*mfactor,nsmpl))
        memory = memory + i0*j0*k0*mfactor*nsmpl
        ALLOCATE(g(i0,j0,k0*mfactor,nsmpl))
        memory = memory + i0*j0*k0*mfactor*nsmpl
!       formerly the solution vector - linear system, now a temporary helper vector
        ALLOCATE(x(i0,j0,k0*mfactor,nsmpl))
        memory = memory + i0*j0*k0*mfactor*nsmpl
!       only for ilu-precond. (shadow vectors)
        ALLOCATE(ud(i0,j0,k0*mfactor,nsmpl))
        memory = memory + i0*j0*k0*mfactor*nsmpl
!       rhs
        ALLOCATE(w(i0,j0,k0*mfactor,nsmpl))
        memory = memory + i0*j0*k0*mfactor*nsmpl
!       r0_hat for BiCGstab
        ALLOCATE(r(i0,j0,k0*mfactor))
        memory = memory + i0*j0*k0*mfactor
!       Dirichlet mask for special diagonal preconditioning
        ALLOCATE(bc_mask(i0*j0*k0*mfactor,nsmpl))
        memory = memory + i0*j0*k0*mfactor*nsmpl

!       state variables
        ALLOCATE(head(i0,j0,k0,nsmpl))
        memory = memory + i0*j0*k0*nsmpl
        ALLOCATE(temp(i0,j0,k0,nsmpl))
        memory = memory + i0*j0*k0*nsmpl
        ALLOCATE(conc(i0,j0,k0,max(ntrans,1),nsmpl))
        memory = memory + i0*j0*k0*max(ntrans,1)*nsmpl
        ALLOCATE(pres(i0,j0,k0,nsmpl))
        memory = memory + i0*j0*k0*nsmpl
        ALLOCATE(tsal(i0,j0,k0,nsmpl))
        memory = memory + i0*j0*k0*nsmpl

!       uindex(allocated in read_model)
!       memory = memory + i0*j0*k0

! inverse
        ALLOCATE(headold(i0*j0*k0,ncgen,nsmpl))
        memory = memory + i0*j0*k0*ncgen*nsmpl
        ALLOCATE(tempold(i0*j0*k0,ncgen,nsmpl))
        memory = memory + i0*j0*k0*ncgen*nsmpl
        ALLOCATE(concold(i0*j0*k0,max(ntrans,1),ncgen,nsmpl))
        memory = memory + i0*j0*k0*ncgen*nsmpl*max(ntrans,1)
        ALLOCATE(presold(i0*j0*k0,ncgen,nsmpl))
        memory = memory + i0*j0*k0*ncgen*nsmpl

        ALLOCATE(delx(i0))
        memory = memory + i0
        ALLOCATE(dely(j0))
        memory = memory + j0
        ALLOCATE(delz(k0))
        memory = memory + k0
        ALLOCATE(delxa(i0))
        memory = memory + i0
        ALLOCATE(delya(j0))
        memory = memory + j0
        ALLOCATE(delza(k0))
        memory = memory + k0

!     time periods - dummy -
        ALLOCATE(bcperiod(ngsmax,3,max(nbctp,1),nsmpl))
        memory = memory + ngsmax*3*max(nbctp,1)*nsmpl
        ALLOCATE(ibcperiod(max(nbctp,1)))
        memory = memory + max(nbctp,1)
        ALLOCATE(lbcperiod(ngsmax,max(nbctp,1)))
        memory = memory + ngsmax*max(nbctp,1)
        ALLOCATE(outt(1))
        memory = memory + 1

        ALLOCATE(smon_idx(nsmpl))
        memory = memory + nsmpl
        CALL set_ival(nsmpl,0,smon_idx)

        ! dummy allocation
        allocate(delta_time(1))
        memory = memory + 1

        ALLOCATE(opti_tp(3,ngsmax*mopti_tp*max(nbctp,1)))
        memory = memory + 3*ngsmax*mopti_tp*max(nbctp,1)
        ALLOCATE(a_bcperiod(ngsmax,2,max(nbctp,1)))
        ALLOCATE(d_bcperiod(ngsmax,2,max(nbctp,1)))
        ALLOCATE(e_bcperiod(ngsmax,2,max(nbctp,1)))
        ALLOCATE(simtime(nsmpl))
        memory = memory + nsmpl

        allocate(tr_switch(nsmpl))
        memory = memory + nsmpl
        do i = 1, nsmpl
          tr_switch(i) = .true.
        end do

        ALLOCATE(fh_table(c_fhandler,tlevel_0))
        memory = memory + c_fhandler*tlevel_0

!     reset file handler table
        CALL omp_new_file_handler(l2,0)

        ALLOCATE(diff_c(max(ntrans,1)))
        memory = memory + max(ntrans,1)
        ALLOCATE(mmas_c(max(ntrans,1)))
        memory = memory + max(ntrans,1)
        ALLOCATE(beta_c(max(ntrans,1)))
        memory = memory + max(ntrans,1)

!     boundary structures
        ALLOCATE(ibc_data(max(nbc_data,1),nibc))
        memory = memory + max(nbc_data,1)*nibc
        ALLOCATE(dbc_data(max(nbc_data,1),ndbc,nsmpl))
        memory = memory + max(nbc_data,1)*ndbc*nsmpl
        ALLOCATE(dbc_dataold(max(nbc_data,1)))
        memory = memory + max(nbc_data,1)

!       borehole logs
        ALLOCATE(ibh_pos(2,nbh_logs))
        memory = memory + 2*nbh_logs
        ALLOCATE(cbh_name(nbh_logs))
        memory = memory + nbh_logs*64

!     - convergency history buffer -
        ALLOCATE(conv_history(conv_hlen,conv_hmax,nsmpl))
        memory = memory + conv_hlen*conv_hmax*nsmpl
        ALLOCATE(conv_chlen(conv_hmax,nsmpl))
        ALLOCATE(conv_ipos(conv_hmax,nsmpl))
        ALLOCATE(lcon(conv_hmax,nsmpl))
        memory = memory + 3*conv_hmax*nsmpl

!     OpenMP specific REDUCTION stuff
        ALLOCATE(omp_dglobal(tlevel_1,9,nsmpl))
        ALLOCATE(omp_iglobal(tlevel_1,3,nsmpl))
        memory = memory + 12*tlevel_1*nsmpl

#ifdef BENCH
        CALL sys_cputime(tend)
        WRITE(*,'(1A,1F8.2,1A)') &
          '  [I] : memory allocation time (serial) =', tend - trun, &
          ' sec'
        trun = tend
#endif

!     initialisation of this kind, is needfull for NUMA architectures !!!
! ----------- NUMA -------------
        IF (nested_build) THEN
!$OMP   parallel default(none) private(l2) shared(nsmpl,Tlevel_0,mfactor) &
!$OMP    num_threads(Tlevel_0)
!           Tlevel_0 = nsmpl !!!
!         openmp-critical to avoid ScaleMP performance bug during first initialisation
!$OMP     critical
          l2 = omp_get_his_thread_num() + 1
          CALL numa_init(mfactor,l2)
!$OMP     end critical
!$OMP   end parallel
        ELSE
          DO l2 = 1, nsmpl
!         Tlevel_0 = 1 !!!
            CALL numa_init(mfactor,l2)
          END DO
        END IF
! ----------- NUMA -------------

#ifdef BENCH
        CALL sys_cputime(tend)
        WRITE(*,'(1A,1F8.2,1A)') &
          '  [I] : memory allocation time (parallel) =', tend - trun, &
          ' sec'
#endif


!     Dummy allocation
!       main_input (allocated in forward_init)
!       memory = memory + max(mpara,1)*nsmpl
!       main_output (allocated in forward_init)
!       memory = memory + max(ndata,1)*nsmpl
        ALLOCATE(sdata(1,nsmpl))
        memory = memory + nsmpl

        RETURN
      END SUBROUTINE alloc_arrays

!>    @brief initialisation of this kind, is needfull for NUMA architectures !!!
!>    @param[in] mfactor multiply factor for multi-phase methods
!>    @param[in] ismpl local sample index
      SUBROUTINE numa_init(mfactor,ismpl)
        use arrays
        use mod_genrl
        use mod_data
        use mod_conc
        use mod_time
        use mod_OMP_TOOLS
        use mod_blocking_size
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k, l

        INCLUDE 'OMP_TOOLS.inc'
        INTRINSIC dble
!     thread stuff
        INTEGER tpos, tanz, l2, loc_mem, t_id, mfactor


! ----------- NUMA -------------
!$OMP parallel default(none) private(tpos,tanz,i,j,k,l,l2,loc_mem,t_id)&
!$OMP  num_threads(Tlevel_1)&
!$OMP  shared(a,b,c,d,e,f,g,w,x,head,temp,pres,r,I0,J0,K0)&
!$OMP  shared(ismpl,conc,concold,tsal,ntrans)&
!$OMP  shared(headold,tempold,presold)&
!$OMP  shared(ud,ntrac,Tlevel_1,mfactor)&
!$OMP  shared(lss_bound_block,lss_dnrm,lss_tmp,lss_lUD,lss_lx,lss_lb)&
!$OMP  shared(lss_lMA,lss_lMB,lss_lMC,lss_lMD,lss_lME,lss_lMF,lss_lMG)&
!$OMP  shared(lss_llocTMP,lss_ldnrm,lss_ud_block,lss_lr0_hat)&
!$OMP  shared(max_blocks,block_i,block_j,block_k,bdim_i,bdim_j,bdim_k)
!$        call omp_binding(ismpl)
!       init
        t_id = omp_get_his_thread_num() + 1
        CALL omp_part(i0*j0*k0,tpos,tanz)
        CALL ijk_m(tpos,i,j,k)
        loc_mem = max_blocks*block_i*block_j*block_k
!
        l2 = (block_i*block_j+block_i*block_k+block_j*block_k)*bdim_i*bdim_j*bdim_k*2
        CALL set_dval(l2,0.D0,lss_bound_block(1,1,1,1,1,ismpl))
!
        CALL set_dval(tanz,0.D0,lss_tmp(tpos,ismpl))
!
        CALL set_dval(loc_mem,0.D0,lss_lma(1,t_id,ismpl))
        CALL set_dval(loc_mem,0.D0,lss_lmb(1,t_id,ismpl))
        CALL set_dval(loc_mem,0.D0,lss_lmc(1,t_id,ismpl))
        CALL set_dval(loc_mem,0.D0,lss_lmd(1,t_id,ismpl))
        CALL set_dval(loc_mem,0.D0,lss_lme(1,t_id,ismpl))
        CALL set_dval(loc_mem,0.D0,lss_lmf(1,t_id,ismpl))
        CALL set_dval(loc_mem,0.D0,lss_lmg(1,t_id,ismpl))
        CALL set_dval(loc_mem,0.D0,lss_lud(1,t_id,ismpl))
        CALL set_dval(loc_mem,0.D0,lss_lx(1,t_id,ismpl))
        CALL set_dval(loc_mem,0.D0,lss_lb(1,t_id,ismpl))
        CALL set_dval(loc_mem*max_loctmp,0.D0,lss_lloctmp(1,1,t_id,ismpl))
        CALL set_dval(loc_mem,0.D0,lss_ldnrm(1,t_id,ismpl))
        CALL set_dval(loc_mem,0.D0,lss_lr0_hat(1,t_id,ismpl))
        l2 = (block_i*block_j+block_i*block_k+block_j*block_k)*max_blocks
        CALL set_dval(l2,0.D0,lss_ud_block(1,1,t_id,ismpl))
!
        DO l = 0, mfactor-1
          CALL set_dval(tanz,0.D0,lss_dnrm(tpos+l*I0*J0*K0,ismpl))
          CALL set_dval(tanz,0.D0,a(i,j,k+l*K0,ismpl))
          CALL set_dval(tanz,0.D0,b(i,j,k+l*K0,ismpl))
          CALL set_dval(tanz,0.D0,c(i,j,k+l*K0,ismpl))
          CALL set_dval(tanz,0.D0,d(i,j,k+l*K0,ismpl))
          CALL set_dval(tanz,0.D0,e(i,j,k+l*K0,ismpl))
          CALL set_dval(tanz,0.D0,f(i,j,k+l*K0,ismpl))
          CALL set_dval(tanz,0.D0,g(i,j,k+l*K0,ismpl))
          CALL set_dval(tanz,0.D0,w(i,j,k+l*K0,ismpl))
          CALL set_dval(tanz,0.D0,x(i,j,k+l*K0,ismpl))
          CALL set_dval(tanz,0.D0,ud(i,j,k+l*K0,ismpl))
        END DO
!
        CALL set_dval(tanz,0.D0,head(i,j,k,ismpl))
        CALL set_dval(tanz,10.0D0,temp(i,j,k,ismpl))
        CALL set_dval(tanz,0.D0,pres(i,j,k,ismpl))
        DO l = 1, max(ntrans,1)
          CALL set_dval(tanz,0.D0,conc(i,j,k,l,ismpl))
        END DO
        CALL set_dval(tanz,0.D0,tsal(i,j,k,ismpl))

        DO l = 1, ncgen
          CALL set_dval(tanz,100.0D0,headold(tpos,l,ismpl))
        END DO
        DO l = 1, ncgen
          CALL set_dval(tanz,20.0D0,tempold(tpos,l,ismpl))
        END DO
        DO l = 1, ncgen
          CALL set_dval(tanz,1.5D7,presold(tpos,l,ismpl))
        END DO
        DO l = 1, ncgen
          DO l2 = 1, max(ntrans,1)
            CALL set_dval(tanz,0.D0,concold(tpos,l2,l,ismpl))
          END DO
        END DO
!       independent of nsmpl
        IF (ismpl==1) THEN
!         [flow,temp,conc*]
          DO l = 1, mfactor
            CALL set_dval(tanz,0.D0,r(i,j,k+(l-1)*k0))
          END DO
!$OMP     barrier
!$OMP     do
          DO k = 1, k0
            DO j = 1, j0
              DO i = 1, i0
!               random vector needed for BiCGStab lineare solver
!                 later can be initialised from a real random-number-generator
                r(i,j,k) = 1.0D3/dble(0.25D0+i0/2-i+j0/2-j+k0/2-k)
!               [flow,temp,conc*], [pres]
                DO l = 2, mfactor
                  r(i,j,k+(l-1)*k0) = r(i,j,k+(l-2)*k0)*0.1D0
                END DO
              END DO
            END DO
          END DO
!$OMP     end do
        END IF
!$OMP end parallel
! ----------- NUMA -------------
        RETURN
      END SUBROUTINE numa_init
