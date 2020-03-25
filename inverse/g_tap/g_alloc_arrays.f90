!>    @brief allocate derivative objects (used for AD)
!>    @param[in] ismpl local sample index
      SUBROUTINE g_alloc_arrays(ismpl)
        use arrays
        use g_arrays
        use mod_genrl
        use mod_genrlc
        use mod_inverse
        use mod_data
        use mod_conc
        use mod_time
        use g_mod_time
        use mod_linfos
        use mod_OMP_TOOLS
        IMPLICIT NONE
        integer :: ismpl
        INCLUDE 'OMP_TOOLS.inc'
        INTRINSIC max

!     thread stuff
        INTEGER tpos, tanz, l2

        IF (linfos(1)>=2) WRITE(*,*) ' [I] : ... g_alloc_arrays'

! allocating g_prad
        ALLOCATE(g_prad(i0,j0,k0,nsmpl))
        memory = memory+i0*j0*k0*nsmpl
! allocating g_conc_conv - Temporary conversions list
        ALLOCATE(g_conc_conv(ntrac,nsmpl))
        memory = memory+ntrac*nsmpl
! allocating g_omp_dglobal
        ALLOCATE(g_omp_dglobal(tlevel_1,9,nsmpl))
        memory = memory+tlevel_1*9*nsmpl

! hydraulic potential,pressure
        ALLOCATE(g_head(i0,j0,k0,nsmpl))
        memory = memory + i0*j0*k0*nsmpl
! temperature
        ALLOCATE(g_temp(i0,j0,k0,nsmpl))
        memory = memory + i0*j0*k0*nsmpl

        ALLOCATE(g_pres(i0,j0,k0,nsmpl))
        memory = memory + i0*j0*k0*nsmpl
 
        ALLOCATE(g_conc(i0,j0,k0,max(ntrans,1),nsmpl))
        memory = memory + i0*j0*k0*max(ntrans,1)*nsmpl

        ALLOCATE(g_tsal(i0,j0,k0,nsmpl))
        memory = memory + i0*j0*k0*nsmpl

!     inversion
        ALLOCATE(g_propunit(nunits,nprop,nsmpl))
        memory = memory + nunits*nprop*nsmpl

! coefficients for equation system solver
        ALLOCATE(g_a(i0,j0,k0,nsmpl))
        memory = memory + i0*j0*k0*nsmpl
        ALLOCATE(g_b(i0,j0,k0,nsmpl))
        memory = memory + i0*j0*k0*nsmpl
        ALLOCATE(g_c(i0,j0,k0,nsmpl))
        memory = memory + i0*j0*k0*nsmpl
        ALLOCATE(g_d(i0,j0,k0,nsmpl))
        memory = memory + i0*j0*k0*nsmpl
        ALLOCATE(g_e(i0,j0,k0,nsmpl))
        memory = memory + i0*j0*k0*nsmpl
        ALLOCATE(g_f(i0,j0,k0,nsmpl))
        memory = memory + i0*j0*k0*nsmpl
        ALLOCATE(g_g(i0,j0,k0,nsmpl))
        memory = memory + i0*j0*k0*nsmpl
        ALLOCATE(g_w(i0,j0,k0,nsmpl))
        memory = memory + i0*j0*k0*nsmpl
        ALLOCATE(g_x(i0,j0,k0,nsmpl))
        memory = memory + i0*j0*k0*nsmpl

! storing old temp,head,conc for iteration
!  4 == index of cgen_opti (needed)
        ALLOCATE(g_headold(i0*j0*k0,4,nsmpl))
        ALLOCATE(g_tempold(i0*j0*k0,4,nsmpl))
        ALLOCATE(g_presold(i0*j0*k0,4,nsmpl))
        memory = memory + i0*j0*k0*4*3*nsmpl
        ALLOCATE(g_concold(i0*j0*k0,max(ntrans,1),4,nsmpl))
        memory = memory + i0*j0*k0*max(ntrans,1)*4*nsmpl

!     boundary conditions
!      value
        ALLOCATE(g_dbc_data(nbc_data,ndbc,nsmpl))
        memory = memory + nbc_data*ndbc*nsmpl
        ALLOCATE(g_sdata(max(ndata,1),nsmpl))
        memory = memory + max(ndata,1)*nsmpl

        ALLOCATE(g_diff_c(max(ntrans,1)))
        memory = memory + max(ntrans,1)
        ALLOCATE(g_mmas_c(max(ntrans,1)))
        memory = memory + max(ntrans,1)

        ALLOCATE(g_bcperiod(ngsmax,3,max(nbctp,1),nsmpl))
        memory = memory + ngsmax*3*max(nbctp,1)*nsmpl
        ALLOCATE(g_delta_time(ntimestep))
        memory = memory + ntimestep
        ALLOCATE(g_simtime(nsmpl))
        memory = memory + nsmpl

!     initialisation of this kind, is needfull for NUMA architectures !!!
! ----------- NUMA -------------
        IF (nested_build) THEN
!$OMP   parallel default(none) private(l2) shared(nsmpl,Tlevel_0) &
!$OMP     num_threads(Tlevel_0)
!           Tlevel_0 = nsmpl !!!
!         ompenmp-critical to avoid ScaleMP performance bug during first initialisation
!$OMP     critical
          l2 = omp_get_his_thread_num() + 1
          CALL g_numa_init(l2)
!$OMP     end critical
!$OMP   end parallel
        ELSE
          DO l2 = 1, nsmpl
!         Tlevel_0 = 1 !!!
            CALL g_numa_init(l2)
          END DO
        END IF
! ----------- NUMA -------------

        g_dbc_data = 0.D0
        g_sdata = 0.D0
!
        g_diff_c = 0.D0
        g_mmas_c = 0.D0
!
        g_delta_time = 0.0d0
!
! inverse
        CALL set_dval(nunits*nprop*nsmpl,0.D0,g_propunit)
!
        g_thetaf = 0.0D0
        g_thetat = 0.0D0
!
        RETURN
      END

!>    @brief initialisation of this kind, is needfull for NUMA architectures !!!
!>    @param[in] ismpl local sample index
      SUBROUTINE g_numa_init(ismpl)
        use arrays
        use g_arrays
        use mod_genrl
        use mod_data
        use mod_conc
        use mod_time
        use mod_OMP_TOOLS
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k, l
        INCLUDE 'OMP_TOOLS.inc'

        INTRINSIC dble
!     thread stuff
        INTEGER tpos, tanz, l3


! ----------- NUMA -------------
!$OMP parallel private(tpos,tanz,i,j,k,l,l3) default(none) &
!$OMP   num_threads(Tlevel_1) shared(Tlevel_1) &
!$OMP   shared(g_head,g_temp,g_pres,g_conc,g_tsal) &
!$OMP   shared(g_a,g_b,g_c,g_d,g_e,g_f,g_g,g_w,g_x) &
!$OMP   shared(g_tempold,g_headold,g_presold) &
!$OMP   shared(g_concold,I0,J0,K0,ntrans,ismpl)
!$      call omp_binding(ismpl)

        CALL omp_part(i0*j0*k0,tpos,tanz)
        CALL ijk_m(tpos,i,j,k)
        CALL set_dval(tanz,0.D0,g_a(i,j,k,ismpl))
        CALL set_dval(tanz,0.D0,g_b(i,j,k,ismpl))
        CALL set_dval(tanz,0.D0,g_c(i,j,k,ismpl))
        CALL set_dval(tanz,0.D0,g_d(i,j,k,ismpl))
        CALL set_dval(tanz,0.D0,g_e(i,j,k,ismpl))
        CALL set_dval(tanz,0.D0,g_f(i,j,k,ismpl))
        CALL set_dval(tanz,0.D0,g_g(i,j,k,ismpl))
        CALL set_dval(tanz,0.D0,g_w(i,j,k,ismpl))
        CALL set_dval(tanz,0.D0,g_x(i,j,k,ismpl))
        CALL set_dval(tanz,0.D0,g_head(i,j,k,ismpl))
        CALL set_dval(tanz,0.D0,g_temp(i,j,k,ismpl))
        CALL set_dval(tanz,0.D0,g_pres(i,j,k,ismpl))
        DO l3 = 1, ntrans
          CALL set_dval(tanz,0.D0,g_conc(i,j,k,l3,ismpl))
        END DO
        CALL set_dval(tanz,0.D0,g_tsal(i,j,k,ismpl))
        DO l = 1, 4
          CALL set_dval(tanz,0.D0,g_tempold(tpos,l,ismpl))
          CALL set_dval(tanz,0.D0,g_headold(tpos,l,ismpl))
          CALL set_dval(tanz,0.D0,g_presold(tpos,l,ismpl))
          DO l3 = 1, ntrans
            CALL set_dval(tanz,0.D0,g_concold(tpos,l3,l,ismpl))
          END DO
        END DO
!$OMP end parallel
! ----------- NUMA -------------
        RETURN
      END SUBROUTINE g_numa_init
