!>    @brief derivative version to solve a linear equation system
!>    @param[in] icode specify the physical value
!>    @param[in] species transport species index
!>    @param[out] x_solution the solution vector of the solved linear system
!>    @param[in] errorc error
!>    @param[in] apar setup explicit or implicit solver (disabled), default 1.0
!>    @param[in] ctrl solver code
!>    @param[in] ismpl local sample index
!>    @details
!>  solve system equations by preconditioned  krylov solvers (bicgstab) or SIP (NAG) and others\n
      SUBROUTINE g_solve(icode,species,x_solution,g_x_solution,errorc, &
          apar,ctrl,ismpl)
        use arrays
        use g_arrays
        use mod_conc
        use mod_flow
        use mod_genrl
        use mod_linfos
        use mod_temp
        use mod_OMP_TOOLS
        IMPLICIT NONE
        integer :: ismpl
        INCLUDE 'OMP_TOOLS.inc'
        DOUBLE PRECISION apar, errorc
        INTEGER icode, itst, mxit, species
        INTEGER ctrl, criteria, precondition, solvername
        DOUBLE PRECISION x_solution(i0*j0*k0)
!     only for tests with 'abbruch'
        DOUBLE PRECISION moderrorc, enough
        DOUBLE PRECISION trun, tend
        INTRINSIC dsqrt, dble, abs
!     -------------
        DOUBLE PRECISION d1
        DOUBLE PRECISION g_x_solution(i0*j0*k0)


        moderrorc = errorc
        IF (ctrl<0) THEN
          WRITE(*,*) 'error: old solver compatibility no &
            &longer        supported'
          STOP
        END IF
! **************************************************************************
! solvername =
!           *0 : BiCGStab [:parallel]
!            1 : NAG [:serial]
!            2 : CG [:parallel] ( prove symmetry, if not then BiCGStab )
!            3 : PLU [:serial] (LAPACK) and math tests (stability)
!         (4-7 : not in use !)
! criteria = < switch to set when should break >
!          0 : relative stopping crit. : ||[res]|| < depsilon*||[res0]||
!          1 : absolute stopping crit. : ||[res]|| < depsilon
!          2 : maximum  stopping crit. : max(abs([res])) < depsilon
!          3 : abs. and rel. stopping crit. : ( ||[res]|| < depsilon ) and ( ||[
!res]|| < 0.99d0*||[res0]|| )
!              0.99d0 is a constant for testing only, it is named 'minRel' in 'a
!bbruch.f'
!         *4 : like '3', but with auto detected range for depsilon, default deps
!ilon used when detection fails
!       (5-7 : not in use !)
! precondition =
!             *0 : ILU [:parallel]
!              1 : SSOR [:serial]
!              2 : Diagonal [:parallel]
!              3 : None [:]
!           (4-7 : not in use !)
! * : recommended (ctrl = 64)
! and (ctrl = 67) : for testing
! [ ctrl = solvername + 16*criteria + 256*precondition ]
!     extract [ ctrl ] :
!      solvername = mod(ctrl,16)
!      ctrl = ctrl/16
!      criteria = mod(ctrl,16)
!      ctrl = ctrl/16
!      precondition = ctrl
        CALL decntrl3(ctrl,solvername,criteria,precondition)
! **************************************************************************
!     -------------
        IF (abs(icode)==pv_pres) THEN
          mxit = lmaxitf
          IF (criteria>=4) THEN
            criteria = 1
            moderrorc = nltolf*1.0D-2
          END IF
        END IF
        IF (abs(icode)==pv_conc) THEN
          mxit = lmaxitc
          IF (criteria>=4) THEN
            criteria = 1
            moderrorc = nltolc*1.0D-2
          END IF
        END IF
        IF (abs(icode)==pv_temp) THEN
          mxit = lmaxitt
          IF (criteria>=4) THEN
            criteria = 1
            moderrorc = nltolt*1.0D-2
          END IF
        END IF
        IF (abs(icode)==pv_head) THEN
          mxit = lmaxitf
          IF (criteria>=4) THEN
            criteria = 1
            moderrorc = nltolf*1.0D-2
          END IF
        END IF
!     -------------
!     -------------

#ifdef fOMP
!$OMP parallel num_threads(Tlevel_1)
!$      call omp_binding(ismpl)
#endif
!       [lss_tmp(1,ismpl)] := [g_A]x[x]
        CALL omp_mvp(i0,j0,k0,x_solution,lss_tmp(1,ismpl), &
          g_a(1,1,1,ismpl),g_b(1,1,1,ismpl),g_c(1,1,1,ismpl), &
          g_d(1,1,1,ismpl),g_e(1,1,1,ismpl),g_f(1,1,1,ismpl), &
          g_g(1,1,1,ismpl))
#ifdef fOMP
!$OMP end parallel
#endif
!       [lss_tmp(1,ismpl)] : [g_b]-[g_A]x[x]=[g_b]-[lss_tmp(1,ismpl)]
        CALL dscal(i0*j0*k0,-1.D0,lss_tmp(1,ismpl),1)
        CALL daxpy(i0*j0*k0,1.D0,g_w(1,1,1,ismpl),1,lss_tmp(1,ismpl), &
          1)

!    ##################################################################
!    #                             solver                             #
!    ##################################################################
#ifdef BENCH
        CALL sys_cputime(trun)
#endif

        CALL solve_type(i0,j0,k0,g_x_solution,lss_tmp(1,ismpl), &
          moderrorc,bc_mask(1,ismpl),solvername,precondition,mxit, &
          criteria,a(1,1,1,ismpl),b(1,1,1,ismpl),c(1,1,1,ismpl), &
          d(1,1,1,ismpl),e(1,1,1,ismpl),f(1,1,1,ismpl),g(1,1,1,ismpl), &
          r,apar,ud(1,1,1,ismpl),ismpl)

#ifdef BENCH
        CALL sys_cputime(tend)
        WRITE(*,*) ' linear system solver time:', tend - trun, 'sec'
#endif
!    ##################################################################
!     -------------
!     "icode" < 0 indicates a default input/output in the [x] vector
        IF (icode<0) THEN
          CALL dcopy(i0*j0*k0,g_x_solution,1,g_x(1,1,1,ismpl),1)
        END IF

        RETURN
      END
