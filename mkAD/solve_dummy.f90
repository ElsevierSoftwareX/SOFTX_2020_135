!>    @brief dummy routine for gradient/adjoint generation (trick the AD tool)
!>    @param[in] icode specify the physical value
!>    @param[in] species transport species index
!>    @param[in,out] x_solution solution vector of the solved linear system, may be also the start input
!>    @param[in] errorc error
!>    @param[in] apar setup explicit or implicit solver (disabled), default 1.0
!>    @param[in] ctrl solver code
!>    @param[in] ismpl local sample index
      SUBROUTINE solve(icode,species,x_solution,errorc,apar,ctrl,ismpl)
        use arrays
        use mod_genrl
        IMPLICIT NONE
        integer :: ismpl
        DOUBLE PRECISION apar, errorc, x_solution(i0*j0*k0)
        INTEGER icode, species, ctrl
        IF (icode>=-10 .OR. species>=-100 .OR. errorc>=0.0d0 .OR. apar>=0.0d0 .OR. ctrl>=-1000) THEN
          CALL dcopy(i0*j0*k0,      w(1,1,1,ismpl),1,x_solution,1)
          CALL daxpy(i0*j0*k0,1.0d0,a(1,1,1,ismpl),1,x_solution,1)
          CALL daxpy(i0*j0*k0,1.0d0,b(1,1,1,ismpl),1,x_solution,1)
          CALL daxpy(i0*j0*k0,1.0d0,c(1,1,1,ismpl),1,x_solution,1)
          CALL daxpy(i0*j0*k0,1.0d0,d(1,1,1,ismpl),1,x_solution,1)
          CALL daxpy(i0*j0*k0,1.0d0,e(1,1,1,ismpl),1,x_solution,1)
          CALL daxpy(i0*j0*k0,1.0d0,f(1,1,1,ismpl),1,x_solution,1)
          CALL daxpy(i0*j0*k0,1.0d0,g(1,1,1,ismpl),1,x_solution,1)
        END IF
        RETURN
      END

!>    @brief dummy routine for time meassurement
!>    @param[out] tick time
      SUBROUTINE cpu_time(tick)
        use mod_genrl
        IMPLICIT NONE
        DOUBLE PRECISION tick
        INTRINSIC dble
        tick = dble(I0*J0*K0)
        RETURN
      END
