!>    @brief AD wrapper routine to compute the derivative pressure (but not the original pressure itself!)
!>    @param[in] init flag: 0-init, 1-normal setup
!>    @param[in] ismpl local sample index
!>    @details
!>    parallelisation wrapper for pressure computation
SUBROUTINE g_head2pres(init, ismpl)
  USE arrays
  USE mod_genrl
  USE mod_linfos
  USE mod_OMP_TOOLS
  IMPLICIT NONE
  integer :: ismpl
  INCLUDE 'OMP_TOOLS.inc'
  INTEGER init

  IF (linfos(3) .GE. 2) WRITE(*, '(A,I1,A)') '  ... g_pressure (init=',init, ')'
#ifdef head_base
! save original pressure
  CALL dcopy(i0*j0*k0,pres(1,1,1,ismpl),1,x(1,1,1,ismpl),1)
!
#ifdef fOMP
!$OMP parallel num_threads(Tlevel_1)
!$      call omp_binding(ismpl)
#endif
!
  CALL g_omp_head2pres(init, ismpl)
!
#ifdef fOMP
!$OMP end parallel
#endif
!
! restore original pressure
  CALL dcopy(i0*j0*k0,x(1,1,1,ismpl),1,pres(1,1,1,ismpl),1)
#endif
  RETURN
END SUBROUTINE g_head2pres
