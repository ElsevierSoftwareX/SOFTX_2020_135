!        Generated by TAPENADE     (INRIA, Ecuador team)
!  tapenade 3.x
!
!  Differentiation of omp_set_tsal in reverse (adjoint) mode (with options noISIZE i8):
!   gradient     of useful results: *tsal *conc
!   with respect to varying inputs: *tsal *conc
!   Plus diff mem management of: tsal:in conc:in
!>    @brief calculate total salinity
!>    @param[in] ismpl local sample index
!>    @details
!>calculate total salinity\n
SUBROUTINE OMP_SET_TSAL_AD(ismpl)
  use arrays

  USE ARRAYS_AD

  USE MOD_CONC
  USE MOD_GENRL
  USE MOD_GENRLC
  USE MOD_LINFOS
  IMPLICIT NONE
  INTEGER :: ismpl
  INTEGER :: i, j, k
  INTEGER :: species
  DOUBLE PRECISION :: summ, fac
  DOUBLE PRECISION :: summ_ad
  DO k=k0,1,-1
    DO j=j0,1,-1
      DO i=i0,1,-1
        summ_ad = tsal_ad(i, j, k, ismpl)
        tsal_ad(i, j, k, ismpl) = 0.D0
        DO species=ntrans,1,-1
          fac = mmas_c(species)/mmas_nacl
          conc_ad(i, j, k, species, ismpl) = conc_ad(i, j, k, species, &
&           ismpl) + fac*summ_ad
        END DO
      END DO
    END DO
  END DO
END SUBROUTINE OMP_SET_TSAL_AD
