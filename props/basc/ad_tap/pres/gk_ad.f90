!        Generated by TAPENADE     (INRIA, Ecuador team)
!  tapenade 3.x
!
!  Differentiation of gk in reverse (adjoint) mode (with options noISIZE i8):
!   gradient     of useful results: *temp *propunit *tsal *pres
!                gk
!   with respect to varying inputs: *temp *propunit *tsal *pres
!   Plus diff mem management of: temp:in propunit:in tsal:in pres:in
!>    @brief average  conductivities on cell faces in z direction
!>    @param[in] i grid indices
!>    @param[in] j grid indices
!>    @param[in] k grid indices
!>    @param[in] ismpl local sample index
!>    @return z  conductivity (m/(Pa s))
SUBROUTINE GK_AD0(i, j, k, ismpl, gk_ad)
  use arrays

  USE ARRAYS_AD

  USE MOD_GENRL
  USE MOD_FLOW
  IMPLICIT NONE
  INTEGER :: ismpl
  INTEGER :: i, j, k
  EXTERNAL KZ, RHOF, VISF
  EXTERNAL KZ_AD, VISF_AD
  DOUBLE PRECISION :: f1, f2, prod, summ, KZ, RHOF, VISF
  DOUBLE PRECISION :: f1_ad, f2_ad, prod_ad, summ_ad
  DOUBLE PRECISION :: result1
  DOUBLE PRECISION :: result1_ad
  DOUBLE PRECISION :: result2
  DOUBLE PRECISION :: result2_ad
  INTEGER :: arg1
  INTEGER :: arg2
  DOUBLE PRECISION :: temporary_ad
  DOUBLE PRECISION :: gk
  DOUBLE PRECISION :: gk_ad
  IF (k0 .GT. 1 .AND. k .LT. k0) THEN
    result1 = KZ(i, j, k, ismpl)
    result2 = VISF(i, j, k, ismpl)
    f1 = result1/result2
    arg1 = k + 1
    CALL PUSHREAL8(result1)
    result1 = KZ(i, j, arg1, ismpl)
    arg2 = k + 1
    CALL PUSHREAL8(result2)
    result2 = VISF(i, j, arg2, ismpl)
    f2 = result1/result2
    prod = f1*f2
    summ = f1*delz(k+1) + f2*delz(k)
    IF (summ .GT. 0.d0) THEN
      temporary_ad = 2.d0*gk_ad/summ
      prod_ad = temporary_ad
      summ_ad = -(prod*temporary_ad/summ)
    ELSE
      prod_ad = 0.D0
      summ_ad = 0.D0
    END IF
    f1_ad = delz(k+1)*summ_ad + f2*prod_ad
    f2_ad = delz(k)*summ_ad + f1*prod_ad
    result1_ad = f2_ad/result2
    result2_ad = -(result1*f2_ad/result2**2)
    CALL POPREAL8(result2)
    CALL VISF_AD(i, j, arg2, ismpl, result2_ad)
    CALL POPREAL8(result1)
    CALL KZ_AD(i, j, arg1, ismpl, result1_ad)
    result1_ad = f1_ad/result2
    result2_ad = -(result1*f1_ad/result2**2)
    CALL VISF_AD(i, j, k, ismpl, result2_ad)
    CALL KZ_AD(i, j, k, ismpl, result1_ad)
  END IF
END SUBROUTINE GK_AD0

