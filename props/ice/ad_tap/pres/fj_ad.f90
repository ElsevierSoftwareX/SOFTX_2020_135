!        Generated by TAPENADE     (INRIA, Ecuador team)
!  tapenade 3.x
!
!  Differentiation of fj in reverse (adjoint) mode (with options noISIZE i8):
!   gradient     of useful results: *temp *propunit *pres fj
!   with respect to varying inputs: *temp *propunit *pres
!   Plus diff mem management of: temp:in propunit:in pres:in
!>    @brief average  conductivities on cell faces in y direction
!>    @param[in] i grid indices
!>    @param[in] j grid indices
!>    @param[in] k grid indices
!>    @param[in] ismpl local sample index
!>    @return y  conductivity (m/(Pa s))
SUBROUTINE FJ_AD(i, j, k, ismpl, fj_adv)
  use arrays

  USE ARRAYS_AD

  USE MOD_GENRL
  USE MOD_FLOW
  IMPLICIT NONE
  double precision :: fj_adv
  INTEGER :: ismpl
  INTEGER :: i, j, k
  EXTERNAL KY, RHOF, VISF
  EXTERNAL KY_AD, VISF_AD
  DOUBLE PRECISION :: f1, f2, prod, summ, KY, RHOF, VISF
  DOUBLE PRECISION :: f1_ad, f2_ad, prod_ad, summ_ad
  DOUBLE PRECISION :: result1
  DOUBLE PRECISION :: result1_ad
  DOUBLE PRECISION :: result2
  DOUBLE PRECISION :: result2_ad
  INTEGER :: arg1
  INTEGER :: arg2
  DOUBLE PRECISION :: temporary_ad
  DOUBLE PRECISION :: fj
  IF (j0 .GT. 1 .AND. j .LT. j0) THEN
    result1 = KY(i, j, k, ismpl)
    result2 = VISF(i, j, k, ismpl)
    f1 = result1/result2
    arg1 = j + 1
    CALL PUSHREAL8(result1)
    result1 = KY(i, arg1, k, ismpl)
    arg2 = j + 1
    CALL PUSHREAL8(result2)
    result2 = VISF(i, arg2, k, ismpl)
    f2 = result1/result2
    prod = f1*f2
    summ = f1*dely(j+1) + f2*dely(j)
    IF (summ .GT. 0.d0) THEN
      temporary_ad = 2.d0*fj_adv/summ
      prod_ad = temporary_ad
      summ_ad = -(prod*temporary_ad/summ)
    ELSE
      prod_ad = 0.D0
      summ_ad = 0.D0
    END IF
    f1_ad = dely(j+1)*summ_ad + f2*prod_ad
    f2_ad = dely(j)*summ_ad + f1*prod_ad
    result1_ad = f2_ad/result2
    result2_ad = -(result1*f2_ad/result2**2)
    CALL POPREAL8(result2)
    CALL VISF_AD(i, arg2, k, ismpl, result2_ad)
    CALL POPREAL8(result1)
    CALL KY_AD(i, arg1, k, ismpl, result1_ad)
    result1_ad = f1_ad/result2
    result2_ad = -(result1*f1_ad/result2**2)
    CALL VISF_AD(i, j, k, ismpl, result2_ad)
    CALL KY_AD(i, j, k, ismpl, result1_ad)
  END IF
END SUBROUTINE FJ_AD

