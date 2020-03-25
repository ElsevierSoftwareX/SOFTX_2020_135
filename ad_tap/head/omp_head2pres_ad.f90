!        Generated by TAPENADE     (INRIA, Ecuador team)
!  tapenade 3.x
!
!  Differentiation of omp_head2pres in reverse (adjoint) mode (with options noISIZE i8):
!   gradient     of useful results: *temp *head *dbc_data *bcperiod
!                *propunit *pres
!   with respect to varying inputs: *temp *head *dbc_data *bcperiod
!                *propunit *pres
!   Plus diff mem management of: temp:in head:in dbc_data:in bcperiod:in
!                propunit:in pres:in simtime:in
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
!>    @brief pressure in Pa from hydraulic potential und heigth above hz=0.0d0
!>    @param[in] init flag: 0-init, 1-normal setup
!>    @param[in] ismpl local sample index
!>    @details
!>    pressure in MPa from hydraulic potential und height above hz=0.0 \n
!>    p(surface) = 0.1 MPa \n
!>    OUTPUT in Pa\n
SUBROUTINE OMP_HEAD2PRES_AD(init, ismpl)
  use arrays

  USE ARRAYS_AD

  USE MOD_GENRL
  USE MOD_GENRLC
  USE MOD_FLOW
  IMPLICIT NONE
  INTEGER :: ismpl
  INTEGER :: i, j, k
  INTEGER :: init
  DOUBLE PRECISION :: psurf
  DOUBLE PRECISION :: dif
  DOUBLE PRECISION :: dif_ad
  DOUBLE PRECISION :: zero
  PARAMETER (zero=0.0d0)
  INTEGER :: branch
!     presetting for dirichlet boundary conditions to avoid site effects
  IF (ALLOCATED(head)) THEN
    CALL PUSHREAL8ARRAY(head, SIZE(head, 1)*SIZE(head, 2)*SIZE(head, 3)*&
&                 SIZE(head, 4))
    CALL PUSHCONTROL1B(1)
  ELSE
    CALL PUSHCONTROL1B(0)
  END IF
  CALL SET_DHBC(ismpl)
!
!
  IF (init .EQ. 0) THEN
!         initialize
    DO k=1,k0
      DO j=1,j0
        DO i=1,i0
          dif = head(i, j, k, ismpl) - delza(k)
          IF (dif .GT. zero) THEN
            CALL PUSHCONTROL1B(1)
          ELSE
            CALL PUSHCONTROL1B(0)
          END IF
        END DO
      END DO
    END DO
    DO k=k0,1,-1
      DO j=j0,1,-1
        DO i=i0,1,-1
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            dif_ad = 0.D0
          ELSE
            dif_ad = rref*grav*pres_ad(i, j, k, ismpl)
            pres_ad(i, j, k, ismpl) = 0.D0
          END IF
          pres_ad(i, j, k, ismpl) = 0.D0
          head_ad(i, j, k, ismpl) = head_ad(i, j, k, ismpl) + dif_ad
        END DO
      END DO
    END DO
  ELSE
!        pres and temp have already appropriate values
    DO k=1,k0
      DO j=1,j0
        DO i=1,i0
          dif = head(i, j, k, ismpl) - delza(k)
! jbr: pres = (h-z)*rref*grav
          IF (dif .GT. zero) THEN
            CALL PUSHCONTROL1B(1)
          ELSE
            CALL PUSHCONTROL1B(0)
          END IF
        END DO
      END DO
    END DO
    DO k=k0,1,-1
      DO j=j0,1,-1
        DO i=i0,1,-1
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            dif_ad = 0.D0
          ELSE
            dif_ad = rref*grav*pres_ad(i, j, k, ismpl)
            pres_ad(i, j, k, ismpl) = 0.D0
          END IF
          pres_ad(i, j, k, ismpl) = 0.D0
          head_ad(i, j, k, ismpl) = head_ad(i, j, k, ismpl) + dif_ad
        END DO
      END DO
    END DO
  END IF
  CALL SET_DTBC_AD(ismpl)
  CALL POPCONTROL1B(branch)
  IF (branch .EQ. 1) CALL POPREAL8ARRAY(head, SIZE(head, 1)*SIZE(head, 2&
&                                 )*SIZE(head, 3)*SIZE(head, 4))
  CALL SET_DHBC_AD(ismpl)
END SUBROUTINE OMP_HEAD2PRES_AD
