!        Generated by TAPENADE     (INRIA, Ecuador team)
!  tapenade 3.x
!
!  Differentiation of lz in forward (tangent) mode:
!   variations   of useful results: lz
!   with respect to varying inputs: *temp *propunit
!   Plus diff mem management of: temp:in propunit:in
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
!> @brief calculates effective thermal conductivity of the cell
!> @param[in] i cell index, direction I0
!> @param[in] j cell index, direction J0
!> @param[in] k cell index, direction K0
!> @param[in] ismpl local sample index
!> @return  thermal conductivity                lz[W/(m*K)]
!> @details
!> calculates effective thermal conductivity of the two phase system
!> matrix-porosity, z-direction.\n\n
!>
!> input:\n
!> porosity                            porlocal [-]\n
!> temperature                         tlocal in [degC]\n
DOUBLE PRECISION FUNCTION g_LZ(i, j, k, ismpl, lz)
  USE ARRAYS

  USE g_ARRAYS

  USE MOD_TEMP
  IMPLICIT NONE
! Location indices
  INTEGER, INTENT(IN) :: i
  INTEGER, INTENT(IN) :: j
  INTEGER, INTENT(IN) :: k
! Sample index
  INTEGER :: ismpl
! Local uindex
  INTEGER :: ui
! Local temperature [degC]
  DOUBLE PRECISION :: tlocal
  DOUBLE PRECISION :: g_tlocal
! Local porosity [-]
  DOUBLE PRECISION :: porlocal
  DOUBLE PRECISION :: g_porlocal
! Reference matrix thermal conductivity [W/(m*K)]
  DOUBLE PRECISION :: lammref
  DOUBLE PRECISION :: g_lammref
! Local fluid thermal conductivity [W/(m*K)]
  DOUBLE PRECISION :: lamfluid
  DOUBLE PRECISION :: g_lamfluid
  DOUBLE PRECISION, EXTERNAL :: LAMF
  DOUBLE PRECISION, EXTERNAL :: g_LAMF
! Local matrix thermal conductivity  [W/(m*K)]
  DOUBLE PRECISION, EXTERNAL :: LAMM
  DOUBLE PRECISION, EXTERNAL :: g_LAMM
  DOUBLE PRECISION :: pwy1
  DOUBLE PRECISION :: g_pwy1
  DOUBLE PRECISION :: pwr1
  DOUBLE PRECISION :: g_pwr1
  DOUBLE PRECISION :: pwr2
  DOUBLE PRECISION :: g_pwr2
  DOUBLE PRECISION :: temp0
  DOUBLE PRECISION :: lz
! Local Temperature in degC
  g_tlocal = g_temp(i, j, k, ismpl)
  tlocal = temp(i, j, k, ismpl)
! Local fluid thermal conductivity [W/(m*K)]
  g_lamfluid = g_LAMF(i, j, k, ismpl, lamfluid)
! Local unit index
  ui = uindex(i, j, k)
! Local porosity
  g_porlocal = g_propunit(ui, idx_por, ismpl)
  porlocal = propunit(ui, idx_por, ismpl)
! Reference matrix thermal conductivity [W/(m*K)]
  g_lammref = g_propunit(ui, idx_lz, ismpl)
  lammref = propunit(ui, idx_lz, ismpl)
! Local matrix thermal conductivity  [W/(m*K)]
  g_lz = g_LAMM(lammref, g_lammref, tlocal, g_tlocal, tref, &
&   ismpl, lz)
  IF (lz .LE. 0.d0 .OR. lamfluid .LE. 0.d0) THEN
    WRITE(*, *) 'Error: "lz" computes bad math !', lz, lamfluid, tlocal
    STOP
  ELSE
    g_pwy1 = -g_porlocal
    pwy1 = 1.d0 - porlocal
    temp0 = lz**pwy1
    IF (lz .LE. 0.0 .AND. (pwy1 .EQ. 0.0 .OR. pwy1 .NE. INT(pwy1))) THEN
      g_pwr1 = 0.D0
    ELSE IF (lz .LE. 0.0) THEN
      g_pwr1 = pwy1*lz**(pwy1-1)*g_lz
    ELSE
      g_pwr1 = pwy1*lz**(pwy1-1)*g_lz + temp0*LOG(lz)*g_pwy1
    END IF
    pwr1 = temp0
    temp0 = lamfluid**porlocal
    IF (lamfluid .LE. 0.0 .AND. (porlocal .EQ. 0.0 .OR. porlocal .NE. &
&       INT(porlocal))) THEN
      g_pwr2 = 0.D0
    ELSE IF (lamfluid .LE. 0.0) THEN
      g_pwr2 = porlocal*lamfluid**(porlocal-1)*g_lamfluid
    ELSE
      g_pwr2 = porlocal*lamfluid**(porlocal-1)*g_lamfluid + temp0*&
&       LOG(lamfluid)*g_porlocal
    END IF
    pwr2 = temp0
    g_lz = pwr2*g_pwr1 + pwr1*g_pwr2
    lz = pwr1*pwr2
    RETURN
  END IF
END FUNCTION g_LZ

