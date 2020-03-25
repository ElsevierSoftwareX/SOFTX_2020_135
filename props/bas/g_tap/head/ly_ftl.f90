!        Generated by TAPENADE     (INRIA, Ecuador team)
!  tapenade 3.x
!
!  Differentiation of ly in forward (tangent) mode:
!   variations   of useful results: ly
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
!> @return  thermal conductivity                ly[W/(m*K)]
!> @details
!> calculates effective thermal conductivity of the two phase system
!> matrix-porosity, y-direction.\n\n
!>
!> input:\n
!> porosity                            porlocal [-]\n
!> temperature                         tlocal in [degC]\n
DOUBLE PRECISION FUNCTION g_LY(i, j, k, ismpl, ly)
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
  DOUBLE PRECISION :: ly
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
  g_lammref = propunit(ui, idx_an_ly, ismpl)*g_propunit(ui, idx_lz, &
&   ismpl) + propunit(ui, idx_lz, ismpl)*g_propunit(ui, idx_an_ly, &
&   ismpl)
  lammref = propunit(ui, idx_lz, ismpl)*propunit(ui, idx_an_ly, ismpl)
! Local matrix thermal conductivity  [W/(m*K)]
  g_ly = g_LAMM(lammref, g_lammref, tlocal, g_tlocal, tref, &
&   ismpl, ly)
  IF (ly .LE. 0.d0 .OR. lamfluid .LE. 0.d0) THEN
    WRITE(*, *) 'Error: "ly" computes bad math !', ly, lamfluid, tlocal
    STOP
  ELSE
    g_pwy1 = -g_porlocal
    pwy1 = 1.d0 - porlocal
    temp0 = ly**pwy1
    IF (ly .LE. 0.0 .AND. (pwy1 .EQ. 0.0 .OR. pwy1 .NE. INT(pwy1))) THEN
      g_pwr1 = 0.D0
    ELSE IF (ly .LE. 0.0) THEN
      g_pwr1 = pwy1*ly**(pwy1-1)*g_ly
    ELSE
      g_pwr1 = pwy1*ly**(pwy1-1)*g_ly + temp0*LOG(ly)*g_pwy1
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
    g_ly = pwr2*g_pwr1 + pwr1*g_pwr2
    ly = pwr1*pwr2
    RETURN
  END IF
END FUNCTION g_LY

