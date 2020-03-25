!        Generated by TAPENADE     (INRIA, Ecuador team)
!  tapenade 3.x
!
!  Differentiation of cpf in reverse (adjoint) mode (with options noISIZE i8):
!   gradient     of useful results: *temp *tsal *pres cpf
!   with respect to varying inputs: *temp *tsal *pres
!   Plus diff mem management of: temp:in tsal:in pres:in
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
!> @brief cpf(i,j,k,ismpl) calculates the isobaric heat capacity of the fluid [J/kg/K]
!> @param[in] i cell index, direction I0
!> @param[in] j cell index, direction J0
!> @param[in] k cell index, direction K0
!> @param[in] ismpl local sample index
!> @return cpf  [J/kg/K]
!> @details
!> cpf(i,j,k,ismpl) calculates the isobaric heat capacity in (in
!> J/kg/K) of brine, given temperature (t, in degC) pressure (p,in MPa),
!> and salinity (s, in mol/L) at node(i,j,k)\\n \n
!>
!> Source:\n
!> Driesner & Heinrich, The system H2O-NaCl. Part I: Correlation
!>    formulae for phase relations in temperature-pressure-composition
!>    space from 0 to 1000 C, 0 to 5000 bar, and 0 to 1 XNaCl
!>    Geochimica et Cosmochimica Acta 71 (2007) 4880-4901\n\n
!>
!> Driesner, The system H2O-NaCl. Part II: Correlations for molar
!>    volume, enthalpy, and isobaric heat capacity from 0 to 1000 C, 1
!>    to 5000 bar, and 0 to 1 XNaCl Geochimica et Cosmochimica Acta 71
!>    (2007) 4902-4919\n \n
SUBROUTINE CPF_AD(i, j, k, ismpl, cpf_adv)
  use arrays

  USE ARRAYS_AD

  USE MOD_FLOW
  IMPLICIT NONE
  double precision :: cpf_adv
! Location indices
  INTEGER, INTENT(IN) :: i
  INTEGER, INTENT(IN) :: j
  INTEGER, INTENT(IN) :: k
! Sample index
  INTEGER :: ismpl
! Local Temperature, [degC]
  DOUBLE PRECISION :: t
  DOUBLE PRECISION :: t_ad
! Local Pressure [MPa]
  DOUBLE PRECISION :: p
  DOUBLE PRECISION :: p_ad
! Local salinity [mol/kg]
  DOUBLE PRECISION :: s
  DOUBLE PRECISION :: s_ad
! Monomials of pressure [bar]
  DOUBLE PRECISION :: pb, pb2
  DOUBLE PRECISION :: pb_ad, pb2_ad
! Pure water compressibility [1/Pa]
  DOUBLE PRECISION :: cw
  DOUBLE PRECISION, EXTERNAL :: COMPW
! Pure water density [kg/m3]
  DOUBLE PRECISION :: rw
  DOUBLE PRECISION :: rw_ad
  DOUBLE PRECISION, EXTERNAL :: RHOW
! Pure water thermal conductivity [J/kg/K]
  DOUBLE PRECISION, EXTERNAL :: CPW
! Molar mass of NaCl [g/mol]
  DOUBLE PRECISION, PARAMETER :: mmnacl=58.44277d0
! Mass fraction of NaCl in solution [-]
  DOUBLE PRECISION :: fracnacl
  DOUBLE PRECISION :: fracnacl_ad
! Molar mass of water, H2O [g/mol]
  DOUBLE PRECISION, PARAMETER :: mmwater=18.01528d0
! Mole fraction of NaCl [-]
  DOUBLE PRECISION :: xnacl
  DOUBLE PRECISION :: xnacl_ad
! Temperature at which pure water has the same specific
! enthalpy as the solution [K]
  DOUBLE PRECISION :: tv
  DOUBLE PRECISION :: tv_ad
! Driesner2007: Coefficients
! Equation (23)
  DOUBLE PRECISION :: q1, q10, q11, q12
  DOUBLE PRECISION :: q1_ad, q10_ad, q11_ad, q12_ad
! Equation (24)
  DOUBLE PRECISION :: q2, q20, q21, q22, q23
  DOUBLE PRECISION :: q2_ad, q20_ad, q21_ad, q22_ad, q23_ad
! Equation (25), (26)
  DOUBLE PRECISION :: q1x, q2x
  DOUBLE PRECISION :: q1x_ad, q2x_ad
  INTRINSIC SQRT
  DOUBLE PRECISION :: result1
  DOUBLE PRECISION :: result1_ad
  DOUBLE PRECISION :: temporary_ad
  DOUBLE PRECISION :: temp0
  DOUBLE PRECISION :: temporary_ad0
  DOUBLE PRECISION :: cpf
! Local temperature [degC]
  t = temp(i, j, k, ismpl)
! Local pressure [MPa]
  p = pres(i, j, k, ismpl)*pa_conv1
! Local salinity [mol/L]
  s = tsal(i, j, k, ismpl)
  IF (s .LE. 0.0d0) THEN
    CALL CPW_AD0(p, p_ad, t, t_ad, cpf_adv)
    s_ad = 0.D0
  ELSE
! pure water density [g/L = kg/m3]
    rw = RHOW(p, t)
! Pressure monomials [bar]
    pb = p*10.0d0
    pb2 = pb*pb
! Mass fraction of NaCl, mol/L > (g/L) / (g/L) mass fraction
    fracnacl = s*mmnacl/(rw+s*mmnacl)
! Mole fraction of NaCl, (mol/L) / (mol/L)
    xnacl = fracnacl/mmnacl/(fracnacl/mmnacl+(1-fracnacl)/mmwater)
! Driesner2007 parameters, Table 5
    q11 = -32.1724d0 + 0.0621255d0*pb
    q21 = -1.69513d0 - 4.52781d-4*pb - 6.04279d-8*pb2
    q22 = 0.0612567d0 + 1.88082d-5*pb
! Driesner2007: q1 and q2 parameters for liquid NaCl
    q1x = 47.9048d0 - 9.36994d-3*pb + 6.51059d-6*pb2
    q2x = 0.241022d0 + 3.45087d-5*pb - 4.28356d-9*pb2
! Driesner2007: Set xnacl=1 in (23) => q1 = q1x
    q10 = q1x
! Driesner2007: Set xnacl=0 in (23) => q1 = 0.0d0
    q12 = -q10 - q11
! Driesner2007: Set xnacl=0 in (24) => q2 = 1.0d0
    q20 = 1.0d0 - q21*SQRT(q22)
! Driesner2007: Set xnacl=1 in (10) => q2 = q2x
    q23 = q2x - q20 - q21*SQRT(1.0d0+q22)
! Driesner2007: Equation (23)
    q1 = q10 + q11*(1.0d0-xnacl) + q12*(1.0d0-xnacl)**2
! Driesner2007: Equation (24)
    q2 = q20 + q21*SQRT(xnacl+q22) + q23*xnacl
! Temperature at which pure water has the same specific
! enthalpy as the solution
    tv = q1 + q2*t
! Isobaric heat capacity of solution
    result1 = CPW(p, tv)
    result1_ad = q2*cpf_adv
    CALL CPW_AD0(p, p_ad, tv, tv_ad, result1_ad)
    q2_ad = result1*cpf_adv + t*tv_ad
    q1_ad = tv_ad
    t_ad = q2*tv_ad
    temp0 = SQRT(xnacl + q22)
    q21_ad = temp0*q2_ad
    IF (xnacl + q22 .EQ. 0.0) THEN
      temporary_ad = 0.D0
    ELSE
      temporary_ad = q21*q2_ad/(2.0*temp0)
    END IF
    q23_ad = xnacl*q2_ad
    q20_ad = q2_ad - q23_ad
    xnacl_ad = q23*q2_ad + temporary_ad - (q11+2*(1.0d0-xnacl)*q12)*q1_ad
    q12_ad = (1.0d0-xnacl)**2*q1_ad
    q10_ad = q1_ad - q12_ad
    q11_ad = (1.0d0-xnacl)*q1_ad - q12_ad
    temp0 = SQRT(q22 + 1.0d0)
    IF (q22 + 1.0d0 .EQ. 0.0) THEN
      q22_ad = temporary_ad
    ELSE
      q22_ad = temporary_ad - q21*q23_ad/(2.0*temp0)
    END IF
    q2x_ad = q23_ad
    q21_ad = q21_ad - temp0*q23_ad
    temp0 = SQRT(q22)
    q21_ad = q21_ad - temp0*q20_ad
    IF (.NOT.q22 .EQ. 0.0) q22_ad = q22_ad - q21*q20_ad/(2.0*temp0)
    q1x_ad = q10_ad
    pb2_ad = 6.51059d-6*q1x_ad - 4.28356d-9*q2x_ad - 6.04279d-8*q21_ad
    pb_ad = 3.45087d-5*q2x_ad + 1.88082d-5*q22_ad - 9.36994d-3*q1x_ad + &
&     0.0621255d0*q11_ad - 4.52781d-4*q21_ad + 2*pb*pb2_ad
    temp0 = mmnacl*(fracnacl/mmnacl+(-fracnacl+1)/mmwater)
    temporary_ad = -(mmnacl*fracnacl*xnacl_ad/temp0**2)
    fracnacl_ad = xnacl_ad/temp0 + (1.0/mmnacl-1.0/mmwater)*temporary_ad
    temporary_ad = mmnacl*fracnacl_ad/(rw+mmnacl*s)
    temporary_ad0 = -(s*temporary_ad/(rw+mmnacl*s))
    s_ad = temporary_ad + mmnacl*temporary_ad0
    rw_ad = temporary_ad0
    p_ad = p_ad + 10.0d0*pb_ad
    CALL RHOW_AD0(p, p_ad, t, t_ad, rw_ad)
  END IF
  tsal_ad(i, j, k, ismpl) = tsal_ad(i, j, k, ismpl) + s_ad
  pres_ad(i, j, k, ismpl) = pres_ad(i, j, k, ismpl) + pa_conv1*p_ad
  temp_ad(i, j, k, ismpl) = temp_ad(i, j, k, ismpl) + t_ad
END SUBROUTINE CPF_AD

