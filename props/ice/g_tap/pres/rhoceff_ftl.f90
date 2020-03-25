!        Generated by TAPENADE     (INRIA, Ecuador team)
!  tapenade 3.x
!
!  Differentiation of rhoceff in forward (tangent) mode:
!   variations   of useful results: rhoceff
!   with respect to varying inputs: *temp *propunit *pres
!   Plus diff mem management of: temp:in propunit:in pres:in
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
!>    @brief calculates effective thermal conductivity of the two phase
!>    @param[in] i cell index, direction I0
!>    @param[in] j cell index, direction J0
!>    @param[in] k cell index, direction K0
!>    @param[in] ismpl local sample index
!>    @return  thermal conductivity                lz[W/(m*K)]
!>    @details
!>    calculates effective thermal conductivity of the two phase\n
!>    system matrix-porosity .\n
!>    input:\n
!>      porosity                            porlocal [-]\n
!>      pressure                            plocal [Mpa]\n
!>      temperature                         tlocal in [C]\n
DOUBLE PRECISION FUNCTION g_RHOCEFF(i, j, k, ismpl, rhoceff)
  USE ARRAYS

  USE g_ARRAYS

  USE ICE
  USE MOD_TEMP
  IMPLICIT NONE
  INTEGER :: i, j, k, ui, ismpl
  DOUBLE PRECISION :: tlocal, rcsolid, rcfluid, rcice, porlocal, fm, fi&
& , ff
  DOUBLE PRECISION :: g_tlocal, g_rcsolid, g_rcfluid, g_rcice, &
& g_porlocal, g_fm, g_fi, g_ff
  DOUBLE PRECISION :: t0, theta, dtheta, w0
  DOUBLE PRECISION :: g_theta, g_dtheta
  EXTERNAL RHOCF, RHOCM, RHOCI, POR, RHOI, &
&     RHOF
  EXTERNAL g_RHOCF, g_RHOCM, g_RHOCI, g_POR, g_RHOF
  DOUBLE PRECISION :: RHOCM, RHOCF, RHOCI, &
& POR, RHOI, RHOF
  DOUBLE PRECISION :: g_RHOCM, g_RHOCF, g_RHOCI, g_POR, g_RHOF
  INTRINSIC ABS
  DOUBLE PRECISION :: abs0
  DOUBLE PRECISION :: result1
  DOUBLE PRECISION :: g_result1
  DOUBLE PRECISION :: rhoceff
  ui = uindex(i, j, k)
  g_tlocal = g_temp(i, j, k, ismpl)
  tlocal = temp(i, j, k, ismpl)
  t0 = liq(i, j, k)
  IF (liq(i, j, k) - sol(i, j, k) .GE. 0.) THEN
    abs0 = liq(i, j, k) - sol(i, j, k)
  ELSE
    abs0 = -(liq(i, j, k)-sol(i, j, k))
  END IF
  w0 = abs0/2.d0
  CALL g_FTHETA(tlocal, g_tlocal, t0, w0, theta, g_theta, dtheta, &
&           g_dtheta, ismpl)
  g_rcfluid = g_RHOCF(i, j, k, ismpl, rcfluid)
  g_rcice = g_RHOCI(i, j, k, ismpl, rcice)
  g_rcsolid = g_RHOCM(i, j, k, ismpl, rcsolid)
  g_porlocal = g_POR(i, j, k, ismpl, porlocal)
  g_fm = -g_porlocal
  fm = 1.d0 - porlocal
  g_ff = theta*g_porlocal + porlocal*g_theta
  ff = porlocal*theta
  g_fi = g_porlocal - g_ff
  fi = porlocal - ff
!DM   Korrektur, 2008/02/21
  g_result1 = g_RHOF(i, j, k, ismpl, result1)
  g_rhoceff = fm*g_rcsolid + rcsolid*g_fm + ff*g_rcfluid + &
&   rcfluid*g_ff + fi*g_rcice + rcice*g_fi + lth*(dtheta*(porlocal&
&   *g_result1+result1*g_porlocal)+result1*porlocal*g_dtheta)
  rhoceff = rcsolid*fm + rcfluid*ff + rcice*fi + result1*porlocal*lth*&
&   dtheta
  RETURN
END FUNCTION g_RHOCEFF

