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
      DOUBLE PRECISION FUNCTION rhoceff(i,j,k,ismpl)


        use arrays
                use ice
        use mod_temp
        IMPLICIT NONE


        INTEGER i, j, k, ui, ismpl
        DOUBLE PRECISION tlocal, rcsolid, rcfluid, rcice, porlocal, &
          fm, fi, ff
        DOUBLE PRECISION t0, theta, dtheta, w0
        DOUBLE PRECISION rhocm, rhocf, rhoci, por, rhoi, rhof
        EXTERNAL rhocf, rhocm, rhoci, por, rhoi, rhof

        ui = uindex(i,j,k)
        tlocal = temp(i,j,k,ismpl)
        t0 = liq(i,j,k)
        w0 = abs(liq(i,j,k)-sol(i,j,k))/2.D0
        CALL ftheta(tlocal,t0,w0,theta,dtheta,ismpl)

        rcfluid = rhocf(i,j,k,ismpl)
        rcice = rhoci(i,j,k,ismpl)
        rcsolid = rhocm(i,j,k,ismpl)

        porlocal = por(i,j,k,ismpl)

        fm = 1.D0 - porlocal
        ff = porlocal*theta
        fi = porlocal - ff

!DM   Korrektur, 2008/02/21
        rhoceff = rcsolid*fm + rcfluid*ff + rcice*fi + &
          rhof(i,j,k,ismpl)*porlocal*lth*dtheta

        RETURN
      END
