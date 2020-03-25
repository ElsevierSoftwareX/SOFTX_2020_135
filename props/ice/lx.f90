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
!>    @return  thermal conductivity                lx[W/(m*K)]
!>    @details
!>    calculates effective thermal conductivity of the two phase\n
!>    system matrix-porosity .\n
!>    input:\n
!>      porosity                            porlocal [-]\n
!>      pressure                            plocal [Mpa]\n
!>      temperature                         tlocal in [C]\n
      DOUBLE PRECISION FUNCTION lx(i,j,k,ismpl)


        use arrays
                use ice
        use mod_temp
        IMPLICIT NONE


        INTEGER i, j, k, ui, ismpl
        DOUBLE PRECISION tlocal, lsolid, lfluid, porlocal, fm, fi, ff
        DOUBLE PRECISION t0, theta, dtheta, w0, lice
        DOUBLE PRECISION lamm, lamf, lami, por
        EXTERNAL lamf, lamm, lami, por

        ui = uindex(i,j,k)
        tlocal = temp(i,j,k,ismpl)

        t0 = liq(i,j,k)
        w0 = abs(liq(i,j,k)-sol(i,j,k))/2.D0
        CALL ftheta(tlocal,t0,w0,theta,dtheta,ismpl)

        lfluid = lamf(i,j,k,ismpl)
        lice = lami(i,j,k,ismpl)

        lsolid = propunit(ui,idx_lz,ismpl)*propunit(ui,idx_an_lx,ismpl &
          )
        lsolid = lamm(lsolid,tlocal,tref,ismpl)

        porlocal = por(i,j,k,ismpl)
        fm = 1.D0 - porlocal
        ff = porlocal*theta
        fi = porlocal - ff



        lx = lsolid**fm*lfluid**ff*lice**fi

        RETURN
      END
