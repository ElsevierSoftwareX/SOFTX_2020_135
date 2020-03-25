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
!>      pressure                            plocal [pa]\n
!>      temperature                         tlocal in [C]\n
      DOUBLE PRECISION FUNCTION lz(i,j,k,ismpl)
        use arrays
        use mod_temp
        IMPLICIT NONE


        INTEGER i, j, k, ui, ismpl
        DOUBLE PRECISION plocal, tlocal, solid, fluid, lamunit, lamf, &
          porlocal, tkelvin, lamm
        EXTERNAL lamf, lamm


!      ploCal = pres(i,j,k,ismpl)*Pa_Conv1
        tlocal = temp(i,j,k,ismpl)
        fluid = lamf(i,j,k,ismpl)
        ui = uindex(i,j,k)
        porlocal = propunit(ui,idx_por,ismpl)
        lamunit = propunit(ui,idx_lz,ismpl)
!         lz=
!     *    (1.d0-porlocal)*lamm(lamunit,tlocal,tref,ismpl)+porlocal*fluid
        lz = lamm(lamunit,tlocal,tref,ismpl)
        IF (lz<=0.D0 .OR. fluid<=0.D0) THEN
          WRITE(*,*) 'warning: "lz" computes bad math !'
        ELSE
          lz = lz**(1.D0-porlocal)*fluid**porlocal
        END IF

        RETURN
      END
