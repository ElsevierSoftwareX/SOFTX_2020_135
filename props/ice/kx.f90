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


!>    @brief assign permeability in x direction to cell
!>    @param[in] i grid indices
!>    @param[in] j grid indices
!>    @param[in] k grid indices
!>    @param[in] ismpl local sample index
!>    @return  permeability                        (m^2)
      DOUBLE PRECISION FUNCTION kx(i,j,k,ismpl)


        use arrays
                use ice
        use mod_temp
        IMPLICIT NONE

        INTEGER i, j, k, ismpl
        DOUBLE PRECISION t0, theta, dtheta, w0, tlocal, permi, ff, fi, &
          porlocal, por
        EXTERNAL por

        tlocal = temp(i,j,k,ismpl)
        t0 = liq(i,j,k)
        w0 = abs(liq(i,j,k)-sol(i,j,k))/2.D0
        CALL ftheta(tlocal,t0,w0,theta,dtheta,ismpl)
        porlocal = por(i,j,k,ismpl)
        ff = porlocal*theta
        fi = porlocal - ff
        permi = 10.D0**(-100.D0*fi/(porlocal+1.D-12))

        IF (permi<1.0D-6) THEN

          permi = 1.0D-6

        END IF

!       permeff=PERM(I,J,K)*permi
        kx = propunit(uindex(i,j,k),idx_kz,ismpl)* &
          propunit(uindex(i,j,k),idx_an_kx,ismpl)*permi

        RETURN
      END
