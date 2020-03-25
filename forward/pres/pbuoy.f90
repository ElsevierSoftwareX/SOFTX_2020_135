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

!>    @brief calculate buoyancy for pressure equation
!>    @param[in] i cell index, direction I0
!>    @param[in] j cell index, direction J0
!>    @param[in] k cell index, direction K0
!>    @param[in] ismpl local sample index
!>    @return buoyancy for pressure
!>    @details
!> calculate buoyancy for pressure equation\n
!> sign convention: negative for positive buoyancy\n
!> if porosity is 1.E-20, the pressure is set to\
!> atmospheric pressure      
      DOUBLE PRECISION FUNCTION buoy(i,j,k,ismpl)
        use arrays
        use mod_genrl
        use mod_flow
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k

        DOUBLE PRECISION rhav, hh, h0, h1, prod, summ
        DOUBLE PRECISION rhof, kz, visf, por, rh1, rh2
        EXTERNAL rhof, kz, visf, por

        buoy = 0.D0
        IF (por(i,j,k,ismpl).LT.1.E-19) THEN 
                rh1 = 1.29E0
        ELSE
                rh1 = rhof(i,j,k,ismpl)
        END IF
        IF (por(i,j,k+1,ismpl).LT.1.E-19) THEN
                 rh2 = 1.29E0
         ELSE
                rh2 = rhof(i,j,k+1,ismpl)
         END IF
        rhav = 0.5D0*(rh1+rh2)

        hh = 0.D0
        h0 = kz(i,j,k,ismpl)/visf(i,j,k,ismpl)
        h1 = kz(i,j,k+1,ismpl)/visf(i,j,k+1,ismpl)
        summ = h0 + h1
        prod = h0*h1
        IF (summ>0.D0) hh = 2.0D0*prod/summ

        buoy = hh*grav*rhav

        RETURN
      END

!>    @brief calculate buoyancy for pressure equation
!>    @param[in] i cell index, direction I0
!>    @param[in] j cell index, direction J0
!>    @param[in] k cell index, direction K0
!>    @param[in] ismpl local sample index
!>    @return buoyancy for pressure
!>    @details
!> calculate buoyancy for pressure equation\n
!> sign convention: negative for positive buoyancy\n
      DOUBLE PRECISION FUNCTION vbuoy(i,j,k,ismpl)
        use arrays
        use mod_genrl
        use mod_flow
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k

        DOUBLE PRECISION h0, h1, prod, summ
        DOUBLE PRECISION rhof, kz, visf, rhav, por
        EXTERNAL rhof, kz, visf, por

        rhav = 0.5D0*(rhof(i,j,k+1,ismpl)+rhof(i,j,k,ismpl))
        vbuoy = 0.D0
        h0 = kz(i,j,k,ismpl)*rhof(i,j,k,ismpl)/visf(i,j,k,ismpl)
        h1 = kz(i,j,k+1,ismpl)*rhof(i,j,k+1,ismpl)/visf(i,j,k+1,ismpl)
        summ = h0 + h1
        prod = h0*h1
        IF (summ>0.D0) vbuoy = 2.0D0*prod/summ
        IF (por(i,j,k,ismpl).LT.1.E-19) vbuoy = 1.29E-8/998.

        vbuoy = vbuoy*grav
!        IF (por(i,j,k,ismpl).LT.1.E-19) vbuoy = 0.

        RETURN
      END
