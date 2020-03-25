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


!>    @brief calculate the thermal conductivity lamf in W/(m*K) of the formation water
!>    @param[in] i cell index, direction I0
!>    @param[in] j cell index, direction J0
!>    @param[in] k cell index, direction K0
!>    @param[in] ismpl local sample index
!>    @return  thermal conductivity                lamf[W/(m*K)]
!>    @details
!> calculate the thermal conductivity lamf in W/(m*K) of\n
!>formation water, given temperature in C, and salinity in mass fraction\n
!>(g/g)of NaCl. Thermal conductivity of freshwater, lamfw is calculated using\n
!>the Phillips (1981) formulation. \n
!> Phillips, S. L. "A technical databook for geothermal energy utilization."(1981).
!>Range of validity:  20 to 330°C and up to 4 molal NaCl\n
!>    input:\n
!>      pressure                            plocal [Mpa]\n
!>      temperature                         tlocal in [C]\n
      DOUBLE PRECISION FUNCTION lamf(i,j,k,ismpl)


        use arrays
        IMPLICIT NONE

        INTEGER i, j, k, ismpl
        DOUBLE PRECISION tlocal, tr, tr2, tr3, tr4


        tlocal = temp(i,j,k,ismpl)
        IF (tlocal<0.D0) tlocal = 0.D0
        IF (tlocal>300.D0) tlocal = 300.D0

        tr = (tlocal+273.15D0)/273.15D0
        tr2 = tr*tr
        tr3 = tr2*tr
        tr4 = tr3*tr
        lamf = (-0.92247D0+2.8395D0*tr-1.8007D0*tr2+0.52577D0*tr3- &
          0.07344D0*tr4)

        RETURN
      END
