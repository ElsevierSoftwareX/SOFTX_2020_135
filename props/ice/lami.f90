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

!>    @brief calculate ice thermal conductivity [W/(m*K)]
!>    @param[in] i cell index, direction I0
!>    @param[in] j cell index, direction J0
!>    @param[in] k cell index, direction K0
!>    @param[in] ismpl local sample index
!#>   return thermal conductivity of ice [W/(m*K)]
!>    @details
!>    Ling, F. & Zhang, T. (2004):
!>    A Numerical Modelfor surface energy balance and the thermal 
!>    regime of the active layer and permafrost containing 
!>    unfrozen water or brine
!>    Cold Regions Science & Technology, 38, 1-15
!>    Ling and Zhang cite Osterkamp, T. E. (1987):
!>    Freezing and Thawing of Soils and Permafrost
!>    Containing Unfrozen Water or Brine
!>    Water Resources Research, Vol. 23, No. 12, pages 2279-2285
!>    Alternative: 
!>    Fukusako, S.:
!>    Thermophysical Properties of Ice, Snow,and Sea Ice
!>    International Journal of Thermophysics, 1990, 11, 353-372
      DOUBLE PRECISION FUNCTION lami(i,j,k,ismpl)


        use arrays
        IMPLICIT NONE

        INTEGER i, j, k, ismpl
        DOUBLE PRECISION plocal, tlocal

        tlocal = temp(i,j,k,ismpl)
        IF (tlocal>0.D0) tlocal = 0.D0
        lami = 0.4685D0 + 488.19D0/(tlocal+273.16D0)
        
! vr (Fukosako1990) 
! vr lami=1.16d0*(1.91d0-8.66d-3*tlocal+2.97d-5*tlocal*tlocal);


        RETURN
      END
