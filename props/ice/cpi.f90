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


!>    @brief calculate ice isobaric heat capacity [J/(kg*K)]
!>    @param[in] i cell index, direction I0
!>    @param[in] j cell index, direction J0
!>    @param[in] k cell index, direction K0
!>    @param[in] ismpl local sample index
!>    @details
!>    After:
!>    Ling, F. & Zhang (2004):
!>    A Numerical Modelfor surface energy balance and the thermal 
!>    regime of the active layer and permafrost containing 
!>    unfrozen water or brine
!>    Cold Regions Science & Technology, 38, 1-15
!>    Ling and Zhang cite Osterkamp, T. E. (1987):
!>    Freezing and Thawing of Soils and Permafrost
!>    Containing Unfrozen Water or Brine
!>    Water Resources Research, Vol. 23, No. 12, pages 2279-2285
!>    Osterkamp cites Dorsey, N. E. (1940):
!>    Properties of Ordinary Water-Substance, Reinhold, New York
!>    Dorsey cites Dickinson, H., & Osborne, N. (1915):
!>    The specific heat and heat of fusion of ice. Journal of
!>    the Washington Academy of Sciences, 5(10), 338-340.
!>
!>    Alternative: 
!>    Fukusako, S.:
!>    Thermophysical Properties of Ice, Snow,and Sea Ice
!>    International Journal of Thermophysics, 1990, 11, 353-372
      DOUBLE PRECISION FUNCTION cpi(i,j,k,ismpl)


        use arrays
        IMPLICIT NONE

        INTEGER i, j, k, ui, ismpl
        DOUBLE PRECISION plocal, tlocal

        tlocal = temp(i,j,k,ismpl)
        IF (tlocal>0.D0) tlocal = 0.D0

!vr (Fukosako 1990)        cpi = 185.D0 + 6.89D0*(273.15D0+tlocal)
        cpi=2110.D0+7.7D0*tlocal

        RETURN
      END
