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

!>    @brief calculate temperature dependent thermal conductivity
!>    @param[in,out] solid thermal conductivity at reference temperature given in input file [W / (m K)]
!>    @param[in,out] tlocal temperature at this sample index [degree Celsius]
!>    @param[in,out] tref reference temperature [degree Celsius]
!>    @param[in] ismpl local sample index
!>    @return thermal conductivity [W / (m K)]
!>    @details
!>    calculate temperature dependent thermal conductivity of the stony matrix\n
!>    (zoth & haenel, 1988)\n
      DOUBLE PRECISION FUNCTION lamm(solid,tlocal,tref,ismpl)


        IMPLICIT NONE

        INTEGER ismpl
        DOUBLE PRECISION solid, tlocal, tref, tlimit
        DOUBLE PRECISION cddz, cddz0, cgt0, wgt, cgt
        PARAMETER (tlimit=800.D0)

        IF (tlocal>tlimit) THEN
          lamm = 770.0D0/(350.0D0+tlocal) + 0.7D0
        ELSE
          cddz = 770.0D0/(350.0D0+tlocal) + 0.7D0
          cddz0 = 770.0D0/(350.0D0+tref) + 0.7D0
          cgt0 = solid/cddz0
          wgt = (tlocal-tref)/(tlimit-tref)
          cgt = cgt0 - (cgt0-1.0D0)*wgt
          lamm = cgt*cddz
        END IF

        RETURN
      END
