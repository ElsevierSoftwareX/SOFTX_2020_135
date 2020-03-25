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

!>    @brief calculates the bulk storativity
!>    @param[in] i cell index, direction I0
!>    @param[in] j cell index, direction J0
!>    @param[in] k cell index, direction K0
!>    @param[in] ismpl local sample index
!>    @return bulk storativity
!>    @details
!>    storb(i,j,k,ismpl) calculates the bulk storativity \n
!>    at node(i,j,k).\n
      DOUBLE PRECISION FUNCTION hstor(i,j,k,ismpl)
        use arrays
        use mod_flow
        IMPLICIT NONE

        INTEGER i, j, k, ismpl
        DOUBLE PRECISION rhof, compm, compf, por
        EXTERNAL rhof, compm, compf, por

        hstor = grav*rhof(i,j,k,ismpl)*(compm(i,j,k,ismpl)+por(i,j,k, &
          ismpl)*compf(i,j,k,ismpl))
        RETURN
      END
