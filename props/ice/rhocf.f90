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


!>    @brief calculates heat capacity times density of water.
!>    @param[in] i cell index, direction I0
!>    @param[in] j cell index, direction J0
!>    @param[in] k cell index, direction K0
!>    @param[in] ismpl local sample index
!>    @return  rhocf                                rhocf[J / (m^3 * K)]
!>    @details
!>    input:\n
!>      heat capacity of water              cpf [J / kg / K]\n
!>      density                             rhof [kg/m^3]\n
!>      pressure                            plocal [Pa]\n
!>      temperature                         tlocal [C]\n
      DOUBLE PRECISION FUNCTION rhocf(i,j,k,ismpl)


        use arrays
        IMPLICIT NONE

        INTEGER i, j, k, ismpl
        DOUBLE PRECISION plocal, tlocal, rfluid, cfluid, cpf, rhof
        EXTERNAL cpf, rhof

!     water density [kg/m**3]
        rfluid = rhof(i,j,k,ismpl)
!     water isobariC heat CapaCity [J/(kg*K)]
        cfluid = cpf(i,j,k,ismpl)
        rhocf = rfluid*cfluid

        RETURN
      END
