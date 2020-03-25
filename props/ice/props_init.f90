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

!>    @brief init routine
!>    @param[in] ismpl local sample index
!>    @details
!>    Set the following module variables \n
!>    - icefraq \n
!>    - liq liquidus temperature in degree Celsius \n
!>    - sol solidus temperature in degree Celsius \n
!>    - lth latent heat in J / kg \n
!>    - liqmin minimal value of fluid/ice partition function theta in [-] \n
!>
!>    Use the following module variables \n
!>    - i0 number of cells in x direction \n
!>    - j0 number of cells in y direction \n
!>    - k0 number of cells in z direction \n
      SUBROUTINE props_init(ismpl)
        use arrays
                use ice
        use mod_genrl
        use mod_time
        use mod_linfos
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k


        IF (linfos(3)>=2) WRITE(*,*) &
          ' ... permafrost properties used'

        CALL ice_allocate

        CALL read_props(ismpl)
        DO k = 1, k0
          DO j = 1, j0
            DO i = 1, i0
              icefrac(i,j,k) = 0.D0
              liq(i,j,k) = 0.D0
              sol(i,j,k) = -2.D0
            END DO
          END DO
        END DO

        lth = 333600.0D0
        liqmin = 0.03D0

        CALL check_props(ismpl)

        RETURN
      END

!>    @brief dummy
!>    @param[in] ismpl local sample index
!>    calls ice_deallocate
      SUBROUTINE props_end(ismpl)
        IMPLICIT NONE
        INTEGER ismpl

        CALL ice_deallocate
        RETURN
      END
