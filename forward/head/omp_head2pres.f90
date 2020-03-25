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

!>    @brief pressure in Pa from hydraulic potential und heigth above hz=0.0d0
!>    @param[in] init flag: 0-init, 1-normal setup
!>    @param[in] ismpl local sample index
!>    @details
!>    pressure in MPa from hydraulic potential und height above hz=0.0 \n
!>    p(surface) = 0.1 MPa \n
!>    OUTPUT in Pa\n
      SUBROUTINE omp_head2pres(init,ismpl)
        use arrays
        use mod_genrl
        use mod_genrlc
        use mod_flow
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k
        INTEGER init
        DOUBLE PRECISION psurf
        DOUBLE PRECISION dif
        DOUBLE PRECISION zero
        PARAMETER (zero=0.0D0)
        

!     presetting for dirichlet boundary conditions to avoid site effects
        CALL set_dhbc(ismpl)
        CALL set_dtbc(ismpl)
!$OMP barrier
!
        psurf = 1.0D5
!
        IF (init==0) THEN
!         initialize
!$OMP     do schedule(static) collapse(3)
          DO k = 1, k0
            DO j = 1, j0
              DO i = 1, i0
                dif = head(i,j,k,ismpl) - delza(k)
                pres(i,j,k,ismpl) = psurf
                IF (dif>zero) pres(i,j,k,ismpl) = psurf + dif*rref*grav
              END DO
            END DO
          END DO
!$OMP     end do nowait
        ELSE
!        pres and temp have already appropriate values
!$OMP     do schedule(static) collapse(3)
          DO k = 1, k0
            DO j = 1, j0
              DO i = 1, i0
                dif = head(i,j,k,ismpl) - delza(k)
                pres(i,j,k,ismpl) = psurf
! jbr: pres = (h-z)*rref*grav
                IF (dif>zero) pres(i,j,k,ismpl) = psurf + &
                  rref*dif*grav
              END DO
            END DO
          END DO
!$OMP     end do nowait
        END IF
!
        RETURN
      END
