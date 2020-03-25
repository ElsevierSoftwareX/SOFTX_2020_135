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

!>    @brief modify coefficents for the temperature equation
!>    @param[in] ismpl local sample index
!>    @details
!> modify coefficents for the temperature equation according to the prescribed sources and sinks\n
!> rhs stored in w.\n
      SUBROUTINE set_tq(ismpl)
        use arrays
        use mod_genrl
        use mod_time
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k

        DOUBLE PRECISION deltat, deltt, qt
        EXTERNAL deltat, qt
! rhs: sources
        IF (transient .AND. tr_switch(ismpl)) THEN
          deltt = deltat(simtime(ismpl),ismpl)
!$OMP     do schedule(static) collapse(3)
          DO k = 1, k0
            DO j = 1, j0
              DO i = 1, i0
                w(i,j,k,ismpl) = w(i,j,k,ismpl) - qt(i,j,k,ismpl)
              END DO
            END DO
          END DO
!$OMP     end do nowait
        ELSE
!$OMP     do schedule(static) collapse(3)
          DO k = 1, k0
            DO j = 1, j0
              DO i = 1, i0
                w(i,j,k,ismpl) = w(i,j,k,ismpl) - qt(i,j,k,ismpl)
              END DO
            END DO
          END DO
!$OMP     end do nowait
        END IF

! Heat sources associated with flow dirichlet nodes - - - - - - -

!      do ib=first_flow,last_flow
!          i=ibc_data(ib,cbc_i)
!          j=ibc_data(ib,cbc_j)
!          k=ibc_data(ib,cbc_k)
!          bctype=ibc_data(ib,cbc_bt)
!          if (bctype.eq.bt_diri) then
!             w(i,j,k)=w(i,j,k)-qheadbcd(i,j,k,ismpl) ??? Zeit-Abhaengikeit ???
!             write(99,'(3i6,4x,e15.5)') i,j,k,qheadbcd(i,j,k,ismpl)
!          endif
!      end do

        RETURN
      END
