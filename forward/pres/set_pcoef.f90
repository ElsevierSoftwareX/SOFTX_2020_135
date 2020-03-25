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

!>    @brief calculate coefficents for the head equation
!>    @param[in] ismpl local sample index
!>    @details
!> calculate coefficents for the head equation\n
!> coefficients are stored as vectors in the diagonals a-g (d center) and rhs in w.\n
      SUBROUTINE set_pcoef(ismpl)
        use arrays
        use mod_genrl
        use mod_genrlc
        use mod_flow
        use mod_time
        use mod_linfos
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k



        DOUBLE PRECISION fi, fj, fk
        EXTERNAL fi, fj, fk


!$OMP master
        IF (linfos(3)>=2) WRITE(*,*) ' ... fcoef'
!$OMP end master

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!$OMP do schedule(static) collapse(3)
        DO k = 1, k0
          DO j = 1, j0
            DO i = 1, i0

              IF (i0>1) THEN
                IF (i<i0) THEN
                  e(i,j,k,ismpl) = fi(i,j,k,ismpl)/delx(i)
                END IF
                IF (i>1) THEN
                  c(i,j,k,ismpl) = fi(i-1,j,k,ismpl)/delx(i)
                END IF
              END IF

              IF (j0>1) THEN
                IF (j<j0) THEN
                  f(i,j,k,ismpl) = fj(i,j,k,ismpl)/dely(j)
                END IF
                IF (j>1) THEN
                  b(i,j,k,ismpl) = fj(i,j-1,k,ismpl)/dely(j)
                END IF
              END IF

              IF (k0>1) THEN
                IF (k<k0) THEN
                  g(i,j,k,ismpl) = fk(i,j,k,ismpl)/delz(k)
                END IF
                IF (k>1) THEN
                  a(i,j,k,ismpl) = fk(i,j,k-1,ismpl)/delz(k)
                END IF
              END IF

              d(i,j,k,ismpl) = -(e(i,j,k,ismpl)+c(i,j,k,ismpl)+f(i,j,k &
                ,ismpl)+b(i,j,k,ismpl)+g(i,j,k,ismpl)+a(i,j,k,ismpl))
            END DO
          END DO
        END DO
!$OMP end do nowait

        RETURN
      END

!>    @brief calculate right hand side for the head equation
!>    @param[in] ismpl local sample index
!>    @details
!> calculate right hand side for the head equation\n
!> coefficients are stored as vectors in the diagonals a-g (d center) and rhs in w.\n
      SUBROUTINE set_pcoefrs(ismpl)
        use arrays
        use mod_genrl
        use mod_genrlc
        use mod_flow
        use mod_time
        use mod_linfos
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k


        DOUBLE PRECISION src, deltf, sijk
        DOUBLE PRECISION buoy, compf, compm, rhof, por, deltat, pstor, &
          visf
        EXTERNAL buoy, compf, compm, rhof, por, deltat, pstor, visf

        deltf = deltat(simtime(ismpl),ismpl)
! rhs: sources
        IF (transient .AND. tr_switch(ismpl)) THEN
! - - - - - - - -  transient - - - - - - - - - - -
          CALL omp_mvp(i0,j0,k0,presold(1,cgen_time,ismpl), &
            x(1,1,1,ismpl),a(1,1,1,ismpl),b(1,1,1,ismpl), &
            c(1,1,1,ismpl),d(1,1,1,ismpl),e(1,1,1,ismpl), &
            f(1,1,1,ismpl),g(1,1,1,ismpl))

!$OMP   do schedule(static) collapse(3)
          DO k = 1, k0
            DO j = 1, j0
              DO i = 1, i0
                src = 0.0D0
!             buoyancy
                IF (k<k0) src = src + buoy(i,j,k,ismpl)/delz(k)
                !IF (k.eq.k0) src = src + buoy(i,j,k-1,ismpl)/delz(k)
                IF (k>1) src = src - buoy(i,j,k-1,ismpl)/delz(k)
                !IF (k.eq.1) src = src - buoy(i,j,k,ismpl)/delz(k)

                sijk = pstor(i,j,k,ismpl)
                d(i,j,k,ismpl) = d(i,j,k,ismpl) - sijk/(deltf*thetaf)
                w(i,j,k,ismpl) = w(i,j,k,ismpl) - &
                  (1.0D0-thetaf)*x(i,j,k,ismpl) - &
                  sijk*presold(i+(j-1)*i0+(k-1)*i0*j0,cgen_time,ismpl) &
                  /deltf - src
                w(i,j,k,ismpl) = w(i,j,k,ismpl)/thetaf
              END DO
            END DO
          END DO
!$OMP   end do nowait

        ELSE
! - - - - - - - - steady state  - - - - - - - - - - - - - - - - - - - - -
!$OMP   do schedule(static) collapse(3)
          DO k = 1, k0
            DO j = 1, j0
              DO i = 1, i0
                src = 0.0D0
! buoyancy  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                IF (k<k0) src = src + buoy(i,j,k,ismpl)/delz(k)
                !IF (k.eq.k0) src = src + buoy(i,j,k-1,ismpl)/delz(k)
                IF (k>1) src = src - buoy(i,j,k-1,ismpl)/delz(k)
                !IF (k.eq.1) src = src - buoy(i,j,k,ismpl)/delz(k)
                w(i,j,k,ismpl) = w(i,j,k,ismpl) - src
              END DO
            END DO
          END DO
!$OMP   end do nowait
        END IF

        RETURN
      END
