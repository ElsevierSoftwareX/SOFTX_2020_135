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

!>    @brief calculate coefficents for the transport equation
!>    @param[in] spec species index
!>    @param[in] ismpl local sample index
!>    @details
!> calculate coefficents for the transport equation\n
!> coefficients are stored as vectors in the diagonals a-g (d center) and rhs in w.\n
      SUBROUTINE set_ccoef(spec,ismpl)
        use arrays
        use mod_genrl
        use mod_genrlc
        use mod_flow
        use mod_conc
        use mod_time
        use mod_linfos
        IMPLICIT NONE
        integer :: i, j, k
        integer :: ismpl
        DOUBLE PRECISION di, dj, dk, por, vx, vy, vz, alfa
        EXTERNAL di, dj, dk, por, vx, vy, vz, alfa
        DOUBLE PRECISION v2, de, alf, p2
        INTEGER spec
!debug      write(99,'(a,3i4)') 'ShemSUITE i0,j0,k0:',i0,j0,k0


!$OMP master
        IF (linfos(3)>=2) WRITE(*,*) ' ... ccoef'
!$OMP end master

! inner points of grid - - - - - - - - - - - - - - - - - - - - - - - - -

!$OMP do schedule(static) collapse(3)
        DO k = 1, k0
          DO j = 1, j0
            DO i = 1, i0

              IF (i0>1) THEN
                IF (i<i0) THEN
                  de = di(i,j,k,spec,ismpl)
                  v2 = 0.5D0*vx(i,j,k,ismpl)
                  IF (de>0.D0) THEN
                    p2 = v2/de
                    alf = alfa(p2)
                  ELSE
                    alf = 0.D0
                    IF (v2<0.D0) alf = -1.D0
                    IF (v2>0.D0) alf = 1.D0
                  END IF
                  e(i,j,k,ismpl) = (de-(1.D0-alf)*v2)/delx(i)
                  d(i,j,k,ismpl) = d(i,j,k,ismpl) - &
                    (de+(1.D0+alf)*v2)/delx(i)
                END IF

                IF (i>1) THEN
                  de = di(i-1,j,k,spec,ismpl)
                  v2 = 0.5*vx(i-1,j,k,ismpl)
                  alf = 0.D0
                  IF (v2==0.D0) THEN
                    alf = 0.D0
                  ELSE
                    IF (de>0.D0) THEN
                      p2 = v2/de
                      alf = alfa(p2)
                    ELSE
                      IF (v2<0.D0) alf = -1.D0
                      IF (v2>0.D0) alf = 1.D0
                    END IF
                  END IF

                  c(i,j,k,ismpl) = (de+(1.D0+alf)*v2)/delx(i)
                  d(i,j,k,ismpl) = d(i,j,k,ismpl) - &
                    (de-(1.D0-alf)*v2)/delx(i)
                END IF
              END IF



              IF (j0>1) THEN

                IF (j<j0) THEN

                  de = dj(i,j,k,spec,ismpl)
                  v2 = 0.5*vy(i,j,k,ismpl)
                  alf = 0.D0
                  IF (v2==0.D0) THEN
                    alf = 0.D0
                  ELSE
                    IF (de>0.D0) THEN
                      p2 = v2/de
                      alf = alfa(p2)
                    ELSE
                      IF (v2<0.D0) alf = -1.D0
                      IF (v2>0.D0) alf = 1.D0
                    END IF
                  END IF

                  f(i,j,k,ismpl) = (de-(1.D0-alf)*v2)/dely(j)
                  d(i,j,k,ismpl) = d(i,j,k,ismpl) - &
                    (de+(1.D0+alf)*v2)/dely(j)
                END IF


                IF (j>1) THEN

                  de = dj(i,j-1,k,spec,ismpl)
                  v2 = 0.5*vy(i,j-1,k,ismpl)
                  alf = 0.D0
                  IF (v2==0.D0) THEN
                    alf = 0.D0
                  ELSE
                    IF (de>0.D0) THEN
                      p2 = v2/de
                      alf = alfa(p2)
                    ELSE
                      IF (v2<0.D0) alf = -1.D0
                      IF (v2>0.D0) alf = 1.D0
                    END IF
                  END IF

                  b(i,j,k,ismpl) = (de+(1.D0+alf)*v2)/dely(j)
                  d(i,j,k,ismpl) = d(i,j,k,ismpl) - &
                    (de-(1.D0-alf)*v2)/dely(j)
                END IF
              END IF

              IF (k0>1) THEN

                IF (k<k0) THEN

                  de = dk(i,j,k,spec,ismpl)
                  v2 = 0.5D0*vz(i,j,k,ismpl)
                  alf = 0.D0
                  IF (v2==0.D0) THEN
                    alf = 0.D0
                  ELSE
                    IF (de>0.D0) THEN
                      p2 = v2/de
                      alf = alfa(p2)
                    ELSE
                      IF (v2<0.D0) alf = -1.D0
                      IF (v2>0.D0) alf = 1.D0
                    END IF
                  END IF

                  g(i,j,k,ismpl) = (de-(1.D0-alf)*v2)/delz(k)
                  d(i,j,k,ismpl) = d(i,j,k,ismpl) - &
                    (de+(1.D0+alf)*v2)/delz(k)
                END IF

                IF (k>1) THEN

                  de = dk(i,j,k-1,spec,ismpl)
                  v2 = 0.5D0*vz(i,j,k-1,ismpl)
                  alf = 0.D0
                  IF (v2==0.D0) THEN
                    alf = 0.D0
                  ELSE
                    IF (de>0.D0) THEN
                      p2 = v2/de
                      alf = alfa(p2)
                    ELSE
                      IF (v2<0.D0) alf = -1.D0
                      IF (v2>0.D0) alf = 1.D0
                    END IF
                  END IF

                  a(i,j,k,ismpl) = (de+(1.D0+alf)*v2)/delz(k)
                  d(i,j,k,ismpl) = d(i,j,k,ismpl) - &
                    (de-(1.D0-alf)*v2)/delz(k)

                END IF

              END IF

            END DO
          END DO
        END DO
!$OMP end do nowait

        RETURN
      END

!>    @brief coefficents for the transport equation (here only the right side)
!>    @param[in] spec species index
!>    @param[in] ismpl local sample index
!>    @details
!> calculate coefficents for the transport equation\n
!> coefficients are stored as vectors in the diagonals a-g (d center) and rhs in w.\n
      SUBROUTINE set_ccoefrs(spec,ismpl)
        use arrays
        use mod_genrl
        use mod_genrlc
        use mod_conc
        use mod_time
        use mod_linfos
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k


        ! DOUBLE PRECISION rhoceff
        DOUBLE PRECISION deltt
        ! INTEGER c1, c2, c3, c4
        INTEGER spec
        DOUBLE PRECISION deltat, por
        EXTERNAL deltat, por


        deltt = deltat(simtime(ismpl),ismpl)

! - - rhs: sources, also terms for transient calculations - - - - - - - - -

        IF (transient .AND. tr_switch(ismpl)) THEN

! - - - - - - - transient - - - - - - - - - - - - - - - - - - - - - - - - -
          CALL omp_mvp(i0,j0,k0,concold(1,spec,cgen_time,ismpl), &
            x(1,1,1,ismpl),a(1,1,1,ismpl),b(1,1,1,ismpl), &
            c(1,1,1,ismpl),d(1,1,1,ismpl),e(1,1,1,ismpl), &
            f(1,1,1,ismpl),g(1,1,1,ismpl))

!$OMP    do schedule(static) collapse(3)
          DO k = 1, k0
            DO j = 1, j0
              DO i = 1, i0

                d(i,j,k,ismpl) = d(i,j,k,ismpl) - &
                  por(i,j,k,ismpl)/(deltt*thetac)
                w(i,j,k,ismpl) = w(i,j,k,ismpl) - &
                  por(i,j,k,ismpl)*concold(i+(j-1)*i0+(k-1)*i0*j0, &
                  spec,cgen_time,ismpl)/deltt - &
                  (1.D0-thetac)*x(i,j,k,ismpl)
                w(i,j,k,ismpl) = w(i,j,k,ismpl)/thetac
              END DO
            END DO
          END DO
!$OMP    end do nowait

        ELSE

! - - - - - - - - steady state - - - - - - - - - - - - - - - - - - - - - 
! C$OMP    do schedule(static) collapse(3)
!          do k=1,k0
!             do j=1,j0
!                do i=1,i0
!                   w(i,j,k,ismpl) = w(i,j,k,ismpl)
!                     - hs(i,j,k,spec,ismpl)
!                end do
!             end do
!          end do
! C$OMP    end do nowait
        END IF

        RETURN
      END
