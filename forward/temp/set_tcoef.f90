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

!>    @brief calculate coefficents for the heat equation
!>    @param[in] ismpl local sample index
!>    @details
!> calculate coefficents for the heat equation\n
!> coefficients are stored as vectors in the diagonals a-g (d center) and rhs in w.\n
      SUBROUTINE set_tcoef(ismpl)
        use arrays
        use mod_genrl
        use mod_genrlc
        use mod_temp
        use mod_time
        use mod_linfos
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k



        DOUBLE PRECISION li, lj, lk, rhocf, vx, vy, vz, alfa, amean
        EXTERNAL li, lj, lk, rhocf, vx, vy, vz, alfa, amean

        DOUBLE PRECISION rijk, ra, rc, va, la, alf, p2


!$OMP master
        IF (linfos(3)>=2) WRITE(*,*) ' ... tcoef'
!$OMP end master

! initialize coefficients for sparse solvers

! inner points of grid - - - - - - - - - - - - - - - - - - - - - - - - -

!$OMP do schedule(static) collapse(3)
        DO k = 1, k0
          DO j = 1, j0
            DO i = 1, i0

              rijk = rhocf(i,j,k,ismpl)

              IF (i0>1) THEN
                IF (i<i0) THEN
                  la = li(i,j,k,ismpl)
                  ra = amean(rhocf(i+1,j,k,ismpl),rijk)
                  va = vx(i,j,k,ismpl)
                  rc = 0.5*ra*va
                  IF (la>0.D0) THEN
                    p2 = rc/la
                    alf = alfa(p2)
                  ELSE
                    alf = 0.D0
                    IF (va<0.D0) alf = -1.D0
                    IF (va>0.D0) alf = 1.D0
                  END IF
                  e(i,j,k,ismpl) = (la-(1.D0-alf)*rc)/delx(i)
                  d(i,j,k,ismpl) = d(i,j,k,ismpl) - &
                    (la+(1.D0+alf)*rc)/delx(i)
                END IF

                IF (i>1) THEN
                  la = li(i-1,j,k,ismpl)
                  ra = amean(rhocf(i-1,j,k,ismpl),rijk)
                  va = vx(i-1,j,k,ismpl)
                  rc = 0.5*ra*va
                  alf = 0.D0
                  IF (va==0.D0) THEN
                    alf = 0.D0
                  ELSE
                    IF (la>0.D0) THEN
                      p2 = rc/la
                      alf = alfa(p2)
                    ELSE
                      IF (va<0.D0) alf = -1.D0
                      IF (va>0.D0) alf = 1.D0
                    END IF
                  END IF
                  c(i,j,k,ismpl) = (la+(1.D0+alf)*rc)/delx(i)
                  d(i,j,k,ismpl) = d(i,j,k,ismpl) - &
                    (la-(1.D0-alf)*rc)/delx(i)
                END IF
              END IF

              IF (j0>1) THEN

                IF (j<j0) THEN
                  la = lj(i,j,k,ismpl)
                  ra = amean(rhocf(i,j+1,k,ismpl),rijk)
                  va = vy(i,j,k,ismpl)
                  rc = 0.5D0*ra*va
                  alf = 0.D0
                  IF (va==0.D0) THEN
                    alf = 0.D0
                  ELSE
                    IF (la>0.D0) THEN
                      p2 = rc/la
                      alf = alfa(p2)
                    ELSE
                      IF (va<0.D0) alf = -1.D0
                      IF (va>0.D0) alf = 1.D0
                    END IF
                  END IF
                  f(i,j,k,ismpl) = (la-(1.D0-alf)*rc)/dely(j)
                  d(i,j,k,ismpl) = d(i,j,k,ismpl) - &
                    (la+(1.D0+alf)*rc)/dely(j)
                END IF

                IF (j>1) THEN
                  la = lj(i,j-1,k,ismpl)
                  ra = amean(rhocf(i,j-1,k,ismpl),rijk)
                  va = vy(i,j-1,k,ismpl)
                  rc = 0.5*ra*va
                  alf = 0.D0
                  IF (va==0.D0) THEN
                    alf = 0.D0
                  ELSE
                    IF (la>0.D0) THEN
                      p2 = rc/la
                      alf = alfa(p2)
                    ELSE
                      IF (va<0.D0) alf = -1.D0
                      IF (va>0.D0) alf = 1.D0
                    END IF
                  END IF
                  b(i,j,k,ismpl) = (la+(1.D0+alf)*rc)/dely(j)
                  d(i,j,k,ismpl) = d(i,j,k,ismpl) - &
                    (la-(1.D0-alf)*rc)/dely(j)
                END IF
              END IF

              IF (k0>1) THEN

                IF (k<k0) THEN
                  la = lk(i,j,k,ismpl)
                  ra = amean(rhocf(i,j,k+1,ismpl),rijk)
                  va = vz(i,j,k,ismpl)
                  rc = 0.5*ra*va
                  alf = 0.D0
                  IF (va==0.D0) THEN
                    alf = 0.D0
                  ELSE
                    IF (la>0.D0) THEN
                      p2 = rc/la
                      alf = alfa(p2)
                    ELSE
                      IF (va<0.D0) alf = -1.D0
                      IF (va>0.D0) alf = 1.D0
                    END IF
                  END IF
                  g(i,j,k,ismpl) = (la-(1.D0-alf)*rc)/delz(k)
                  d(i,j,k,ismpl) = d(i,j,k,ismpl) - &
                    (la+(1.D0+alf)*rc)/delz(k)
                END IF

                IF (k>1) THEN
                  la = lk(i,j,k-1,ismpl)
                  ra = amean(rhocf(i,j,k-1,ismpl),rijk)
                  va = vz(i,j,k-1,ismpl)
                  rc = 0.5*ra*va
                  alf = 0.D0
                  IF (va==0.D0) THEN
                    alf = 0.D0
                  ELSE
                    IF (la>0.D0) THEN
                      p2 = rc/la
                      alf = alfa(p2)
                    ELSE
                      IF (va<0.D0) alf = -1.D0
                      IF (va>0.D0) alf = 1.D0
                    END IF
                  END IF
                  a(i,j,k,ismpl) = (la+(1.D0+alf)*rc)/delz(k)
                  d(i,j,k,ismpl) = d(i,j,k,ismpl) - &
                    (la-(1.D0-alf)*rc)/delz(k)
                END IF
              END IF

            END DO
          END DO
        END DO
!$OMP end do nowait

        RETURN
      END

!>    @brief coefficents for the heat equation (here right hand side)
!>    @param[in] ismpl local sample index
!>    @details
!> calculate coefficents for the heat equation\n
!> coefficients are stored as vectors in the diagonals a-g (d center) and rhs in w.\n
      SUBROUTINE set_tcoefrs(ismpl)
        use arrays
        use mod_genrl
        use mod_genrlc
        use mod_temp
        use mod_time
        use mod_linfos
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k



        DOUBLE PRECISION deltt, rce
        ! INTEGER c1, c2, c3, c4

        DOUBLE PRECISION rhoceff, por, qt, deltat
        EXTERNAL rhoceff, por, qt, deltat


        deltt = deltat(simtime(ismpl),ismpl)

! - - rhs: sources, also terms for transient caolculations - - - - - - - - -

        IF (transient .AND. tr_switch(ismpl)) THEN
! - - - - - - - transient  - - - - - - - - - - - - - - - - - - - - - - - - - -

          CALL omp_mvp(i0,j0,k0,tempold(1,cgen_time,ismpl), &
            x(1,1,1,ismpl),a(1,1,1,ismpl),b(1,1,1,ismpl), &
            c(1,1,1,ismpl),d(1,1,1,ismpl),e(1,1,1,ismpl), &
            f(1,1,1,ismpl),g(1,1,1,ismpl))

!$OMP do schedule(static) collapse(3)
          DO k = 1, k0
            DO j = 1, j0
              DO i = 1, i0
!                 add  solid and fluid heat production, boundary terms
                rce = rhoceff(i,j,k,ismpl)
                d(i,j,k,ismpl) = d(i,j,k,ismpl) - rce/(deltt*thetat)
                w(i,j,k,ismpl) = w(i,j,k,ismpl) - &
                  rce*tempold(i+(j-1)*i0+(k-1)*i0*j0,cgen_time,ismpl)/ &
                  deltt - (1.D0-thetat)*x(i,j,k,ismpl) - &
                  hpf*por(i,j,k,ismpl)
                w(i,j,k,ismpl) = w(i,j,k,ismpl)/thetat
              END DO
            END DO
          END DO
!$OMP    end do nowait
        ELSE
! - - - - - - - - steady state - - - - - - - - - - - - - - - - - - - - - 
!$OMP do schedule(static) collapse(3)
          DO k = 1, k0
            DO j = 1, j0
              DO i = 1, i0
!                 solid and fluid heat production
                w(i,j,k,ismpl) = w(i,j,k,ismpl) - hpf*por(i,j,k,ismpl)
              END DO
            END DO
          END DO
!$OMP    end do nowait
        END IF

        RETURN
      END
