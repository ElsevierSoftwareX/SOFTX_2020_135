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

!>    @brief calculate x mass flux at cell faces
!>    @param[in] i grid indices
!>    @param[in] j grid indices
!>    @param[in] k grid indices
!>    @param[in] l species
!>    @param[in] ismpl local sample index
!>    @return x mass flux
      DOUBLE PRECISION FUNCTION sx(i,j,k,l,ismpl)
        use arrays
        use mod_genrl
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k, l
        DOUBLE PRECISION dif, di
        EXTERNAL di

        sx = 0.D0
        IF (i0>1 .AND. i<i0) THEN
          dif = conc(i+1,j,k,l,ismpl) - conc(i,j,k,l,ismpl)
          sx = -di(i,j,k,l,ismpl)*dif
        END IF

        RETURN
      END

!>    @brief calculate y mass flux at cell faces
!>    @param[in] i grid indices
!>    @param[in] j grid indices
!>    @param[in] k grid indices
!>    @param[in] l species
!>    @param[in] ismpl local sample index
!>    @return y mass flux
      DOUBLE PRECISION FUNCTION sy(i,j,k,l,ismpl)
        use arrays
        use mod_genrl
        use mod_temp
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k, l

        DOUBLE PRECISION dif, dj
        EXTERNAL dj

        sy = 0.D0
        IF (j0>1 .AND. j<j0) THEN
          dif = conc(i,j+1,k,l,ismpl) - conc(i,j,k,l,ismpl)
          sy = -dj(i,j,k,l,ismpl)*dif
        END IF

        RETURN
      END

!>    @brief calculate z mass flux at cell  faces
!>    @param[in] i grid indices
!>    @param[in] j grid indices
!>    @param[in] k grid indices
!>    @param[in] l species
!>    @param[in] ismpl local sample index
!>    @return z mass flux
      DOUBLE PRECISION FUNCTION sz(i,j,k,l,ismpl)
        use arrays
        use mod_genrl
        use mod_temp
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k, l

        DOUBLE PRECISION dif, dk
        EXTERNAL dk

        sz = 0.D0
        IF (k0>1 .AND. k<k0) THEN
          dif = conc(i,j,k+1,l,ismpl) - conc(i,j,k,l,ismpl)
          sz = -dk(i,j,k,l,ismpl)*dif
        END IF

        RETURN
      END

!>    @brief calculate x mass flux at cell centers
!>    @param[in] i grid indices
!>    @param[in] j grid indices
!>    @param[in] k grid indices
!>    @param[in] l species
!>    @param[in] ismpl local sample index
!>    @return x mass flux
      DOUBLE PRECISION FUNCTION sxc(i,j,k,l,ismpl)
        use arrays
        use mod_genrl
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k, l
        DOUBLE PRECISION d1, d2, di, amean
        EXTERNAL di, amean

        sxc = 0.D0
        IF (i0<=1) RETURN
        IF (i>1 .AND. i<i0) THEN
          d1 = conc(i+1,j,k,l,ismpl) - conc(i,j,k,l,ismpl)
          d2 = conc(i,j,k,l,ismpl) - conc(i-1,j,k,l,ismpl)
          sxc = amean(-di(i,j,k,l,ismpl)*d1,-di(i-1,j,k,l,ismpl)*d2)
        ELSE IF (i==1) THEN
          sxc = -di(i,j,k,l,ismpl)*(conc(i+1,j,k,l,ismpl)-conc(i,j,k,l &
            ,ismpl))
        ELSE IF (i==i0) THEN
          sxc = -di(i-1,j,k,l,ismpl)*(conc(i,j,k,l,ismpl)-conc(i-1,j,k &
            ,l,ismpl))
        END IF

        RETURN
      END

!>    @brief calculate y mass fluxat cell center
!>    @param[in] i grid indices
!>    @param[in] j grid indices
!>    @param[in] k grid indices
!>    @param[in] l species
!>    @param[in] ismpl local sample index
!>    @return y mass flux
      DOUBLE PRECISION FUNCTION syc(i,j,k,l,ismpl)
        use arrays
        use mod_genrl
        use mod_temp
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k, l

        DOUBLE PRECISION d1, d2, dj, amean
        EXTERNAL dj, amean

        syc = 0.D0
        IF (j0<=1) RETURN
        IF (j>1 .AND. j<j0) THEN
          d1 = conc(i,j+1,k,l,ismpl) - conc(i,j,k,l,ismpl)
          d2 = conc(i,j,k,l,ismpl) - conc(i,j-1,k,l,ismpl)
          syc = amean(-dj(i,j,k,l,ismpl)*d1,-dj(i,j-1,k,l,ismpl)*d2)
        ELSE IF (j==1) THEN
          syc = -dj(i,j,k,l,ismpl)*(conc(i,j+1,k,l,ismpl)-conc(i,j,k,l &
            ,ismpl))
        ELSE IF (j==j0) THEN
          syc = -dj(i,j-1,k,l,ismpl)*(conc(i,j,k,l,ismpl)-conc(i,j-1,k &
            ,l,ismpl))
        END IF

        RETURN
      END

!>    @brief calculate z mass flux at cell center
!>    @param[in] i grid indices
!>    @param[in] j grid indices
!>    @param[in] k grid indices
!>    @param[in] l species
!>    @param[in] ismpl local sample index
!>    @return z mass flux
      DOUBLE PRECISION FUNCTION szc(i,j,k,l,ismpl)
        use arrays
        use mod_genrl
        use mod_temp
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k, l

        DOUBLE PRECISION d1, d2, dk, amean
        EXTERNAL dk, amean

        szc = 0.D0
        IF (k0<=1) RETURN
        IF (k>1 .AND. k<k0) THEN
          d1 = conc(i,j,k+1,l,ismpl) - conc(i,j,k,l,ismpl)
          d2 = conc(i,j,k,l,ismpl) - conc(i,j,k-1,l,ismpl)
          szc = amean(-dk(i,j,k,l,ismpl)*d1,-dk(i,j,k-1,l,ismpl)*d2)
        ELSE IF (k==1) THEN
          szc = -dk(i,j,k,l,ismpl)*(conc(i,j,k+1,l,ismpl)-conc(i,j,k,l &
            ,ismpl))
        ELSE IF (k==k0) THEN
          szc = -dk(i,j,k-1,l,ismpl)*(conc(i,j,k,l,ismpl)-conc(i,j,k-1 &
            ,l,ismpl))
        END IF

        RETURN
      END

!>    @brief average effective diffusivities on cell faces in x direction
!>    @param[in] i grid indices
!>    @param[in] j grid indices
!>    @param[in] k grid indices
!>    @param[in] spec species
!>    @param[in] ismpl local sample index
!>    @return effective diffusivities (J/mK)
      DOUBLE PRECISION FUNCTION di(i,j,k,spec,ismpl)
        use arrays
        use mod_genrl
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k
        INTEGER spec
        DOUBLE PRECISION f1, f2, prod, summ, betx, bety, betz, bet
        DOUBLE PRECISION por, disp, vx, vy, vz
        EXTERNAL por, disp, vx, vy, vz

        di = 0.D0
        betx = 0.D0
        bety = 0.D0
        betz = 0.D0
        IF (k0>1 .AND. k<k0) THEN
               betz = vz(i,j,k,ismpl)
               betz = betz*betz
        END IF
        IF (j0>1 .AND. j<j0) THEN
               bety = vy(i,j,k,ismpl)
               bety = bety*bety
        END IF
        IF (i0>1 .AND. i<i0) THEN
          betx = vx(i,j,k,ismpl)
          betx = betx*betx
          bet = SQRT(betx + bety + betz)
          f1 = por(i,j,k,ismpl)*diff_c(spec) + disp(i,j,k,ismpl)*bet
          f2 = por(i+1,j,k,ismpl)*diff_c(spec) + &
            disp(i+1,j,k,ismpl)*bet
          prod = f1*f2
          summ = f1*delx(i+1) + f2*delx(i)
          IF (summ>0.D0) di = 2.D0*prod/summ
        END IF

        RETURN
      END

!>    @brief average effective diffusivities on cell faces in y direction
!>    @param[in] i grid indices
!>    @param[in] j grid indices
!>    @param[in] k grid indices
!>    @param[in] spec species
!>    @param[in] ismpl local sample index
!>    @return effective diffusivities (J/mK)
      DOUBLE PRECISION FUNCTION dj(i,j,k,spec,ismpl)
        use arrays
        use mod_genrl
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k
        INTEGER spec
        DOUBLE PRECISION f1, f2, prod, summ, betx, bety, betz, bet
        DOUBLE PRECISION por, disp, vy, vx, vz
        EXTERNAL por, disp, vy, vx, vz

        dj = 0.D0
        betx = 0.D0
        bety = 0.D0
        betz = 0.D0
        IF (k0>1 .AND. k<k0) THEN
               betz = vz(i,j,k,ismpl)
               betz = betz*betz
        END IF
        IF (i0>1 .AND. i<i0) THEN
               betx = vx(i,j,k,ismpl)
               betx = betx*betx
        END IF
        IF (j0>1 .AND. j<j0) THEN
          bety = vy(i,j,k,ismpl)
          bety = bety*bety
          bet = SQRT(betx + bety + betz)
          f1 = por(i,j,k,ismpl)*diff_c(spec) + disp(i,j,k,ismpl)*bet
          f2 = por(i,j+1,k,ismpl)*diff_c(spec) + &
            disp(i,j+1,k,ismpl)*bet
          prod = f1*f2
          summ = f1*dely(j+1) + f2*dely(j)
          IF (summ>0.D0) dj = 2.D0*prod/summ
        END IF

        RETURN
      END

!>    @brief average effective diffusivities on cell faces in z direction
!>    @param[in] i grid indices
!>    @param[in] j grid indices
!>    @param[in] k grid indices
!>    @param[in] spec species
!>    @param[in] ismpl local sample index
!>    @return effective diffusivities (J/mK)
      DOUBLE PRECISION FUNCTION dk(i,j,k,spec,ismpl)
        use arrays
        use mod_genrl
        IMPLICIT NONE
        integer :: i, j, k
        integer :: ismpl
        INTEGER spec
        DOUBLE PRECISION f1, f2, prod, summ, betx, bety, betz, bet
        DOUBLE PRECISION por, disp, vz, vy, vx
        EXTERNAL por, disp, vz, vy, vx

        dk = 0.D0
        betx = 0.D0
        bety = 0.D0
        betz = 0.D0
        IF (j0>1 .AND. j<j0) THEN
               bety = vy(i,j,k,ismpl)
               bety = bety*bety
        END IF
        IF (i0>1 .AND. i<i0) THEN
               betx = vx(i,j,k,ismpl)
               betx = betx*betx
        END IF
        IF (k0>1 .AND. k<k0) THEN
          betz = vz(i,j,k,ismpl)
          betz = betz*betz
          bet = SQRT(betx + bety + betz)
          f1 = por(i,j,k,ismpl)*diff_c(spec) + disp(i,j,k,ismpl)*bet
          f2 = por(i,j,k+1,ismpl)*diff_c(spec) + &
            disp(i,j,k+1,ismpl)*bet
          prod = f1*f2
          summ = f1*delz(k+1) + f2*delz(k)
          IF (summ>0.D0) dk = 2.D0*prod/summ
        END IF

        RETURN
      END
