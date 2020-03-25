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

!>    @brief calculate velocities at cell faces
!>    @param[in] i grid indices
!>    @param[in] j grid indices
!>    @param[in] k grid indices
!>    @param[in] ismpl local sample index
!>    @return x velocity (m/(Pa s))
      DOUBLE PRECISION FUNCTION vx(i,j,k,ismpl)
        use arrays
        use mod_genrl
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k
        DOUBLE PRECISION f1, f2, dif, gi
        EXTERNAL gi

        vx = 0.D0
        IF (i0>1 .AND. i<i0) THEN
          dif = pres(i+1,j,k,ismpl) - pres(i,j,k,ismpl)
          vx = -gi(i,j,k,ismpl)*dif
        END IF
        RETURN
      END

!>    @brief calculate velocities at cell faces
!>    @param[in] i grid indices
!>    @param[in] j grid indices
!>    @param[in] k grid indices
!>    @param[in] ismpl local sample index
!>    @return y velocity (m/(Pa s))
      DOUBLE PRECISION FUNCTION vy(i,j,k,ismpl)
        use arrays
        use mod_genrl
        use mod_flow
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k

        DOUBLE PRECISION f1, f2, dif, gj
        EXTERNAL gj

        vy = 0.D0
        IF (j0>1 .AND. j<j0) THEN
          dif = pres(i,j+1,k,ismpl) - pres(i,j,k,ismpl)
          vy = -gj(i,j,k,ismpl)*dif
        END IF
        RETURN
      END

!>    @brief calculate velocities at cell faces
!>    @param[in] i grid indices
!>    @param[in] j grid indices
!>    @param[in] k grid indices
!>    @param[in] ismpl local sample index
!>    @return z velocity (m/(Pa s))
      DOUBLE PRECISION FUNCTION vz(i,j,k,ismpl)
        use arrays
        use mod_genrl
        use mod_flow
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k

        DOUBLE PRECISION f1, f2, dif, gk, vbuoy, por, rhof, rhav
        EXTERNAL gk, rhof, por

        vz = 0.D0
        IF (k0>1 .AND. k<k0) THEN
!          rhav = 0.5D0*(rhof(i,j,k+1,ismpl)+rhof(i,j,k,ismpl))&
!                *(delza(k+1 )-delza(k))
          dif = pres(i,j,k+1,ismpl) - pres(i,j,k,ismpl)
          IF (por(i,j,k,ismpl).GT.1.E-19) THEN
                   rhav = rhof(i,j,k,ismpl)
          ELSE
                   rhav = 1.29E0
          END IF
          IF (por(i,j,k+1,ismpl).GT.1.E-19) THEN
                   rhav = rhav + rhof(i,j,k+1,ismpl)
          ELSE
                   rhav = rhav + 1.29E0
          END IF
          rhav = 0.5*rhav*(delza(k+1 )-delza(k))
          dif =  dif + rhav*grav
          vz = -gk(i,j,k,ismpl)*dif! - vbuoy(i,j,k,ismpl)
        END IF

        RETURN
      END

!>    @brief calculate x-velocities at cell centers
!>    @param[in] i grid indices
!>    @param[in] j grid indices
!>    @param[in] k grid indices
!>    @param[in] ismpl local sample index
!>    @return x velocity (m/(Pa s))
      DOUBLE PRECISION FUNCTION vxc(i,j,k,ismpl)
        use arrays
        use mod_genrl
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k
        DOUBLE PRECISION vx, amean
        EXTERNAL vx, amean

        vxc = 0.D0
        IF (i0<=1) RETURN
        IF (i>1 .AND. i<i0) THEN
          vxc = amean(vx(i,j,k,ismpl),vx(i-1,j,k,ismpl))
        ELSE IF (i==1) THEN
          vxc = vx(i,j,k,ismpl)
        ELSE IF (i==i0) THEN
          vxc = vx(i-1,j,k,ismpl)
        END IF
        RETURN
      END

!>    @brief calculate y-velocities at cell centers
!>    @param[in] i grid indices
!>    @param[in] j grid indices
!>    @param[in] k grid indices
!>    @param[in] ismpl local sample index
!>    @return y velocity (m/(Pa s))
      DOUBLE PRECISION FUNCTION vyc(i,j,k,ismpl)
        use arrays
        use mod_genrl
        use mod_flow
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k

        DOUBLE PRECISION vy, amean
        EXTERNAL vy, amean

        vyc = 0.D0
        IF (j0<=1) RETURN
        IF (j>1 .AND. j<j0) THEN
          vyc = amean(vy(i,j,k,ismpl),vy(i,j-1,k,ismpl))
        ELSE IF (j==1) THEN
          vyc = vy(i,j,k,ismpl)
        ELSE IF (j==j0) THEN
          vyc = vy(i,j-1,k,ismpl)
        END IF
        RETURN
      END

!>    @brief calculate z-velocities at cell centers
!>    @param[in] i grid indices
!>    @param[in] j grid indices
!>    @param[in] k grid indices
!>    @param[in] ismpl local sample index
!>    @return z velocity (m/(Pa s))
      DOUBLE PRECISION FUNCTION vzc(i,j,k,ismpl)
        use arrays
        use mod_genrl
        use mod_flow
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k

        DOUBLE PRECISION vz, amean
        EXTERNAL vz, amean

        vzc = 0.D0
        IF (k0<=1) RETURN
        IF (k>1 .AND. k<k0) THEN
          vzc = amean(vz(i,j,k,ismpl),vz(i,j,k-1,ismpl))
        ELSE IF (k==1) THEN
          vzc = vz(i,j,k,ismpl)
        ELSE IF (k==k0) THEN
          vzc = vz(i,j,k-1,ismpl)
        END IF
        RETURN
      END

!>    @brief average  conductivities on cell faces in x direction
!>    @param[in] i grid indices
!>    @param[in] j grid indices
!>    @param[in] k grid indices
!>    @param[in] ismpl local sample index
!>    @return x  conductivity (m/(Pa s))
      DOUBLE PRECISION FUNCTION fi(i,j,k,ismpl)
        use arrays
        use mod_genrl
        use mod_flow
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k

        DOUBLE PRECISION f1, f2, prod, summ, kx, rhof, visf
        EXTERNAL kx, rhof, visf

        fi = 0.D0
        IF (i0>1 .AND. i<i0) THEN
          f1 = kx(i,j,k,ismpl)/visf(i,j,k,ismpl)
          f2 = kx(i+1,j,k,ismpl)/ &
            visf(i+1,j,k,ismpl)
          prod = f1*f2
          summ = f1*delx(i+1) + f2*delx(i)
          IF (summ>0.D0) fi = 2.D0*prod/summ
        END IF
        RETURN
      END

!>    @brief average  conductivities on cell faces in y direction
!>    @param[in] i grid indices
!>    @param[in] j grid indices
!>    @param[in] k grid indices
!>    @param[in] ismpl local sample index
!>    @return y  conductivity (m/(Pa s))
      DOUBLE PRECISION FUNCTION fj(i,j,k,ismpl)
        use arrays
        use mod_genrl
        use mod_flow
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k

        DOUBLE PRECISION f1, f2, prod, summ, ky, rhof, visf
        EXTERNAL ky, rhof, visf

        fj = 0.D0
        IF (j0>1 .AND. j<j0) THEN
          f1 = ky(i,j,k,ismpl)/visf(i,j,k,ismpl)
          f2 = ky(i,j+1,k,ismpl)/ &
            visf(i,j+1,k,ismpl)
          prod = f1*f2
          summ = f1*dely(j+1) + f2*dely(j)
          IF (summ>0.D0) fj = 2.D0*prod/summ
        END IF
        RETURN
      END

!>    @brief average  conductivities on cell faces in z direction
!>    @param[in] i grid indices
!>    @param[in] j grid indices
!>    @param[in] k grid indices
!>    @param[in] ismpl local sample index
!>    @return z  conductivity (m/(Pa s))
      DOUBLE PRECISION FUNCTION fk(i,j,k,ismpl)
        use arrays
        use mod_genrl
        use mod_flow
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k

        DOUBLE PRECISION f1, f2, prod, summ, kz, rhof, visf
        EXTERNAL kz, rhof, visf

        fk = 0.D0
        IF (k0>1 .AND. k<k0) THEN
          f1 = kz(i,j,k,ismpl)/visf(i,j,k,ismpl)
          f2 = kz(i,j,k+1,ismpl)/ &
            visf(i,j,k+1,ismpl)
          prod = f1*f2
          summ = f1*delz(k+1) + f2*delz(k)
          IF (summ>0.D0) fk = 2.D0*prod/summ
        END IF
        RETURN
      END

!>    @brief average  conductivities on cell faces in x direction
!>    @param[in] i grid indices
!>    @param[in] j grid indices
!>    @param[in] k grid indices
!>    @param[in] ismpl local sample index
!>    @return x  conductivity (m/(Pa s))
      DOUBLE PRECISION FUNCTION gi(i,j,k,ismpl)
        use arrays
        use mod_genrl
        use mod_flow
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k

        DOUBLE PRECISION f1, f2, prod, summ, kx, rhof, visf
        EXTERNAL kx, rhof, visf

        gi = 0.D0
        IF (i0>1 .AND. i<i0) THEN
          f1 = kx(i,j,k,ismpl)/visf(i,j,k,ismpl)
          f2 = kx(i+1,j,k,ismpl)/visf(i+1,j,k,ismpl)
          prod = f1*f2
          summ = f1*delx(i+1) + f2*delx(i)
          IF (summ>0.D0) gi = 2.D0*prod/summ
        END IF
        RETURN
      END

!>    @brief average  conductivities on cell faces in y direction
!>    @param[in] i grid indices
!>    @param[in] j grid indices
!>    @param[in] k grid indices
!>    @param[in] ismpl local sample index
!>    @return y  conductivity (m/(Pa s))
      DOUBLE PRECISION FUNCTION gj(i,j,k,ismpl)
        use arrays
        use mod_genrl
        use mod_flow
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k

        DOUBLE PRECISION f1, f2, prod, summ, ky, rhof, visf
        EXTERNAL ky, rhof, visf

        gj = 0.D0
        IF (j0>1 .AND. j<j0) THEN
          f1 = ky(i,j,k,ismpl)/visf(i,j,k,ismpl)
          f2 = ky(i,j+1,k,ismpl)/visf(i,j+1,k,ismpl)
          prod = f1*f2
          summ = f1*dely(j+1) + f2*dely(j)
          IF (summ>0.D0) gj = 2.D0*prod/summ
        END IF
        RETURN
      END

!>    @brief average  conductivities on cell faces in z direction
!>    @param[in] i grid indices
!>    @param[in] j grid indices
!>    @param[in] k grid indices
!>    @param[in] ismpl local sample index
!>    @return z  conductivity (m/(Pa s))
      DOUBLE PRECISION FUNCTION gk(i,j,k,ismpl)
        use arrays
        use mod_genrl
        use mod_flow
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k

        DOUBLE PRECISION f1, f2, prod, summ, kz, rhof, visf
        EXTERNAL kz, rhof, visf

        gk = 0.D0
        IF (k0>1 .AND. k<k0) THEN
          f1 = kz(i,j,k,ismpl)/visf(i,j,k,ismpl)
          f2 = kz(i,j,k+1,ismpl)/visf(i,j,k+1,ismpl)
          prod = f1*f2
          summ = f1*delz(k+1) + f2*delz(k)
          IF (summ>0.D0) gk = 2.D0*prod/summ
        END IF
        RETURN
      END
