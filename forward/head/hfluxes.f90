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

!>    @brief calculate velocities at cell  faces
!>    @param[in] i grid indices
!>    @param[in] j grid indices
!>    @param[in] k grid indices
!>    @param[in] ismpl local sample index
!>    @return x velocity (m/s)
      DOUBLE PRECISION FUNCTION vx(i,j,k,ismpl)
        use arrays
        use mod_genrl
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k

        DOUBLE PRECISION dif, ti
        EXTERNAL ti

        vx = 0.D0
        if (.not. head_active .and. vdefaultswitch) then
           vx = vdefault(1,ismpl)
        end if
        IF (i0>1 .AND. i<i0 .AND. head_active) THEN
#ifdef head_base
           dif = head(i+1,j,k,ismpl) - head(i,j,k,ismpl)
           vx = -ti(i,j,k,ismpl)*dif
#endif
        END IF
        RETURN
      END

!>    @brief calculate velocities at cell  faces
!>    @param[in] i grid indices
!>    @param[in] j grid indices
!>    @param[in] k grid indices
!>    @param[in] ismpl local sample index
!>    @return y velocity (m/s)
      DOUBLE PRECISION FUNCTION vy(i,j,k,ismpl)
        use arrays
        use mod_genrl
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k

        DOUBLE PRECISION dif, tj
        EXTERNAL tj

        vy = 0.D0
        if (.not. head_active .and. vdefaultswitch) then
           vy = vdefault(2,ismpl)
        end if
        IF (j0>1 .AND. j<j0 .AND. head_active) THEN
#ifdef head_base
           dif = head(i,j+1,k,ismpl) - head(i,j,k,ismpl)
           vy = -tj(i,j,k,ismpl)*dif
#endif
        END IF
        RETURN
      END

!>    @brief calculate velocities at cell  faces
!>    @param[in] i grid indices
!>    @param[in] j grid indices
!>    @param[in] k grid indices
!>    @param[in] ismpl local sample index
!>    @return z velocity (m/s)
      DOUBLE PRECISION FUNCTION vz(i,j,k,ismpl)
        use arrays
        use mod_genrl
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k

        DOUBLE PRECISION dif, tk, buoy
        EXTERNAL tk, buoy

        vz = 0.D0
        if (.not. head_active .and. vdefaultswitch) then
           vz = vdefault(3,ismpl)
        end if
        IF (k0>1 .AND. k<k0 .AND. head_active) THEN
#ifdef head_base
           dif = head(i,j,k+1,ismpl) - head(i,j,k,ismpl)
           vz = -tk(i,j,k,ismpl)*dif - buoy(i,j,k,ismpl)
#endif
        END IF
        RETURN
      END

!>    @brief calculate velocities at cell centers
!>    @param[in] i grid indices
!>    @param[in] j grid indices
!>    @param[in] k grid indices
!>    @param[in] ismpl local sample index
!>    @return x velocity (m/s)
      DOUBLE PRECISION FUNCTION vxc(i,j,k,ismpl)
        use arrays
        use mod_genrl
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k
        DOUBLE PRECISION vx, amean
        EXTERNAL vx, amean

        vxc = 0.D0
        IF (i0<=1 .OR. .NOT.head_active) RETURN
#ifdef head_base
        IF (i>1 .AND. i<i0) THEN
           vxc = amean(vx(i,j,k,ismpl),vx(i-1,j,k,ismpl))
        ELSE IF (i==1) THEN
           vxc = vx(i,j,k,ismpl)
        ELSE IF (i==i0) THEN
           vxc = vx(i-1,j,k,ismpl)
        END IF
#endif
        RETURN
      END

!>    @brief calculate y-velocities at cell center
!>    @param[in] i grid indices
!>    @param[in] j grid indices
!>    @param[in] k grid indices
!>    @param[in] ismpl local sample index
!>    @return y velocity (m/s)
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
        IF (j0<=1 .OR. .NOT.head_active) RETURN
#ifdef head_base
        IF (j>1 .AND. j<j0) THEN
           vyc = amean(vy(i,j,k,ismpl),vy(i,j-1,k,ismpl))
        ELSE IF (j==1) THEN
           vyc = vy(i,j,k,ismpl)
        ELSE IF (j==j0) THEN
           vyc = vy(i,j-1,k,ismpl)
        END IF
#endif
        RETURN
      END

!>    @brief calculate z-velocities at cell center
!>    @param[in] i grid indices
!>    @param[in] j grid indices
!>    @param[in] k grid indices
!>    @param[in] ismpl local sample index
!>    @return z velocity (m/s)
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
        IF (k0<=1 .OR. .NOT.head_active) RETURN
#ifdef head_base
        IF (k>1 .AND. k<k0) THEN
           vzc = amean(vz(i,j,k,ismpl),vz(i,j,k-1,ismpl))
        ELSE IF (k==1) THEN
           vzc = vz(i,j,k,ismpl)
        ELSE IF (k==k0) THEN
           vzc = vz(i,j,k-1,ismpl)
        END IF
#endif
        RETURN
      END

!>    @brief harmonic mean Kx on cell faces in x direction over delx*
!>    @param[in] i grid indices
!>    @param[in] j grid indices
!>    @param[in] k grid indices
!>    @param[in] ismpl local sample index
!>    @return x hydraulic conductivity over delx*  (1/s)
!>    @details
!>    Compute the harmonic mean of Kx on the cell face in positive
!>    x-direction from the current node (i, j, k) divided by the
!>    x-distance of the current node (i, j, k) to the neighboring node
!>    (i+1, j, k).
!>
!>    delx* = 0.5 ( delx(i) + delx(i+1) )
!>
!>    Kx / delx* = [ 0.5*delx(i+1)/K(i+1) + 0.5*delx(i)/K(i) ]**-1
!>
!>               = [ ( 0.5* K(i)* delx(i+1) + 0.5*K(i+1)* delx(i) ) / ( K(i)*K(i+1) ) ]**-1
!>
!>               = [ ( 0.5 * summ ) / (prod) ]**-1 = 2.0*prod/summ
      DOUBLE PRECISION FUNCTION ti(i,j,k,ismpl)
        use arrays
        use mod_genrl
        use mod_flow
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k

        DOUBLE PRECISION f1, f2, prod, summ, kx, rhof, visf
        EXTERNAL kx, rhof, visf

        ti = 0.D0
        IF (i0>1 .AND. i<i0) THEN
          f1 = kx(i,j,k,ismpl)*rhof(i,j,k,ismpl)*grav/ &
            visf(i,j,k,ismpl)
          f2 = kx(i+1,j,k,ismpl)*rhof(i+1,j,k,ismpl)*grav/ &
            visf(i+1,j,k,ismpl)
          prod = f1*f2
          summ = f1*delx(i+1) + f2*delx(i)
          IF (summ>0.D0) ti = 2.D0*prod/summ
        END IF
        RETURN
      END

!>    @brief harmonic mean Ky on cell faces in y direction over dely*
!>    @param[in] i grid indices
!>    @param[in] j grid indices
!>    @param[in] k grid indices
!>    @param[in] ismpl local sample index
!>    @return y hydraulic conductivity over dely*  (1/s)
!>    @details
!>    Compute the harmonic mean of Ky on the cell face in positive
!>    y-direction from the current node (i, j, k) divided by the
!>    y-distance of the current node (i, j, k) to the neighboring node
!>    (i, j+1, k).
!>
!>    dely* = 0.5 ( dely(j) + dely(j+1) )
!>
!>    Ky / dely* = [ 0.5*dely(j+1)/K(j+1) + 0.5*dely(j)/K(j) ]**-1
!>
!>               = [ ( 0.5* K(j)* dely(j+1) + 0.5*K(j+1)* dely(j) ) / ( K(j)*K(j+1) ) ]**-1
!>
!>               = [ ( 0.5 * summ ) / (prod) ]**-1 = 2.0*prod/summ
      DOUBLE PRECISION FUNCTION tj(i,j,k,ismpl)
        use arrays
        use mod_genrl
        use mod_flow
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k

        DOUBLE PRECISION f1, f2, prod, summ, ky, rhof, visf
        EXTERNAL ky, rhof, visf

        tj = 0.D0
        IF (j0>1 .AND. j<j0) THEN
          f1 = ky(i,j,k,ismpl)*rhof(i,j,k,ismpl)*grav/ &
            visf(i,j,k,ismpl)
          f2 = ky(i,j+1,k,ismpl)*rhof(i,j+1,k,ismpl)*grav/ &
            visf(i,j+1,k,ismpl)
          prod = f1*f2
          summ = f1*dely(j+1) + f2*dely(j)
          IF (summ>0.D0) tj = 2.D0*prod/summ
        END IF
        RETURN
      END

!>    @brief harmonic mean Kz on cell faces in z direction over delz*
!>    @param[in] i grid indices
!>    @param[in] j grid indices
!>    @param[in] k grid indices
!>    @param[in] ismpl local sample index
!>    @return z hydraulic conductivity over delz*  (1/s)
!>    @details
!>    Compute the harmonic mean of Kz on the cell face in positive
!>    z-direction from the current node (i, j, k) divided by the
!>    z-distance of the current node (i, j, k) to the neighboring node
!>    (i, j, k+1).
!>
!>    delz* = 0.5 ( delz(k) + delz(k+1) )
!>
!>    Kz / delz* = [ 0.5*delz(k+1)/K(k+1) + 0.5*delz(k)/K(k) ]**-1
!>
!>               = [ ( 0.5* K(k)* delz(k+1) + 0.5*K(k+1)* delz(k) ) / ( K(k)*K(k+1) ) ]**-1
!>
!>               = [ ( 0.5 * summ ) / (prod) ]**-1 = 2.0*prod/summ
      DOUBLE PRECISION FUNCTION tk(i,j,k,ismpl)
        use arrays
        use mod_genrl
        use mod_flow
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k

        DOUBLE PRECISION f1, f2, prod, summ, kz, rhof, visf
        EXTERNAL kz, rhof, visf

        tk = 0.D0
        IF (k0>1 .AND. k<k0) THEN
          f1 = kz(i,j,k,ismpl)*rhof(i,j,k,ismpl)*grav/ &
            visf(i,j,k,ismpl)
          f2 = kz(i,j,k+1,ismpl)*rhof(i,j,k+1,ismpl)*grav/ &
            visf(i,j,k+1,ismpl)
          prod = f1*f2
          summ = f1*delz(k+1) + f2*delz(k)
          IF (summ>0.D0) tk = 2.D0*prod/summ
        END IF
        RETURN
      END
