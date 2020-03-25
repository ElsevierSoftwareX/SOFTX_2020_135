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

!>    @brief calculate x heat flux at cell  faces
!>    @param[in] i grid indices
!>    @param[in] j grid indices
!>    @param[in] k grid indices
!>    @param[in] ismpl local sample index
!>    @return x heat flux (W/m^2)
      DOUBLE PRECISION FUNCTION qx(i,j,k,ismpl)
        use arrays
        use mod_genrl
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k
        DOUBLE PRECISION dif, li
        EXTERNAL li

        qx = 0.D0
        IF (i0>1 .AND. i<i0) THEN
          dif = temp(i+1,j,k,ismpl) - temp(i,j,k,ismpl)
          qx = -li(i,j,k,ismpl)*dif
        END IF
        RETURN
      END

!>    @brief calculate y heat flux at cell  faces
!>    @param[in] i grid indices
!>    @param[in] j grid indices
!>    @param[in] k grid indices
!>    @param[in] ismpl local sample index
!>    @return y heat flux (W/m^2)
      DOUBLE PRECISION FUNCTION qy(i,j,k,ismpl)
        use arrays
        use mod_genrl
        use mod_temp
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k

        DOUBLE PRECISION dif, lj
        EXTERNAL lj

        qy = 0.D0
        IF (j0>1 .AND. j<j0) THEN
          dif = temp(i,j+1,k,ismpl) - temp(i,j,k,ismpl)
          qy = -lj(i,j,k,ismpl)*dif
        END IF

        RETURN
      END

!>    @brief calculate z heat flux at cell  faces
!>    @param[in] i grid indices
!>    @param[in] j grid indices
!>    @param[in] k grid indices
!>    @param[in] ismpl local sample index
!>    @return z heat flux (W/m^2)
      DOUBLE PRECISION FUNCTION qz(i,j,k,ismpl)
        use arrays
        use mod_genrl
        use mod_temp
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k

        DOUBLE PRECISION dif, lk
        EXTERNAL lk

        qz = 0.D0
        IF (k0>1 .AND. k<k0) THEN
          dif = temp(i,j,k+1,ismpl) - temp(i,j,k,ismpl)
          qz = -lk(i,j,k,ismpl)*dif
        END IF
        RETURN
      END

!>    @brief calculate x heat flux at cell centers
!>    @param[in] i grid indices
!>    @param[in] j grid indices
!>    @param[in] k grid indices
!>    @param[in] ismpl local sample index
!>    @return x heat flux (W/m^2)
      DOUBLE PRECISION FUNCTION qxc(i,j,k,ismpl)
        use arrays
        use mod_genrl
        IMPLICIT NONE
        DOUBLE PRECISION d1, d2, li, amean
        integer :: ismpl
        integer :: i, j, k
        EXTERNAL li, amean

        qxc = 0.D0
        IF (i0<=1) RETURN
        IF (i>1 .AND. i<i0) THEN
          d1 = temp(i+1,j,k,ismpl) - temp(i,j,k,ismpl)
          d2 = temp(i,j,k,ismpl) - temp(i-1,j,k,ismpl)
          qxc = amean(-li(i,j,k,ismpl)*d1,-li(i-1,j,k,ismpl)*d2)
        ELSE IF (i==1) THEN
          qxc = -li(i,j,k,ismpl)*(temp(i+1,j,k,ismpl)-temp(i,j,k,ismpl &
            ))
        ELSE IF (i==i0) THEN
          qxc = -li(i-1,j,k,ismpl)*(temp(i,j,k,ismpl)-temp(i-1,j,k, &
            ismpl))
        END IF
        RETURN
      END

!>    @brief calculate y heat fluxat cell center
!>    @param[in] i grid indices
!>    @param[in] j grid indices
!>    @param[in] k grid indices
!>    @param[in] ismpl local sample index
!>    @return y heat flux (W/m^2)
      DOUBLE PRECISION FUNCTION qyc(i,j,k,ismpl)
        use arrays
        use mod_genrl
        use mod_temp
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k

        DOUBLE PRECISION d1, d2, lj, amean
        EXTERNAL lj, amean

        qyc = 0.D0
        IF (j0<=1) RETURN
        IF (j>1 .AND. j<j0) THEN
          d1 = temp(i,j+1,k,ismpl) - temp(i,j,k,ismpl)
          d2 = temp(i,j,k,ismpl) - temp(i,j-1,k,ismpl)
          qyc = amean(-lj(i,j,k,ismpl)*d1,-lj(i,j-1,k,ismpl)*d2)
        ELSE IF (j==1) THEN
          qyc = -lj(i,j,k,ismpl)*(temp(i,j+1,k,ismpl)-temp(i,j,k,ismpl &
            ))
        ELSE IF (j==j0) THEN
          qyc = -lj(i,j-1,k,ismpl)*(temp(i,j,k,ismpl)-temp(i,j-1,k, &
            ismpl))
        END IF
        RETURN
      END

!>    @brief calculate z heat flux at cell center
!>    @param[in] i grid indices
!>    @param[in] j grid indices
!>    @param[in] k grid indices
!>    @param[in] ismpl local sample index
!>    @return z heat flux (W/m^2)
      DOUBLE PRECISION FUNCTION qzc(i,j,k,ismpl)
        use arrays
        use mod_genrl
        use mod_temp
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k

        DOUBLE PRECISION d1, d2, lk, amean
        EXTERNAL lk, amean

        qzc = 0.D0
        IF (k0<=1) RETURN
        IF (k>1 .AND. k<k0) THEN
          d1 = temp(i,j,k+1,ismpl) - temp(i,j,k,ismpl)
          d2 = temp(i,j,k,ismpl) - temp(i,j,k-1,ismpl)
          qzc = amean(-lk(i,j,k,ismpl)*d1,-lk(i,j,k-1,ismpl)*d2)
        ELSE IF (k==1) THEN
          qzc = -lk(i,j,k,ismpl)*(temp(i,j,k+1,ismpl)-temp(i,j,k,ismpl &
            ))
        ELSE IF (k==k0) THEN
          qzc = -lk(i,j,k-1,ismpl)*(temp(i,j,k,ismpl)-temp(i,j,k-1, &
            ismpl))
        END IF
        RETURN
      END

!>    @brief average thermal conductivities on cell faces in x direction
!>    @param[in] i grid indices
!>    @param[in] j grid indices
!>    @param[in] k grid indices
!>    @param[in] ismpl local sample index
!>    @return x thermal conductivity (J/mK)
      DOUBLE PRECISION FUNCTION li(i,j,k,ismpl)
        use arrays
        use mod_genrl
        use mod_temp
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k

        DOUBLE PRECISION f1, f2, prod, summ, lx
        EXTERNAL lx

        li = 0.D0
        IF (i0>1 .AND. i<i0) THEN
          f1 = lx(i,j,k,ismpl)
          f2 = lx(i+1,j,k,ismpl)
          prod = f1*f2
          summ = f1*delx(i+1) + f2*delx(i)
          IF (summ>0.D0) li = 2.D0*prod/summ
        END IF
!       write(99,'(a,5G14.5)') 'in li ',f1,f2,prod,summ,li
        RETURN
      END

!>    @brief average thermal conductivities on cell faces in y direction
!>    @param[in] i grid indices
!>    @param[in] j grid indices
!>    @param[in] k grid indices
!>    @param[in] ismpl local sample index
!>    @return y thermal conductivity (J/mK)
      DOUBLE PRECISION FUNCTION lj(i,j,k,ismpl)
        use arrays
        use mod_genrl
        use mod_temp
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k

        DOUBLE PRECISION f1, f2, prod, summ, ly
        EXTERNAL ly

        lj = 0.D0
        IF (j0>1 .AND. j<j0) THEN
          f1 = ly(i,j,k,ismpl)
          f2 = ly(i,j+1,k,ismpl)
          prod = f1*f2
          summ = f1*dely(j+1) + f2*dely(j)
          IF (summ>0.D0) lj = 2.D0*prod/summ
        END IF
        RETURN
      END

!>    @brief average thermal conductivities on cell faces in z direction
!>    @param[in] i grid indices
!>    @param[in] j grid indices
!>    @param[in] k grid indices
!>    @param[in] ismpl local sample index
!>    @return z thermal conductivity (J/mK,ismpl)
      DOUBLE PRECISION FUNCTION lk(i,j,k,ismpl)
        use arrays
        use mod_genrl
        use mod_temp
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k

        DOUBLE PRECISION f1, f2, prod, summ, lz
        EXTERNAL lz

        lk = 0.D0
        IF (k0>1 .AND. k<k0) THEN
          f1 = lz(i,j,k,ismpl)
          f2 = lz(i,j,k+1,ismpl)
          prod = f1*f2
          summ = f1*delz(k+1) + f2*delz(k)
          IF (summ>0.D0) lk = 2.D0*prod/summ
        END IF
        RETURN
      END

!>    @brief calculate advective heat flux at cell faces
!>    @param[in] i grid indices
!>    @param[in] j grid indices
!>    @param[in] k grid indices
!>    @param[in] ismpl local sample index
!>    @return x advective heat flux (W/m^2)
      DOUBLE PRECISION FUNCTION qvx(i,j,k,ismpl)
        use arrays
        use mod_genrl
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k
        DOUBLE PRECISION tcf, vx, rhof, cpf
        EXTERNAL vx, rhof, cpf

        qvx = 0.D0
        IF (i0>1 .AND. i<i0) THEN
          tcf = (delx(i+1)*temp(i+1,j,k,ismpl)*rhof(i+1,j,k,ismpl)* &
            cpf(i+1,j,k,ismpl)+delx(i)*temp(i,j,k,ismpl)*rhof(i,j,k, &
            ismpl)*cpf(i,j,k,ismpl))/(delx(i+1)+delx(i))
          qvx = delz(k)*dely(j)*tcf*vx(i,j,k,ismpl)
        END IF
        RETURN
      END

!>    @brief calculate advective heat flux at cell faces
!>    @param[in] i grid indices
!>    @param[in] j grid indices
!>    @param[in] k grid indices
!>    @param[in] ismpl local sample index
!>    @return y advective heat flux (W/m^2)
      DOUBLE PRECISION FUNCTION qvy(i,j,k,ismpl)
        use arrays
        use mod_genrl
        use mod_temp
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k

        DOUBLE PRECISION tcf, vy, rhof, cpf
        EXTERNAL vy, rhof, cpf

        qvy = 0.D0
        IF (j0>1 .AND. j<j0) THEN
          tcf = (dely(j+1)*temp(i,j+1,k,ismpl)*rhof(i,j+1,k,ismpl)* &
            cpf(i,j+1,k,ismpl)+dely(j)*temp(i,j,k,ismpl)*rhof(i,j,k, &
            ismpl)*cpf(i,j,k,ismpl))/(dely(j+1)+dely(j))
          qvy = delx(i)*delz(k)*tcf*vy(i,j,k,ismpl)
        END IF
        RETURN
      END

!>    @brief calculate advective heat flux at cell faces
!>    @param[in] i grid indices
!>    @param[in] j grid indices
!>    @param[in] k grid indices
!>    @param[in] ismpl local sample index
!>    @return z advective heat flux (W/m^2)
      DOUBLE PRECISION FUNCTION qvz(i,j,k,ismpl)
        use arrays
        use mod_genrl
        use mod_temp
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k

        DOUBLE PRECISION tcf, vz, rhof, cpf
        EXTERNAL vz, rhof, cpf

        qvz = 0.D0
        IF (k0>1 .AND. k<k0) THEN
          tcf = (delz(k+1)*temp(i,j,k+1,ismpl)*rhof(i,j,k+1,ismpl)* &
            cpf(i,j,k+1,ismpl)+delz(k)*temp(i,j,k,ismpl)*rhof(i,j,k, &
            ismpl)*cpf(i,j,k,ismpl))/(delz(k+1)+delz(k))
          qvz = delx(i)*dely(j)*tcf*vz(i,j,k,ismpl)
        END IF
        RETURN
      END
