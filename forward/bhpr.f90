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

!> @brief calculate well pressure according to Shu,2005
!> @param[in] ismpl local sample index
!> @param[in] ii i cell-index
!> @param[in] jj i cell-index
!> @param[in] kk i cell-index
!> @details
!> Reference: J. Shu\n COMPARISON OF VARIOUS TECHINQUES FOR COMPUTING
!> WELL INDEX\n Master Thesis, Stanford 2005\n\n
!>
!> modify coefficents for the head equation according to the boundary
!> conditions, coefficients are stored as vectors in the diagonals a-g
!> (d center), and rhs in w.\n
      DOUBLE PRECISION FUNCTION bhpr(ii,jj,kk,ismpl)

        use arrays
        use mod_genrl
        use mod_genrlc
        use mod_temp
        use mod_flow
        use mod_time
        use mod_linfos

        IMPLICIT NONE

        integer :: ismpl
        integer :: i, j, k
        integer :: ib
!
        INTEGER ii, jj, kk, bcu, tpbcu, bctype, i_dir
        DOUBLE PRECISION val, malfa, mbeta
        DOUBLE PRECISION d_x,d_y,d_z,k_x,k_y,k_z,wi_x,wi_y,wi_z,wi_pj
        ! DOUBLE PRECISION l_x,l_y
        DOUBLE PRECISION l_z
        DOUBLE PRECISION r_b,r_w,skin
        ! DOUBLE PRECISION xz,zx,yz,zy
        DOUBLE PRECISION xy,yx 
        DOUBLE PRECISION kx,ky,kz,visf
!
        PARAMETER (skin=0.d0,r_w=0.15)
!
        EXTERNAL kx,ky,kz,visf
        INTRINSIC max


! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! dirichlet nodes - - - - - - - - - - - - - - - - - - - - - - - - - - -
! neumann  nodes    - - - - - - - - - - - - - - - - - - - - - - - - - -
!
        ! default well pressure
        bhpr = 0.0d0

        DO ib = first_flow, last_flow

          i = ibc_data(ib,cbc_i)
          j = ibc_data(ib,cbc_j)
          k = ibc_data(ib,cbc_k)
          bctype = ibc_data(ib,cbc_bt)

          ! neumann bc?, skip otherwise
          ! data point? skip otherwise
          IF (ii==i.AND.jj==j .AND. kk==k .AND. bctype==bt_neuw) THEN

            bcu = ibc_data(ib,cbc_bcu)
            tpbcu = max(ibc_data(ib,cbc_bctp),0)
            i_dir = ibc_data(ib,cbc_dir)

            ! WELLMODEL
            ! discrete values
            val = dbc_data(ib,1,ismpl)

            IF ((tpbcu>0) .AND. nbctp>0) THEN

              ! time-dependent bc:  val=ac*val+bc
              ! get Alfa and Beta modificators
              CALL get_tpbcalbe(malfa,mbeta,tpbcu,ismpl)

              ! update time dependend modification of the bc-value
              val = malfa + mbeta*val

            END IF
!
            IF (tpbcu>=0) THEN

              k_x = kx(i,j,k,ismpl)
              k_y = ky(i,j,k,ismpl)
              k_z = kz(i,j,k,ismpl)

              d_x = delx(i)
              d_y = dely(j)
              d_z = delz(k)

!             l_x=0.d0
!             yz=k_y/k_z
!             zy=k_z/k_y
!             r_b=0.28d0*sqrt((sqrt(yz)*d_z**2+sqrt(zy)*d_y**2))/ &
!                                       (yz**0.25+zy**0.25)
!             wi_x = 2.0d0*pi*sqrt(k_y*k_z)*l_x/(log(r_b/r_w)+skin)
              wi_x = 0.0d0

!             l_y=0.d0
!             xz=k_x/k_z
!             zx=k_z/k_x
!             r_b=0.28d0*sqrt((sqrt(xz)*d_z**2+sqrt(zx)*d_x**2))/ &
!                                       (xz**0.25+zx**0.25)
!             wi_y = 2.0d0*pi*sqrt(k_x*k_z)*l_y/(log(r_b/r_w)+skin)
              wi_y = 0.0d0

              l_z = d_z
              xy = k_x/k_y
              yx = k_y/k_x
              r_b = 0.28d0*sqrt((sqrt(xy)*d_y**2+sqrt(yx)*d_x**2))/ &
                   (xy**0.25+yx**0.25)
              wi_z = 2.0d0*pi*sqrt(k_x*k_y)*l_z/(log(r_b/r_w)+skin)

              wi_pj = sqrt(wi_x**2 + wi_y**2 + wi_z**2)

              bhpr = pres(i,j,k,ismpl)+val*visf(i,j,k,ismpl)/wi_pj

              dbc_data(ib,3,ismpl) = bhpr

            END IF
          END IF
!
        END DO
!
        RETURN
      END
