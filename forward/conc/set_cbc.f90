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

!>    @brief modify coefficents for the species equation according to the boundary
!>    @param[in] spec species index
!>    @param[in] ismpl local sample index
!>    @details
!> modify coefficents for the species equation according to the boundary conditions\n
!> coefficients are stored as vectors in the diagonals a-g (d center) and rhs in w.\n
      SUBROUTINE set_cbc(spec,ismpl)
        use arrays
        use mod_genrl
        use mod_genrlc
        use mod_conc
        use mod_time
        use mod_linfos
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k
        integer :: ib
        INTEGER bcu, tpbcu, bctype, i_dir, spec, i_spec
        DOUBLE PRECISION val, malfa, mbeta, vx, vy, vz, dv, ds, vv
        EXTERNAL vx, vy, vz
        INTRINSIC max


! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  flow due to flow neumann nodes / wellars  - - - - - - - - - - - - - -

        DO ib = first_flow, last_flow
          i = ibc_data(ib,cbc_i)
          j = ibc_data(ib,cbc_j)
          k = ibc_data(ib,cbc_k)
          bcu = ibc_data(ib,cbc_bcu)
          tpbcu = max(ibc_data(ib,cbc_bctp),0)
          bctype = ibc_data(ib,cbc_bt)
!        "neumann"?, skip otherwise
          IF (bctype==bt_neum.OR.bctype==bt_neuw) THEN
!           discrete values
            IF (bcu<=0) THEN
              val = dbc_data(ib,1,ismpl)
            ELSE
              val = propunit(bcu,idx_hbc,ismpl)
            END IF

            IF ((tpbcu>0) .AND. nbctp>0) THEN
!               time-dependent bc:  val=ac*val+bc
!               get Alfa and Beta modificators
              CALL get_tpbcalbe(malfa,mbeta,tpbcu,ismpl)
!               update time dependend modification of the bc-value
              val = malfa + mbeta*val
            END IF

!           wellar test
            IF (val<0.0D0 .AND. tpbcu>=0) THEN
              ds = 0.D0
              dv = 0.D0
              IF (i>1) THEN
                vv = abs(vx(i-1,j,k,ismpl))
                ds = ds + vv
                dv = dv + conc(i-1,j,k,spec,ismpl)*vv
              END IF
              IF (i<i0) THEN
                vv = abs(vx(i,j,k,ismpl))
                ds = ds + vv
                dv = dv + conc(i+1,j,k,spec,ismpl)*vv
              END IF
              IF (j>1) THEN
                vv = abs(vy(i,j-1,k,ismpl))
                ds = ds + vv
                dv = dv + conc(i,j-1,k,spec,ismpl)*vv
              END IF
              IF (j<j0) THEN
                vv = abs(vy(i,j,k,ismpl))
                ds = ds + vv
                dv = dv + conc(i,j+1,k,spec,ismpl)*vv
              END IF
              IF (k>1) THEN
                vv = abs(vz(i,j,k-1,ismpl))
                ds = ds + vv
                dv = dv + conc(i,j,k-1,spec,ismpl)*vv
              END IF
              IF (k<k0) THEN
                vv = abs(vz(i,j,k,ismpl))
                ds = ds + vv
                dv = dv + conc(i,j,k+1,spec,ismpl)*vv
              END IF
              dv = dv/ds
!              apply dirichlet update [dv]
#ifdef BCMY
!              D = D+my
              d(i,j,k,ismpl) = d(i,j,k,ismpl) - dbc_data(ib,2,ismpl)
              w(i,j,k,ismpl) = w(i,j,k,ismpl) - &
                dbc_data(ib,2,ismpl)*dv
#else
!              standard boundary condition handling
              a(i,j,k,ismpl) = 0.0D0
              b(i,j,k,ismpl) = 0.0D0
              c(i,j,k,ismpl) = 0.0D0
              e(i,j,k,ismpl) = 0.0D0
              f(i,j,k,ismpl) = 0.0D0
              g(i,j,k,ismpl) = 0.0D0
              d(i,j,k,ismpl) = 1.0D0
              w(i,j,k,ismpl) = dv
              conc(i,j,k,spec,ismpl) = dv
!              mark as boundary for normalising the lin. system
              bc_mask(i+(j-1)*i0+(k-1)*i0*j0,ismpl) = '0'
#endif
            END IF
          END IF
        END DO

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! dirichlet nodes - - - - - - - - - - - - - - - - - - - - - - - - - - -

        DO ib = first_conc, last_conc
          i = ibc_data(ib,cbc_i)
          j = ibc_data(ib,cbc_j)
          k = ibc_data(ib,cbc_k)
          bcu = ibc_data(ib,cbc_bcu)
          tpbcu = max(ibc_data(ib,cbc_bctp),0)
          bctype = ibc_data(ib,cbc_bt)
          i_spec = ibc_data(ib,cbc_si)
!        "dirichlet" and "spec"?, skip otherwise
          IF (bctype==bt_diri .AND. i_spec==spec) THEN
!           discrete values
            IF (bcu<=0) THEN
              val = dbc_data(ib,1,ismpl)
            ELSE
              val = propunit(bcu,idx_cbc,ismpl)
            END IF

            IF ((tpbcu>0) .AND. nbctp>0) THEN
!               time-dependent bc:  val=ac*val+bc
!               get Alfa and Beta modificators
              CALL get_tpbcalbe(malfa,mbeta,tpbcu,ismpl)
!               update time dependend modification of the bc-value
              val = malfa + mbeta*val
            END IF

            IF (tpbcu>=0) THEN
#ifdef BCMY
!              D = D+my
              d(i,j,k,ismpl) = d(i,j,k,ismpl) - dbc_data(ib,2,ismpl)
              w(i,j,k,ismpl) = w(i,j,k,ismpl) - &
                dbc_data(ib,2,ismpl)*val
#else
!              standard boundary condition handling
              a(i,j,k,ismpl) = 0.0D0
              b(i,j,k,ismpl) = 0.0D0
              c(i,j,k,ismpl) = 0.0D0
              e(i,j,k,ismpl) = 0.0D0
              f(i,j,k,ismpl) = 0.0D0
              g(i,j,k,ismpl) = 0.0D0
              d(i,j,k,ismpl) = 1.0D0
              w(i,j,k,ismpl) = val
              conc(i,j,k,spec,ismpl) = val
!              mark as boundary for normalising the lin. system
              bc_mask(i+(j-1)*i0+(k-1)*i0*j0,ismpl) = '0'
#endif
            END IF
          END IF
        END DO

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! neumann  nodes    - - - - - - - - - - - - - - - - - - - - - - - - - -

        DO ib = first_conc, last_conc
          i = ibc_data(ib,cbc_i)
          j = ibc_data(ib,cbc_j)
          k = ibc_data(ib,cbc_k)
          bcu = ibc_data(ib,cbc_bcu)
          tpbcu = max(ibc_data(ib,cbc_bctp),0)
          bctype = ibc_data(ib,cbc_bt)
          i_spec = ibc_data(ib,cbc_si)
          i_dir = ibc_data(ib,cbc_dir)
!        "neumann" and "spec"?, skip otherwise
          IF (bctype==bt_neum .AND. i_spec==spec) THEN
!           discrete values
            IF (bcu<=0) THEN
              val = dbc_data(ib,1,ismpl)
            ELSE
              val = propunit(bcu,idx_cbc,ismpl)
            END IF

            IF ((tpbcu>0) .AND. nbctp>0) THEN
!               time-dependent bc:  val=ac*val+bc
!               get Alfa and Beta modificators
              CALL get_tpbcalbe(malfa,mbeta,tpbcu,ismpl)
!               update time dependend modification of the bc-value
              val = malfa + mbeta*val
            END IF

            IF (tpbcu>=0) THEN
              IF ((i_dir==0)) val = val/(delx(i)*dely(j)*delz(k))
              IF ((i_dir==1) .OR. (i_dir==2)) val = val/delx(i)
              IF ((i_dir==3) .OR. (i_dir==4)) val = val/dely(j)
              IF ((i_dir==5) .OR. (i_dir==6)) val = val/delz(k)

              w(i,j,k,ismpl) = w(i,j,k,ismpl) - val
            END IF
          END IF
        END DO

        RETURN
      END
