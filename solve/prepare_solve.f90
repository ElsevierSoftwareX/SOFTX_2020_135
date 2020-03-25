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

!>    @brief setup ILU(0) preconditioner diagonal [UD]
!>    @param[in] NI lengths of I dimension of local matrix [M]
!>    @param[in] NJ lengths of J dimension of local matrix [M]
!>    @param[in] NK lengths of K dimension of local matrix [M]
!>    @param[in] LA 1. diagonal of the system matrix [M]
!>    @param[in] LB 2. diagonal of the system matrix [M]
!>    @param[in] LC 3. diagonal of the system matrix [M]
!>    @param[in] LD 4. diagonal of the system matrix [M]
!>    @param[in] LE 5. diagonal of the system matrix [M]
!>    @param[in] LF 6. diagonal of the system matrix [M]
!>    @param[in] LG 7. diagonal of the system matrix [M]
!>    @param[out] LUD helper diagonal elements for preconditioning, vector [UD]
      SUBROUTINE prepare_ilu(ni,nj,nk,la,lb,lc,ld,le,lf,lg,lud)
        use mod_linfos
        IMPLICIT NONE
        INTEGER ni, nj, nk
        INTEGER i, j, k
        DOUBLE PRECISION la(ni,nj,nk), lb(ni,nj,nk), lc(ni,nj,nk)
        DOUBLE PRECISION ld(ni,nj,nk), le(ni,nj,nk), lf(ni,nj,nk)
        DOUBLE PRECISION lg(ni,nj,nk)
        DOUBLE PRECISION lud(ni,nj,nk), tmp, tmp_ud, tmp_d, tmp_bug
        INTEGER tnul, tnul_
        LOGICAL test_null, lintel_dummy
        EXTERNAL test_null


        tmp_bug = 0.D0
!$OMP master
        CALL dcopy(ni*nj*nk,ld,1,lud,1)

        DO k = 1, nk
          DO j = 1, nj
            DO i = 1, ni
              tmp_ud = lud(i,j,k)
              tnul = 0
              CALL test_zero(tmp_ud,1,tnul)
              tmp_d = ld(i,j,k)
              tnul_ = 0
              CALL test_zero(tmp_d,1,tnul_)

              IF ((tnul_/=1) .AND. (tnul/=1)) THEN
                tmp = 1.0D0/tmp_ud

                IF (i<ni) THEN
!buggy          lUD(i+1,j,k) = lUD(i+1,j,k) -lC(i+1,j,k)*tmp*lE(i,j,k)
                  tmp_bug = lc(i+1,j,k)*tmp
                  lud(i+1,j,k) = lud(i+1,j,k) - tmp_bug*le(i,j,k)
                END IF
                IF (j<nj) THEN
!buggy          lUD(i,j+1,k) = lUD(i,j+1,k) -lB(i,j+1,k)*tmp*lF(i,j,k)
                  tmp_bug = lb(i,j+1,k)*tmp
                  lud(i,j+1,k) = lud(i,j+1,k) - tmp_bug*lf(i,j,k)
                END IF
                IF (k<nk) THEN
!buggy          lUD(i,j,k+1) = lUD(i,j,k+1) -lA(i,j,k+1)*tmp*lG(i,j,k)
                  tmp_bug = la(i,j,k+1)*tmp
                  lud(i,j,k+1) = lud(i,j,k+1) - tmp_bug*lg(i,j,k)
                END IF
              ELSE
                WRITE(*,'(A,3I5,A,2e15.8)') &
                  'error in prepare_solve.f90: main diagonal element equal to zero at ', i, &
                  j, k, '=', ld(i,j,k), lud(i,j,k)
                STOP
              END IF
            END DO
          END DO
        END DO
!$OMP end master
!$OMP barrier

!     break intel-compiler optimisation,
!       to avoid numerical instabilities for
!         lUD(i+1,j,k) = lUD(i+1,j,k) - [ lC(i+1,j,k)*tmp ] *lE(i,j,k)
!         lUD(i,j+1,k) = lUD(i,j+1,k) - [ lB(i,j+1,k)*tmp ] *lF(i,j,k)
!         lUD(i,j,k+1) = lUD(i,j,k+1) - [ lA(i,j,k+1)*tmp ] *lG(i,j,k)
        lintel_dummy = test_null(tmp_bug)

!     faster ...
!$OMP do schedule(static)
        DO k = 1, nk
          DO j = 1, nj
            DO i = 1, ni
!           for UD, an near NULL criteria is saver
              tmp_ud = lud(i,j,k)
              tnul = 0
              CALL test_zero(tmp_ud,1,tnul)
              IF ((ld(i,j,k)/=0.0D0) .AND. (tnul/=1)) THEN
                lud(i,j,k) = 1.0D0/tmp_ud
              END IF
            END DO
          END DO
        END DO
!$OMP end do nowait

        RETURN
      END

!>    @brief copy "private", distribute the "shared" [full] vector into the [local] "private" block area
!>    @param[in] NI lengths of I dimension of local matrix [M]
!>    @param[in] NJ lengths of J dimension of local matrix [M]
!>    @param[in] NK lengths of K dimension of local matrix [M]
!>    @param[in] full global "shared" vector [full]
!>    @param[out] local "private" block vector [local]
!>    @param[in] lxyz_block number of blocks
!>    @param[in] xyz_block block dimensions
      SUBROUTINE lcopy_ilu(ni,nj,nk,full,local,lxyz_block,xyz_block)
!        use mod_genrl
        use mod_blocking_size
        IMPLICIT NONE
        integer :: i, j, k
        INTEGER ni,nj,nk,i1, j1, k1, i2, j2, k2, i_max, j_max, k_max, ii
        DOUBLE PRECISION full(ni,nj,nk)
        DOUBLE PRECISION local(block_i,block_j,block_k,*)
        INTEGER lxyz_block, xyz_block(3,lxyz_block)

!     search all blocks - identify private
        DO ii = 1, lxyz_block
          i = xyz_block(1,ii)
          j = xyz_block(2,ii)
          k = xyz_block(3,ii)
!        global block offset (-1)
          i2 = block_i*(i-1)
          j2 = block_j*(j-1)
          k2 = block_k*(k-1)
!        block length
          k_max = min(block_k,nk-k2)
          j_max = min(block_j,nj-j2)
          i_max = min(block_i,ni-i2)
!        clear all (dummy-)values
          DO k1 = k_max + 1, block_k
            DO j1 = 1, block_j
              DO i1 = 1, block_i
!           | + + | | * + |
!           | # # | | # # |
                local(i1,j1,k1,ii) = 0.0D0
              END DO
            END DO
          END DO
          DO k1 = 1, k_max
            DO j1 = j_max + 1, block_j
              DO i1 = 1, block_i
!           | + # | | * # |
!           | # # | | # # |
                local(i1,j1,k1,ii) = 0.0D0
              END DO
            END DO
          END DO
          DO k1 = 1, k_max
            DO j1 = 1, j_max
              DO i1 = i_max + 1, block_i
!           | # # | | * # |
!           | # # | | # # |
                local(i1,j1,k1,ii) = 0.0D0
              END DO
            END DO
          END DO
!        copy private blocks in "local"
          DO k1 = 1, k_max
            DO j1 = 1, j_max
              DO i1 = 1, i_max
                local(i1,j1,k1,ii) = full(i2+i1,j2+j1,k2+k1)
              END DO
            END DO
          END DO
        END DO
!
        RETURN
      END

!>    @brief copy "shared", collect the [local] "private" blocks into the "shared" [full] vector
!>    @param[in] NI lengths of I dimension of local matrix [M]
!>    @param[in] NJ lengths of J dimension of local matrix [M]
!>    @param[in] NK lengths of K dimension of local matrix [M]
!>    @param[out] full global "shared" vector [full]
!>    @param[in] local "private" block vector [local]
!>    @param[in] lxyz_block number of blocks
!>    @param[in] xyz_block block dimensions
      SUBROUTINE lcopy_bak_ilu(ni,nj,nk,full,local,lxyz_block,xyz_block)
!        use mod_genrl
        use mod_blocking_size
        IMPLICIT NONE
        integer :: i, j, k
        INTEGER ni,nj,nk,i1, j1, k1, i2, j2, k2, i_max, j_max, k_max, ii
        DOUBLE PRECISION full(ni,nj,nk)
        DOUBLE PRECISION local(block_i,block_j,block_k,*)
        INTEGER lxyz_block, xyz_block(3,lxyz_block)

!     search all blocks - identify private
        DO ii = 1, lxyz_block
          i = xyz_block(1,ii)
          j = xyz_block(2,ii)
          k = xyz_block(3,ii)
!        global block offset (-1)
          i2 = block_i*(i-1)
          j2 = block_j*(j-1)
          k2 = block_k*(k-1)
!        block length
          k_max = min(block_k,nk-k2)
          j_max = min(block_j,nj-j2)
          i_max = min(block_i,ni-i2)
!        copy bak private blocks in "full"
          DO k1 = 1, k_max
            DO j1 = 1, j_max
              DO i1 = 1, i_max
                full(i2+i1,j2+j1,k2+k1) = local(i1,j1,k1,ii)
              END DO
            END DO
          END DO
        END DO
!$OMP barrier
!
        RETURN
      END

!>    @brief copy surface/boundary (position-1) - make private
!>    @param[in] NI lengths of I dimension of local matrix [M]
!>    @param[in] NJ lengths of J dimension of local matrix [M]
!>    @param[in] NK lengths of K dimension of local matrix [M]
!>    @param[in] ii block index number
!>    @param[in] full global "shared" vector
!>    @param[out] local "private" boundary buffer
!>    @param[in] lxyz_block number of blocks
!>    @param[in] xyz_block block dimensions
      SUBROUTINE lsurf_ilu(ni,nj,nk,ii,full,local,lxyz_block,xyz_block)
!        use mod_genrl
        use mod_blocking_size
        IMPLICIT NONE
        integer :: i, j, k
        INTEGER lxyz_block, xyz_block(3,lxyz_block)
        INTEGER ni,nj,nk,i1, j1, k1, i2, j2, k2, i_max, j_max, k_max, ii, offs
        DOUBLE PRECISION full(ni,nj,nk)
        DOUBLE PRECISION local(block_i*block_j+block_i*block_k+ &
          block_j*block_k)

!     search all blocks - identify private
        i = xyz_block(1,ii)
        j = xyz_block(2,ii)
        k = xyz_block(3,ii)
!     global block offset (-1)
        i2 = block_i*(i-1)
        j2 = block_j*(j-1)
        k2 = block_k*(k-1)
!     block length
        i_max = min(block_i,ni-i2)
        j_max = min(block_j,nj-j2)
        k_max = min(block_k,nk-k2)

!     I-dim
        offs = 0
        IF (i2>=1) THEN
          DO k1 = 1, k_max
            DO j1 = 1, j_max
              local(offs+j1+block_j*(k1-1)) = full(i2,j2+j1,k2+k1)
            END DO
          END DO
        ELSE
          DO k1 = 1, k_max
            DO j1 = 1, j_max
              local(offs+j1+block_j*(k1-1)) = 0.0D0
            END DO
          END DO
        END IF
!     J-dim
        offs = offs + block_j*block_k
        IF (j2>=1) THEN
          DO k1 = 1, k_max
            DO i1 = 1, i_max
              local(offs+i1+block_i*(k1-1)) = full(i2+i1,j2,k2+k1)
            END DO
          END DO
        ELSE
          DO k1 = 1, k_max
            DO i1 = 1, i_max
              local(offs+i1+block_i*(k1-1)) = 0.0D0
            END DO
          END DO
        END IF
!     K-dim
        offs = offs + block_i*block_k
        IF (k2>=1) THEN
          DO j1 = 1, j_max
            DO i1 = 1, i_max
              local(offs+i1+block_i*(j1-1)) = full(i2+i1,j2+j1,k2)
            END DO
          END DO
        ELSE
          DO j1 = 1, j_max
            DO i1 = 1, i_max
              local(offs+i1+block_i*(j1-1)) = 0.0D0
            END DO
          END DO
        END IF
!
        RETURN
      END
