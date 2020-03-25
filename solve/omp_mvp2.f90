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

!>    @brief apply thread local 7point-star matrix multiply [as]:=[M]x[s], (OpenMP version)
!>    @param[in] N_I lengths of I dimension of local matrix [M]
!>    @param[in] N_J lengths of J dimension of local matrix [M]
!>    @param[in] N_K lengths of K dimension of local matrix [M]
!>    @param[in] lxyz_block number of blocks
!>    @param[in] xloc thread local blocks of the vector [s]
!>    @param[out] Mxloc thread local blocks of the vector [as]
!>    @param[in] LMA thread local blocks of the 1. diagonal of the system matrix [M]
!>    @param[in] LMB thread local blocks of the 2. diagonal of the system matrix [M]
!>    @param[in] LMC thread local blocks of the 3. diagonal of the system matrix [M]
!>    @param[in] LMD thread local blocks of the 4. (main) diagonal of the system matrix [M]
!>    @param[in] LME thread local blocks of the 5. diagonal of the system matrix [M]
!>    @param[in] LMF thread local blocks of the 6. diagonal of the system matrix [M]
!>    @param[in] LMG thread local blocks of the 7. diagonal of the system matrix [M]
!>    @param[in] xyz_block block dimensions
!>    @param[in] bound_block boundary exchange buffer for each block, between the threads
!>    @details
!>    OpenMP parallelised with special blocking\n
      SUBROUTINE omp_mvp2(n_i,n_j,n_k,lxyz_block,xloc,mxloc,lma,lmb, &
          lmc,lmd,lme,lmf,lmg,bound_block,xyz_block)
        use mod_OMP_TOOLS
        use mod_blocking_size
        IMPLICIT NONE
        INCLUDE 'OMP_TOOLS.inc'

!     N : length of all vectors r,z
        INTEGER n_i, n_j, n_k
        INTEGER k_e, j_e, i_e
        INTEGER ii
        INTEGER xi, yi, zi
        INTEGER xip1, yip1, zip1, xim1, yim1, zim1

!      use mod_blocking_size

        INTEGER lxyz_block, xyz_block(3,lxyz_block)
        DOUBLE PRECISION lme(block_i*block_j*block_k,lxyz_block)
        DOUBLE PRECISION lmf(block_i*block_j*block_k,lxyz_block)
        DOUBLE PRECISION lmg(block_i*block_j*block_k,lxyz_block)
        DOUBLE PRECISION lmb(block_i*block_j*block_k,lxyz_block)
        DOUBLE PRECISION lmc(block_i*block_j*block_k,lxyz_block)
        DOUBLE PRECISION lmd(block_i*block_j*block_k,lxyz_block)
        DOUBLE PRECISION lma(block_i*block_j*block_k,lxyz_block)
        DOUBLE PRECISION xloc(block_i*block_j*block_k,lxyz_block)
        DOUBLE PRECISION mxloc(block_i*block_j*block_k,lxyz_block)
!     global boundary buffer [x]
        DOUBLE PRECISION bound_block(block_i*block_j+block_i*block_k+ &
          block_j*block_k,bdim_i,bdim_j,bdim_k,2)

        INTEGER isurf, ijsurf, ijksurf


!     surface offset: 1=I-dim, isurf=J-dim, ijsurf=K-dim
        isurf = block_j*block_k + 1
        ijsurf = isurf + block_i*block_k
        ijksurf = block_j*block_k + block_i*block_k + block_i*block_j
!     -----------------------------------------------------------------

!     exchange boundaries of all blocks
        DO ii = 1, lxyz_block
          xi = xyz_block(1,ii)
          yi = xyz_block(2,ii)
          zi = xyz_block(3,ii)
!        compute begin & end
          i_e = min(block_i,n_i-((xi-1)*block_i))
          j_e = min(block_j,n_j-((yi-1)*block_j))
          k_e = min(block_k,n_k-((zi-1)*block_k))
!        block index dim-1
          xim1 = max(xi-1,1)
          yim1 = max(yi-1,1)
          zim1 = max(zi-1,1)
!        block index dim+1
          xip1 = min(xi+1,bdim_i)
          yip1 = min(yi+1,bdim_j)
          zip1 = min(zi+1,bdim_k)

          CALL exchange_xloc(block_i,block_j,block_k,i_e,j_e,k_e, &
            xloc(1,ii),bound_block(1,xim1,yi,zi,2), &
            bound_block(isurf,xi,yim1,zi,2),bound_block(ijsurf,xi,yi, &
            zim1,2),bound_block(1,xip1,yi,zi,1), &
            bound_block(isurf,xi,yip1,zi,1),bound_block(ijsurf,xi,yi, &
            zip1,1),bound_block(1,xi,yi,zi,1), &
            bound_block(isurf,xi,yi,zi,1),bound_block(ijsurf,xi,yi,zi, &
            1),bound_block(1,xi,yi,zi,2),bound_block(isurf,xi,yi,zi,2) &
            ,bound_block(ijsurf,xi,yi,zi,2),xim1-xi+1,yim1-yi+1, &
            zim1-zi+1,xip1-xi-1,yip1-yi-1,zip1-zi-1)
        END DO

!$OMP barrier
!     MVP of all blocks
        DO ii = 1, lxyz_block
          xi = xyz_block(1,ii)
          yi = xyz_block(2,ii)
          zi = xyz_block(3,ii)
!        compute begin & end
          i_e = min(block_i,n_i-((xi-1)*block_i))
          j_e = min(block_j,n_j-((yi-1)*block_j))
          k_e = min(block_k,n_k-((zi-1)*block_k))

          CALL mvp_xloc(block_i,block_j,block_k,i_e,j_e,k_e, &
            xloc(1,ii),mxloc(1,ii),lma(1,ii),lmb(1,ii),lmc(1,ii), &
            lmd(1,ii),lme(1,ii),lmf(1,ii),lmg(1,ii), &
            bound_block(1,xi,yi,zi,1),bound_block(isurf,xi,yi,zi,1), &
            bound_block(ijsurf,xi,yi,zi,1),bound_block(1,xi,yi,zi,2), &
            bound_block(isurf,xi,yi,zi,2),bound_block(ijsurf,xi,yi,zi, &
            2))
        END DO

        RETURN
      END

!>    @brief exchange boundary for the [xloc] block
!>    @param[in] ldI leading dimensions in I direction
!>    @param[in] ldJ leading dimensions in J direction
!>    @param[in] ldK leading dimensions in K direction
!>    @param[in] i_e block size in I direction
!>    @param[in] j_e block size in J direction
!>    @param[in] k_e block size in K direction
!>    @param[out] xloc_im1 global boundary buffer to the -1 neighbour in I direction
!>    @param[out] xloc_jm1 global boundary buffer to the -1 neighbour in J direction
!>    @param[out] xloc_km1 global boundary buffer to the -1 neighbour in K direction
!>    @param[out] xloc_ip1 global boundary buffer to the +1 neighbour in I direction
!>    @param[out] xloc_jp1 global boundary buffer to the +1 neighbour in J direction
!>    @param[out] xloc_kp1 global boundary buffer to the +1 neighbour in K direction
!>    @param[out] sxloc_im1 (self) own boundary buffer from the -1 neighbour in I direction
!>    @param[out] sxloc_jm1 (self) own boundary buffer from the -1 neighbour in J direction
!>    @param[out] sxloc_km1 (self) own boundary buffer from the -1 neighbour in K direction
!>    @param[out] sxloc_ip1 (self) own boundary buffer from the +1 neighbour in I direction
!>    @param[out] sxloc_jp1 (self) own boundary buffer from the +1 neighbour in J direction
!>    @param[out] sxloc_kp1 (self) own boundary buffer from the +1 neighbour in K direction
!>    @param[in] cmi 0: copy boundary to the -1 neighbour in I direction, <>0: clear local buffer
!>    @param[in] cmj 0: copy boundary to the -1 neighbour in J direction, <>0: clear local buffer
!>    @param[in] cmk 0: copy boundary to the -1 neighbour in K direction, <>0: clear local buffer
!>    @param[in] cpi 0: copy boundary to the +1 neighbour in I direction, <>0: clear local buffer
!>    @param[in] cpj 0: copy boundary to the +1 neighbour in J direction, <>0: clear local buffer
!>    @param[in] cpk 0: copy boundary to the +1 neighbour in K direction, <>0: clear local buffer
!>    @param[in] xloc local block (elements)
      SUBROUTINE exchange_xloc(ldi,ldj,ldk,i_e,j_e,k_e,xloc,xloc_im1, &
          xloc_jm1,xloc_km1,xloc_ip1,xloc_jp1,xloc_kp1,sxloc_im1, &
          sxloc_jm1,sxloc_km1,sxloc_ip1,sxloc_jp1,sxloc_kp1,cmi,cmj, &
          cmk,cpi,cpj,cpk)
        use mod_OMP_TOOLS
        use mod_blocking_size
        IMPLICIT NONE
        INCLUDE 'OMP_TOOLS.inc'
!     N : length of all vectors r,z
        INTEGER i_e, j_e, k_e, i, j, k, ldi, ldj, ldk
        INTEGER cmi, cmj, cmk, cpi, cpj, cpk
!     vectors [r], [z]
        DOUBLE PRECISION xloc(ldi,ldj,ldk)

!      use mod_blocking_size
!     global boundary buffer
        DOUBLE PRECISION xloc_im1(block_j,block_k)
        DOUBLE PRECISION xloc_jm1(block_i,block_k)
        DOUBLE PRECISION xloc_km1(block_i,block_j)
        DOUBLE PRECISION xloc_ip1(block_j,block_k)
        DOUBLE PRECISION xloc_jp1(block_i,block_k)
        DOUBLE PRECISION xloc_kp1(block_i,block_j)
!     (self) owner of this bloCk
        DOUBLE PRECISION sxloc_im1(block_j,block_k)
        DOUBLE PRECISION sxloc_jm1(block_i,block_k)
        DOUBLE PRECISION sxloc_km1(block_i,block_j)
        DOUBLE PRECISION sxloc_ip1(block_j,block_k)
        DOUBLE PRECISION sxloc_jp1(block_i,block_k)
        DOUBLE PRECISION sxloc_kp1(block_i,block_j)

!     boundary exchange
!     [2D-I]
        IF (cmi==0) THEN
          DO k = 1, k_e
            DO j = 1, j_e
              xloc_im1(j,k) = xloc(1,j,k)
            END DO
          END DO
        ELSE
!        own buffer
          DO k = 1, k_e
            DO j = 1, j_e
              sxloc_im1(j,k) = 0.0D0
            END DO
          END DO
        END IF
        IF (cpi==0) THEN
          DO k = 1, k_e
            DO j = 1, j_e
              xloc_ip1(j,k) = xloc(i_e,j,k)
            END DO
          END DO
        ELSE
!        own buffer
          DO k = 1, k_e
            DO j = 1, j_e
              sxloc_ip1(j,k) = 0.0D0
            END DO
          END DO
        END IF

!     [2D-J]
        IF (cmj==0) THEN
          DO k = 1, k_e
            DO i = 1, i_e
              xloc_jm1(i,k) = xloc(i,1,k)
            END DO
          END DO
        ELSE
!        own buffer
          DO k = 1, k_e
            DO i = 1, i_e
              sxloc_jm1(i,k) = 0.0D0
            END DO
          END DO
        END IF
        IF (cpj==0) THEN
          DO k = 1, k_e
            DO i = 1, i_e
              xloc_jp1(i,k) = xloc(i,j_e,k)
            END DO
          END DO
        ELSE
!        own buffer
          DO k = 1, k_e
            DO i = 1, i_e
              sxloc_jp1(i,k) = 0.0D0
            END DO
          END DO
        END IF

!     [2D-K]
        IF (cmk==0) THEN
          DO j = 1, j_e
            DO i = 1, i_e
              xloc_km1(i,j) = xloc(i,j,1)
            END DO
          END DO
        ELSE
!        own buffer
          DO j = 1, j_e
            DO i = 1, i_e
              sxloc_km1(i,j) = 0.0D0
            END DO
          END DO
        END IF
        IF (cpk==0) THEN
          DO j = 1, j_e
            DO i = 1, i_e
              xloc_kp1(i,j) = xloc(i,j,k_e)
            END DO
          END DO
        ELSE
!        own buffer
          DO j = 1, j_e
            DO i = 1, i_e
              sxloc_kp1(i,j) = 0.0D0
            END DO
          END DO
        END IF

        RETURN
      END

!>    @brief thread local block matrix vector product [as]:=[M]x[s], one single block
!>    @param[in] ldI leading dimensions in I direction
!>    @param[in] ldJ leading dimensions in J direction
!>    @param[in] ldK leading dimensions in K direction
!>    @param[in] i_e block size in I direction
!>    @param[in] j_e block size in J direction
!>    @param[in] k_e block size in K direction
!>    @param[in] xloc local block (elements) - thread local part of vector [s]
!>    @param[out] Mxloc result vector - thread local part of vector [as]
!>    @param[in] MA thread local blocks of the 1. diagonal of the system matrix [M]
!>    @param[in] MB thread local blocks of the 2. diagonal of the system matrix [M]
!>    @param[in] MC thread local blocks of the 3. diagonal of the system matrix [M]
!>    @param[in] MD thread local blocks of the 4. (main) diagonal of the system matrix [M]
!>    @param[in] ME thread local blocks of the 5. diagonal of the system matrix [M]
!>    @param[in] MF thread local blocks of the 6. diagonal of the system matrix [M]
!>    @param[in] MG thread local blocks of the 7. diagonal of the system matrix [M]
!>    @param[in] xloc_im1 global boundary buffer from the -1 neighbour in I direction
!>    @param[in] xloc_jm1 global boundary buffer from the -1 neighbour in J direction
!>    @param[in] xloc_km1 global boundary buffer from the -1 neighbour in K direction
!>    @param[in] xloc_ip1 global boundary buffer from the +1 neighbour in I direction
!>    @param[in] xloc_jp1 global boundary buffer from the +1 neighbour in J direction
!>    @param[in] xloc_kp1 global boundary buffer from the +1 neighbour in K direction
      SUBROUTINE mvp_xloc(ldi,ldj,ldk,i_e,j_e,k_e,xloc,mxloc,ma,mb,mc, &
          md,me,mf,mg,xloc_im1,xloc_jm1,xloc_km1,xloc_ip1,xloc_jp1, &
          xloc_kp1)
        use mod_OMP_TOOLS
        use mod_blocking_size
        IMPLICIT NONE
        INCLUDE 'OMP_TOOLS.inc'

!     N : length of all vectors r,z
        INTEGER i_e, j_e, k_e, i, j, k, ldi, ldj, ldk
!     vectors [r], [z]
        DOUBLE PRECISION xloc(ldi,ldj,ldk), mxloc(ldi,ldj,ldk)
        DOUBLE PRECISION ma(ldi,ldj,ldk), mb(ldi,ldj,ldk)
        DOUBLE PRECISION mc(ldi,ldj,ldk), md(ldi,ldj,ldk)
        DOUBLE PRECISION me(ldi,ldj,ldk), mf(ldi,ldj,ldk)
        DOUBLE PRECISION mg(ldi,ldj,ldk)

!      use mod_blocking_size
!     global boundary buffer
        DOUBLE PRECISION xloc_im1(block_j,block_k)
        DOUBLE PRECISION xloc_jm1(block_i,block_k)
        DOUBLE PRECISION xloc_km1(block_i,block_j)
        DOUBLE PRECISION xloc_ip1(block_j,block_k)
        DOUBLE PRECISION xloc_jp1(block_i,block_k)
        DOUBLE PRECISION xloc_kp1(block_i,block_j)

!     [I] direction
        i = 1
        DO k = 1, k_e
          DO j = 1, j_e
            mxloc(i,j,k) = md(i,j,k)*xloc(i,j,k) + &
              mc(i,j,k)*xloc_im1(j,k)
          END DO
        END DO
        DO k = 1, k_e
          DO j = 1, j_e
            DO i = 2, i_e
              mxloc(i,j,k) = md(i,j,k)*xloc(i,j,k) + &
                mc(i,j,k)*xloc(i-1,j,k)
            END DO
          END DO
        END DO
        i = i_e
        DO k = 1, k_e
          DO j = 1, j_e
            mxloc(i,j,k) = mxloc(i,j,k) + me(i,j,k)*xloc_ip1(j,k)
          END DO
        END DO
        DO k = 1, k_e
          DO j = 1, j_e
            DO i = 1, i_e - 1
              mxloc(i,j,k) = mxloc(i,j,k) + me(i,j,k)*xloc(i+1,j,k)
            END DO
          END DO
        END DO

!     [J] direction
        j = 1
        DO k = 1, k_e
          DO i = 1, i_e
            mxloc(i,j,k) = mxloc(i,j,k) + mb(i,j,k)*xloc_jm1(i,k)
          END DO
        END DO
        DO k = 1, k_e
          DO j = 2, j_e
            DO i = 1, i_e
              mxloc(i,j,k) = mxloc(i,j,k) + mb(i,j,k)*xloc(i,j-1,k)
            END DO
          END DO
        END DO
        j = j_e
        DO k = 1, k_e
          DO i = 1, i_e
            mxloc(i,j,k) = mxloc(i,j,k) + mf(i,j,k)*xloc_jp1(i,k)
          END DO
        END DO
        DO k = 1, k_e
          DO j = 1, j_e - 1
            DO i = 1, i_e
              mxloc(i,j,k) = mxloc(i,j,k) + mf(i,j,k)*xloc(i,j+1,k)
            END DO
          END DO
        END DO

!     [K] direction
        k = 1
        DO j = 1, j_e
          DO i = 1, i_e
            mxloc(i,j,k) = mxloc(i,j,k) + ma(i,j,k)*xloc_km1(i,j)
          END DO
        END DO
        DO k = 2, k_e
          DO j = 1, j_e
            DO i = 1, i_e
              mxloc(i,j,k) = mxloc(i,j,k) + ma(i,j,k)*xloc(i,j,k-1)
            END DO
          END DO
        END DO
        k = k_e
        DO j = 1, j_e
          DO i = 1, i_e
            mxloc(i,j,k) = mxloc(i,j,k) + mg(i,j,k)*xloc_kp1(i,j)
          END DO
        END DO
        DO k = 1, k_e - 1
          DO j = 1, j_e
            DO i = 1, i_e
              mxloc(i,j,k) = mxloc(i,j,k) + mg(i,j,k)*xloc(i,j,k+1)
            END DO
          END DO
        END DO

        RETURN
      END
