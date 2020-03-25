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

!>    @brief ILU(0) preconditioner, solving block-based [L][U] x [x] = [b], (OpenMP version)
!>    @param[in] N_I lengths of I dimension of local matrix [M]~[L][U]
!>    @param[in] N_J lengths of J dimension of local matrix [M]~[L][U]
!>    @param[in] N_K lengths of K dimension of local matrix [M]~[L][U]
!>    @param[in] lxyz_block number of blocks
!>    @param[in] bloc thread local blocks of vector [b]
!>    @param[out] xloc thread local blocks of vector [x]
!>    @param[in] lMA thread local blocks of the 1. diagonal of [L]
!>    @param[in] lMB thread local blocks of the 2. diagonal of [L]
!>    @param[in] lMC thread local blocks of the 3. diagonal of [L]
!>    @param[in] lUD thread local blocks of the helper diagonal [UD] of [L] and [U]
!>    @param[in] lME thread local blocks of the 2. diagonal of [U]
!>    @param[in] lMF thread local blocks of the 3. diagonal of [U]
!>    @param[in] lMG thread local blocks of the 4. diagonal of [U]
!>    @param[out] ud_block block buffer for helper diagonal [UD]
!>    @param[out] bound_block boundary exchange buffer for each block, between the threads
!>    @param[in] xyz_block block dimensions
!>    @param[in,out] ProzA_lock locks to mark alreads computed blocks
      SUBROUTINE omp_lu_solve2(n_i,n_j,n_k,lxyz_block,bloc,xloc,lma, &
          lmb,lmc,lud,lme,lmf,lmg,ud_block,bound_block,xyz_block, &
          proza_lock)
        use arrays
        use mod_OMP_TOOLS
        use mod_blocking_size
        IMPLICIT NONE
        INCLUDE 'OMP_TOOLS.inc'
!     N : length of all vectors r,z
        INTEGER n_i, n_j, n_k
        INTEGER k_e, j_e, i_e
        INTEGER ii, jj, kk, ll
        INTEGER xi, yi, zi
        INTEGER xip1, yip1, zip1, xim1, yim1, zim1
!      use mod_blocking_size
        INTEGER lxyz_block, xyz_block(3,lxyz_block)
        DOUBLE PRECISION lme(block_i*block_j*block_k,lxyz_block)
        DOUBLE PRECISION lmf(block_i*block_j*block_k,lxyz_block)
        DOUBLE PRECISION lmg(block_i*block_j*block_k,lxyz_block)
        DOUBLE PRECISION lmb(block_i*block_j*block_k,lxyz_block)
        DOUBLE PRECISION lmc(block_i*block_j*block_k,lxyz_block)
        DOUBLE PRECISION lud(block_i*block_j*block_k,lxyz_block)
        DOUBLE PRECISION lma(block_i*block_j*block_k,lxyz_block)
        DOUBLE PRECISION bloc(block_i*block_j*block_k,lxyz_block)
        DOUBLE PRECISION xloc(block_i*block_j*block_k,lxyz_block)
!     [UD] surface on position-1
        DOUBLE PRECISION ud_block(block_i*block_j+block_i*block_k+ &
          block_j*block_k,lxyz_block)
!     global boundary buffer [x]
        DOUBLE PRECISION bound_block(block_i*block_j+block_i*block_k+ &
          block_j*block_k,bdim_i,bdim_j,bdim_k)

        INTEGER proza_lock(bdim_i*bdim_j*bdim_k)
        INTEGER par_name, isurf, ijsurf, ijksurf
        EXTERNAL par_name


!     surface offset: 1=I-dim, isurf=J-dim, ijsurf=K-dim
        isurf = block_j*block_k + 1
        ijsurf = isurf + block_i*block_k
        ijksurf = block_j*block_k + block_i*block_k + block_i*block_j
!     -----------------------------------------------------------------

!     over all local blocks
        DO ii = 1, lxyz_block
!        'jj' is the linear block index corresponding [xi,yi,zi]
100       DO jj = 1, lxyz_block
            kk = 0
            xi = xyz_block(1,jj)
            yi = xyz_block(2,jj)
            zi = xyz_block(3,jj)
!           [i,j,k]-index
            ll = xi + bdim_i*(yi-1) + bdim_i*bdim_j*(zi-1)
!           search for locked blocks (needs to be computed)
            IF (par_name(proza_lock(ll))==0) THEN
!              search for unlocked blocks (already computed)
!              test [i-1,j,k]
              IF (xi>1) THEN
!                 [i-1,j,k]-index
                xip1 = ll - 1
                IF (par_name(proza_lock(xip1))==1) kk = kk + 1
              ELSE
                kk = kk + 1
              END IF
!              test [i,j-1,k]
              IF (yi>1) THEN
!                 [i,j-1,k]-index
                xip1 = ll - bdim_i
                IF (par_name(proza_lock(xip1))==1) kk = kk + 1
              ELSE
                kk = kk + 1
              END IF
!              test [i,j,k-1]
              IF (zi>1) THEN
!                 [i,j,k-1]-index
                xip1 = ll - bdim_i*bdim_j
                IF (par_name(proza_lock(xip1))==1) kk = kk + 1
              ELSE
                kk = kk + 1
              END IF
!              can the block be computed ?
              IF (kk==3) THEN
!                 now: [xi,yi,zi] is the current computable block 'll', private 'jj'
                GO TO 101
              END IF
            END IF
          END DO
!        polling ... needs to found a free block
          GO TO 100
!        found a computable block
101       CONTINUE

!        compute begin & end
          i_e = min(block_i,n_i-((xi-1)*block_i))
          j_e = min(block_j,n_j-((yi-1)*block_j))
          k_e = min(block_k,n_k-((zi-1)*block_k))
!        block index dim+1
          xip1 = min(xi+1,bdim_i)
          yip1 = min(yi+1,bdim_j)
          zip1 = min(zi+1,bdim_k)

!        solve [L]x[t]=[b]
          CALL flu_left2(block_i,block_j,block_k,i_e,j_e,k_e, &
            bloc(1,jj),xloc(1,jj),lma(1,jj),lmb(1,jj),lmc(1,jj), &
            lud(1,jj),bound_block(1,xi,yi,zi), &
            bound_block(isurf,xi,yi,zi),bound_block(ijsurf,xi,yi,zi), &
            bound_block(1,xip1,yi,zi),bound_block(isurf,xi,yip1,zi), &
            bound_block(ijsurf,xi,yi,zip1),ud_block(1,jj), &
            ud_block(isurf,jj),ud_block(ijsurf,jj),xip1-xi-1, &
            yip1-yi-1,zip1-zi-1)

          CALL par_disab(proza_lock(ll))

          IF ((xi==bdim_i) .AND. (yi==bdim_j) .AND. (zi==bdim_k)) THEN
            CALL par_reset(proza_lock)
            CALL set_dval(ijksurf,0.D0,bound_block(1,xi,yi,zi))
          END IF

        END DO

!     -----------------------------------------------------------------

!     waiting for finishing L-step
!$OMP barrier

!     -----------------------------------------------------------------

!     over all local blocks
        DO ii = 1, lxyz_block
!        'jj' is the linear block index corresponding [xi,yi,zi]
200       DO jj = lxyz_block, 1, -1
            kk = 0
            xi = xyz_block(1,jj)
            yi = xyz_block(2,jj)
            zi = xyz_block(3,jj)
!           [i,j,k]-index
            ll = xi + bdim_i*(yi-1) + bdim_i*bdim_j*(zi-1)
!           search for locked blocks (needs to be computed)
            IF (par_name(proza_lock(ll))==0) THEN
!              search for unlocked blocks (already computed)
!              test [i+1,j,k]
              IF (xi<bdim_i) THEN
!                 [i+1,j,k]-index
                xip1 = ll + 1
                IF (par_name(proza_lock(xip1))==1) kk = kk + 1
              ELSE
                kk = kk + 1
              END IF
!              test [i,j+1,k]
              IF (yi<bdim_j) THEN
!                 [i,j+1,k]-index
                xip1 = ll + bdim_i
                IF (par_name(proza_lock(xip1))==1) kk = kk + 1
              ELSE
                kk = kk + 1
              END IF
!              test [i,j,k+1]
              IF (zi<bdim_k) THEN
!                 [i,j,k+1]-index
                xip1 = ll + bdim_i*bdim_j
                IF (par_name(proza_lock(xip1))==1) kk = kk + 1
              ELSE
                kk = kk + 1
              END IF
!              can the block be computed ?
              IF (kk==3) THEN
!                 now: [xi,yi,zi] is the current computable block 'll', private 'jj'
                GO TO 201
              END IF
            END IF
          END DO
!        polling ... needs to found a free block
          GO TO 200
!        found a computable block
201       CONTINUE

!        compute begin & end
          i_e = min(block_i,n_i-((xi-1)*block_i))
          j_e = min(block_j,n_j-((yi-1)*block_j))
          k_e = min(block_k,n_k-((zi-1)*block_k))
!        block index dim-1
          xim1 = max(xi-1,1)
          yim1 = max(yi-1,1)
          zim1 = max(zi-1,1)

!        solve [U]x[t]=[b]
          CALL flu_right2(block_i,block_j,block_k,i_e,j_e,k_e, &
            xloc(1,jj),lud(1,jj),lme(1,jj),lmf(1,jj),lmg(1,jj), &
            bound_block(1,xim1,yi,zi),bound_block(isurf,xi,yim1,zi), &
            bound_block(ijsurf,xi,yi,zim1),bound_block(1,xi,yi,zi), &
            bound_block(isurf,xi,yi,zi),bound_block(ijsurf,xi,yi,zi), &
            xim1-xi+1,yim1-yi+1,zim1-zi+1)

          CALL par_disab(proza_lock(ll))

          IF ((xi==1) .AND. (yi==1) .AND. (zi==1)) THEN
            CALL par_reset(proza_lock)
            CALL set_dval(ijksurf,0.D0,bound_block(1,xi,yi,zi))
          END IF

        END DO

!     -----------------------------------------------------------------

        RETURN
      END

!>    @brief compute [L]-part by solving [L] x [t] = [b]
!>    @param[in] bloc vector [b]
!>    @param[in,out] xloc result vector [t]
!>    @param[in] UA thread local blocks of the 1. diagonal of [L]
!>    @param[in] UB thread local blocks of the 2. diagonal of [L]
!>    @param[in] UC thread local blocks of the 3. diagonal of [L]
!>    @param[in] UD thread local blocks of the helper diagonal [UD] of [L]
!>    @param[in] ldI leading dimensions in I direction
!>    @param[in] ldJ leading dimensions in J direction
!>    @param[in] ldK leading dimensions in K direction
!>    @param[in] i_e block size in I direction
!>    @param[in] j_e block size in J direction
!>    @param[in] k_e block size in K direction
!>    @param[in] xloc_im1 global boundary buffer from the -1 neighbour in I direction
!>    @param[in] xloc_jm1 global boundary buffer from the -1 neighbour in J direction
!>    @param[in] xloc_km1 global boundary buffer from the -1 neighbour in K direction
!>    @param[out] xloc_ip1 global boundary buffer to the +1 neighbour in I direction
!>    @param[out] xloc_jp1 global boundary buffer to the +1 neighbour in J direction
!>    @param[out] xloc_kp1 global boundary buffer to the +1 neighbour in K direction
!>    @param[in] UD_im1 local UD-buffer from the -1 neighbour in I direction
!>    @param[in] UD_jm1 local UD-buffer from the -1 neighbour in J direction
!>    @param[in] UD_km1 local UD-buffer from the -1 neighbour in K direction
!>    @param[in] ci 0: copy boundary values in [xloc_ip1], <>0: do nothing
!>    @param[in] cj 0: copy boundary values in [xloc_jp1], <>0: do nothing
!>    @param[in] ck 0: copy boundary values in [xloc_kp1], <>0: do nothing
      SUBROUTINE flu_left2(ldi,ldj,ldk,i_e,j_e,k_e,bloc,xloc,ua,ub,uc, &
          ud,xloc_im1,xloc_jm1,xloc_km1,xloc_ip1,xloc_jp1,xloc_kp1, &
          ud_im1,ud_jm1,ud_km1,ci,cj,ck)
        use mod_OMP_TOOLS
        use mod_blocking_size
        IMPLICIT NONE
        INCLUDE 'OMP_TOOLS.inc'

!     N : length of all vectors r,z
        INTEGER i_e, j_e, k_e, i, j, k, ldi, ldj, ldk, ci, cj, ck
!     vectors [r], [z]
        DOUBLE PRECISION bloc(ldi,ldj,ldk), xloc(ldi,ldj,ldk)
        DOUBLE PRECISION ua(ldi,ldj,ldk), ub(ldi,ldj,ldk)
        DOUBLE PRECISION uc(ldi,ldj,ldk), ud(ldi,ldj,ldk)

!      use mod_blocking_size
!     global boundary buffer
        DOUBLE PRECISION xloc_im1(block_j,block_k)
        DOUBLE PRECISION xloc_jm1(block_i,block_k)
        DOUBLE PRECISION xloc_km1(block_i,block_j)
        DOUBLE PRECISION xloc_ip1(block_j,block_k)
        DOUBLE PRECISION xloc_jp1(block_i,block_k)
        DOUBLE PRECISION xloc_kp1(block_i,block_j)
!     local UD-buffer (-1)
        DOUBLE PRECISION ud_im1(block_j,block_k)
        DOUBLE PRECISION ud_jm1(block_i,block_k)
        DOUBLE PRECISION ud_km1(block_i,block_j)

!$OMP flush(xloc_im1, xloc_jm1, xloc_km1)

!     solve [L]x[t]=[b]
!     [0D]
        xloc(1,1,1) = bloc(1,1,1) - uc(1,1,1)*ud_im1(1,1)*xloc_im1(1,1 &
          ) - ub(1,1,1)*ud_jm1(1,1)*xloc_jm1(1,1) - &
          ua(1,1,1)*ud_km1(1,1)*xloc_km1(1,1)


!     [1D-I]
        DO i = 2, i_e
          xloc(i,1,1) = bloc(i,1,1) - uc(i,1,1)*ud(i-1,1,1)*xloc(i-1,1 &
            ,1) - ub(i,1,1)*ud_jm1(i,1)*xloc_jm1(i,1) - &
            ua(i,1,1)*ud_km1(i,1)*xloc_km1(i,1)
        END DO

!     [1D-J]
        DO j = 2, j_e
          xloc(1,j,1) = bloc(1,j,1) - uc(1,j,1)*ud_im1(j,1)*xloc_im1(j &
            ,1) - ub(1,j,1)*ud(1,j-1,1)*xloc(1,j-1,1) - &
            ua(1,j,1)*ud_km1(1,j)*xloc_km1(1,j)
        END DO

!     [1D-K]
        DO k = 2, k_e
          xloc(1,1,k) = bloc(1,1,k) - uc(1,1,k)*ud_im1(1,k)*xloc_im1(1 &
            ,k) - ub(1,1,k)*ud_jm1(1,k)*xloc_jm1(1,k) - &
            ua(1,1,k)*ud(1,1,k-1)*xloc(1,1,k-1)
        END DO


!     [2D-I]
        DO k = 2, k_e
          DO j = 2, j_e
            xloc(1,j,k) = bloc(1,j,k) - uc(1,j,k)*ud_im1(j,k)*xloc_im1 &
              (j,k) - ub(1,j,k)*ud(1,j-1,k)*xloc(1,j-1,k) - &
              ua(1,j,k)*ud(1,j,k-1)*xloc(1,j,k-1)
          END DO
        END DO

!     [2D-J]
        DO k = 2, k_e
          DO i = 2, i_e
            xloc(i,1,k) = bloc(i,1,k) - uc(i,1,k)*ud(i-1,1,k)*xloc(i-1 &
              ,1,k) - ub(i,1,k)*ud_jm1(i,k)*xloc_jm1(i,k) - &
              ua(i,1,k)*ud(i,1,k-1)*xloc(i,1,k-1)
          END DO
        END DO

!     [2D-K]
        DO j = 2, j_e
          DO i = 2, i_e
            xloc(i,j,1) = bloc(i,j,1) - uc(i,j,1)*ud(i-1,j,1)*xloc(i-1 &
              ,j,1) - ub(i,j,1)*ud(i,j-1,1)*xloc(i,j-1,1) - &
              ua(i,j,1)*ud_km1(i,j)*xloc_km1(i,j)
          END DO
        END DO


!     [3D]
        DO k = 2, k_e
          DO j = 2, j_e
            DO i = 2, i_e
              xloc(i,j,k) = bloc(i,j,k) - uc(i,j,k)*ud(i-1,j,k)*xloc(i &
                -1,j,k) - ub(i,j,k)*ud(i,j-1,k)*xloc(i,j-1,k) - &
                ua(i,j,k)*ud(i,j,k-1)*xloc(i,j,k-1)
            END DO
          END DO
        END DO


!     boundary exchange
!     [2D-I]
        IF (ci==0) THEN
          DO k = 1, k_e
            DO j = 1, j_e
              xloc_ip1(j,k) = xloc(i_e,j,k)
            END DO
          END DO
        END IF

!     [2D-J]
        IF (cj==0) THEN
          DO k = 1, k_e
            DO i = 1, i_e
              xloc_jp1(i,k) = xloc(i,j_e,k)
            END DO
          END DO
        END IF

!     [2D-K]
        IF (ck==0) THEN
          DO j = 1, j_e
            DO i = 1, i_e
              xloc_kp1(i,j) = xloc(i,j,k_e)
            END DO
          END DO
        END IF

!$OMP flush(xloc_ip1, xloc_jp1, xloc_kp1)

        RETURN
      END

!>    @brief compute [U]-part by solving [U] x [x] = [t]
!>    @param[in] ldI leading dimensions in I direction
!>    @param[in] ldJ leading dimensions in J direction
!>    @param[in] ldK leading dimensions in K direction
!>    @param[in] i_e block size in I direction
!>    @param[in] j_e block size in J direction
!>    @param[in] k_e block size in K direction
!>    @param[in,out] xloc result vector [x], [t] as input
!>    @param[in] UD thread local blocks of the helper diagonal [UD] of [U]
!>    @param[in] UE thread local blocks of the 2. diagonal of [U]
!>    @param[in] UF thread local blocks of the 3. diagonal of [U]
!>    @param[in] UG thread local blocks of the 4. diagonal of [U]
!>    @param[out] xloc_im1 global boundary buffer to the -1 neighbour in I direction
!>    @param[out] xloc_jm1 global boundary buffer to the -1 neighbour in J direction
!>    @param[out] xloc_km1 global boundary buffer to the -1 neighbour in K direction
!>    @param[in] xloc_ip1 global boundary buffer from the +1 neighbour in I direction
!>    @param[in] xloc_jp1 global boundary buffer from the +1 neighbour in J direction
!>    @param[in] xloc_kp1 global boundary buffer from the +1 neighbour in K direction
!>    @param[in] ci 0: copy boundary values in [xloc_im1], <>0: do nothing
!>    @param[in] cj 0: copy boundary values in [xloc_jm1], <>0: do nothing
!>    @param[in] ck 0: copy boundary values in [xloc_km1], <>0: do nothing
      SUBROUTINE flu_right2(ldi,ldj,ldk,i_e,j_e,k_e,xloc,ud,ue,uf,ug, &
          xloc_im1,xloc_jm1,xloc_km1,xloc_ip1,xloc_jp1,xloc_kp1,ci,cj, &
          ck)
        use mod_OMP_TOOLS
        use mod_blocking_size
        IMPLICIT NONE
        INCLUDE 'OMP_TOOLS.inc'

!     N : length of all vectors r,z
        INTEGER k_e, j_e, i_e, i, j, k, ldi, ldj, ldk, ci, cj, ck

!     vectors [z]
        DOUBLE PRECISION xloc(ldi,ldj,ldk)
        DOUBLE PRECISION ud(ldi,ldj,ldk), ue(ldi,ldj,ldk)
        DOUBLE PRECISION uf(ldi,ldj,ldk), ug(ldi,ldj,ldk)

!      use mod_blocking_size
!     global boundary buffer
        DOUBLE PRECISION xloc_im1(block_j,block_k)
        DOUBLE PRECISION xloc_jm1(block_i,block_k)
        DOUBLE PRECISION xloc_km1(block_i,block_j)
        DOUBLE PRECISION xloc_ip1(block_j,block_k)
        DOUBLE PRECISION xloc_jp1(block_i,block_k)
        DOUBLE PRECISION xloc_kp1(block_i,block_j)

!$OMP flush(xloc_ip1, xloc_jp1, xloc_kp1)

!     solve [U]x[x]=[t]
!     [0D]
        xloc(i_e,j_e,k_e) = (xloc(i_e,j_e,k_e)-ue(i_e,j_e,k_e)* &
          xloc_ip1(j_e,k_e)-uf(i_e,j_e,k_e)*xloc_jp1(i_e,k_e)- &
          ug(i_e,j_e,k_e)*xloc_kp1(i_e,j_e))*ud(i_e,j_e,k_e)


!     [1D-I]
        DO i = i_e - 1, 1, -1
          xloc(i,j_e,k_e) = (xloc(i,j_e,k_e)-ue(i,j_e,k_e)*xloc(i+1, &
            j_e,k_e)-uf(i,j_e,k_e)*xloc_jp1(i,k_e)- &
            ug(i,j_e,k_e)*xloc_kp1(i,j_e))*ud(i,j_e,k_e)
        END DO

!     [1D-J]
        DO j = j_e - 1, 1, -1
          xloc(i_e,j,k_e) = (xloc(i_e,j,k_e)-ue(i_e,j,k_e)*xloc_ip1(j, &
            k_e)-uf(i_e,j,k_e)*xloc(i_e,j+1,k_e)- &
            ug(i_e,j,k_e)*xloc_kp1(i_e,j))*ud(i_e,j,k_e)
        END DO

!     [1D-K]
        DO k = k_e - 1, 1, -1
          xloc(i_e,j_e,k) = (xloc(i_e,j_e,k)-ue(i_e,j_e,k)*xloc_ip1( &
            j_e,k)-uf(i_e,j_e,k)*xloc_jp1(i_e,k)- &
            ug(i_e,j_e,k)*xloc(i_e,j_e,k+1))*ud(i_e,j_e,k)
        END DO


!     [2D-I]
        DO k = k_e - 1, 1, -1
          DO j = j_e - 1, 1, -1
            xloc(i_e,j,k) = (xloc(i_e,j,k)-ue(i_e,j,k)*xloc_ip1(j,k)- &
              uf(i_e,j,k)*xloc(i_e,j+1,k)-ug(i_e,j,k)*xloc(i_e,j,k+1)) &
              *ud(i_e,j,k)
          END DO
        END DO

!     [2D-J]
        DO k = k_e - 1, 1, -1
          DO i = i_e - 1, 1, -1
            xloc(i,j_e,k) = (xloc(i,j_e,k)-ue(i,j_e,k)*xloc(i+1,j_e,k) &
              -uf(i,j_e,k)*xloc_jp1(i,k)-ug(i,j_e,k)*xloc(i,j_e,k+1))* &
              ud(i,j_e,k)
          END DO
        END DO

!     [2D-K]
        DO j = j_e - 1, 1, -1
          DO i = i_e - 1, 1, -1
            xloc(i,j,k_e) = (xloc(i,j,k_e)-ue(i,j,k_e)*xloc(i+1,j,k_e) &
              -uf(i,j,k_e)*xloc(i,j+1,k_e)-ug(i,j,k_e)*xloc_kp1(i,j))* &
              ud(i,j,k_e)
          END DO
        END DO


!     [3D]
        DO k = k_e - 1, 1, -1
          DO j = j_e - 1, 1, -1
            DO i = i_e - 1, 1, -1
              xloc(i,j,k) = (xloc(i,j,k)-ue(i,j,k)*xloc(i+1,j,k)-uf(i, &
                j,k)*xloc(i,j+1,k)-ug(i,j,k)*xloc(i,j,k+1))*ud(i,j,k)
            END DO
          END DO
        END DO


!     boundary exchange
!     [2D-I]
        IF (ci==0) THEN
          DO k = 1, k_e
            DO j = 1, j_e
              xloc_im1(j,k) = xloc(1,j,k)
            END DO
          END DO
        END IF

!     [2D-J]
        IF (cj==0) THEN
          DO k = 1, k_e
            DO i = 1, i_e
              xloc_jm1(i,k) = xloc(i,1,k)
            END DO
          END DO
        END IF

!     [2D-K]
        IF (ck==0) THEN
          DO j = 1, j_e
            DO i = 1, i_e
              xloc_km1(i,j) = xloc(i,j,1)
            END DO
          END DO
        END IF

!$OMP flush(xloc_im1, xloc_jm1, xloc_km1)

        RETURN
      END

!>    @brief compute a processor grid distribution
!>    @param[in] lbl_size cache (block) size
!>    @param[in] P number of threads
!>    @param[out] bi block size in each dimension
!>    @param[out] bj block size in each dimension
!>    @param[out] bk block size in each dimension
!>    @param[in] ni grid dimension in I0 direction
!>    @param[in] nj grid dimension in I0 direction
!>    @param[in] nk grid dimension in I0 direction
      SUBROUTINE proz_grid(lbl_size,p,bi,bj,bk,ni,nj,nk)
        use mod_genrl
        use mod_linfos
        use mod_blocking_size
        IMPLICIT NONE
        integer :: i, j, k, l, m
        INTEGER bi, bj, bk, p, ii, jj, kk,ni,nj,nk
        INTEGER qwurz, lbl_size
        INTRINSIC sqrt, max, dble, int
        INTEGER vec_anz, ct_cm
!     vec_anz : [A B C + D + W + X] =#6
        PARAMETER (vec_anz=6)
        PARAMETER (ct_cm=55)
       INTEGER cache2d_size                   
!       8192 Double-Precisions (64KByte)       
        PARAMETER (cache2d_size = 8192)        
        DOUBLE PRECISION dii, djj, dkk, dti, dtj, dtk
!       ! manuell block size for benchmarking !
        IF (linfos(4)==-1000) THEN
          bi = block_i
          bj = block_j
          bk = block_k
          RETURN
        END IF
! --------- NEW METHOD
!       compute block size for big models
        IF (ni*nj*nk>=4*cache2d_size) THEN
!         - preset: I0 dimension
          dii = -1.0d0
          bi = -1
          DO ii = min(ni,20), min(ni,70)
            i = (ni-1) /ii
            i = ni -i*ii
            dti = dble(i)/dble(ii)
            IF (dti>=dii) THEN
              dii = dti -1.0D-10
              bi = ii
            END IF
          END DO
!         - preset: J0 dimension
          djj = -1.0d0
          bj = -1
          DO jj = min(nj,20), min(nj,40)
            j = (nj-1) /jj
            j = nj -j*jj
            dtj = dble(j)/dble(jj)
            IF (dtj>=djj) THEN
              djj = dtj -1.0D-10
              bj = jj
            END IF
          END DO
!         - preset: K0 dimension
          dkk = -1.0d0
          bk = -1
          DO kk = nk, min(nk,4), -1
            k = (nk-1) /kk
            k = nk -k*kk
            dtk = dble(k)/dble(kk)
            IF (dtk>=dkk) THEN
              dkk = dtk -1.0D-10
              bk = kk
            END IF
          END DO
!         search in the I0,J0,K0 combination for a 64K-cache-fit
          l = ni*nj*nk
          DO kk = nk, min(nk,4), -1
            k = (nk-1) /kk
            k = nk -k*kk
            dtk = dble(k)/dble(kk)
            IF (dtk>=dkk) THEN
              DO jj = min(nj,20), min(nj,40)
                j = (nj-1) /jj
                j = nj -j*jj
                dtj = dble(j)/dble(jj)
                IF (dtj>=djj) THEN
                  DO ii = min(ni,20), min(ni,70)
                    i = (ni-1) /ii
                    i = ni -i*ii
                    dti = dble(i)/dble(ii)
                    IF (dti>=dii) THEN
!                     proof block size [m] close to cache block size [cache2d_size]
                      m = abs(ii*jj*kk -cache2d_size)
                      IF (m<=l) THEN
                        dii = dti -1.0D-10
                        bi = ii
                        djj = dtj -1.0D-10
                        bj = jj
                        dkk = dtk -1.0D-10
                        bk = kk
                        l = m
                      END IF
                    END IF
                  END DO
                END IF
              END DO
            END IF
          END DO
        ELSE
! --------- OLD METHOD
!         ! automatic block size computation !
!         minimal block length in I0 direction
          IF (ni>ct_cm) THEN
!          2D or 3D partitioning
!          a little bit more is better
!          0.0-1.99: 1; 2.0-2.99: 2; ...
            bi = max(ni/ct_cm,1)
            bi = (ni-1)/bi + 1
!          number of blocks in J0*K0
!            ik_bl = bi*J0*K0*veC_anz/lbl_size
!          size of J0*K0 blocks matching the cache size
!            J0*K0/ik_bl
!          ideal size for J0 and K0
!          qwurz = int(sqrt(dble(J0*K0/ik_bl)))
            qwurz = max(int(sqrt(dble(lbl_size/(bi*vec_anz)))),1)
!          J0 dimension
!          a little bit more is better
!          0.0-1.99: 1; 2.0-2.99: 2; ...
            bj = (nj-1)/qwurz + 1
            bj = (nj-1)/bj + 1
!          rest for K0 dimension
!          a little bit more is better
!          0.0-1.99: 1; 2.0-2.99: 2; ...
            bk = vec_anz*(nk-1)*bj*bi/lbl_size + 1
            bk = (nk-1)/bk + 1
          ELSE IF (ni*nj>ct_cm) THEN
!          simple 2D partitioning, when enough elements
!          I0 dimension
            bi = ni
!          number of cache blocks in J0*K0
!            (J0 *K0 *vec_anz *bi) /lbl_size
!          ideal size for J0 and K0
            qwurz = max(int(sqrt(dble((nj*nk*vec_anz*bi)/lbl_size))),1)
!          J0 dimension
!          a little bit more is better
!          0.0-1.99: 1; 2.0-2.99: 2; ...
            bj = max(ni*nj/ct_cm,1)
            bj = (nj-1)/bj + 1
            bj = max(bj,(nj-1)/qwurz+1)
!          rest for K0 dimension
!          a little bit more is better
!          0.0-1.99: 1; 2.0-2.99: 2; ...
            bk = vec_anz*(nk-1)*bj*bi/lbl_size + 1
            bk = (nk-1)/bk + 1
          ELSE
!         !  in case of a 2D model and a too small I0, no parallelization !
            bi = ni
            bj = nj
            bk = nk
          END IF
        END IF
! ---------
        RETURN
      END

!>    @param[in] ni grid dimension in I0 direction
!>    @param[in] nj grid dimension in J0 direction
!>    @param[in] nk grid dimension in K0 direction
!>    @brief initialise [proza]-thread association for parallel computation
      SUBROUTINE par_init2(ni,nj,nk)
        use arrays
        use mod_genrl
        use mod_linfos
        use mod_OMP_TOOLS
        use mod_blocking_size
        IMPLICIT NONE
        integer :: i, j, k
        INCLUDE 'OMP_TOOLS.inc'
        INTEGER bi, bj, bk, p, pa, ni,nj,nk
        INTEGER li, lj, lk, maxlevel
        INTEGER, ALLOCATABLE :: speed(:)
        INTEGER, ALLOCATABLE :: nbl_proz(:), sbl_proz(:), lev_idx(:,:)
        INTEGER ilev_b, ilev_e, pro_b, pro_e, t_max
        INTEGER subbl, parmax, parlen
        PARAMETER (subbl=16)
        DOUBLE PRECISION speedup, speedup_vec
!      use mod_blocking_size
        INTRINSIC max, min, log10, dble, int
#ifdef BENCH        
        character (len=40) :: bench_format
#endif

!     number of processors
        p = tlevel_1

!     minimal block length for each dimension
        i = subbl + 1
123     i = i - 1
        CALL proz_grid(i*bl_size/subbl,p,bi,bj,bk,ni,nj,nk)
        li = (ni-1)/bi + 1
        lj = (nj-1)/bj + 1
        lk = (nk-1)/bk + 1
!       max. useable threads
        parmax = max(max(min(li,lj),min(li,lk)),min(lj,lk))
!       # of max. thread runs
        parlen = max(max(li,lj),lk) - min(parmax,p)
!debug        write(*,*) '[parmax,parlen,i]:',parmax,parlen,i
        IF ((parmax<p .OR. parlen<p) .AND. i>1) GO TO 123

        block_i = bi
        block_j = bj
        block_k = bk
        bdim_i = li
        bdim_j = lj
        bdim_k = lk

        ALLOCATE(proza(li,lj,lk))

!     compute levels, all blocks for the same level are independent
        maxlevel = 0
        DO k = 1, lk
          DO j = 1, lj
            DO i = 1, li
              proza(i,j,k) = 0
!           set to maximum of all previous neighboars
              IF (i>1) proza(i,j,k) = max(proza(i,j,k),proza(i-1,j,k))
              IF (j>1) proza(i,j,k) = max(proza(i,j,k),proza(i,j-1,k))
              IF (k>1) proza(i,j,k) = max(proza(i,j,k),proza(i,j,k-1))
!           increase level
              proza(i,j,k) = proza(i,j,k) + 1
              maxlevel = max(maxlevel,proza(i,j,k))
            END DO
          END DO
        END DO

        ALLOCATE(speed(maxlevel))
        DO i = 1, maxlevel
          speed(i) = 0
        END DO
!     estimate max number of threads per level
        DO k = 1, lk
          DO j = 1, lj
            DO i = 1, li
!           ProzA from 1..(P)
              speed(proza(i,j,k)) = speed(proza(i,j,k)) + 1
            END DO
          END DO
        END DO
!     PA: maximum number of independent threads over all levels
!     speedup : approx. speedup
        pa = 0
        speedup = 0.0D0
        DO i = 1, maxlevel
          pa = max(speed(i),pa)
          speedup = speedup + dble(speed(i))/dble(min(speed(i),p))
        END DO
        speedup = dble(li*lj*lk)/speedup


!     max number of blocks for each thread
        ALLOCATE(nbl_proz(p))
!     used (setup) number of blocks for each thread during each level
        ALLOCATE(sbl_proz(p))
!     index of blocks for this level
        ALLOCATE(lev_idx(3,pa))

!     number of blocks equal for each thread
        i = (li*lj*lk)/p
        DO j = 1, p
          nbl_proz(j) = i
        END DO
!     overhead of blocks distributed over all threads
        i = li*lj*lk - i*p
        DO j = 1, i
          nbl_proz(j) = nbl_proz(j) + 1
        END DO
        speedup_vec = dble(li*lj*lk)/dble(nbl_proz(1))
!     max number of blocks (the first thread has the highest number)
        max_blocks = nbl_proz(1)

!     for each level distribute the blocks of each thread (bi-directional)
        ilev_b = 1
        ilev_e = maxlevel
1000    CONTINUE
!     first: forward direction
!       init
        pro_b = 0
        t_max = 0
        DO i = 1, p
          sbl_proz(i) = 0
        END DO
!       search for the 'ilev_b' level
        DO k = 1, lk
          DO j = 1, lj
            DO i = 1, li
!              right [i,j,k]-index for this level
              IF (proza(i,j,k)==ilev_b) THEN
                t_max = t_max + 1
!                 next thread with free blocks
1001            pro_b = mod(pro_b,p) + 1
!                 skip this thread:
!                   1. no blocks
!                   2. already one block
!                      and not enough for later levels
                IF (nbl_proz(pro_b)<=0 .OR. (sbl_proz(pro_b)>0 .AND. ( &
                  maxlevel-2*ilev_b+1)>=nbl_proz(pro_b))) GO TO 1001
!                 mark the block for this level
                nbl_proz(pro_b) = nbl_proz(pro_b) - 1
                sbl_proz(pro_b) = sbl_proz(pro_b) + 1
                lev_idx(1,t_max) = i
                lev_idx(2,t_max) = j
                lev_idx(3,t_max) = k
              END IF
            END DO
          END DO
        END DO
!       distribute the blocks
        pro_b = 1
        CALL distribute_bl(t_max,p,pro_b,sbl_proz,lev_idx)
!     second: backward direction
        IF (ilev_b==ilev_e) GO TO 1020
!       init
        pro_e = 1
        t_max = 0
        DO i = 1, p
          sbl_proz(i) = 0
        END DO
!       search for the 'ilev_e' level
        DO k = lk, 1, -1
          DO j = lj, 1, -1
            DO i = li, 1, -1
!              right [i,j,k]-index for this level
              IF (proza(i,j,k)==ilev_e) THEN
                t_max = t_max + 1
!                 next thread with free blocks
1011            pro_e = mod(pro_e+p-2,p) + 1
!                 skip this thread:
!                   1. no blocks
!                   2. already one block
!                      and not enough for later levels
                IF (nbl_proz(pro_e)<=0 .OR. (sbl_proz(pro_e)>0 .AND. ( &
                  maxlevel-2*ilev_b)>=nbl_proz(pro_e))) GO TO 1011
!                 mark the block for this level
                nbl_proz(pro_e) = nbl_proz(pro_e) - 1
                sbl_proz(pro_e) = sbl_proz(pro_e) + 1
                lev_idx(1,t_max) = i
                lev_idx(2,t_max) = j
                lev_idx(3,t_max) = k
              END IF
            END DO
          END DO
        END DO
!       distribute the blocks
        pro_e = p
        CALL distribute_bl(t_max,p,pro_e,sbl_proz,lev_idx)

!     next levels
1020    CONTINUE
        ilev_b = ilev_b + 1
        ilev_e = ilev_e - 1
        IF (ilev_b<=ilev_e) GO TO 1000


        DEALLOCATE(lev_idx)
        DEALLOCATE(sbl_proz)
        DEALLOCATE(nbl_proz)
        DEALLOCATE(speed)

!     setup return values
        DO k = 1, lk
          DO j = 1, lj
            DO i = 1, li
!           ProzA : [0..(PA-1)]
              proza(i,j,k) = -proza(i,j,k) - 1
            END DO
          END DO
        END DO

#ifdef BENCH
        IF (p>1) THEN
          WRITE(*,'(A)') '  [I] : ILU block - thread assignment:'
          WRITE(bench_format,'(A,I2,A,I2,A,I2,A)') '(7X,', lk, '(', &
            li, 'I', int(log10(dble(pa))) + 2, ',2X))'
          DO j = 1, lj
            WRITE(*,bench_format) ((proza(i,j,k),i=1,li),k=1,lk)
          END DO
        END IF
!     'speed' should now be the number of blocks for each thread
!aw      write(*,'(A)') '  [I] : Number of bloCks for eaCh thread:'
!aw      write (benCh_format,'(A,I2,A,I2,A)') '(A,',P,'I',
!aw     &  int(log10(dble(speed(1))))+2,',A)'
!aw      write(*,benCh_format) '        [',(i-1,i=1,P),'] thread'
!aw      write(*,benCh_format) '        [',(speed(i),i=1,P),'] bloCks'
        WRITE(*,'(A,3(I5,A),3(I3,A))') '  [I] : ILU: block size =[', &
          bi, ',', bj, ',', bk, '], dim=[', li, ',', lj, ',', lk, ']'
        WRITE(*,'(A,F7.2)') &
          '        theoretical speedup for solver only :', speedup_vec
        WRITE(*,'(A,F7.2)') &
          '        theoretical speedup for ILU only    :', speedup
        CALL write_grid()
#endif
        IF (linfos(4)>=1) THEN
          WRITE(*,'(A,F7.2)') &
            '  [I] : approx. speedup for the full solver :', &
            (2.0D0*speedup+1.0D0*speedup_vec)/3.0D0
        END IF

        RETURN
      END

!>    @brief search for nearest neighbours and distributes the blocks, writes output in [proza]
!>    @param[in] t_max number of blocks
!>    @param[in] P number of threads (processor cores)
!>    @param[in] pro_be initial thread index, switch with 1: forward counting, [P]: backward counting
!>    @param[in,out] sbl_proz blocks per thread counter
!>    @param[in] lev_idx i,j,k index for each block
      SUBROUTINE distribute_bl(t_max,p,pro_be,sbl_proz,lev_idx)
        use arrays
        use mod_OMP_TOOLS
        use mod_genrl
        use mod_linfos
        use mod_blocking_size
        IMPLICIT NONE
        integer :: i, j, k
        INCLUDE 'OMP_TOOLS.inc'
        INTEGER p, t_max, pro_be
        INTEGER sbl_proz(p), lev_idx(3,t_max)
        INTEGER tidx, ii, jj, kk, ll, mm

        tidx = pro_be
        DO ii = 1, t_max
!         thread for the next block
2001      IF (sbl_proz(tidx)<=0) THEN
            IF (pro_be==p) THEN
!             backward
              tidx = mod(tidx+p-2,p) + 1
            ELSE
!             forward
              tidx = mod(tidx,p) + 1
            END IF
            GO TO 2001
          END IF
          sbl_proz(tidx) = sbl_proz(tidx) - 1
!         search the best fit (kk: most neighboars)
          ll = -1
          mm = 1
          DO jj = 1, t_max
            i = lev_idx(1,jj)
            j = lev_idx(2,jj)
            k = lev_idx(3,jj)
            kk = 0
            IF (i>1) THEN
              IF (proza(i-1,j,k)==-tidx) kk = kk + 1
            END IF
            IF (j>1) THEN
              IF (proza(i,j-1,k)==-tidx) kk = kk + 1
            END IF
            IF (k>1) THEN
              IF (proza(i,j,k-1)==-tidx) kk = kk + 1
            END IF
            IF (i<bdim_i) THEN
              IF (proza(i+1,j,k)==-tidx) kk = kk + 1
            END IF
            IF (j<bdim_j) THEN
              IF (proza(i,j+1,k)==-tidx) kk = kk + 1
            END IF
            IF (k<bdim_k) THEN
              IF (proza(i,j,k+1)==-tidx) kk = kk + 1
            END IF
            IF (kk>ll .AND. proza(i,j,k)>=0) THEN
!              better fit
              ll = kk
              mm = jj
            END IF
          END DO
          proza(lev_idx(1,mm),lev_idx(2,mm),lev_idx(3,mm)) = -tidx
!         old strait-forward distribution
!aw          ProzA(lev_idx(1,i),lev_idx(2,i),lev_idx(3,i)) = -tidx
        END DO
        RETURN
      END

!>    @brief function wrapper for array [proza]
!>    @param[in] i cell index, direction I0
!>    @param[in] j cell index, direction J0
!>    @param[in] k cell index, direction K0
!>    @return value of [proza(i,j,k)]
      INTEGER FUNCTION fkt_proza(i,j,k)
        use arrays
        IMPLICIT NONE
        INTEGER i, j, k

        fkt_proza = proza(i,j,k)

        RETURN
      END

!>    @brief writes ILU thread grid in tecplot-format
!>    @details
!> writes ILU thread grid in tecplot-format\n
!> for each block an eight-node cube is created\n
      SUBROUTINE write_grid()
        use arrays
        use mod_blocking_size
        IMPLICIT NONE
!      use mod_blocking_size
        INTEGER i, j, k


        OPEN(76,file='thread_grid.plt',status='unknown',blank='null')

        WRITE(76,'(1A)') 'variables = "i", "j", "k", "thread"'
        WRITE(76,'(3(1A,I5),1A)') 'zone i=', bdim_i*2, ', j=', &
          bdim_j*2, ', k=', bdim_k*2, ', f=point'

        DO k = 1, bdim_k
          DO j = 1, bdim_j
            DO i = 1, bdim_i
              WRITE(76,'(3F6.2,1I6)') dble(i) - 0.99D0, &
                dble(j) - 0.99D0, dble(k) - 0.99D0, proza(i,j,k)
              WRITE(76,'(3F6.2,1I6)') dble(i), dble(j) - 0.99D0, &
                dble(k) - 0.99D0, proza(i,j,k)
            END DO
            DO i = 1, bdim_i
              WRITE(76,'(3F6.2,1I6)') dble(i) - 0.99D0, dble(j), &
                dble(k) - 0.99D0, proza(i,j,k)
              WRITE(76,'(3F6.2,1I6)') dble(i), dble(j), &
                dble(k) - 0.99D0, proza(i,j,k)
            END DO
          END DO
          DO j = 1, bdim_j
            DO i = 1, bdim_i
              WRITE(76,'(3F6.2,1I6)') dble(i) - 0.99D0, &
                dble(j) - 0.99D0, dble(k), proza(i,j,k)
              WRITE(76,'(3F6.2,1I6)') dble(i), dble(j) - 0.99D0, &
                dble(k), proza(i,j,k)
            END DO
            DO i = 1, bdim_i
              WRITE(76,'(3F6.2,1I6)') dble(i) - 0.99D0, dble(j), &
                dble(k), proza(i,j,k)
              WRITE(76,'(3F6.2,1I6)') dble(i), dble(j), dble(k), &
                proza(i,j,k)
            END DO
          END DO
        END DO

        CLOSE(76)
        RETURN
      END
