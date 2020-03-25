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

!>    @brief solve of : [M] x [x] = [b], BICGSTAB algorithm based with ILU preconditioning
!>    @param[in] N_I lengths of I dimension of local matrix [M]
!>    @param[in] N_J lengths of J dimension of local matrix [M]
!>    @param[in] N_K lengths of K dimension of local matrix [M]
!>    @param[in,out] x solution vector [x], on start = start vector
!>    @param[in] b right side, vector [b]
!>    @param[in] r0_hat random vector [r0_hat] ^= [r0^]
!>    @param[in] depsilon precision criteria to break iterations
!>    @param[in] mbc_mask boundary condition pattern (mask)
!>    @param[in] max_It max iteration number
!>    @param[in] criteria precision criteria mode to break iterations\n
!>    - 0 : relative stopping crit.: ||[res]||       < depsilon*||[res0]||\n
!>    - 1 : absolute stopping crit.: ||[res]||       < depsilon\n
!>    - 2 : maximum  stopping crit.: max(abs([res])) < depsilon\n
!>    first [res]^=[r], later (if precise enough): [res]^=([M]x[x]-[b])
!>    @param[in] MA 1. diagonal of the system matrix [M]
!>    @param[in] MB 2. diagonal of the system matrix [M]
!>    @param[in] MC 3. diagonal of the system matrix [M]
!>    @param[in] MD 4. diagonal of the system matrix [M]
!>    @param[in] ME 5. diagonal of the system matrix [M]
!>    @param[in] MF 6. diagonal of the system matrix [M]
!>    @param[in] MG 7. diagonal of the system matrix [M]
!>    @param[in] UD helper diagonal elements for preconditioning
!>    @param[out] bound_block boundary exchange buffer for each block, between the threads
!>    @param[out] dnrm normalisation vector, temporary use
!>    @param[out] lMA temporary thread local elements of the 1. diagonal of [M]
!>    @param[out] lMB temporary thread local elements of the 2. diagonal of [M]
!>    @param[out] lMC temporary thread local elements of the 3. diagonal of [M]
!>    @param[out] lMD temporary thread local elements of the 4. diagonal of [M]
!>    @param[out] lME temporary thread local elements of the 5. diagonal of [M]
!>    @param[out] lMF temporary thread local elements of the 6. diagonal of [M]
!>    @param[out] lMG temporary thread local elements of the 7. diagonal of [M]
!>    @param[out] lUD temporary thread local elements of the helper diagonal [UD]
!>    @param[out] lx temporary thread local elements of the solution vector [x]
!>    @param[out] lb temporary thread local elements of the right side [b]
!>    @param[out] lr0_hat temporary thread local elements of the random vector [r0_hat]
!>    @param[out] ldnrm temporary thread local elements of the normalisation vector
!>    @param[out] llocTMP temporary thread local elements of the local temporary vectors
!>    @param[out] ud_block block buffer for helper diagonal [UD]
!>    @param[in] ismpl local sample index
!>    @details
!>    solve of : [M] x [x] = [b]\n
!>    [M] is general Matrix, only used in 'omp_MVP'\n
!>    Technics :\n
!>               - use reverse communication technics.\n
!>                 each vector should be dense full without any hole,\n
!>                 ( you can copy your elements from your structure to a \n
!>                 temporary dense full vector, befor you use this algorithm \n
!>                 and give the correct number of elements in 'N' ).\n
!>                 if you have setup all vectors by a specific composition,\n
!>                 each vector (x,b,r,...) on the same thread should use\n
!>                 the same composition (same structure for all vectors on\n
!>                 one thread).\n
      SUBROUTINE omp_gen_solve_ilu(n_i,n_j,n_k,x,b,r0_hat,depsilon, &
          mbc_mask,max_it,criteria,ma,mb,mc,md,me,mf,mg,ud, &
          bound_block,dnrm,lma,lmb,lmc,lmd,lme,lmf,lmg,lud,lx,lb, &
          lr0_hat,ldnrm,lloctmp,ud_block,ismpl)
        use mod_OMP_TOOLS
        use mod_blocking_size
        IMPLICIT NONE
        INCLUDE 'OMP_TOOLS.inc'

!     N_I*N_J*N_K : length of all vector r,z,s,t,v,p,y,t_pc,s_pc
        INTEGER n_i, n_j, n_k, max_it, ismpl

!     vector x and b for [M]x[x]=[b]
!     res0 ^= ||res0||, start residuel, given for 'criteria=0'
        DOUBLE PRECISION x(n_i*n_j*n_k), b(n_i*n_j*n_k)
        DOUBLE PRECISION r0_hat(n_i*n_j*n_k), res0, ud(n_i*n_j*n_k)
        DOUBLE PRECISION ma(n_i*n_j*n_k), mb(n_i*n_j*n_k)
        DOUBLE PRECISION mc(n_i*n_j*n_k), md(n_i*n_j*n_k)
        DOUBLE PRECISION me(n_i*n_j*n_k), mf(n_i*n_j*n_k)
        DOUBLE PRECISION mg(n_i*n_j*n_k)
        CHARACTER mbc_mask(n_i*n_j*n_k)

!     definitions of 'work' and 'locTMP'
        INCLUDE 'pre_bicgstab.inc'
!       locTMP   : space for local vectors, using to exchange data with
!                  'matrix-vector-product'(MVP) and 'pre-conditioners'(L/R),
!                  for definitions see more in 'pre_bicgstab.inc'

!     global buffer for boundary exchange
        DOUBLE PRECISION dnrm(n_i*n_j*n_k)
        DOUBLE PRECISION bound_block(block_i*block_j+block_i*block_k+block_j*block_k,bdim_i,bdim_j,bdim_k,2)
!     private copy for preconditioning
        DOUBLE PRECISION lma(max_blocks*block_i*block_j*block_k, tlevel_1)
        DOUBLE PRECISION lmb(max_blocks*block_i*block_j*block_k, tlevel_1)
        DOUBLE PRECISION lmc(max_blocks*block_i*block_j*block_k, tlevel_1)
        DOUBLE PRECISION lmd(max_blocks*block_i*block_j*block_k, tlevel_1)
        DOUBLE PRECISION lme(max_blocks*block_i*block_j*block_k, tlevel_1)
        DOUBLE PRECISION lmf(max_blocks*block_i*block_j*block_k, tlevel_1)
        DOUBLE PRECISION lmg(max_blocks*block_i*block_j*block_k, tlevel_1)
        DOUBLE PRECISION lud(max_blocks*block_i*block_j*block_k, tlevel_1)
        DOUBLE PRECISION lx(max_blocks*block_i*block_j*block_k, tlevel_1)
        DOUBLE PRECISION lb(max_blocks*block_i*block_j*block_k, tlevel_1)
        DOUBLE PRECISION lr0_hat(max_blocks*block_i*block_j*block_k, tlevel_1)
        DOUBLE PRECISION ldnrm(max_blocks*block_i*block_j*block_k, tlevel_1)
        DOUBLE PRECISION lloctmp(max_blocks*block_i*block_j*block_k, max_loctmp,tlevel_1)
        DOUBLE PRECISION ud_block(block_i*block_j+block_i*block_k+block_j*block_k,max_blocks,tlevel_1)

!     x,y,z-grid index for each block position
        INTEGER, ALLOCATABLE :: xyz_block(:,:)
        INTEGER lxyz_block, tid

!     Pre_BiCGStab stuff
!     work     : control variable : what is to do in this subroutine,
!                see more discription in 'pre_bicgstab.inc',
!                on startup should set to 'work=START'
        INTEGER work

!     break with enough precision
        DOUBLE PRECISION depsilon, summ

        INTEGER xi, yi, zi, loc_mem, loc_memm1
        INTEGER criteria, i, j, k, fkt_proza
        EXTERNAL fkt_proza
!     openmp-private variables
        INTEGER tpos, tanz
!     openmp-shared variables
        INTEGER ipar(5), iii
        DOUBLE PRECISION, ALLOCATABLE :: rpar(:)
        INTEGER, ALLOCATABLE :: proza_lock(:,:,:)
        LOGICAL lpar(1)


!     start values
        work = start
        res0 = 1.D+99

        ALLOCATE(proza_lock(bdim_i,bdim_j,bdim_k))
!     init
        CALL set_ival(5,0,ipar)
        CALL par_reset(proza_lock)
        lpar(1) = .FALSE.

!$OMP  parallel &
!$OMP num_threads(Tlevel_1) &
!$OMP default(none) shared(Tlevel_1,ismpl) &
!$OMP shared(block_i,block_j,block_k, bdim_i,bdim_j,bdim_k,ipar,rpar) &
!$OMP shared(N_I,N_J,N_K, depsilon, max_It, criteria, res0, work,lpar) &
!$OMP shared(bound_block,MA,MB,MC,MD,ME,MF,MG,UD,x,b,r0_hat,ProzA_lock) &
!$OMP shared(mbc_mask, dnrm, ud_block, ldnrm, max_blocks) &
!$OMP shared(lMA,lMB,lMC,lMD,lUD,lME,lMF,lMG,lx,lb,lr0_hat,llocTMP) &
!$OMP private(i,j,k, tid, loc_mem,loc_memm1, lxyz_block, xyz_block) &
!$OMP private(xi,yi,zi,summ,iii, tpos, tanz)
!$      call omp_binding(ismpl)

        CALL omp_part(n_i*n_j*n_k,tpos,tanz)
        tid = omp_get_his_thread_num() + 1

!$OMP master
!       how many entries ?
        iii = 5 + 4*omp_get_num_of_threads()
        ALLOCATE(rpar(iii))
        CALL set_dval(iii,0.D0,rpar)
!$OMP end master
!     normalise the linear system, use [dnrm] to normalise the system
        CALL norm_linsys(tanz,mbc_mask(tpos),b(tpos),x(tpos),ma(tpos), &
          mb(tpos),mc(tpos),md(tpos),me(tpos),mf(tpos),mg(tpos), &
          dnrm(tpos))
!$OMP barrier
!     init ILU(0)-preconditioner, other start values
        CALL prepare_ilu(n_i,n_j,n_k,ma,mb,mc,md,me,mf,mg,ud)

!     max number of private blocks = bdim_i*bdim_j*bdim_k
        ALLOCATE(xyz_block(3,bdim_i*bdim_j*bdim_k))

        lxyz_block = 0
        DO k = 1, bdim_k
          DO j = 1, bdim_j
            DO i = 1, bdim_i
              IF (fkt_proza(i,j,k)==tid-1) THEN
!           compute number of private blocks
                lxyz_block = lxyz_block + 1
!           setup index information for each private block
                xyz_block(1,lxyz_block) = i
                xyz_block(2,lxyz_block) = j
                xyz_block(3,lxyz_block) = k
              END IF
            END DO
          END DO
        END DO
!     memory requirements
        loc_mem = lxyz_block*block_i*block_j*block_k
        loc_memm1 = max_blocks*block_i*block_j*block_k

!     make private copies from the global arrays
!$OMP barrier
        CALL lcopy_ilu(n_i,n_j,n_k,ma,lma(1,tid),lxyz_block,xyz_block)
        CALL lcopy_ilu(n_i,n_j,n_k,mb,lmb(1,tid),lxyz_block,xyz_block)
        CALL lcopy_ilu(n_i,n_j,n_k,mc,lmc(1,tid),lxyz_block,xyz_block)
        CALL lcopy_ilu(n_i,n_j,n_k,md,lmd(1,tid),lxyz_block,xyz_block)
        CALL lcopy_ilu(n_i,n_j,n_k,ud,lud(1,tid),lxyz_block,xyz_block)
        CALL lcopy_ilu(n_i,n_j,n_k,me,lme(1,tid),lxyz_block,xyz_block)
        CALL lcopy_ilu(n_i,n_j,n_k,mf,lmf(1,tid),lxyz_block,xyz_block)
        CALL lcopy_ilu(n_i,n_j,n_k,mg,lmg(1,tid),lxyz_block,xyz_block)
        CALL lcopy_ilu(n_i,n_j,n_k,x,lx(1,tid),lxyz_block,xyz_block)
        CALL lcopy_ilu(n_i,n_j,n_k,b,lb(1,tid),lxyz_block,xyz_block)
        CALL lcopy_ilu(n_i,n_j,n_k,dnrm,ldnrm(1,tid),lxyz_block,xyz_block)
        CALL lcopy_ilu(n_i,n_j,n_k,r0_hat,lr0_hat(1,tid),lxyz_block,xyz_block)
!     lcopy_ilu of [locTMP] not needed, but of the cleanup
!aw-test      Call set_dval(loC_memm1*max_loCTMP,0.d0,lloCTMP(1,1,tid))

!     copy private [UD] surface (position-1)
        DO i = 1, lxyz_block
          CALL lsurf_ilu(n_i,n_j,n_k,i,ud,ud_block(1,i,tid),lxyz_block,xyz_block)
        END DO


!     INIT
!     preload ([M]x[x]) in [z]
        CALL omp_mvp2(n_i,n_j,n_k,lxyz_block,lx(1,tid), &
          lloctmp(1,z,tid),lma(1,tid),lmb(1,tid),lmc(1,tid), &
          lmd(1,tid),lme(1,tid),lmf(1,tid),lmg(1,tid),bound_block, &
          xyz_block)

!**************************************************************

10      CONTINUE

!     BiCGStab routine
        CALL pre_bicgstab(loc_mem,lx(1,tid),lb(1,tid),lr0_hat(1,tid), &
          loc_memm1,lloctmp(1,1,tid),depsilon,ldnrm(1,tid),max_it, &
          criteria,res0,work,ipar,rpar,lpar)
!     implicit barrier here

!     preconditioner [y]:=[K^-1]x[p],
!     matrix-vector product [v]:=[M]x[y]
        IF ((work==do_y_p_v) .OR. (work==more_y_p_v)) THEN
!       left precond.
          CALL omp_lu_solve2(n_i,n_j,n_k,lxyz_block,lloctmp(1,p,tid), &
            lloctmp(1,t_pc,tid),lma(1,tid),lmb(1,tid),lmc(1,tid), &
            lud(1,tid),lme(1,tid),lmf(1,tid),lmg(1,tid), &
            ud_block(1,1,tid),bound_block,xyz_block,proza_lock)
!       needs barrier here, [LU] and [MVP] modify "bound_block"
!$OMP   barrier
!       right precond.
!       : call myPrCo(N,locTMP(1,t_pc),locTMP(1,y))
          CALL dcopy(loc_mem,lloctmp(1,t_pc,tid),1,lloctmp(1,y,tid),1)
!       MVP
          CALL omp_mvp2(n_i,n_j,n_k,lxyz_block,lloctmp(1,y,tid), &
            lloctmp(1,v,tid),lma(1,tid),lmb(1,tid),lmc(1,tid), &
            lmd(1,tid),lme(1,tid),lmf(1,tid),lmg(1,tid),bound_block, &
            xyz_block)
        END IF

!     [z]:=[M]x[x], for advanced precision
        IF (work==more_y_p_v) THEN
!       needs barrier here, both [MVP] calls modify "bound_block"
!$OMP   barrier
          CALL omp_mvp2(n_i,n_j,n_k,lxyz_block,lx(1,tid), &
            lloctmp(1,z,tid),lma(1,tid),lmb(1,tid),lmc(1,tid), &
            lmd(1,tid),lme(1,tid),lmf(1,tid),lmg(1,tid),bound_block, &
            xyz_block)
        END IF

!     preconditioner [z]:=[K^-1]x[s],
!     preconditioner [s_pc]:=[L_K^-1]x[s],
!     matrix-vector product [t]:=[M]x[z],
!     preconditioner [t_pc]:=[L_K^-1]x[t]
        IF (work==do_z_s_t) THEN
!       left precond.
          CALL omp_lu_solve2(n_i,n_j,n_k,lxyz_block,lloctmp(1,s,tid), &
            lloctmp(1,s_pc,tid),lma(1,tid),lmb(1,tid),lmc(1,tid), &
            lud(1,tid),lme(1,tid),lmf(1,tid),lmg(1,tid), &
            ud_block(1,1,tid),bound_block,xyz_block,proza_lock)
!       needs barrier here, [LU] and [MVP] modify "bound_block"
!$OMP   barrier
!       right precond.
!       : call myPrCo(N,locTMP(1,s_pc),locTMP(1,z))
          CALL dcopy(loc_mem,lloctmp(1,s_pc,tid),1,lloctmp(1,z,tid),1)
!       MVP
          CALL omp_mvp2(n_i,n_j,n_k,lxyz_block,lloctmp(1,z,tid), &
            lloctmp(1,t,tid),lma(1,tid),lmb(1,tid),lmc(1,tid), &
            lmd(1,tid),lme(1,tid),lmf(1,tid),lmg(1,tid),bound_block, &
            xyz_block)
!       needs barrier here, [MVP] and [LU] modify "bound_block"
!$OMP   barrier
!       left precond.
          CALL omp_lu_solve2(n_i,n_j,n_k,lxyz_block,lloctmp(1,t,tid), &
            lloctmp(1,t_pc,tid),lma(1,tid),lmb(1,tid),lmc(1,tid), &
            lud(1,tid),lme(1,tid),lmf(1,tid),lmg(1,tid), &
            ud_block(1,1,tid),bound_block,xyz_block,proza_lock)
        END IF


!     precision not enough ?
        IF ((work/=fine) .AND. (work/=abort)) GO TO 10
!     at "work=ABORT", we can startup with a new [r^]

!**************************************************************


!     get the global [x] from private
        CALL lcopy_bak_ilu(n_i,n_j,n_k,x,lx(1,tid),lxyz_block,xyz_block)

        DEALLOCATE(xyz_block)

!$OMP end parallel
        DEALLOCATE(rpar)
        DEALLOCATE(proza_lock)

        RETURN
      END
