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

!>    @brief solve of : [M] x [x] = [b], BICGSTAB algorithm based with SSOR preconditioning
!>    @param[in] N_I lengths of I dimension of local matrix [M]
!>    @param[in] N_J lengths of J dimension of local matrix [M]
!>    @param[in] N_K lengths of K dimension of local matrix [M]
!>    @param[in,out] x solution vector [x], on start = start vector
!>    @param[in] b right side, vector [b]
!>    @param[in] r0_hat random vector [r0_hat] ^= [r0^]
!>    @param[in] depsilon precision criteria to break iterations
!>    @param[in] max_It max iteration number
!>    @param[in] criteria precision criteria mode to break iterations\n
!>    - 0 : relative stopping crit.: ||[res]||       < depsilon*||[res0]||\n
!>    - 1 : absolute stopping crit.: ||[res]||       < depsilon\n
!>    - 2 : maximum  stopping crit.: max(abs([res])) < depsilon\n
!>    first [res]^=[r], later (if precise enough): [res]^=([M]x[x]-[b])
!>    @param[in] mbc_mask boundary condition pattern (mask)
!>    @param[in] MA 1. diagonal of the system matrix [M]
!>    @param[in] MB 2. diagonal of the system matrix [M]
!>    @param[in] MC 3. diagonal of the system matrix [M]
!>    @param[in] MD 4. diagonal of the system matrix [M]
!>    @param[in] ME 5. diagonal of the system matrix [M]
!>    @param[in] MF 6. diagonal of the system matrix [M]
!>    @param[in] MG 7. diagonal of the system matrix [M]
!>    @param[out] locTMP local temporary vectors
!>    @param[out] dnrm normalisation vector, temporary use
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
      SUBROUTINE omp_gen_solve_ssor(n_i,n_j,n_k,x,b,r0_hat,depsilon, &
          mbc_mask,max_it,criteria,ma,mb,mc,md,me,mf,mg,loctmp,dnrm, &
          ismpl)
        use mod_OMP_TOOLS
        use mod_blocking_size
        IMPLICIT NONE
        INCLUDE 'OMP_TOOLS.inc'

!     N     : length of all vector r,z,s,t,v,p,y,t_pc,s_pc
        INTEGER n, n_i, n_j, n_k, max_it, ismpl

!     thread stuff
        INTEGER tpos, tanz

!     vector x and b for [M]x[x]=[b]
!     res0 ^= ||res0||, start residuel, given for 'criteria=0'
        DOUBLE PRECISION x(n_i*n_j*n_k), b(n_i*n_j*n_k)
        DOUBLE PRECISION r0_hat(n_i*n_j*n_k), res0
        CHARACTER mbc_mask(n_i*n_j*n_k)
        DOUBLE PRECISION ma(n_i*n_j*n_k), mb(n_i*n_j*n_k), &
          mc(n_i*n_j*n_k)
        DOUBLE PRECISION md(n_i*n_j*n_k), me(n_i*n_j*n_k), &
          mf(n_i*n_j*n_k)
        DOUBLE PRECISION mg(n_i*n_j*n_k)

!     definitions of 'work' and 'locTMP'
        INCLUDE 'pre_bicgstab.inc'
!       locTMP   : space for local vectors, using to exchange data with
!                  'matrix-vector-product'(MVP) and 'pre-conditioners'(L/R),
!                  for definitions see more in 'pre_bicgstab.inc'
        DOUBLE PRECISION loctmp(n_i*n_j*n_k,max_loctmp)
        DOUBLE PRECISION dnrm(n_i*n_j*n_k)

!     Pre_BiCGStab stuff
!     work     : control variable : what is to do in this subroutine,
!                see more discription in 'pre_bicgstab.inc',
!                on startup should set to 'work=START'
        INTEGER work

!     break with enough precision
        DOUBLE PRECISION depsilon

        INTEGER criteria
!     openmp-shared variables
        INTEGER ipar(5), iii
        DOUBLE PRECISION, ALLOCATABLE :: rpar(:)
        LOGICAL lpar(1)


!     full number of elements
        n = n_i*n_j*n_k
        res0 = 1.D+99


!**************************************************************

!     start values
        work = start

!$OMP  parallel &
!$OMP num_threads(Tlevel_1) &
!$OMP default(shared) &
!$OMP private(tanz,tpos,iii)
!$      call omp_binding(ismpl)

        CALL omp_part(n,tpos,tanz)

!$OMP master
        iii = 5 + 4*omp_get_num_of_threads()
        ALLOCATE(rpar(iii))
        CALL set_dval(iii,0.D0,rpar)
!      allocate(rpar(5 +4 *OMP_GET_NUM_of_THREADS()))
!$OMP end master
!     normalise the linear system, use [dnrm] to normalise the system
        CALL norm_linsys(tanz,mbc_mask(tpos),b(tpos),x(tpos),ma(tpos), &
          mb(tpos),mc(tpos),md(tpos),me(tpos),mf(tpos),mg(tpos), &
          dnrm(tpos))
!$OMP barrier

!     ###################### not parallel Code !!! #########################
!     prepare [b^] from [b]
!     [b^] = D*(D+L)^(-1) * [b]
        CALL ddl(n_i,n_j,n_k,b,loctmp(1,b_hat),ma,mb,mc,md)
!     implicit barrier here

!     prepare [r0~] from [r0^]
!     [r0~] = D*(D+L)^(-1) * [r0^]
        CALL ddl(n_i,n_j,n_k,r0_hat,loctmp(1,r0_tilde),ma,mb,mc,md)
!     implicit barrier here
!     ################### above not parallel Code !!! ######################

!     preload ([M]x[x]) in [z]
        CALL ssor_mvp_single(n_i,n_j,n_k,x,loctmp(1,z),loctmp(1,mt), &
          ma,mb,mc,md,me,mf,mg)
!     implicit barrier here


10      CONTINUE


!     BiCGStab routine
        CALL pre_bicgstab(tanz,x(tpos),loctmp(tpos,b_hat), &
          loctmp(tpos,r0_tilde),n,loctmp(tpos,1),depsilon,dnrm(tpos), &
          max_it,criteria,res0,work,ipar,rpar,lpar)
!     implicit barrier here

!     preconditioner [y]:=[K^-1]x[p],
!     matrix-vector product [v]:=[M]x[y]
        IF ((work==do_y_p_v) .OR. (work==more_y_p_v)) THEN
!       Left precond.
          CALL myprco(n,loctmp(1,p),loctmp(1,t_pc))
!       Right precond.
          CALL myprco(n,loctmp(1,t_pc),loctmp(1,y))
!$OMP   barrier

          CALL ssor_mvp_single(n_i,n_j,n_k,loctmp(1,y),loctmp(1,v), &
            loctmp(1,mt),ma,mb,mc,md,me,mf,mg)
!       impliciete barrier here
        END IF

!     [z]:=[M]x[x], for advanced precision
        IF (work==more_y_p_v) CALL ssor_mvp_single(n_i,n_j,n_k,x, &
          loctmp(1,z),loctmp(1,mt),ma,mb,mc,md,me,mf,mg)
!     impliciete barrier here

!     preconditioner [z]:=[K^-1]x[s],
!     preconditioner [s_pc]:=[L_K^-1]x[s],
!     matrix-vector product [t]:=[M]x[z],
!     preconditioner [t_pc]:=[L_K^-1]x[t]
        IF (work==do_z_s_t) THEN
!       Left precond.
          CALL myprco(n,loctmp(1,s),loctmp(1,s_pc))
!       Right precond.
          CALL myprco(n,loctmp(1,s_pc),loctmp(1,z))
!$OMP   barrier

          CALL ssor_mvp_single(n_i,n_j,n_k,loctmp(1,z),loctmp(1,t), &
            loctmp(1,mt),ma,mb,mc,md,me,mf,mg)
!       implicit barrier here

!       Left precond.
          CALL myprco(n,loctmp(1,t),loctmp(1,t_pc))
!$OMP   barrier
        END IF


!     precision not enough ?
        IF ((work/=fine) .AND. (work/=abort)) GO TO 10
!     at "work=ABORT", we can startup with a new [r^]

!$OMP end parallel

!     ###################### not parallel Code !!! #########################
!     compute [x]
!     [x] = (D+U)^(-1) * [x^]
        CALL du(n_i,n_j,n_k,x,x,md,me,mf,mg)
!     ################### above not parallel Code !!! ######################

        DEALLOCATE(rpar)

!**************************************************************


        RETURN
      END
