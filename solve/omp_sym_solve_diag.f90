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

!>    @brief solve of : [M] x [x] = [b], CG algorithm based with Diagonal preconditioning
!>    @param[in] N_I lengths of I dimension of local matrix [M]
!>    @param[in] N_J lengths of J dimension of local matrix [M]
!>    @param[in] N_K lengths of K dimension of local matrix [M]
!>    @param[in,out] x solution vector [x], on start = start vector
!>    @param[in] b right side, vector [b]
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
!>    [M] is s.p.d. Matrix, only used in 'omp_MVP'\n
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
      SUBROUTINE omp_sym_solve_diag(n_i,n_j,n_k,x,b,depsilon,mbc_mask, &
          max_it,criteria,ma,mb,mc,md,me,mf,mg,loctmp,dnrm,ismpl)
        use mod_OMP_TOOLS
        use mod_blocking_size
        IMPLICIT NONE
        INCLUDE 'OMP_TOOLS.inc'

!     N    : length of all vector x, r, z, s, b, p, q
!     ld_N : leading dimension of [locTMP]
        INTEGER n, n_i, n_j, n_k, max_it, ismpl

!     thread stuff
        INTEGER tanz, tpos

!     vector x and b for [M]x[x]=[b]
!     res0 ^= ||res0||, start residuel, given for 'criteria=0'
        DOUBLE PRECISION x(n_i*n_j*n_k), b(n_i*n_j*n_k)
        DOUBLE PRECISION res0
        CHARACTER mbc_mask(n_i*n_j*n_k)
        DOUBLE PRECISION ma(n_i*n_j*n_k), mb(n_i*n_j*n_k), &
          mc(n_i*n_j*n_k)
        DOUBLE PRECISION md(n_i*n_j*n_k), me(n_i*n_j*n_k), &
          mf(n_i*n_j*n_k)
        DOUBLE PRECISION mg(n_i*n_j*n_k)

!     definitions of 'work' and 'locTMP'
        INCLUDE 'pre_cg.inc'
!       locTMP   : space for local vectors, using to exchange data with
!                  'matrix-vector-product'(MVP) and 'pre-conditioners'(L/R),
!                  for definitions see more in 'pre_bicgstab.inc'
        DOUBLE PRECISION loctmp(n_i*n_j*n_k,max_loctmp)
        DOUBLE PRECISION dnrm(n_i*n_j*n_k)

!     Pre_CG stuff
!     work     : control variable : what is to do out of this subroutine,
!                see more discription in 'pre_cg.inc',
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
        iii = 6 + 4*omp_get_num_of_threads()
        ALLOCATE(rpar(iii))
        CALL set_dval(iii,0.D0,rpar)
!      allocate(rpar(6 +4 *OMP_GET_NUM_of_THREADS()))
!$OMP end master
!     normalise the linear system, use [dnrm] to normalise the system
        CALL norm_linsys(tanz,mbc_mask(tpos),b(tpos),x(tpos),ma(tpos), &
          mb(tpos),mc(tpos),md(tpos),me(tpos),mf(tpos),mg(tpos), &
          dnrm(tpos))
!$OMP barrier

!     preload ([M]x[x]) in [s]
        CALL omp_mvp(n_i,n_j,n_k,x,loctmp(1,s),ma,mb,mc,md,me,mf,mg)
!     implicit barrier here

10      CONTINUE


!     cg-routines without Pre-Cond.
!     step : is the return-lable for pre_cg()
        CALL pre_cg(tanz,x(tpos),b(tpos),n,loctmp(tpos,1),depsilon, &
          dnrm(tpos),max_it,criteria,res0,work,ipar,rpar,lpar)
!     implicit barrier here

        IF ((work==mvp) .OR. (work==mvpx)) THEN
!       solve: [M]x[z]=[r] -> [z]:=[A^-1]x[r]
!       here [M] can be a substituted Matrix [T] -> solve: [T]x[z]=[r]
          CALL diagprco(n,md,loctmp(1,r),loctmp(1,z))
!$OMP   barrier

!       [s]:=[M]x[z]
          CALL omp_mvp(n_i,n_j,n_k,loctmp(1,z),loctmp(1,s),ma,mb,mc, &
            md,me,mf,mg)
!       implicit barrier here
        END IF

!     [v]:=[M]x[x], for advanced precision
        IF (work==mvpx) CALL omp_mvp(n_i,n_j,n_k,x,loctmp(1,v),ma,mb, &
          mc,md,me,mf,mg)
!     implicit barrier here


!     precision not enough ?
        IF ((work/=fine) .AND. (work/=abort)) GO TO 10

!$OMP  end parallel

        DEALLOCATE(rpar)

!**************************************************************


        RETURN
      END
