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

!>    @brief solve of : [Bayes] x [x] = [b], BICGSTAB algorithm based
!>    @param[in] N dimension of the matrix [Bayes]
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
!>    @param[out] locTMP local temporary vectors
!>    @param[out] DATTMP additional data temporary vectors
!>    @param[in,out] dnrm normalisation vector, temporary use
!>    @param[in] ismpl local sample index
!>    @details
!>    solve of : [Bayes] x [x] = [b]\n
!>    [Bayes] is a indirect defined Matrix\n
!>    Technics :\n
!>      - use reverse communication technics.\n
!>        each vector should be dense full without any hole,\n
!>        ( you can copy your elements from your structure to a \n
!>        temporary dense full vector, befor you use this algorithm \n
!>        and give the correct number of elements in 'N' ).\n
!>        if you have setup all vectors by a specific composition,\n
!>        each vector (x,b,r,...) on the same thread should use\n
!>        the same composition (same structure for all vectors on\n
!>        one thread).\n
      SUBROUTINE omp_bayes_solve(N,x,b,r0_hat,depsilon,max_it,criteria,loctmp,dattmp,dnrm,ismpl)
        use mod_OMP_TOOLS
        use mod_blocking_size
        use mod_linfos
        IMPLICIT NONE
        INCLUDE 'OMP_TOOLS.inc'
!     N     : length of all vector r,z,s,t,v,p,y,t_pc,s_pc
        INTEGER N, max_it, ismpl
!     thread stuff
        INTEGER tpos, tanz
!     vector x and b for [M]x[x]=[b]
!     res0 ^= ||res0||, start residuel, given for 'criteria=0'
        DOUBLE PRECISION x(N), b(N)
        DOUBLE PRECISION r0_hat(N), res0
!     definitions of 'work' and 'locTMP'
        INCLUDE 'pre_bicgstab.inc'
!       locTMP   : space for local vectors, using to exchange data with
!                  'matrix-vector-product'(MVP) and 'pre-conditioners'(L/R),
!                  for definitions see more in 'pre_bicgstab.inc'
        DOUBLE PRECISION loctmp(N,max_loctmp)
        DOUBLE PRECISION dattmp(*), dnrm(N)
!     Pre_BiCGStab stuff
!     work     : control variable : what is to do in this subroutine,
!                see more discription in 'pre_bicgstab.inc',
!                on startup should set to 'work=START'
        INTEGER work
!     break with enough precision
        DOUBLE PRECISION depsilon
!
        INTEGER criteria, linfos_old
        INTEGER l, wieviele
        DOUBLE PRECISION dmaxi, dmini, durchs
!     openmp-shared variables
        INTEGER ipar(5), iii
        DOUBLE PRECISION, ALLOCATABLE :: rpar(:)
        LOGICAL lpar(1)

!
!     full number of elements
        res0 = 1.D+99

!**************************************************************

!     start values
        work = start

!$OMP parallel &
!$OMP num_threads(Tlevel_0)&
!$OMP default(shared)&
!$OMP private(tanz,tpos,iii,linfos_old)
!$      call omp_binding(ismpl)

        CALL omp_part(n,tpos,tanz)

!$OMP master
        iii = 5 + 4*omp_get_num_of_threads()
        ALLOCATE(rpar(iii))
        CALL set_dval(iii,0.D0,rpar)
!$OMP end master
        linfos_old = linfos(4)
!$OMP barrier

!     preload ([M]x[x]) in [z]
!$OMP master
        CALL step_bayes_tmpm(x,dattmp,loctmp(1,z),ismpl)
!$OMP end master
!$OMP barrier

! Reverse communication loop for bicgstab
10      CONTINUE


!     BiCGStab routine
        linfos(4) = 1
        CALL pre_bicgstab(tanz,x(tpos),b(tpos),r0_hat(tpos),n, &
          loctmp(tpos,1),depsilon,dnrm(tpos),max_it,criteria,res0, &
          work,ipar,rpar,lpar)
!     implicite barrier here
        linfos(4) = linfos_old
!$OMP barrier

!     preconditioner [y]:=[K^-1]x[p],
!     matrix-vector product [v]:=[M]x[y]
        IF ((work==do_y_p_v) .OR. (work==more_y_p_v)) THEN
!       left precond.
          CALL myprco(n,loctmp(1,p),loctmp(1,t_pc))
!       right precond.
          CALL myprco(n,loctmp(1,t_pc),loctmp(1,y))
!$OMP   barrier

!$OMP   master
          CALL step_bayes_tmpm(loctmp(1,y),dattmp,loctmp(1,v),ismpl)
!$OMP   end master
!$OMP   barrier
        END IF

!     [z]:=[M]x[x], for advanced precision
        IF (work==more_y_p_v) THEN
!$OMP     master
            CALL step_bayes_tmpm(x,dattmp,loctmp(1,z),ismpl)
!$OMP     end master
!$OMP     barrier
        END IF

!     preconditioner [z]:=[K^-1]x[s],
!     preconditioner [s_pc]:=[L_K^-1]x[s],
!     matrix-vector product [t]:=[M]x[z],
!     preconditioner [t_pc]:=[L_K^-1]x[t]
        IF (work==do_z_s_t) THEN
!       left precond.
          CALL myprco(n,loctmp(1,s),loctmp(1,s_pc))
!       right precond.
          CALL myprco(n,loctmp(1,s_pc),loctmp(1,z))
!$OMP     barrier

!$OMP     master
            CALL step_bayes_tmpm(loctmp(1,z),dattmp,loctmp(1,t),ismpl)
!$OMP     end master
!$OMP     barrier

!       left precond.
          CALL myprco(n,loctmp(1,t),loctmp(1,t_pc))
!$OMP   barrier
        END IF


!     precision not enough ?
        IF ((work/=fine) .AND. (work/=abort)) GO TO 10
!     at "work=ABORT", we can startup with a new [r^]

!$OMP end parallel

        DEALLOCATE(rpar)

!**************************************************************


        RETURN
      END
