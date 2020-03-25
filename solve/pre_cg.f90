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

!>    @brief compute the solution for [A]x[x]=[b], [A] sym.pos.def. matrix
!>    @param[in] N number of elements, vector length, equal for all vectors
!>    @param[in,out] x starting-vector [x0], on exit = solution-vector [x]
!>    @param[in] b right side, vector [b]
!>    @param[in] ld_lv leading Dimension of 'locTMP'
!>    @param[in,out] locTMP space for local vectors, using to exchange data with\n
!>      'matrix-vector-product'(MVP) and 'pre-conditioners'(L/R),\n
!>      for definitions see more in 'pre_bicgstab.inc'
!>    @param[in] depsilon precision criteria to break iterations
!>    @param[in] max_It maximum of iterations, counted with 'iter'
!>    @param[in] criteria  precision criteria mode to break iterations\n
!>    - 0 : relative stopping crit.: ||[res]||       < depsilon*||[res0]||\n
!>    - 1 : absolute stopping crit.: ||[res]||       < depsilon\n
!>    - 2 : maximum  stopping crit.: max(abs([res])) < depsilon\n
!>    first [res]^=[r], later (if precise enough): [res]^=([M]x[x]-[b])
!>    @param[in] res0 res0 ^= ||[res0]||, start residuel, given for 'criteria=0'
!>    @param[in,out] work control variable : what is to do out of this subroutine,\n
!>    see more discription in 'pre_bicgstab.inc',\n
!>    on startup should set to 'work=START'
!>    @param[in,out] dnrm normalisation vector, temporary use
!>    @param[out] ipar integer type help vector - OpenMP "shared"
!>    @param[out] rpar floating type point help vector - OpenMP "shared"
!>    @param[out] lpar logical type help vector - OpenMP "shared"
!>    @details
!>    compute the solution for [A]x[x]=[b], [A] sym.pos.def. matrix\n
!>    make use of [K] as a preconditioner for [A], (function : myPrCo)\n
!>    Technics :20;\n
!>               - use reverse communication technics.\n
!>                 each vector should be dense full without any hole,\n
!>                 ( you can copy your elements from your structure to a \n
!>                 temporary dense full vector, befor you use this algorithm \n
!>                 and give the correct number of elements in 'N' ).\n
!>                 if you have setup all vectors by a specific composition,\n
!>                 each vector (x,b,r,...) on the same thread should use\n
!>                 the same composition (same structure for all vectors on\n
!>                 one thread).
!>
!>    CG algorithm :\n
!>    before begin : result of [A]x[x] should given in [s]-vector (see locTMP)\n
!>      100 -----------   (... begin ...)\n
!>          -             [r]:=[b]-[s]\n
!>          -             abbruch(N,r,depsilon, ...) ? -> goto 400\n
!>          -             [p]:=[q]:=0\n
!>          -             beta:=0\n
!>          -             "solve [K][z]=[r]":\n
!>            -              myPrCo->[z]     [z]:=[K^-1]x[r]\n
!>          -             "myMVP->[s]":\n
!>            -              [s]:=[A]x[z]\n
!>
!>      200 -----------   \n
!>          -             sigma:=([r]^T*[z])\n
!>          -             mu:=([s]^T*[z])\n
!>          -             alpha:=sigma/mu\n
!>          -             [p]:=[z]+beta*[p]\n
!>          -             [q]:=[s]+beta*[q]\n
!>          -             [x]:=[x]+alpha*[p]\n
!>          -             [r]:=[r]-alpha*[q]\n
!>          -             "solve [K][z]=[r]":\n
!>            -              myPrCo->[z]     [z]:=[K^-1]x[r]\n
!>          -             "myMVP->[s]":\n
!>            -              [s]:=[A]x[z]\n
!>
!>      300 -----------   (... iterations ...)\n
!>          -             abbruch(N,r,depsilon, ...) ? -> goto 400\n
!>          -             sigma_new:=([r]^T*[z])\n
!>          -             mu:=([s]^T*[z])\n
!>          -             beta:=sigma_new/sigma\n
!>          -             sigma:=sigma_new\n
!>          -             alpha:=sigma/(mu-sigma*beta/alpha)\n
!>          -             [p]:=[z]+beta*[p]\n
!>          -             [q]:=[s]+beta*[q]\n
!>          -             [x]:=[x]+alpha*[p]\n
!>          -             [r]:=[r]-alpha*[q]\n
!>          -             "solve [K][z]=[r]":\n
!>            -              myPrCo->[z]     [z]:=[K^-1]x[r]\n
!>          -             "myMVP->[s]":\n
!>            -              [s]:=[A]x[z]\n
!>          -             -> goto 300 (... iterations ...)\n
!>
!>      400 -----------   (... end ...)\n
!>          -             -> END\n
      SUBROUTINE pre_cg(n,x,b,ld_lv,loctmp,depsilon,dnrm,max_it, &
          criteria,res0,work,ipar,rpar,lpar)
        use mod_OMP_TOOLS
        use mod_linfos
        use mod_blocking_size
        IMPLICIT NONE
        INCLUDE 'OMP_TOOLS.inc'

!     N     : number of elements of all vector
!     ld_lv : leading Dimension for vectors in 'locTMP'
        INTEGER n, ld_lv
!     vectors [x], [b]
        DOUBLE PRECISION x(max(n,1)), b(max(n,1)), dnrm(max(n,1))

!     temporary local variables
        DOUBLE PRECISION mu, sigma_new, enough, res0

!     break with enough precision
        DOUBLE PRECISION depsilon, p_e_old
        INTEGER criteria, p_e_count

!     used external functions
        LOGICAL omp_abbruch
        EXTERNAL omp_abbruch

!     definitions of 'work' and 'locTMP'
        INCLUDE 'pre_cg.inc'
        DOUBLE PRECISION loctmp(ld_lv,max_loctmp)

!     work : next Work to do
        INTEGER work

!     max_It : max iterations
        INTEGER max_it

!     first times [Ax]
        INTEGER p_first_ax

!     private only for the master thread
        INTEGER p_master

        DOUBLE PRECISION d_one, dzero
        PARAMETER (d_one=1.0D0,dzero=0.0D0)

!     blocking staff
        INTEGER von, bis

!     for definitions of r,z,s,p,q,v and locTMP see in 'pre_cg.inc'

!     p_* : temporary local private copy of *-variables
        DOUBLE PRECISION p_beta, p_sigma, p_alpha
        INTEGER p_iter, p_step, p_divide_zero
        LOGICAL p_need_ax

        INTRINSIC dsqrt

!     needing saved and shared variables
!       iter : count iterations, max_It : max iterations
!       step : jump-lable, next line to continue
!       divide_zero : control error level
!       first_Ax : first times [Ax]
!       need_Ax : switch to compute an extra MVP:([A]x[x]) in [z]
!       sh_help : openmp-shared space
        INTEGER iiter, istep, idivide_zero, ie_count, ifirst_ax
        INTEGER rbeta, rsigma, ralpha, re_old, rb_e_old, rsh_help, &
          rsh_vhelp
        INTEGER lneed_ax
        PARAMETER (iiter=1,istep=2,idivide_zero=3,ie_count=4, &
          ifirst_ax=5)
        PARAMETER (rbeta=1,rsigma=2,ralpha=3,re_old=4,rb_e_old=5)
!     sh_vhelp needs [4*#threads] entries
        PARAMETER (rsh_help=6,rsh_vhelp=7)
        PARAMETER (lneed_ax=1)
        INTEGER ipar(5)
        DOUBLE PRECISION rpar(*)
        LOGICAL lpar(1)
!     save iter, need_Ax, step, divide_zero, beta, sigma, alpha
!     save e_Count, e_old, b_e_old, first_Ax


!      shared variables should be private ...
        p_iter = ipar(iiter)
        p_need_ax = lpar(lneed_ax)
        p_step = ipar(istep)
        p_divide_zero = ipar(idivide_zero)
        p_beta = rpar(rbeta)
        p_sigma = rpar(rsigma)
        p_alpha = rpar(ralpha)
        p_e_count = ipar(ie_count)
        p_e_old = rpar(re_old)
        p_first_ax = ipar(ifirst_ax)
        enough = 2.0D0*rpar(rb_e_old)

!     jump-table
        IF (work==start) GO TO 100
        IF (p_step==200) GO TO 200
        IF (p_step==300) GO TO 300

        WRITE(*,*) &
          'error : no jump lable, startup value should work=START'
        GO TO 400


100     CONTINUE

!      need barrier here for fewer data races, because of the init above ...
!$OMP  barrier
!      initialize some 'save'-variables at the beginning of a new system
!aw    p_iter        = 0
!aw    p_need_Ax     = .false.
        p_divide_zero = 0

        p_e_count = 0
        p_e_old = 0.0D0
        rpar(rb_e_old) = 1.0D99
        enough = 2.0D0*rpar(rb_e_old)
        p_first_ax = 0

!      beta:=0
        p_beta = dzero

!      [r]:=[b]-[s]
        CALL dcopy(n,b(1),1,loctmp(1,r),1)
        CALL daxpy(n,-d_one,loctmp(1,s),1,loctmp(1,r),1)

!      abbruch(iter,r,depsilon, ...) ? -> goto 400
!      temp. use of 'enough' to prove precision
        IF (criteria==2) THEN
          CALL omp_damax(n,loctmp(1,r),enough,rpar(rsh_help))
        ELSE
          CALL omp_ddot(n,loctmp(1,r),loctmp(1,r),enough, &
            rpar(rsh_vhelp))
          res0 = dsqrt(enough)
        END IF
        p_need_ax = .TRUE.
        p_iter = 1
        IF (omp_abbruch(enough,p_iter,max_it,depsilon,p_need_ax, &
          criteria,res0,p_divide_zero,p_e_count,p_e_old)) GO TO 400

!      initialize some 'save'-variables at the beginning of a new system
        p_need_ax = .FALSE.
        p_iter = 0

!      [p]:=[q]:=0
        CALL set_dval(n,0.D0,loctmp(1,p))
        CALL set_dval(n,0.D0,loctmp(1,q))

!      mu:= [b]^T*[b], trivial solution test
        CALL omp_ddot(n,b(1),b(1),mu,rpar(rsh_vhelp))

!      mu = 0 ? aborting, trivial solution
        CALL test_zero(mu,4,p_divide_zero)
        IF (p_divide_zero/=0) GO TO 400

!      solve [K][z]=[r]:
!       myPrCo->[z]     [z]:=[K^-1]x[r]
!      myMVP->[s]      :[s]:=[A]x[z]

!      barrier above ...
!$OMP  master
        work = mvp
!$OMP  end master
        p_step = 200
!      no barrier here ... need barrier later ...

        GO TO 1000


200     CONTINUE

!      sigma:=([r]^T*[z])
!      call omp_ddot(N,locTMP(1,r),locTMP(1,z),sigma,rpar(Rsh_vhelp))

!      mu:=([s]^T*[z])
!      call omp_ddot(N,locTMP(1,s),locTMP(1,z),mu,rpar(Rsh_vhelp))

        CALL omp_2ddot(n,loctmp(1,r),loctmp(1,z),p_sigma,loctmp(1,s), &
          loctmp(1,z),mu,rpar(rsh_vhelp))
!      implicit barrier here ...

!      error prevention
        CALL test_zero(mu,1,p_divide_zero)

!      alpha:=sigma/mu
        p_alpha = p_sigma/mu

!      [p]:=[z]+beta*[p]
        CALL dscal(n,p_beta,loctmp(1,p),1)
        CALL daxpy(n,d_one,loctmp(1,z),1,loctmp(1,p),1)

!      [q]:=[s]+beta*[q]
        CALL dscal(n,p_beta,loctmp(1,q),1)
        CALL daxpy(n,d_one,loctmp(1,s),1,loctmp(1,q),1)

!      [x]:=[x]+alpha*[p]
        CALL daxpy(n,p_alpha,loctmp(1,p),1,x(1),1)

!      [r]:=[r]-alpha*[q]
        CALL daxpy(n,-p_alpha,loctmp(1,q),1,loctmp(1,r),1)

!       solve [K][z]=[r]:
!        myPrCo->[z]     [z]:=[K^-1]x[r]
!       myMVP->[s]      :[s]:=[A]x[z]

!$OMP  barrier
!      need barrier here ...
!$OMP  master
        work = mvp
!$OMP  end master
        p_step = 300
!      no barrier here ... need barrier later ...

        GO TO 1000


300     CONTINUE

!      preparing values for 'omp_abbruch'-function
!      temp. use of 'enough' to compute precision
        IF (p_need_ax) THEN
!        test absulute residuel : ([b]-[A]x[x])
!        [A]x[x] is given in [v]
!        compute ([A]x[x]-[b]) instead above
          CALL daxpy(n,-d_one,b,1,loctmp(1,v),1)
          CALL norm_resid(n,dnrm,x,loctmp(1,v))
          IF (criteria==2) THEN
            CALL omp_damax(n,loctmp(1,v),enough,rpar(rsh_help))
!          implicit barrier here ...

!          sigma_new:=([r]^T*[z])
!          mu:=([s]^T*[z])
            CALL omp_2ddot(n,loctmp(1,r),loctmp(1,z),sigma_new, &
              loctmp(1,s),loctmp(1,z),mu,rpar(rsh_vhelp))
!          implicit barrier here ...
          ELSE
!          call omp_ddot(N,locTMP(1,v),locTMP(1,v),enough,rpar(Rsh_vhelp))

!          sigma_new:=([r]^T*[z])
!          mu:=([s]^T*[z])
            CALL omp_3ddot(n,loctmp(1,v),loctmp(1,v),enough, &
              loctmp(1,r),loctmp(1,z),sigma_new,loctmp(1,s), &
              loctmp(1,z),mu,rpar(rsh_vhelp))
!          implicit barrier here ...
          END IF
        ELSE
!        test [r] instead of residuel
          IF (criteria==2) THEN
            CALL omp_damax(n,loctmp(1,r),enough,rpar(rsh_help))
!          implicit barrier here ...

!          sigma_new:=([r]^T*[z])
!          mu:=([s]^T*[z])
            CALL omp_2ddot(n,loctmp(1,r),loctmp(1,z),sigma_new, &
              loctmp(1,s),loctmp(1,z),mu,rpar(rsh_vhelp))
!          implicit barrier here ...
          ELSE
!          call omp_ddot(N,locTMP(1,r),locTMP(1,r),enough,rpar(Rsh_vhelp))

!          sigma_new:=([r]^T*[z])
!          mu:=([s]^T*[z])
            CALL omp_3ddot(n,loctmp(1,r),loctmp(1,r),enough, &
              loctmp(1,r),loctmp(1,z),sigma_new,loctmp(1,s), &
              loctmp(1,z),mu,rpar(rsh_vhelp))
!          implicit barrier here ...
          END IF
        END IF

!      abbruch(iter,r,depsilon, ...) ? -> goto 400
!      temp. use of 'enough' to prove precision
        IF (omp_abbruch(enough,p_iter,max_it,depsilon,p_need_ax, &
          criteria,res0,p_divide_zero,p_e_count,p_e_old)) GO TO 400

!      sigma_new:=([r]^T*[z])
!      is compute above ...

!      mu:=([s]^T*[z])
!      is compute above ...

!      error prevention
        CALL test_zero(p_sigma,2,p_divide_zero)

!      beta:=sigma_new/sigma
        p_beta = sigma_new/p_sigma

!      sigma:=sigma_new
        p_sigma = sigma_new

!      alpha:=sigma/(mu-sigma*beta/alpha)
        p_alpha = mu - p_sigma*p_beta/p_alpha

!      error prevention
        CALL test_zero(p_alpha,3,p_divide_zero)

        p_alpha = p_sigma/p_alpha

!      iteration counting
        p_iter = p_iter + 1


!     begin blocking
        von = 1
        bis = min(n,int(bl_size/bldiv_cg))
1       CONTINUE

!      [p]:=[z]+beta*[p]
        CALL dscal(bis-von+1,p_beta,loctmp(von,p),1)
        CALL daxpy(bis-von+1,d_one,loctmp(von,z),1,loctmp(von,p),1)

!      [q]:=[s]+beta*[q]
        CALL dscal(bis-von+1,p_beta,loctmp(von,q),1)
        CALL daxpy(bis-von+1,d_one,loctmp(von,s),1,loctmp(von,q),1)

!      [x]:=[x]+alpha*[p]
        CALL daxpy(bis-von+1,p_alpha,loctmp(von,p),1,x(von),1)

!      [r]:=[r]-alpha*[q]
        CALL daxpy(bis-von+1,-p_alpha,loctmp(von,q),1,loctmp(von,r),1)

!     end blocking
        von = bis + 1
        bis = min(n,bis+int(bl_size/bldiv_cg))
        IF (von<=n) GO TO 1


!      solve [K][z]=[r]:
!       myPrCo->[z]     [z]:=[K^-1]x[r]
!      myMVP->[s]      :[s]:=[A]x[z]

!      barrier above ... (in omp_ddot)
!$OMP  master
        work = mvp
!      need a second MVP, compute ([A]x[x]) in [z]
        IF (p_need_ax) work = mvpx
!$OMP  end master
        p_step = 300
!      no barrier here ... need barrier later ...

        GO TO 1000


400     CONTINUE

!      -> END

!     restore best [x] value
        IF (rpar(rb_e_old)<p_e_old .AND. rpar(rb_e_old)<enough) THEN
          CALL dcopy(n,loctmp(1,best),1,x,1)
!        so "b_e_old" will not changed later by the master
          p_e_old = rpar(rb_e_old)
!$OMP  master
          IF (linfos(4)>=1) THEN
            IF (criteria==2) THEN
              WRITE(*,'(1A,1e20.13,1A)',advance='NO') ', {', &
                rpar(rb_e_old), '}'
            ELSE
              WRITE(*,'(1A,1e20.13,1A)',advance='NO') ', {', &
                dsqrt(rpar(rb_e_old)), '}'
            END IF
          END IF
!$OMP  end master
        END IF

!      need barrier before ...
!$OMP  master
        work = fine
!$OMP  end master
        p_step = 0

        IF (p_divide_zero/=0) THEN
!$OMP    master
          work = abort
          IF (linfos(4)>=0) THEN
            IF (p_divide_zero==4) THEN
              WRITE(*,*) ' '
              WRITE(*,'(2A)') ' warning: trivial solution, no boundary'&
                , ' condition different from zero'
            ELSE
              WRITE(*,*) ' '
              WRITE(*,'(A,I2)') ' warning: Division by zero with code:'&
                , p_divide_zero
            END IF
          END IF
!$OMP    end master
        END IF
!      no barrier here ... need it later ...

!$OMP  master
        IF (linfos(4)>=1) THEN
!$       write(*,'(1A,1I4)',ADVANCE='NO') ', ',OMP_GET_NUM_of_THREADS()
!$       write(*,'(1A)',ADVANCE='NO') ' threads used'
          IF (p_iter>max_it) THEN
            WRITE(*,'(A,I5,A)') ', Number of Iterations: ', p_iter, &
              ':MAXIMUM (CG)'
          ELSE
            WRITE(*,'(A,I5,A)') ', Number of Iterations: ', p_iter, &
              ' (CG)'
          END IF
        END IF
!$OMP  end master

1000    CONTINUE

        IF (p_need_ax .AND. p_divide_zero==0 .AND. p_first_ax>=0) &
          p_first_ax = p_first_ax + 1

!     save best [x] value
        p_master = 0
        IF (p_e_old<rpar(rb_e_old) .OR. p_first_ax==3) THEN
!       save [x]
          CALL dcopy(n,x,1,loctmp(1,best),1)
          IF (p_first_ax==3) p_first_ax = -1
!$OMP   master
          p_master = 1
!$OMP   end master
        END IF

!      Danger : here is an important syncronistaion for consistent view
!               of 'work' and all save variables on exit ...
!$OMP  barrier

!$OMP  master
        IF (p_master==1) rpar(rb_e_old) = p_e_old
!      private variables should be shared to save ...
        ipar(iiter) = p_iter
        lpar(lneed_ax) = p_need_ax
        ipar(istep) = p_step
        ipar(idivide_zero) = p_divide_zero
        rpar(rbeta) = p_beta
        rpar(rsigma) = p_sigma
        rpar(ralpha) = p_alpha
        ipar(ie_count) = p_e_count
        rpar(re_old) = p_e_old
        ipar(ifirst_ax) = p_first_ax
!$OMP  end master

        RETURN
      END
