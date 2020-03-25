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

!>    @brief compute the solution for [A]x[x]=[b], [A] general matrix
!>    @param[in] N number of elements, vector length, equal for all vectors
!>    @param[in,out] x starting-vector [x0], on exit = solution-vector [x]
!>    @param[in] b right side, vector [b]
!>    @param[in] r0_hat random vector [r0_hat] ^= [r0^]
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
!>    compute the solution for [A]x[x]=[b], [A] general matrix\n
!>    make use of [K]=[L_K]*[R_K], [K] is preconditioner for [A]\n
!>    -  [L_K]    : left preconditioner (function : L_PrCo)\n
!>    -  [R_K]    : right preconditioner (function : R_PrCo)
!>
!>    Technics :\n
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
!>    BiCGStab algorithm :\n
!>    before begin : result of [A]x[x] should given in [z]-vector (see locTMP)\n
!>    100 -----------   (... begin ...)\n
!>          -            [r]:=[b]-[z] ^= [b]-[A]x[x]\n
!>          -            abbruch(N,r,z,depsilon, ...) ? -> goto 400\n
!>          -            sigma:=[r0_hat]^T*[r]\n
!>          -            sigma = 0 ? -> aborting\n
!>          -            rho:=alpha:=omega=1\n
!>          -            [v]:=[p]:=0\n
!>          -            beta:=(sigma/rho)*(alpha/omega)\n
!>          -            rho:=sigma\n
!>          -            [p]:=[r]+beta*([p]-omega*[v])
!>          -            "solve [K]x[y]=[p]":\n
!>            -             L_PrCo->[t_pc]   [t_pc]:=[L_K^-1]x[p]\n
!>            -             R_PrCo->[y]      [y]:=[R_K^-1]x[t_pc]\n
!>          -            "myMVP->[v]":\n
!>            -             [v]:=[A]x[y]\n
!>
!>    200 -----------   (... iterations ...)\n
!>          -            abbruch(N,r,z,depsilon, ...) ? -> goto 400\n
!>          -            sigma:=[r0_hat]^T*[v]\n
!>          -            alpha:=rho/sigma\n
!>          -            [s]:=[r]-alpha*[v]\n
!>          -            "solve [K]x[z]=[s]":\n
!>            -             L_PrCo->[s_pc]   [s_pc]:=[L_K^-1]x[s]\n
!>            -             R_PrCo->[z]      [z]:=[R_K^-1]x[s_pc]\n
!>          -            "myMVP->[t]":\n
!>            -             [t]:=[A]x[z]\n
!>          -            L_PrCo->[t_pc]   :[t_pc]:=[L_K^-1]x[t]\n
!>
!>    300 -----------   \n
!>          -            omega:=[t_pc]^T*[s_pc]\n
!>          -            sigma:=[t_pc]^T*[t_pc]\n
!>          -            d1:=[r0_hat]^T*[s]\n
!>          -            d2:=[r0_hat]^T*[t]\n
!>          -            omega:=omega/sigma\n
!>          -            [x]:=[x]+alpha*[y]+omega*[z]\n
!>          -            [r]:=[s]-omega*[t]\n
!>          -            sigma:=d1-omega*d2\n
!>          -            sigma = 0 ? -> aborting\n
!>          -            beta:=(sigma/rho)*(alpha/omega)\n
!>          -            rho:=sigma\n
!>          -            [p]:=[r]+beta*([p]-omega*[v])\n
!>          -            "solve [K]x[y]=[p]":\n
!>            -             L_PrCo->[t_pc]   [t_pc]:=[L_K^-1]x[p]\n
!>            -             R_PrCo->[y]      [y]:=[R_K^-1]x[t_pc]\n
!>          -            "myMVP->[v]":\n
!>            -             [v]:=[A]x[y]\n
!>          -            -> goto 200 (... iterations ...)\n
!>
!>    400 -----------   (... end ...)\n
!>          -            -> END\n
      SUBROUTINE pre_bicgstab(n,x,b,r0_hat,ld_lv,loctmp,depsilon,dnrm, &
          max_it,criteria,res0,work,ipar,rpar,lpar)
        use mod_OMP_TOOLS
        use mod_linfos
        use mod_blocking_size
        IMPLICIT NONE
        INCLUDE 'OMP_TOOLS.inc'

!     N     : number of elements of all vector
!     ld_lv : leading Dimension for vectors in 'locTMP'
        INTEGER n, ld_lv
!     vectors [x], [b], [r0_hat]
        DOUBLE PRECISION x(max(n,1)), b(max(n,1))
        DOUBLE PRECISION r0_hat(max(n,1)), dnrm(max(n,1))

!     temporary local variables
        DOUBLE PRECISION sigma, beta, omega, res0, d1, d2, de

!     break with enough precision
        DOUBLE PRECISION depsilon, p_e_old
        INTEGER criteria, p_e_count

!     used external functions
        LOGICAL omp_abbruch
        EXTERNAL omp_abbruch

!     definitions of 'work' and 'locTMP'
        INCLUDE 'pre_bicgstab.inc'
        DOUBLE PRECISION loctmp(ld_lv,max_loctmp)

!     work : next Work to do
        INTEGER work

!     max_It : max iterations
        INTEGER max_it

!     first times [Ax]
        INTEGER p_first_ax

!     private only for the master thread
        INTEGER p_master

        DOUBLE PRECISION d_one
        PARAMETER (d_one=1.0D0)

!     blocking staff
        INTEGER von, bis

!     for definitions of r,z,s,t,v,p,y,t_pc,s_pc and locTMP see in 'pre_bicgstab.inc'

!     p_*  : temporary local private copy of *-variables
        DOUBLE PRECISION p_rho, p_alpha
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
        INTEGER rrho, ralpha, re_old, rb_e_old, rsh_help, rsh_vhelp
        INTEGER lneed_ax
        PARAMETER (iiter=1,istep=2,idivide_zero=3,ie_count=4, &
          ifirst_ax=5)
        PARAMETER (rrho=1,ralpha=2,re_old=3,rb_e_old=4)
!     sh_vhelp needs [4*#threads] entries
        PARAMETER (rsh_help=5,rsh_vhelp=6)
        PARAMETER (lneed_ax=1)
        INTEGER ipar(5)
        DOUBLE PRECISION rpar(*)
        LOGICAL lpar(1)
!     save iter, rho, alpha, need_Ax, step, divide_zero, e_count
!     save e_old, b_e_old, first_Ax


!      shared variables should be private ...
        p_iter = ipar(iiter)
        p_rho = rpar(rrho)
        p_alpha = rpar(ralpha)
        p_need_ax = lpar(lneed_ax)
        p_step = ipar(istep)
        p_divide_zero = ipar(idivide_zero)
        p_e_count = ipar(ie_count)
        p_e_old = rpar(re_old)
        p_first_ax = ipar(ifirst_ax)
        de = 2.0D0*rpar(rb_e_old)

!     jump-table
        IF (work==start) GO TO 100
        IF (p_step==200) GO TO 200
        IF (p_step==300) GO TO 300

        WRITE(*,'(A,I3,A)') ' error : no jump label ', p_step, &
          ', startup value should be "work=START" !'
        GO TO 400


100     CONTINUE

!      need barrier here for fewer data races, because of the init part above
!$OMP  barrier
!      initialize some 'save'-variables at the beginning of a new system
!aw    p_iter        = 0
!aw    p_need_Ax     = .false.
        p_divide_zero = 0

        p_e_count = 0
        p_e_old = 0.0D0
        rpar(rb_e_old) = 1.0D99
        de = 2.0D0*rpar(rb_e_old)
        p_first_ax = 0

!      rho:=alpha:=omega=1
!      p_rho         = d_one
        p_alpha = d_one
        omega = d_one

!      [r]:=[b]-[z] ^= [b]-[A]x[x]
        CALL dcopy(n,b(1),1,loctmp(1,r),1)
        CALL daxpy(n,-d_one,loctmp(1,z),1,loctmp(1,r),1)

!      abbruch(N,r,z,depsilon, ...) ? -> goto 400
!      temp. use of 'de' to prove precision
        IF (criteria==2) THEN
          CALL omp_damax(n,loctmp(1,r),de,rpar(rsh_help))
        ELSE
          CALL omp_ddot(n,loctmp(1,r),loctmp(1,r),de,rpar(rsh_vhelp))
!$OMP    master
          res0 = dsqrt(de)
!$OMP    end master
!        need barrier here, because of the next DDOT (rpar-reuse)
!$OMP    barrier
        END IF
        p_need_ax = .TRUE.
        p_iter = 1

        IF (omp_abbruch(de,p_iter,max_it,depsilon,p_need_ax,criteria, &
          res0,p_divide_zero,p_e_count,p_e_old)) GO TO 400

!      initialize some 'save'-variables at the beginning of a new system
        p_need_ax = .FALSE.
        p_iter = 0

!      sigma:=[r0_hat]^T*[r]
!      d1:=[b]^T*[b]; trivial solution test
        CALL omp_2ddot(n,r0_hat(1),loctmp(1,r),sigma,b(1),b(1),d1, &
          rpar(rsh_vhelp))
!      impliCit barrier here ...

!      sigma = 0 ? -> aborting, error
        CALL test_zero(sigma,1,p_divide_zero)
!      d1 = 0 ? -> aborting, trivial solution
        CALL test_zero(d1,7,p_divide_zero)
        IF (p_divide_zero/=0) GO TO 400

!      beta:=(sigma/rho)*(alpha/omega)
        beta = sigma

!      rho:=sigma
        p_rho = sigma

!      [v]:=[p]:=0
        CALL set_dval(n,0.D0,loctmp(1,v))
        CALL set_dval(n,0.D0,loctmp(1,p))

!      [p]:=[r]+beta*([p]-omega*[v])
        CALL daxpy(n,-omega,loctmp(1,v),1,loctmp(1,p),1)
        CALL dscal(n,beta,loctmp(1,p),1)
        CALL daxpy(n,d_one,loctmp(1,r),1,loctmp(1,p),1)

!      solve [K]x[y]=[p]:
!       L_PrCo->[t_pc]   [t_pc]:=[L_K^-1]x[p]
!       R_PrCo->[y]      [y]:=[R_K^-1]x[t_pc]
!      myMVP->[v]       :[v]:=[A]x[y]

!      barrier above ...
!$OMP  master
        work = do_y_p_v
!$OMP  end master
        p_step = 200
!      no barrier here ... need barrier later ...

        GO TO 1000


200     CONTINUE

!     prepare values for 'omp_abbruch'-function
!     temp. use of 'de' to compute precision
!      test absulute residuel : ([b]-[A]x[x])
!      [A]x[x] is given in [z]
        IF (p_need_ax) THEN
          CALL daxpy(n,-d_one,b,1,loctmp(1,z),1)
          CALL norm_resid(n,dnrm,x,loctmp(1,z))
        END IF
!      impliCit barrier here ...
        IF (criteria==2) THEN
          IF (p_need_ax) THEN
            CALL omp_damax(n,loctmp(1,z),de,rpar(rsh_help))
          ELSE
            CALL omp_damax(n,loctmp(1,r),de,rpar(rsh_help))
          END IF
!        impliCit barrier here ...
!        sigma:=[r0_hat]^T*[v]
          CALL omp_ddot(n,r0_hat(1),loctmp(1,v),sigma,rpar(rsh_vhelp))
!        impliCit barrier here ...
        ELSE
!        sigma:=[r0_hat]^T*[v]
!        call omp_ddot(N,r0_hat(1),locTMP(1,v),sigma,rpar(Rsh_vhelp))
          IF (p_need_ax) THEN
            CALL omp_2ddot(n,loctmp(1,z),loctmp(1,z),de,r0_hat(1), &
              loctmp(1,v),sigma,rpar(rsh_vhelp))
!          impliCit barrier here ...
          ELSE
            CALL omp_2ddot(n,loctmp(1,r),loctmp(1,r),de,r0_hat(1), &
              loctmp(1,v),sigma,rpar(rsh_vhelp))
!          impliCit barrier here ...
          END IF
        END IF

!      abbruch(N,r,z,depsilon, ...) ? -> goto 400
!      temp. use of 'de' to prove precision
        IF (omp_abbruch(de,p_iter,max_it,depsilon,p_need_ax,criteria, &
          res0,p_divide_zero,p_e_count,p_e_old)) GO TO 400

!      sigma:=[r0_hat]^T*[v]
!      is compute above ...

!      error prevention
        CALL test_zero(sigma,2,p_divide_zero)

!      alpha:=rho/sigma
        p_alpha = p_rho/sigma


!     begin blocking
        von = 1
        bis = min(n,int(bl_size/bldiv_bicg(1)))
1       CONTINUE

!      [s]:=[r]-alpha*[v]
        CALL dcopy(bis-von+1,loctmp(von,r),1,loctmp(von,s),1)
        CALL daxpy(bis-von+1,-p_alpha,loctmp(von,v),1,loctmp(von,s),1)

!     end blocking
        von = bis + 1
        bis = min(n,bis+int(bl_size/bldiv_bicg(1)))
        IF (von<=n) GO TO 1


!      solve [K]x[z]=[s]:
!       L_PrCo->[s_pc]   [s_pc]:=[L_K^-1]x[s]
!       R_PrCo->[z]      [z]:=[R_K^-1]x[s_pc]
!      myMVP->[t]       :[t]:=[A]x[z]
!      L_PrCo->[t_pc]   :[t_pc]:=[L_K^-1]x[t]

!      barrier above ... (in omp_ddot)
!$OMP  master
        work = do_z_s_t
!$OMP  end master
        p_step = 300
!      no barrier here ... need barrier later ...

        GO TO 1000


300     CONTINUE

!      omega:=[t_pc]^T*[s_pc]
!      call omp_ddot(N,locTMP(1,t_pc),locTMP(1,s_pc),omega,rpar(Rsh_vhelp))

!      sigma:=[t_pc]^T*[t_pc]
!      call omp_ddot(N,locTMP(1,t_pc),locTMP(1,t_pc),sigma,rpar(Rsh_vhelp))

!      d1:=[r0_hat]^T*[s]
!      call omp_ddot(N,r0_hat(1),locTMP(1,s),d1,rpar(Rsh_vhelp))

!      d2:=[r0_hat]^T*[t]
!      call omp_ddot(N,r0_hat(1),locTMP(1,t),d2,rpar(Rsh_vhelp))


        CALL omp_4ddot(n,loctmp(1,t_pc),loctmp(1,s_pc),omega, &
          loctmp(1,t_pc),loctmp(1,t_pc),sigma,r0_hat(1),loctmp(1,s), &
          d1,r0_hat(1),loctmp(1,t),d2,rpar(rsh_vhelp))
!      impliCit barrier here ...

!      error prevention
        CALL test_zero(sigma,3,p_divide_zero)

!      omega:=omega/sigma
        omega = omega/sigma

!      change computing order ...

!      sigma:=d1-omega*d2
        sigma = d1 - omega*d2

!      sigma = 0 ? -> aborting
        CALL test_zero(sigma,4,p_divide_zero)

!      error prevention
        CALL test_zero(p_rho,5,p_divide_zero)
        CALL test_zero(omega,6,p_divide_zero)

!      beta:=(sigma/rho)*(alpha/omega)
        beta = (sigma*p_alpha)/(p_rho*omega)


!     begin blocking
        von = 1
        bis = min(n,int(bl_size/bldiv_bicg(2)))
2       CONTINUE

!      [x]:=[x]+alpha*[y]+omega*[z]
        CALL daxpy(bis-von+1,p_alpha,loctmp(von,y),1,x(von),1)
        CALL daxpy(bis-von+1,omega,loctmp(von,z),1,x(von),1)

!      [r]:=[s]-omega*[t]
        CALL dcopy(bis-von+1,loctmp(von,s),1,loctmp(von,r),1)
        CALL daxpy(bis-von+1,-omega,loctmp(von,t),1,loctmp(von,r),1)

!      [p]:=[r]+beta*([p]-omega*[v])
        CALL daxpy(bis-von+1,-omega,loctmp(von,v),1,loctmp(von,p),1)
        CALL dscal(bis-von+1,beta,loctmp(von,p),1)
        CALL daxpy(bis-von+1,d_one,loctmp(von,r),1,loctmp(von,p),1)

!     end blocking
        von = bis + 1
        bis = min(n,bis+int(bl_size/bldiv_bicg(2)))
        IF (von<=n) GO TO 2


!      rho:=sigma
        p_rho = sigma

!      iteration counting (+2 Matrix-Vector-Products)
        p_iter = p_iter + 1

!      if division by zero detected, try to testing for good precision
!      can be disabled !!!
        IF (p_divide_zero/=0) p_need_ax = .TRUE.

!      solve [K]x[y]=[p]:
!       L_PrCo->[t_pc]   [t_pc]:=[L_K^-1]x[p]
!       R_PrCo->[y]      [y]:=[R_K^-1]x[t_pc]
!      myMVP->[v]       :[v]:=[A]x[y]

!$OMP  master
        work = do_y_p_v

!      need a second MVP, compute ([A]x[x]) in [z]
        IF (p_need_ax) work = more_y_p_v
!$OMP  end master
        p_step = 200
!      no barrier here ... need barrier later ...

        GO TO 1000


400     CONTINUE

!      -> END

!     restore best [x] value
        IF (rpar(rb_e_old)<p_e_old .AND. rpar(rb_e_old)<de) THEN
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
          IF (p_divide_zero==1) THEN
            WRITE(*,*) ' '
            WRITE(*,'(A)') ' warning: ([r^]^T*[r]) near Null'
          ELSE IF (p_divide_zero==7) THEN
            IF (linfos(4)>=0) THEN
              WRITE(*,*) ' '
              WRITE(*,'(2A)') ' warning: trivial solution, no boundary'&
                , ' condition different from zero defined'
            END IF
          ELSE
            IF (linfos(4)>=0) THEN
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
              ':MAXIMUM (BiCGStab)'
          ELSE
            WRITE(*,'(A,I5,A)') ', Number of Iterations: ', p_iter, &
              ' (BiCGStab)'
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
        rpar(rrho) = p_rho
        rpar(ralpha) = p_alpha
        lpar(lneed_ax) = p_need_ax
        ipar(istep) = p_step
        ipar(idivide_zero) = p_divide_zero
        ipar(ie_count) = p_e_count
        rpar(re_old) = p_e_old
        ipar(ifirst_ax) = p_first_ax
!$OMP  end master

        RETURN
      END
