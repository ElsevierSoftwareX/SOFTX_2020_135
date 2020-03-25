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

!>    @brief wrapper routine to call each kind of linear system solver
!>    @param[in] ni lengths of I dimension of local matrix [M]
!>    @param[in] nj lengths of J dimension of local matrix [M]
!>    @param[in] nk lengths of K dimension of local matrix [M]
!>    @param[in,out] l_x solution vector [x], on start = start vector
!>    @param[in] l_w right side, vector [b]
!>    @param[in] moderrc_ only for tests with 'abbruch'
!>    @param[in] l_bc_mask boundary condition pattern (mask)
!>    @param[in] solvername solver name, code for CG, BiCGStab...
!>    @param[in] precondition preconditioner code for None, Diagonal, SSOR or ILU(0)
!>    @param[in] MXIT max iteration number
!>    @param[in] criteria_ precision criteria mode to break iterations\n
!>    - 0 : relative stopping crit.: ||[res]||       < depsilon*||[res0]||\n
!>    - 1 : absolute stopping crit.: ||[res]||       < depsilon\n
!>    - 2 : maximum  stopping crit.: max(abs([res])) < depsilon\n
!>    first [res]^=[r], later (if precise enough): [res]^=([M]x[x]-[b])
!>    @param[in] L_A 1. diagonal of the system matrix [M]
!>    @param[in] L_B 2. diagonal of the system matrix [M]
!>    @param[in] L_C 3. diagonal of the system matrix [M]
!>    @param[in] L_D 4. diagonal of the system matrix [M]
!>    @param[in] L_E 5. diagonal of the system matrix [M]
!>    @param[in] L_F 6. diagonal of the system matrix [M]
!>    @param[in] L_G 7. diagonal of the system matrix [M]
!>    @param[in] l_r random vector [r0_hat] ^= [r0^]
!>    @param[in] APAR explicit - implicit solver switch (default 1.0)
!>    @param[in] L_UD helper diagonal elements for preconditioning
!>    @param[in] ismpl local sample index
!>    @details
!>    Function: wrapper for any solver, spec. by criteria\n
!>    Input/Output: see in 'solver.f90'\n
      SUBROUTINE solve_type(ni,nj,nk,l_x,l_w,moderrc_,l_bc_mask, & 
          solvername,precondition,mxit,criteria_,l_a,l_b,l_c,l_d,l_e, & 
          l_f,l_g,l_r,apar,l_ud,ismpl)
        use arrays
        use mod_linfos
        use mod_genrl
        IMPLICIT NONE
        INTEGER ni, nj, nk
        DOUBLE PRECISION moderrc_, l_a(ni,nj,nk), l_b(ni,nj,nk), & 
          l_c(ni,nj,nk), l_d(ni,nj,nk), l_e(ni,nj,nk), l_f(ni,nj,nk), & 
          l_g(ni,nj,nk), l_ud(ni,nj,nk)
        DOUBLE PRECISION apar, l_x(ni,nj,nk), l_w(ni,nj,nk), & 
          l_r(ni,nj,nk)
        CHARACTER l_bc_mask(ni,nj,nk)
        INTEGER criteria_, solvername, mxit, ii
        INTEGER precondition, ismpl
        LOGICAL test_symmetry, ldummy
        EXTERNAL test_symmetry


!     ##################################################################
!     #                             solver                             #
!     ##################################################################

!     ** nag solver not supported anymore **
        IF (solvername==1) THEN
           write(unit = *, fmt = *) '[E] error: nag solver not&
                & supported anymore.'
           stop
        END IF

! ---------------------------------------------------------

!     ** use Cg solver **
        IF (solvername==2) THEN
!       test about symmetry
          IF (test_symmetry(ismpl)) THEN
!          symmetriC solver ilu-preCo
            IF (precondition==0) THEN
              IF (linfos(4)>=2) WRITE(*,*) ' ilu precondition'
              CALL omp_sym_solve_ilu(ni,nj,nk,l_x,l_w,moderrc_, & 
                l_bc_mask,mxit,criteria_,l_a,l_b,l_c,l_d,l_e,l_f,l_g, & 
                l_ud,lss_bound_block(1,1,1,1,1,ismpl), & 
                lss_dnrm(1,ismpl),lss_lma(1,1,ismpl), & 
                lss_lmb(1,1,ismpl),lss_lmc(1,1,ismpl), & 
                lss_lmd(1,1,ismpl),lss_lme(1,1,ismpl), & 
                lss_lmf(1,1,ismpl),lss_lmg(1,1,ismpl), & 
                lss_lud(1,1,ismpl),lss_lx(1,1,ismpl), & 
                lss_lb(1,1,ismpl),lss_ldnrm(1,1,ismpl), & 
                lss_lloctmp(1,1,1,ismpl),lss_ud_block(1,1,1,ismpl), & 
                ismpl)
              RETURN
            END IF

            IF (precondition==1) THEN
              IF (linfos(4)>=2) WRITE(*,*) ' ssor precondition'
              CALL omp_sym_solve_ssor(ni,nj,nk,l_x,l_w,moderrc_, & 
                l_bc_mask,mxit,criteria_,l_a,l_b,l_c,l_d,l_e,l_f,l_g, & 
                lss_lloctmp(1,1,1,ismpl),lss_dnrm(1,ismpl),ismpl)
              RETURN
            END IF

            IF (precondition==2) THEN
              IF (linfos(4)>=2) WRITE(*,*) ' diagonal precondition'
              CALL omp_sym_solve_diag(ni,nj,nk,l_x,l_w,moderrc_, & 
                l_bc_mask,mxit,criteria_,l_a,l_b,l_c,l_d,l_e,l_f,l_g, & 
                lss_lloctmp(1,1,1,ismpl),lss_dnrm(1,ismpl),ismpl)
              RETURN
            END IF

            IF (precondition==3) THEN
              IF (linfos(4)>=2) WRITE(*,*) ' no precondition'
              CALL omp_sym_solve(ni,nj,nk,l_x,l_w,moderrc_,l_bc_mask, & 
                mxit,criteria_,l_a,l_b,l_c,l_d,l_e,l_f,l_g, & 
                lss_lloctmp(1,1,1,ismpl),lss_dnrm(1,ismpl),ismpl)
              RETURN
            END IF
          ELSE
            solvername = 0
          END IF
        END IF

! ---------------------------------------------------------

!     ** use biCgstab solver **
        IF (solvername==0) THEN
!        generic solver ilu-preco
          IF (precondition==0) THEN
            IF (linfos(4)>=2) WRITE(*,*) ' ilu precondition'
            CALL omp_gen_solve_ilu(ni,nj,nk,l_x,l_w,l_r,moderrc_, & 
              l_bc_mask,mxit,criteria_,l_a,l_b,l_c,l_d,l_e,l_f,l_g, & 
              l_ud,lss_bound_block(1,1,1,1,1,ismpl),lss_dnrm(1,ismpl), & 
              lss_lma(1,1,ismpl),lss_lmb(1,1,ismpl), & 
              lss_lmc(1,1,ismpl),lss_lmd(1,1,ismpl), & 
              lss_lme(1,1,ismpl),lss_lmf(1,1,ismpl), & 
              lss_lmg(1,1,ismpl),lss_lud(1,1,ismpl),lss_lx(1,1,ismpl), & 
              lss_lb(1,1,ismpl),lss_lr0_hat(1,1,ismpl), & 
              lss_ldnrm(1,1,ismpl),lss_lloctmp(1,1,1,ismpl), & 
              lss_ud_block(1,1,1,ismpl),ismpl)
            RETURN
          END IF

          IF (precondition==1) THEN
            IF (linfos(4)>=2) WRITE(*,*) ' ssor precondition'
            CALL omp_gen_solve_ssor(ni,nj,nk,l_x,l_w,l_r,moderrc_, & 
              l_bc_mask,mxit,criteria_,l_a,l_b,l_c,l_d,l_e,l_f,l_g, & 
              lss_lloctmp(1,1,1,ismpl),lss_dnrm(1,ismpl),ismpl)
            RETURN
          END IF

          IF (precondition==2) THEN
            IF (linfos(4)>=2) WRITE(*,*) ' diagonal precondition'
            CALL omp_gen_solve_diag(ni,nj,nk,l_x,l_w,l_r,moderrc_, & 
              l_bc_mask,mxit,criteria_,l_a,l_b,l_c,l_d,l_e,l_f,l_g, & 
              lss_lloctmp(1,1,1,ismpl),lss_dnrm(1,ismpl),ismpl)
            RETURN
          END IF

          IF (precondition==3) THEN
            IF (linfos(4)>=2) WRITE(*,*) ' no precondition'
            CALL omp_gen_solve(ni,nj,nk,l_x,l_w,l_r,moderrc_, & 
              l_bc_mask,mxit,criteria_,l_a,l_b,l_c,l_d,l_e,l_f,l_g, & 
              lss_lloctmp(1,1,1,ismpl),lss_dnrm(1,ismpl),ismpl)
            RETURN
          END IF
        END IF

! ---------------------------------------------------------

!     ** use PLU (LAPACK) and make math tests **
        IF (solvername==3) THEN
!       --------------------------------
!       test about symmetry and stability
          ldummy = test_symmetry(ismpl)

          CALL test_matrix(ismpl)
!       --------------------------------
!       memory space enough ?
          ii = ni*nj
          IF (ni==1) ii = ni*nj
          IF (nj==1) ii = ni*nj
          IF (nk==1) ii = ni

          IF (8.0D0*dble(3*ii+1)*dble(ni*nj*nk)<=2.0D9) THEN
            CALL direct_solve(ni,nj,nk,l_x,l_w,l_a,l_b,l_c,l_d,l_e, & 
              l_f,l_g)
          ELSE
            IF (linfos(4)>=2) WRITE(*,'(A,F12.2,A)') & 
              '   choose BICGStab (PLU need :', & 
              (8.0D0*dble(3*ii+1)*dble(ni*nj*nk)/(1024.0D0*1024.0D0)), & 
              ' MByte, which is more than 2GByte)'
!         BiCGStab with ILU(0) precondition
            CALL omp_gen_solve_ilu(ni,nj,nk,l_x,l_w,l_r,moderrc_, & 
              l_bc_mask,mxit,criteria_,l_a,l_b,l_c,l_d,l_e,l_f,l_g, & 
              l_ud,lss_bound_block(1,1,1,1,1,ismpl),lss_dnrm(1,ismpl), & 
              lss_lma(1,1,ismpl),lss_lmb(1,1,ismpl), & 
              lss_lmc(1,1,ismpl),lss_lmd(1,1,ismpl), & 
              lss_lme(1,1,ismpl),lss_lmf(1,1,ismpl), & 
              lss_lmg(1,1,ismpl),lss_lud(1,1,ismpl),lss_lx(1,1,ismpl), & 
              lss_lb(1,1,ismpl),lss_lr0_hat(1,1,ismpl), & 
              lss_ldnrm(1,1,ismpl),lss_lloctmp(1,1,1,ismpl), & 
              lss_ud_block(1,1,1,ismpl),ismpl)
          END IF
        END IF

! ---------------------------------------------------------

!       multi-phase newton bicgstb solver ilu-preco
        IF (solvername==4) THEN
           write(unit = *, fmt = *) '[E] error: multi-phase newton solver not'// &
                'supported anymore.'
           stop
        END IF


! ---------------------------------------------------------

        RETURN
      END
