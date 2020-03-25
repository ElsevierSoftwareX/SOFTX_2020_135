!>    @brief derivative version to solve a linear equation system
!>    @param[in] icode specify the physical value
!>    @param[in] species transport species index
!>    @param[out] x_solution the solution vector of the solved linear system
!>    @param[in] errorc error
!>    @param[in] apar setup explicit or implicit solver (disabled), default 1.0
!>    @param[in] ctrl solver code
!>    @param[in] ismpl local sample index
!>    @details
!>  solve system equations by preconditioned  krylov solvers (bicgstab) or SIP (NAG) and others\n
!>  According to Tapenade FAQ:
!>  subroutine SOLVE_AD(A,Ab,y,yb,b,bb,n)
!>  INTEGER n
!>  REAL A(n,n), Ab(n,n)
!>  REAL y(n), yb(n), b(n), bb(n)
!>  INTEGER i,j
!>  REAL AT(n,n), incrbb(n)
!>! Compute A^T
!>  DO i=1,n
!>     DO j=1,n
!>        AT(i,j) =A(j,i)
!>     ENDDO
!>  ENDDO
!>! Solve A^T * incrbb=yb
!>  call SOLVE(AT,incrbb,yb,n)
!>  DO j=1,n
!>     bb(j)=bb(j)+incrbb(j)
!>  ENDDO
!>! Solve A*y=b
!>  call SOLVE(A,y,b,n) 
!>! Ab=Ab-y*incrbb
!>  DO i=1,n
!>     DO j=1,n
!>        Ab(i,j)=Ab(i,j)-y(j)*incrbb(i)
!>     ENDDO
!>  ENDDO
!>  END

      SUBROUTINE solve_ad(icode,species,x_solution,x_solution_ad,errorc, &
          apar,ctrl,ismpl)
        use arrays
        use arrays_ad
        use mod_conc
        use mod_flow
        use mod_genrl
        use mod_linfos
        use mod_temp
        use mod_OMP_TOOLS
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k
        INCLUDE 'OMP_TOOLS.inc'
        DOUBLE PRECISION apar, errorc
        INTEGER icode, itst, mxit, species
        INTEGER ctrl, criteria, precondition, solvername
        DOUBLE PRECISION x_solution(i0*j0*k0)
!     only for tests with 'abbruch'
        DOUBLE PRECISION moderrorc, enough
        DOUBLE PRECISION trun, tend
        INTRINSIC dsqrt, dble, abs
!     -------------
        DOUBLE PRECISION d1
        DOUBLE PRECISION x_solution_ad(i0*j0*k0)
        INTEGER pm, am

!
!
!



!KFP helper arrays for transpose matrix 
!        DOUBLE PRECISION, ALLOCATABLE :: at_h(:,:,:,:)
!        DOUBLE PRECISION, ALLOCATABLE :: bt_h(:,:,:,:)
!        DOUBLE PRECISION, ALLOCATABLE :: ct_h(:,:,:,:)
!        DOUBLE PRECISION, ALLOCATABLE :: dt_h(:,:,:,:)
!        DOUBLE PRECISION, ALLOCATABLE :: et_h(:,:,:,:)
!        DOUBLE PRECISION, ALLOCATABLE :: ft_h(:,:,:,:)
!        DOUBLE PRECISION, ALLOCATABLE :: gt_h(:,:,:,:)
!        DOUBLE PRECISION, ALLOCATABLE :: w_ad_h(:,:,:,:)
!
!        ALLOCATE(at_h(I0,J0,K0,nsmpl))
!        memory = memory + i0*j0*k0*nsmpl
!        ALLOCATE(bt_h(I0,J0,K0,nsmpl))
!        memory = memory + i0*j0*k0*nsmpl
!        ALLOCATE(ct_h(I0,J0,K0,nsmpl))
!        memory = memory + i0*j0*k0*nsmpl
!        ALLOCATE(dt_h(I0,J0,K0,nsmpl))
!        memory = memory + i0*j0*k0*nsmpl
!        ALLOCATE(et_h(I0,J0,K0,nsmpl))
!        memory = memory + i0*j0*k0*nsmpl
!        ALLOCATE(ft_h(I0,J0,K0,nsmpl))
!        memory = memory + i0*j0*k0*nsmpl
!        ALLOCATE(gt_h(I0,J0,K0,nsmpl))
!        memory = memory + i0*j0*k0*nsmpl
!        ALLOCATE(w_ad_h(I0,J0,K0,nsmpl))
!        memory = memory + i0*j0*k0*nsmpl
!end KFP


        moderrorc = errorc
        IF (ctrl<0) THEN
          WRITE(*,*) 'error: old solver compatibility no &
            &longer        supported'
          STOP
        END IF
! **************************************************************************
! solvername =
!           *0 : BiCGStab [:parallel]
!            1 : NAG [:serial]
!            2 : CG [:parallel] ( prove symmetry, if not then BiCGStab )
!            3 : PLU [:serial] (LAPACK) and math tests (stability)
!         (4-7 : not in use !)
! criteria = < switch to set when should break >
!          0 : relative stopping crit. : ||[res]|| < depsilon*||[res0]||
!          1 : absolute stopping crit. : ||[res]|| < depsilon
!          2 : maximum  stopping crit. : max(abs([res])) < depsilon
!          3 : abs. and rel. stopping crit. : ( ||[res]|| < depsilon ) and ( ||[
!res]|| < 0.99d0*||[res0]|| )
!              0.99d0 is a constant for testing only, it is named 'minRel' in 'a
!bbruch.f'
!         *4 : like '3', but with auto detected range for depsilon, default deps
!ilon used when detection fails
!       (5-7 : not in use !)
! precondition =
!             *0 : ILU [:parallel]
!              1 : SSOR [:serial]
!              2 : Diagonal [:parallel]
!              3 : None [:]
!           (4-7 : not in use !)
! * : recommended (ctrl = 64)
! and (ctrl = 67) : for testing
! [ ctrl = solvername + 16*criteria + 256*precondition ]
!     extract [ ctrl ] :
!      solvername = mod(ctrl,16)
!      ctrl = ctrl/16
!      criteria = mod(ctrl,16)
!      ctrl = ctrl/16
!      precondition = ctrl
        CALL decntrl3(ctrl,solvername,criteria,precondition)
! **************************************************************************
!     -------------
        IF (abs(icode)==pv_pres) THEN
          mxit = lmaxitf
          IF (criteria>=4) THEN
            criteria = 1
            moderrorc = nltolf*1.0D-2
          END IF
        END IF
        IF (abs(icode)==pv_conc) THEN
          mxit = lmaxitc
          IF (criteria>=4) THEN
            criteria = 1
            moderrorc = nltolc*1.0D-2
          END IF
        END IF
        IF (abs(icode)==pv_temp) THEN
          mxit = lmaxitt
          IF (criteria>=4) THEN
            criteria = 1
            moderrorc = nltolt*1.0D-2
          END IF
        END IF
        IF (abs(icode)==pv_head) THEN
          mxit = lmaxitf
          IF (criteria>=4) THEN
            criteria = 1
            moderrorc = nltolf*1.0D-2
          END IF
        END IF
!     -------------

!AW - original function call - if needed !!!
!AW        CALL solve_type(i0,j0,k0,x_solution,w(1,1,1,ismpl),moderrorc, &
!AW          bc_mask(1,ismpl),solvername,precondition,mxit,criteria, &
!AW          a(1,1,1,ismpl),b(1,1,1,ismpl),c(1,1,1,ismpl),d(1,1,1,ismpl), &
!AW          e(1,1,1,ismpl),f(1,1,1,ismpl),g(1,1,1,ismpl),r,apar, &
!AW          ud(1,1,1,ismpl),ismpl)

!KFP
!Ax=b: [x_ad]=[A_transpose] [b_ad_update] => [b_ad_update],
!      [A_ad_update]=-[b_ad_update] [x_transpose]
![b_ad]=[b_ad_old]+[b_ad_update], store [b_ad_update] in lss_tmp 
![A_ad]=[A_ad_old]+[A_ad_update]
!
!
!Note: only those entries A_ad(i,j) are needed for which A(i,j)!=0
!      i.e., where A(i,j) is specified by a,b,c,d,e,f,g
!Reason: product A_ad(i,j)*dA(i,j) in all adjoint updates using A_ad 
!        if A(i,j)=0, dA(i,j)=0 => A_ad(i,j) does not enter final result 
!
!1) transpose [A]: a,b,c,d,e,f,g => at_h, bt_h etc.
!2) solve [x_ad]=[A_transpose] [b_ad] (solve_type) => [b_ad]
!3) calculate relevant [A_ad] entries from [b_ad] and [x] => a_ad, b_ad, ...
!
!end KFP

! Populate helper arrays for transpose
        am=i0*j0*k0
!        do i=1,i0
!           do j=1,j0
!              do k=1,k0
!                 pm=i+(j-1)*I0+(k-1)*I0*J0
!                 if (pm>=i0*j0+1) at_h(i,j,k,ismpl)=g(i,j,k-1,ismpl)
!                 if (pm>=i0+1) bt_h(i,j,k,ismpl)=f(i,j-1,k,ismpl)
!                 if (pm>=2) ct_h(i,j,k,ismpl)=e(i-1,j,k,ismpl)
!                 dt_h(i,j,k,ismpl)=d(i,j,k,ismpl)
!                 if (pm<=am-1) et_h(i,j,k,ismpl)=c(i+1,j,k,ismpl)
!                 if (pm<=am-i0) ft_h(i,j,k,ismpl)=b(i,j+1,k,ismpl)
!                 if (pm<=am-i0*j0) gt_h(i,j,k,ismpl)=a(i,j,k+1,ismpl)
!              end do
!           end do
!        end do

!KFP: wrong at borders i=i,j=j0->i=i+1,j=1, etc.
!        if (k0.ge.2) at_h(:,:,2:k0,:)=g(:,:,1:k0-1,:)
!        if (j0.ge.2) bt_h(:,2:j0,:,:)=f(:,1:j0-1,:,:)
!        if (i0.ge.2) ct_h(2:i0,:,:,:)=e(1:i0-1,:,:,:)
!        dt_h(:,:,:,:)=d(:,:,:,:)
!        if (i0.ge.2) et_h(1:i0-1,:,:,:)=c(2:i0,:,:,:)
!        if (j0.ge.2) ft_h(:,1:j0-1,:,:)=b(:,2:j0,:,:)
!        if (k0.ge.2) gt_h(:,:,1:k0-1,:)=a(:,:,2:k0,:)

!AW: faster/safer, but remember to swap back
        if (K0>=2) call dswap(I0*J0*K0-I0*J0,a(1,1,2,ismpl),1,g(1,1,1,ismpl),1)
        if (J0>=2) call dswap(I0*J0*K0-I0,   b(1,2,1,ismpl),1,f(1,1,1,ismpl),1)
        if (I0>=2) call dswap(I0*J0*K0-1,    c(2,1,1,ismpl),1,e(1,1,1,ismpl),1)

! Solve for [b_ad]=w_ad
        CALL solve_type(i0,j0,k0,lss_tmp(1,ismpl),x_solution_ad, &
          moderrorc,bc_mask(1,ismpl),solvername,precondition,mxit, &
          criteria,a(1,1,1,ismpl),b(1,1,1,ismpl),c(1,1,1,ismpl), &
          d(1,1,1,ismpl),e(1,1,1,ismpl),f(1,1,1,ismpl), &
          g(1,1,1,ismpl),r,apar,ud(1,1,1,ismpl),ismpl)
!        CALL solve_type(i0,j0,k0,lss_tmp(1,ismpl),x_solution_ad, &
!          moderrorc,bc_mask(1,ismpl),solvername,precondition,mxit, &
!          criteria,at_h(1,1,1,ismpl),bt_h(1,1,1,ismpl),ct_h(1,1,1,ismpl), &
!          dt_h(1,1,1,ismpl),et_h(1,1,1,ismpl),ft_h(1,1,1,ismpl), &
!          gt_h(1,1,1,ismpl),r,apar,ud(1,1,1,ismpl),ismpl)

! [A_ad((i1,j1,k1),(i2,j2,j2))]=-[b_ad(i1,j1,k1)][x(i2,j2,k2)]
! Select a_ad=(A_ad(i,j,k),(i,j,k-1)), b_ad=(A_ad(i,j,k),(i,j-1,k)),...
        do i=1,i0
           do j=1,j0
              do k=1,k0
                 pm=i+(j-1)*I0+(k-1)*I0*J0
                 w_ad(i,j,k,ismpl)=w_ad(i,j,k,ismpl)+lss_tmp(pm,ismpl)
                 if (pm>=i0*j0+1) a_ad(i,j,k,ismpl)=a_ad(i,j,k,ismpl) &
                      -lss_tmp(pm,ismpl)*x_solution(pm-i0*j0)
                 if (pm>=i0+1) b_ad(i,j,k,ismpl)=b_ad(i,j,k,ismpl) &
                      -lss_tmp(pm,ismpl)*x_solution(pm-i0)
                 if (pm>=2) c_ad(i,j,k,ismpl)=c_ad(i,j,k,ismpl) &
                      -lss_tmp(pm,ismpl)*x_solution(pm-1)
                 d_ad(i,j,k,ismpl)=d_ad(i,j,k,ismpl) &
                      -lss_tmp(pm,ismpl)*x_solution(pm)
                 if (pm<=am-1) e_ad(i,j,k,ismpl)=e_ad(i,j,k,ismpl) &
                      -lss_tmp(pm,ismpl)*x_solution(pm+1)
                 if (pm<=am-i0) f_ad(i,j,k,ismpl)=f_ad(i,j,k,ismpl) &
                      -lss_tmp(pm,ismpl)*x_solution(pm+i0)
                 if (pm<=am-i0*j0) g_ad(i,j,k,ismpl)=g_ad(i,j,k,ismpl) &
                      -lss_tmp(pm,ismpl)*x_solution(pm+i0*j0)
              end do
           end do
        end do

! swap back
        if (K0>=2) call dswap(I0*J0*K0-I0*J0,a(1,1,2,ismpl),1,g(1,1,1,ismpl),1)
        if (J0>=2) call dswap(I0*J0*K0-I0,   b(1,2,1,ismpl),1,f(1,1,1,ismpl),1)
        if (I0>=2) call dswap(I0*J0*K0-1,    c(2,1,1,ismpl),1,e(1,1,1,ismpl),1)

        call set_dval(I0*J0*K0,0.0d0,x_solution_ad)

!        DEALLOCATE(at_h)
!        memory = memory - i0*j0*k0*nsmpl
!        DEALLOCATE(bt_h)
!        memory = memory - i0*j0*k0*nsmpl
!        DEALLOCATE(ct_h)
!        memory = memory - i0*j0*k0*nsmpl
!        DEALLOCATE(dt_h)
!        memory = memory - i0*j0*k0*nsmpl
!        DEALLOCATE(et_h)
!        memory = memory - i0*j0*k0*nsmpl
!        DEALLOCATE(ft_h)
!        memory = memory - i0*j0*k0*nsmpl
!        DEALLOCATE(gt_h)
!        memory = memory - i0*j0*k0*nsmpl
        RETURN
      END
