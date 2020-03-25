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

!>    @brief main routine to solve a linear equation system, calls "solve_type" (wrapper)
!>    @param[in] icode specify the physical value
!>    @param[in] species transport species index
!>    @param[in,out] x_solution solution vector of the solved linear system, may be also the start input
!>    @param[in] errorc error
!>    @param[in] apar setup explicit or implicit solver (disabled), default 1.0
!>    @param[in] ctrl solver code
!>    @param[in] ismpl local sample index
!>    @details
!>  solve system equations by preconditioned  krylov solvers (bicgstab) or SIP (NAG) and others\n
      SUBROUTINE solve(icode,species,x_solution,errorc,apar,ctrl,ismpl)
        use arrays
        use mod_genrl
        use mod_flow
        use mod_temp
        use mod_conc
        use mod_linfos
        IMPLICIT NONE
        integer :: ismpl
        double precision :: apar, errorc
        ! double precision :: dnormalise
        integer :: icode
        ! integer :: itst
        integer :: mxit
        integer :: species
        integer :: mfactor
        integer :: ctrl, criteria, precondition, solvername
        double precision :: x_solution(i0*j0*k0)
!     only for tests with 'abbruch'
        double precision :: moderrorc
        ! double precision :: enough
#ifdef BENCH
        double precision :: trun
        double precision :: tend
#endif
        intrinsic dsqrt, dble, abs

!     -------------
        moderrorc = errorc
        IF (ctrl<0) THEN
          WRITE(*,*) & 
            'error: old solver compatibility no longer supported'
          STOP
        END IF

! **************************************************************************
! solvername =
!           *0 : BiCGStab [:parallel]
!            1 : not supported, formerly: NAG [:serial]
!            2 : CG [:parallel] ( prove symmetry, if not then BiCGStab )
!            3 : PLU [:serial] (LAPACK) and math tests (stability)
!         (4-7 : not in use !)
! criteria = < switch to set when should break >
!          0 : relative stopping crit. : ||[res]|| < depsilon*||[res0]||
!          1 : absolute stopping crit. : ||[res]|| < depsilon
!          2 : maximum  stopping crit. : max(abs([res])) < depsilon
!          3 : abs. and rel. stopping crit. : ( ||[res]|| < depsilon ) and ( ||[res]|| < 0.99d0*||[res0]|| )
!              0.99d0 is a constant for testing only, it is named 'minRel' in 'abbruch.f'
!         *4 : like '3', but with auto detected range for depsilon, default depsilon used when detection fails
!       (5-7 : not in use !)
! precondition =
!             *0 : ILU [:parallel]
!              1 : SSOR [:serial]
!              2 : Diagonal [:parallel]
!              3 : None [:]
!           (4-7 : not in use !)
! * : recommended (ctrl = 64)
! and (ctrl = 67) : for testing
!
! [ ctrl = solvername + 16*criteria + 256*precondition ]
!     extract [ ctrl ] :
!      solvername = mod(ctrl,16)
!      ctrl = ctrl/16
!      criteria = mod(ctrl,16)
!      ctrl = ctrl/16
!      precondition = ctrl
!
        CALL decntrl3(ctrl,solvername,criteria,precondition)

! **************************************************************************

!     -------------
!     defualt multiply factoor
        mfactor = 1
!
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

!     -------------
        IF (tec_out) THEN
!        save current state of the right side vector - may be modified after "solve_type"
          CALL dcopy(i0*j0*k0,w(1,1,1,ismpl),1,lss_tmp(1,ismpl),1)
        END IF

!    ##################################################################
!    #                             solver                             #
!    ##################################################################

#ifdef BENCH
        CALL sys_cputime(trun)
#endif
        CALL solve_type(i0,j0,k0*mfactor,x_solution,w(1,1,1,ismpl),moderrorc, & 
          bc_mask(1,ismpl),solvername,precondition,mxit,criteria, & 
          a(1,1,1,ismpl),b(1,1,1,ismpl),c(1,1,1,ismpl),d(1,1,1,ismpl), & 
          e(1,1,1,ismpl),f(1,1,1,ismpl),g(1,1,1,ismpl),r,apar, & 
          ud(1,1,1,ismpl),ismpl)
#ifdef BENCH
        CALL sys_cputime(tend)
        WRITE(*,*) ' linear system solver time:', tend - trun, 'sec'
#endif

!    ##################################################################

!     -------------
!     "icode" < 0 indicates a default input/output in the [x] vector
        IF (icode<0) THEN
          CALL dcopy(i0*j0*k0,x_solution,1,x(1,1,1,ismpl),1)
        END IF

!     -------------
        IF (tec_out) THEN
!        restore the state of the right side vector
          CALL dcopy(i0*j0*k0,lss_tmp(1,ismpl),1,w(1,1,1,ismpl),1)
        END IF


        RETURN
      END
