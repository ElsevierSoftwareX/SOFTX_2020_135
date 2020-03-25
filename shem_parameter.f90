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

!>    @brief function library version to call Shemat-Suite as a library, needed for EFCOSS Parameter Estimation ( e.g. for SpaceMapping)
!>    @details
!>an inverse aquifer simulation tool for coupled flow, heat transfer, and transport\n
!>
!>based on numerical modules developed by: \n
!> -  Jörn Bartels     (gtn neubrandenburg gmbh)\n
!> -  Michael Kühn     (tu hamburg-harburg)\n
!> -  Roland Wagner    (ag,rwth aachen)\n
!> -  Martin H. Bücker   (sc,rwth aachen)\n
!> -  Christof Clauser (ag,rwth aachen)\n
      SUBROUTINE shemat(pp0,pn0,pparm0,pddata,infilename,fdim,omp_outer,omp_inner)
        use arrays
        use g_arrays
        use mod_genrl
        use mod_genrlc
        use mod_flow
        use mod_temp
        use mod_time
        use mod_data
        use mod_linfos
        use mod_OMP_TOOLS
        use mod_inverse
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k
        double precision :: tsglobal, tfglobal
        double precision :: tflocal
        INCLUDE 'version.inc'
        INCLUDE 'OMP_TOOLS.inc'

        CHARACTER*80 fname,filename,fname_bin, fname_sha,sha1sum_new,sha1sum_old
        character*80, intent(in), optional :: infilename
        DOUBLE PRECISION qval_old, rms_para_old, qval_new
        INTEGER ll, lli, llp, llh, llk,filemode,fstat
        INTEGER lblank, max_it, criteria
        integer, optional :: fdim
        integer startc,endc
        EXTERNAL lblank,chln
        DOUBLE PRECISION trun, tend

        LOGICAL lshem, g_lshem, hu_g_lshem,ad_matvec,g_matvec
        LOGICAL write_bin,use_bin
        PARAMETER ( write_bin=.FALSE.)

        INTEGER h_p_, g_p_
        INTEGER h_pmax_, g_pmax_
        COMMON /mpmax/h_pmax_, g_pmax_
!     OMP Stuff
        INTEGER, INTENT(in),optional :: omp_inner,omp_outer
!     local and direct parameter for each subroutine interface
        INTEGER, PARAMETER :: ncoord = 2
!     pP0:#parameters, pN0:#borehole-cells
        INTEGER pn0, pk0
        DOUBLE PRECISION ptmax, hu_ptmax(h_p_)
        INTEGER pp0
        DOUBLE PRECISION pparm0(pp0), g_pparm0(g_p_,pp0)
!     main values
        DOUBLE PRECISION pddata(pn0), g_pddata(g_p_,pn0), &
          hu_g_pddata(h_p_,g_p_,pn0)
       double precision:: xvec(pp0),yvec(pn0)
        double precision:: xvec_ad(pn0),yvec_ad(pp0)
!AW-      double precision pdhead(pK0,pN0),g_pdhead(g_p_,pK0,pN0),
!AW-     &  hu_g_pdhead(h_p_,g_p_,pK0,pN0)
!AW-      double precision pdtemp(pK0,pN0),g_pdtemp(g_p_,pK0,pN0),
!AW-     &  hu_g_pdtemp(h_p_,g_p_,pK0,pN0)
!AW-      double precision pdconc(pK0,pN0),g_pdconc(g_p_,pK0,pN0),
!AW-     &  hu_g_pdconc(h_p_,g_p_,pK0,pN0)
!       seeding sorting
        INTEGER, ALLOCATABLE :: lin_seed(:)
! ############################################
!     default: shemat, only forward evaluation
        lshem = .TRUE.
        g_lshem = .FALSE.
        hu_g_lshem = .FALSE.
        g_matvec = .FALSE.
        ad_matvec = .FALSE.
        mpara = 0
        g_pmax_ = 0
        h_pmax_ = 0
        GO TO 200

#ifdef matvec
        ENTRY g_shemat_matvec(pp0,pn0,pparm0,xvec,pddata,yvec,infilename,fdim,omp_inner,omp_outer)
        lshem = .TRUE.
         g_lshem = .FALSE.
        hu_g_lshem = .FALSE.
        ad_matvec = .FALSE.
       
        g_matvec = .TRUE.
        g_pmax_=1
         go to 200

        ENTRY g_shemat_mattvec(pp0,pn0,pparm0,xvec_ad,pddata,yvec_ad,infilename,fdim,omp_inner,omp_outer)
        lshem = .TRUE.
        g_lshem = .FALSE.
        hu_g_lshem = .FALSE.
        g_matvec = .FALSE.
        ad_matvec = .TRUE.
        g_pmax_=1
        write(*,*) "Called mattvec interface with xvec=",yvec_ad," yvec=",xvec_ad
        go to 200

#endif



! ############################################
!     default: g_shemat, dependency from [main-parameters]
        ENTRY g_shemat_proc(g_p_,pp0,pn0,pparm0,g_pparm0, &
          pddata,g_pddata,infilename,fdim,omp_inner,omp_outer)
        lshem = .TRUE.
        g_lshem = .TRUE.
        g_matvec=.FALSE.
        ad_matvec= .FALSE.
        hu_g_lshem = .FALSE.
        g_pmax_ = g_p_
        h_pmax_ = 0
        write(*,*) 'g_p_=',g_pmax_, infilename(1:fdim)
        GO TO 200
! ############################################
!     default: hu_g_shemat, dependency from [main-parameters] and [x-coords]
        ENTRY hu_shemat_proc(h_p_,g_p_,pp0,pn0, &
          pparm0,g_pparm0,pddata,g_pddata,infilename,fdim,omp_inner, omp_outer)
        lshem = .TRUE.
        g_lshem = .TRUE.
        hu_g_lshem = .TRUE.
        g_pmax_ = g_p_
        h_pmax_ = h_p_
        IF (h_p_/=1+ncoord) THEN
          WRITE(*,'(1A)') 'error: wrong number of parameters (h_p_)!'
          STOP
        END IF
        GO TO 200

! ############################################


!     begin main routine
200     CONTINUE
        CALL sys_cputime(tsglobal)
        write_disable = .FALSE.
        write_iter_disable = .FALSE.
        nested_build = .TRUE.
        ismpl = 1
        def_binary = 'inverse'
        lib_override = .TRUE.
        filemode=0
        if (.not. present(infilename) .or. infilename=='') then
            fname='shemade.job'
            OPEN(66,file=fname)
         else
            filemode=1
            filename(1:fdim)=infilename
            filename(fdim+1:80)=''
            write(*,*) "Shemat_Parameter_Estimation for EFCOSS ", filename(1:fdim),fdim
        end if

! -----------------------------------------

        write(*,*) "OMP config:",omp_inner,omp_outer
        if (omp_inner>0 .or. omp_outer>0) fl_omp_set=.true.
        if (present(omp_inner)) fl_omp_inner=omp_inner
        if (present(omp_outer)) fl_omp_outer=omp_outer
!     read new model
10      if (filemode .eq. 0) READ(66,'(A)',end=99999) filename
        CALL sys_cputime(tslocal)


        IF (filename(1:1)=='!') GO TO 10
        IF (filename(1:1)=='%') GO TO 10
        IF (filename(1:1)=='#') GO TO 10

        WRITE(*,*) ' '
        WRITE(*,*) ' '
        WRITE(*,*) ' *** NEW MODEL '

        project = filename
        runmode = 0
        transient = .FALSE.
        nperiod = 0
        is_init_flow_trafo_needed = .true.

!     no DD-mode
        update_dd = .FALSE.
        use_bin= .FALSE. 
!     read forward model
        CALL sys_cputime(trun)
        CALL read_model(filename,ismpl)
        CALL read_control(filename,ismpl)
        CALL sys_cputime(tend)
        write(*,*) ' [T] : "read_model" in ',tend-trun,' sec'
!     enable data reading
        runmode = max(runmode,1)
!     read timestep parameter
        CALL sys_cputime(trun)
        CALL read_time(filename,ismpl)
        CALL sys_cputime(tend)
        write(*,*) ' [T] : "read_time" in ',tend-trun,' sec'
!     read data
        CALL sys_cputime(trun)
        IF (runmode>0) CALL read_data(filename_data,ismpl)
        CALL sys_cputime(tend)
        write(*,*) ' [T] : "read_data" in ',tend-trun,' sec'
        CALL sys_cputime(trun)
        IF (runmode>1) CALL read_inverse(filename_inverse,ismpl)
        CALL sys_cputime(tend)
        write(*,*) ' [T] : "read_inverse" in ',tend-trun,' sec'
!
        write(*,*) ' <D> : ignore UNIT splitting when in LIBRARY mode'
!AW-not here -> reading INV-version !       split units
!AW-not here -> reading INV-version         CALL read_split_inv(filename,ismpl)
!
!      stop
!     setup log file name
        status_log = filename(1:lblank(filename)) // '_status.log'
        status_log_inv = filename(1:lblank(filename)) // &
          '_status-inv.log'
        restart_name = filename(1:lblank(filename)) // '_restart'
!
!     initialisation for different solver (optimisation methods)
        CALL init_step()
!
!       initialize
        iter_inv = 1
!        write(*,*) "temp_active=",temp_active
!
! ------
!       add output for new optimization observation
        WRITE(*,'(3A)') '  [W] : "', status_log(1:lblank(status_log)), '"'
        OPEN(76,file=status_log)
        WRITE(76,'(2A)') '% Shemat-Suite version: ', version
        WRITE(76,'(2A)') '%          build: ', datum
        WRITE(76,'(2A)') '%   build command line: ', makecmd
        WRITE(76,'(3A)') '%        Project: "', &
          filename(1:lblank(filename)), '"'
        WRITE(76,'(1A)') '%'
        IF (transient) THEN
          WRITE(76,'(2A)') '% <time step>', ' <cum time>'
        END IF
#ifdef head_base
        WRITE(76,'(4A)') '%    <iteration>',' <delta head>',' <delta temp>',' (<delta conc> ...)'
#endif
#ifdef pres_base
        WRITE(76,'(4A)') '%    <iteration>',' <delta pres>',' <delta temp>',' (<delta conc> ...)'
#endif
        CLOSE(76)
!
! ------
!       sanitity check
        IF (g_pmax_>mpara) THEN
          WRITE(*,'(1A)') 'error: wrong number of parameters (g_p_)!'
          write(*,*) 'g_p_=',g_pmax_
          write(*,*)' mpara=',mpara
          STOP
        END IF
        IF (pp0/=mpara) THEN
          WRITE(*,'(1A)') 'error: wrong number of parameters (P0)!'
          write(*,*)' pp0= ', pp0
          write(*,*)' mpara=',mpara
          STOP
        END IF
!        IF (pk0/=1) THEN
!          WRITE(*,'(1A)') 'error: wrong (>1) output vector length!'
!          STOP
!        END IF
!
!       overwrite parameters
        WRITE(*,'(1A)') '  [I] : modify parameter(s)'
!       init offset
        k = 0
!       for pparm0
        DO ll = 1, pp0
          CALL set_optip(k +ll,pparm0(ll),ismpl)
        END DO
        k = k +pp0
! ------
!
!       forward modeling (initialization)
        CALL forward_init(ismpl)
        write(*,*) "Forward init done"
!
        CALL sys_cputime(tfglobal)
!       compute linear seeding factors
        if (.not. ad_matvec .and. .not. g_matvec) then
           ALLOCATE(lin_seed(g_pmax_))
           DO i = 1, g_pmax_
             lin_seed(i) = 0
           END DO
!          init offset
           k = 0
!          for pparm0
           DO i = 1, g_pmax_
             DO j = 1, pp0
               IF (g_pparm0(i,j)>0.5d0) lin_seed(i) = k +j
             END DO
           END DO
           k = k +pp0
        end if
!          CALL forward_iter(simtime_0,max_simtime,iter_inv,0,ismpl)

! ------ original function (forward model)
        IF (lshem .AND. .NOT. g_lshem .AND. .NOT. hu_g_lshem .AND. .NOT. g_matvec .AND. .NOT. ad_matvec) THEN
!         forward modeling (nonlinear iteration)
          write(*,*) "Launching Forward Iter"
          CALL forward_iter(simtime_0,max_simtime,iter_inv,0,ismpl)
          GO TO 300
        END IF
! ------

! ------ first derivatives (jacobi)
        IF (lshem .AND. g_lshem .AND. .NOT. hu_g_lshem .and. .not. g_matvec .and. .not. ad_matvec) THEN
          write(*,*) "Launching Jacobian Computation"
          CALL prepare_jacobian(ismpl)
!          ***********************************************
          CALL jacobi_compute(ismpl)
!          ***********************************************
          GO TO 299
        END IF
! ------
#ifdef matvec
        if (g_matvec) then       
           call single_init(ismpl)
          CALL prepare_jacobian(ismpl)
           CALL jac_mv(xvec,yvec,ismpl)
        go to 299
        end if

        write(*,*) d(:,1,1,1)

        if (ad_matvec) then
           call single_init(ismpl)
          CALL prepare_jacobian(ismpl)
          call jact_mv(xvec_ad,yvec_ad,ismpl)
          write(*,*) "jactvec=",yvec_ad
        go to 299
        end if
#endif 
! ------ second derivatives (hesse)
        IF (lshem .AND. g_lshem .AND. hu_g_lshem) THEN
          WRITE(*,*) &
            'error: second order derivitive not supported yet!'
          STOP
! -> loop over `?`
!          hu-seeding
!?         do ll = 1, pN0
!?            do llp = 1+(ll-1)*21, ll*21
!?               g_ddata(llp,cdd_time) = hu_pTmax(?)/dble(pN0)*dble(ll)
!?            enddo
!?         enddo
!?         call prepare_jacobian(ismpl)
!?         if (.not.transient) then
!            save state (before), to avoid side effects !
!            steady-state: use state from the last computation (improvement)
!?           call old_save(cgen_opti,ismpl)
!?         else
!            transient   : use state from file init (necessary)
!?           call old_restore(cgen_opti,ismpl)
!?        endif
!          ***********************************************
!?#ifndef JACOBI_FREE
!           full Jacobi matrix computation - one times
!?            CALL jacobi_compute(ismpl)
!?#endif
!          ***********************************************


          GO TO 298
        END IF


! ------

! ------
298     CONTINUE
!?          DO ll = 1, pn0
!?              llp = (ll-1)*3 + 1
!          - head -
!? [dhilf] has included the jacobian values as a 7-point-star
!?            call h_g_one_hole(I0,J0,K0,i,j,k,pxcoord,hu_pxcoord,
!?     &        delxa,delya,delza,
!?     &        phi(1,1,1,ismpl),dhilf,
!?     &        pdhead(ll),g_pdhead(1,ll),hu_g_pdhead(1,1,ll))
!          - temp -
!?            call h_g_one_hole(I0,J0,K0,i,j,k,pxcoord,hu_pxcoord,
!?     &        delxa,delya,delza,
!?     &        temp(1,1,1,ismpl),dhilf,
!?     &        pdtemp(ll),g_pdtemp(1,ll),hu_g_pdtemp(1,1,ll))
!          - conc -
!?            call h_g_one_hole(I0,J0,K0,i,j,k,pxcoord,hu_pxcoord,
!?     &        delxa,delya,delza,
!?     &        conc(1,1,1,1,ismpl),dhilf,
!?     &        pdconc(ll),g_pdconc(1,ll),hu_g_pdconc(1,1,ll))
!?          END DO
!        WRITE(*,*) ' '
!        WRITE(*,'(1A,3e17.8)') '  [I] : ddHEAD[1,1]/dpdt value= ', &
!          hu_g_pdhead(1:3,1,1)
!        WRITE(*,'(1A,3e17.8)') '  [I] : ddTEMP[1,1]/dpdt value= ', &
!          hu_g_pdtemp(1:3,1,1)
!        WRITE(*,'(1A,3e17.8)') '  [I] : ddCONC[1,1]/dpdt value= ', &
!          hu_g_pdconc(1:3,1,1)
!
299     CONTINUE
        WRITE(*,*) ' '
        WRITE(*,'(1A,5I4)') ' [I] : EFCOSS integer parameter: ',g_p_,pp0,pn0
        WRITE(*,'(1A,999I4)') ' [I] : seeding index vector: ',lin_seed
!#ifndef JACOBI_FREE
!        write(*,*) "JAC(",ndata,",",mpara,")="
!        do i=1,ndata
!          write(*,*) ( jac(i,j),j=1,mpara )
!        end do
        if ( .not. ad_matvec .and. .not. g_matvec) then
        DO ll = 1,pn0
!          write(*,*) "llp=",llp
!         - head -
          DO llh = 1, g_p_
            IF (lin_seed(llh)>0) THEN
              g_pddata(llh,ll) = jac(ll,lin_seed(llh))
            ELSE
              g_pddata(llh,ll) = 0.0d0
            END IF
          END DO
        END DO
        end if
        WRITE(*,*) ' '
!#endif
!
300     CONTINUE
        DO ll = 1, pn0
          pddata(ll) = sdata(ll,ismpl)
        END DO
! ------

!     --- free memory
        CALL dealloc_arrays(ismpl)
        CALL props_end(ismpl)
        CALL dealloc_data(ismpl)

        IF (runmode>1) THEN
           if (.not. ad_matvec .and. .not. g_matvec) then
          DEALLOCATE(lin_seed)
       end if
          CALL dealloc_inverse(ismpl)
!          if (.not. ad_matvec) then
          CALL g_dealloc_arrays(ismpl)
!          else
#ifdef matvec
          CALL dealloc_arrays_ad(ismpl)
#endif
!          end if
!AW-noexist-   call g_props_end
        END IF
!     ---
        WRITE(*,'(3A)') '  [W] : "', status_log(1:lblank(status_log)) &
          , '"'
        OPEN(76,file=status_log,status='unknown',position='append')
        CALL sys_cputime(tflocal)
        WRITE(76,'(A,F11.2,A)') '% job finished, total cpu time: ', &
          (tflocal-tslocal)/60.0D0, ' min'
        CLOSE(76)
!     Another file to load?
        if (filemode .eq. 0) then
          GO TO 10
        end if
!
! -----------------------------------------
!
!     finis terrae
99999   CONTINUE
!
        if (filemode .eq. 0) then
        CLOSE(66)
        end if
        CALL sys_cputime(tfglobal)
        tslocal = tfglobal - tsglobal
        i = tslocal/60.D0
        WRITE(*,'(1A,1I4,1A,1F5.2,1A)') ' total cpu time: ', i, ':', &
          tslocal - dble(i)*60.D0, ' min'
        WRITE(*,*) 'RUN O.K.'
!
        RETURN
      END
