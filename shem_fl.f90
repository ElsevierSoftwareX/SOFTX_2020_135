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

!>  @brief Find block in x direction
!>  @details Based on <xin> (coordinates) it finds the corresponding block in
!>   simulation, where <xin> is located at
      function getXBlock(xin)
         use arrays
         use g_arrays
         use mod_genrl
         use mod_genrlc
         implicit none
         double precision, intent(in) :: xin
         integer :: getXBlock,ll
         
         getXBlock=0
         do ll=1,i0
            if (delxa(ll)-0.5D0*delx(ll)<=xin)getXBlock=ll
         end do
         return
      end

!>  @brief Find block in y direction
!>  @details Based on <yin> (coordinates) it finds the corresponding block in
!>   simulation, where <yin> is located at
      function getYBlock(yin)
         use arrays
         use g_arrays
         use mod_genrl
         use mod_genrlc
         implicit none
         double precision, intent(in) :: yin
         integer :: getYBlock,ll
         
         getYBlock=0
         do ll=1,j0
            if (delya(ll)-0.5D0*dely(ll)<=yin)getYBlock=ll
         end do
         return
      end

!>  @brief Find block in x direction
!>  @details Based on <zin> (coordinates) it finds the corresponding block in
!>   simulation, where <zin> is located at
      function getZBlock(zin)
         use arrays
         use g_arrays
         use mod_genrl
         use mod_genrlc
         implicit none
         double precision, intent(in) :: zin
         integer :: getZBlock,ll
         getZBlock=0
         do ll=1,k0
            if (delza(ll)-0.5D0*delz(ll)<=zin)getZBlock=ll
         end do
         return
      end

!>    @brief function library version to call Shemat-Suite as a library, needed for OED- EFCOSS
!>    @details
!>an inverse aquifer simulation tool for coupled flow, heat transfer, and transport\n
!>
!>based on numerical modules developed by: \n
!> -  Jörn Bartels     (gtn neubrandenburg gmbh)\n
!> -  Michael Kühn     (tu hamburg-harburg)\n
!> -  Roland Wagner    (ag,rwth aachen)\n
!> -  Martin H. Bücker   (sc,rwth aachen)\n
!> -  Christof Clauser (ag,rwth aachen)\n
      SUBROUTINE shemat(pp0,pk0,pn0,ptmax,pparm0,pdhead,pdtemp,pdconc, &
          pxcoord,zmin,zmax,omp_inner,omp_outer)
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

        CHARACTER*80 filename
        DOUBLE PRECISION qval_old, rms_para_old, qval_new
        INTEGER ll, lli, llp, llh, llk
        INTEGER lblank, max_it, criteria
        integer zmin,zmax
        DOUBLE PRECISION qfunc
!     time stuff
!     if "restart" option is used
        LOGICAL test_option
        EXTERNAL qfunc, lblank, test_option
!     time of the last saving of the restart data
        DOUBLE PRECISION backuptime, depsilon, trun, tend

        LOGICAL lshem, g_lshem, hu_g_lshem,getGridJacobian,setGridJacobian

!     instead of constants, only possible in Fortran90
        INTEGER h_p_, g_p_
        INTEGER h_pmax_, g_pmax_
        COMMON /mpmax/h_pmax_, g_pmax_
!     OMP Stuff
        INTEGER, INTENT(in),optional :: omp_inner,omp_outer
!     local and direct parameter for each subroutine interface
        INTEGER, PARAMETER :: ncoord = 2
!     pP0:#parameter, pK0:#borehole-number, pN0:#borehole-cell
        INTEGER pn0, pk0
        DOUBLE PRECISION ptmax, hu_ptmax(h_p_)
        INTEGER pp0
        DOUBLE PRECISION pparm0(pp0), g_pparm0(g_p_,pp0)
!     main values
        DOUBLE PRECISION pdhead(pk0,pn0), g_pdhead(g_p_,pk0,pn0), &
          hu_g_pdhead(h_p_,g_p_,pk0,pn0)
        DOUBLE PRECISION pdtemp(pk0,pn0), g_pdtemp(g_p_,pk0,pn0), &
          hu_g_pdtemp(h_p_,g_p_,pk0,pn0)
        DOUBLE PRECISION pdconc(pk0,pn0), g_pdconc(g_p_,pk0,pn0), &
          hu_g_pdconc(h_p_,g_p_,pk0,pn0)
       integer gridX,gridY,gridZ
       DOUBLE PRECISION res(3,gridX,gridY,gridZ)
       DOUBLE PRECISION jacobiFullGrid(3,gridX,gridY,gridZ,g_p_)
!AW-      double precision pdhead(pK0,pN0),g_pdhead(g_p_,pK0,pN0),
!AW-     &  hu_g_pdhead(h_p_,g_p_,pK0,pN0)
!AW-      double precision pdtemp(pK0,pN0),g_pdtemp(g_p_,pK0,pN0),
!AW-     &  hu_g_pdtemp(h_p_,g_p_,pK0,pN0)
!AW-      double precision pdconc(pK0,pN0),g_pdconc(g_p_,pK0,pN0),
!AW-     &  hu_g_pdconc(h_p_,g_p_,pK0,pN0)
!       seeding sorting
        INTEGER, ALLOCATABLE :: lin_seed(:)
!     xyz-coords
        DOUBLE PRECISION pxcoord(pk0,ncoord), hu_pxcoord(h_p_,ncoord)
        integer :: kts(pk0,pn0)
        integer :: getXBlock,getYBlock,getZBlock
! ############################################
!     default: shemat, only forward evaluation
        lshem = .TRUE.
        g_lshem = .FALSE.
        hu_g_lshem = .FALSE.
        getGridJacobian=.False.
        setGridJacobian=.False.
        mpara = 0
        g_pmax_ = 0
        h_pmax_ = 0
        GO TO 200
! ############################################
!     default: g_shemat, dependency from [main-parameters]
        ENTRY g_shemat_proc(g_p_,pp0,pk0,pn0,ptmax,pparm0,g_pparm0, &
          pdhead,g_pdhead,pdtemp,g_pdtemp,pdconc,g_pdconc,pxcoord,zmin,zmax,&
          omp_inner,omp_outer)
        lshem = .TRUE.
        g_lshem = .TRUE.
        hu_g_lshem = .FALSE.
        getGridJacobian = .FALSE.
        setGridJacobian = .FALSE.
        g_pmax_ = g_p_
        h_pmax_ = 0
!        write(*,*) 'g_p_=',g_pmax_
        GO TO 200
! ############################################
!     default: hu_g_shemat, dependency from [main-parameters] and [x-coords]
        ENTRY hu_shemat_proc(h_p_,g_p_,pp0,pk0,pn0,ptmax,hu_ptmax, &
          pparm0,g_pparm0,pdhead,g_pdhead,hu_g_pdhead,pdtemp,g_pdtemp, &
          hu_g_pdtemp,pdconc,g_pdconc,hu_g_pdconc,pxcoord,hu_pxcoord,zmin,zmax,&
          omp_inner, omp_outer)
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
!     default: g_shemat, dependency from [main-parameters]
        ENTRY g_shemat_proc_get_fullgrid(g_p_,pp0,pk0,pn0,ptmax,pparm0,g_pparm0, &
          pdhead,g_pdhead,pdtemp,g_pdtemp,pdconc,g_pdconc,pxcoord,zmin,zmax,&
          omp_inner,omp_outer,gridX,gridY,gridZ,res,jacobiFullGrid)
        lshem = .TRUE.
        g_lshem = .TRUE.
        hu_g_lshem = .FALSE.
        getGridJacobian = .TRUE.
        setGridJacobian = .FALSE.
        g_pmax_ = g_p_
        h_pmax_ = 0
!        write(*,*) 'g_p_=',g_pmax_
        GO TO 200
!     default: g_shemat, dependency from [main-parameters]
        ENTRY g_shemat_proc_set_fullgrid(g_p_,pp0,pk0,pn0,ptmax,pparm0,g_pparm0, &
          pdhead,g_pdhead,pdtemp,g_pdtemp,pdconc,g_pdconc,pxcoord,zmin,zmax,&
          omp_inner,omp_outer,gridX,gridY,gridZ,res,jacobiFullGrid)
        lshem = .TRUE.
        g_lshem = .TRUE.
        hu_g_lshem = .FALSE.
        getGridJacobian=.FALSE.
        setGridJacobian = .TRUE.
        g_pmax_ = g_p_
        h_pmax_ = 0
!        write(*,*) 'g_p_=',g_pmax_
        GO TO 200

! ############################################


!     begin main routine
200     CONTINUE
        CALL sys_cputime(tsglobal)
        backuptime = tsglobal
        write_disable = .FALSE.
        write_iter_disable = .FALSE.
        nested_build = .TRUE.
        ismpl = 1
        def_binary = 'inverse'
        lib_override = .TRUE.

        OPEN(66,file='shemade.job')

! -----------------------------------------

        if (omp_inner>0 .or. omp_outer>0) fl_omp_set=.true.
        if (present(omp_inner)) fl_omp_inner=omp_inner
        if (present(omp_outer)) fl_omp_outer=omp_outer
!     read new model
10     READ(66,'(A)',end=99999) filename
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
!AW-not here -> reading INV-version !     read data
!AW-not here -> reading INV-version         IF (runmode>0) CALL read_data(filename_data,ismpl)
!AW-not here -> reading INV-version !     split units
!AW-not here -> reading INV-version         CALL read_split(filename,ismpl)
!AW-not here -> reading INV-version       call read_split(filename,ismpl)

!     3-physical main values, pN0-measurements in z direction, pK0 number of
!     boreholes
        ndata = 3*pn0*pk0
        CALL alloc_data(ismpl)
        ndata_h = pn0*pk0
        ndata_t = pn0*pk0
        ndata_c = pn0*pk0
        ndata_p = 0
        ndata_s = 0

!        write(*,*) "coords=",pxcoord

      do llk=1,pk0
! ------
!     search the cell index -> [i,j,k]
        i = getXBlock(pxcoord(llk,1))
        j = getYBlock(pxcoord(llk,2))
!        k = 1
        do ll = 1,pn0
           kts(llk,ll) = getZBlock(dble(zmin)+(dble(zmax-zmin)/dble(pn0))*dble(ll-1))
        end do
!        write(*,*) "Inserting ",i," ",j," ",kt(llk,:)
! ------
!AW-not here -> reading INV-version !       modify data-collector
!AW-not here -> reading INV-version         CALL dealloc_data(ismpl)
!
        DO ll = 1, pn0
          llp = ((ll-1)*3 + 1)+(llk-1)*pn0*3
          idata(llp+0,cid_i) = i
          idata(llp+1,cid_i) = i
          idata(llp+2,cid_i) = i
          idata(llp+0,cid_j) = j
          idata(llp+1,cid_j) = j
          idata(llp+2,cid_j) = j
          idata(llp+0,cid_k) = kts(llk,ll)
          idata(llp+1,cid_k) = kts(llk,ll)
          idata(llp+2,cid_k) = kts(llk,ll)
          idata(llp+0,cid_pv) = pv_head
          idata(llp+1,cid_pv) = pv_temp
          idata(llp+2,cid_pv) = pv_conc
          idata(llp+0,cid_si) = 0
          idata(llp+1,cid_si) = 0
          idata(llp+2,cid_si) = 1
          idata(llp+0,cid_obs) = ll
          idata(llp+1,cid_obs) = ll
          idata(llp+2,cid_obs) = ll
          ddata(llp+0,cdd_pv) = 0.D0
          ddata(llp+1,cdd_pv) = 0.D0
          ddata(llp+2,cdd_pv) = 0.D0
          ddata(llp+0,cdd_w) = 1.D0
          ddata(llp+1,cdd_w) = 1.D0
          ddata(llp+2,cdd_w) = 1.D0
          ddata(llp+0,cdd_time) = 1.D0 !ptmax*tunit
          ddata(llp+1,cdd_time) = 1.D0 !ptmax*tunit
          ddata(llp+2,cdd_time) = 1.D0 !ptmax*tunit
          ddata(llp+0,cdd_i) = pxcoord(llk,1)
          ddata(llp+1,cdd_i) = pxcoord(llk,1)
          ddata(llp+2,cdd_i) = pxcoord(llk,1)
          ddata(llp+0,cdd_j) = pxcoord(llk,2)
          ddata(llp+1,cdd_j) = pxcoord(llk,2)
          ddata(llp+2,cdd_j) = pxcoord(llk,2)
          ddata(llp+0,cdd_k) = zmin+(zmax-zmin)/pn0*(ll-1)
          ddata(llp+1,cdd_k) = zmin+(zmax-zmin)/pn0*(ll-1)
          ddata(llp+2,cdd_k) = zmin+(zmax-zmin)/pn0*(ll-1)
          sdata(llp+0,ismpl) = 0.D0
          sdata(llp+1,ismpl) = 0.D0
          sdata(llp+2,ismpl) = 0.D0
        END DO
      end do
!      write(*,*) idata(:,cid_i)
! ------
!
        CALL sys_cputime(trun)
        IF (runmode>1) CALL read_inverse(filename_inverse,ismpl)
        CALL sys_cputime(tend)
        write(*,*) ' [T] : "read_inverse" in ',tend-trun,' sec'
!
        write(*,*) ' <D> : ignore UNIT splitting when in LIBRARY mode'
!AW-not here -> reading INV-version !       split units
!AW-not here -> reading INV-version         CALL read_split_inv(filename,ismpl)
!
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
!
        CALL sys_cputime(tfglobal)
!       compute linear seeding factors
        ALLOCATE(lin_seed(g_pmax_))
        DO i = 1, g_pmax_
          lin_seed(i) = 0
        END DO
!       init offset
        k = 0
!       for pparm0
        DO i = 1, g_pmax_
          DO j = 1, pp0
            IF (g_pparm0(i,j)>0.5d0) lin_seed(i) = k +j
          END DO
        END DO
        k = k +pp0

! ------ original function (forward model)
        IF (lshem .AND. .NOT. g_lshem .AND. .NOT. hu_g_lshem) THEN
!         forward modeling (nonlinear iteration)
          CALL forward_iter(simtime_0,max_simtime,iter_inv,0,ismpl)
          GO TO 300
        END IF
! ------

! ------ first derivatives (jacobi)
        IF (lshem .AND. g_lshem .AND. .NOT. hu_g_lshem) THEN
          CALL prepare_jacobian(ismpl)
!          ***********************************************
#ifndef JACOBI_FREE
!           full Jacobi matrix computation - one times
          if (.NOT. setGridJacobian) then 
             CALL jacobi_compute(ismpl)
          else
            write (*,*) "setting simulation and jacobi with old values"
            do i = 1,gridX
             do j=1,gridY
               do k=1,gridZ
                 head(i,j,k,1)=res(1,i,j,k)
                 temp(i,j,k,1)=res(2,i,j,k)
                 conc(i,j,k,1,1)=res(3,i,j,k)
                 do ll=1,g_p_
!                   write(*,*) i,j,k,ll,jacobiFullGrid(:,i,j,k,ll)
                   g_head(i,j,k,ll)=jacobiFullGrid(1,i,j,k,ll)
                   g_temp(i,j,k,ll)=jacobiFullGrid(2,i,j,k,ll)
                   g_conc(i,j,k,1,ll)=jacobiFullGrid(3,i,j,k,ll)
                 end do
                end do
              end do
             end do
             write(*,*) "saving data"
             do i = 1,g_p_
                call g_save_data(i)
                jac(1,i) = 0.D0
                CALL g_DCOPY(ndata, sdata(1, i), g_sdata(1, i), 1,&
                main_output(i,1), jac(1,i), 1)
            end do

          end if
#endif
!          ***********************************************
          GO TO 299
        END IF
! ------

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
        WRITE(*,'(1A,5I4)') ' [I] : EFCOSS integer parameter: ',g_p_,pp0,pk0,pn0
        WRITE(*,'(1A,999I4)') ' [I] : seeding index vector: ',lin_seed
#ifndef JACOBI_FREE
!        write(*,*) "JAC(",ndata,",",mpara,")="
!        do i=1,ndata
!          write(*,*) ( jac(i,j),j=1,mpara )
!        end do
        if (getGridJacobian) then
           write(*,*) "Getting new results and jacobian"
           do i = 1,gridX
             do j=1,gridY
               do k=1,gridZ
                 res(1,i,j,k)=head(i,j,k,1)
                 res(2,i,j,k)=temp(i,j,k,1)
                 res(3,i,j,k)=conc(i,j,k,1,1)
                 do ll=1,g_p_
                   jacobiFullGrid(1,i,j,k,ll)=g_head(i,j,k,ll)
                   jacobiFullGrid(2,i,j,k,ll)=g_temp(i,j,k,ll)
                   jacobiFullGrid(3,i,j,k,ll)=g_conc(i,j,k,1,ll)
                 end do
                end do
              end do
           end do
        end if
        do llk=1,pk0
        DO ll = 1,pn0
          llp = ((ll-1)*3 + 1)+(llk-1)*pn0*3
!          write(*,*) "llp=",llp
!         - head -
          DO llh = 1, g_p_
            IF (lin_seed(llh)>0) THEN
              g_pdhead(llh,llk,ll) = jac(llp+0,lin_seed(llh))
              g_pdtemp(llh,llk,ll) = jac(llp+1,lin_seed(llh))
              g_pdconc(llh,llk,ll) = jac(llp+2,lin_seed(llh))
            ELSE
              g_pdhead(llh,llk,ll) = 0.0d0
              g_pdtemp(llh,llk,ll) = 0.0d0
              g_pdconc(llh,llk,ll) = 0.0d0
            END IF
          END DO
        END DO
        end do
        WRITE(*,*) ' '
        WRITE(*,'(1A,999e17.8)') '  [I] : dHEAD[1]/dp value= ', &
          g_pdhead(1:g_p_,1,1)
        WRITE(*,'(1A,999e17.8)') '  [I] : dTEMP[1]/dp value= ', &
          g_pdtemp(1:g_p_,1,1)
        WRITE(*,'(1A,999e17.8)') '  [I] : dCONC[1]/dp value= ', &
          g_pdconc(1:g_p_,1,1)
#endif
!
300     CONTINUE
        do llk=1,pk0
        DO ll = 1, pn0
          llp = (ll-1)*3 + 1 + (llk-1)*pn0
!         - head -
          pdhead(llk,ll) = sdata(llp+0,ismpl)
!         - temp -
          pdtemp(llk,ll) = sdata(llp+1,ismpl)
!         - conc -
          pdconc(llk,ll) = sdata(llp+2,ismpl)
        END DO
        end do
        WRITE(*,*) ' '
        WRITE(*,'(1A,1e17.8)') '  [I] : HEAD[1] value= ', pdhead(1,1)
        WRITE(*,'(1A,1e17.8)') '  [I] : TEMP[1] value= ', pdtemp(1,1)
        WRITE(*,'(1A,1e17.8)') '  [I] : CONC[1] value= ', pdconc(1,1)
        WRITE(*,*) ' '
! ------

!     --- free memory
        DEALLOCATE(lin_seed)
        CALL dealloc_arrays(ismpl)
        CALL props_end(ismpl)
        CALL dealloc_data(ismpl)

        IF (runmode>1) THEN
          CALL dealloc_inverse(ismpl)
          CALL g_dealloc_arrays(ismpl)
#ifdef JACOBI_FREE
          CALL dealloc_arrays_ad(ismpl)
#endif
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
        GO TO 10
!
! -----------------------------------------
!
!     finis terrae
99999   CONTINUE
!
        CLOSE(66)
        CALL sys_cputime(tfglobal)
        tslocal = tfglobal - tsglobal
        i = tslocal/60.D0
        WRITE(*,'(1A,1I4,1A,1F5.2,1A)') ' total cpu time: ', i, ':', &
          tslocal - dble(i)*60.D0, ' min'
        WRITE(*,*) 'RUN O.K.'
!
        RETURN
      END
