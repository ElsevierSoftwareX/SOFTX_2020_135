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

!>    @brief nonlinear picard iteration for flow, heat and transport equation
!>    @param[in] iter_time forward iteration counter
!>    @param[in] ismpl local sample index
!>    @details
!> nonlinear iteration loop (convergency) for steady state and\n
!> transient case (one time step)\n
      SUBROUTINE forward_picard(iter_time,ismpl)

        use arrays
        use mod_genrl
        use mod_genrlc
        use mod_flow
        use mod_temp
        use mod_conc
        use mod_time
        use mod_data
        use mod_linfos

        IMPLICIT NONE

        integer :: ismpl
        integer :: i
        integer :: iter_nl

        ! Maximal difference and RMS-difference between old and new
        ! variable array
        double precision difmaxf, difrmsf
        double precision difmaxt, difrmst
        double precision difmaxc, difrmsc
        double precision difmaxs, difrmss
!
        INTEGER ijk, mode, iter_time, lblank, species
        DOUBLE PRECISION neu_max, pec_max, cou_max
        DOUBLE PRECISION difmaxfold, difmaxtold, difmaxcold
        DOUBLE PRECISION relaxold, relaxf, relaxt, relaxc, relaxs

        ! Local non linear tolerance variables for convergence check
        DOUBLE PRECISION loc_nltolf,loc_nltolt,loc_nltolc,loc_nltols
!
        LOGICAL converged, lcon_sum
        EXTERNAL converged, lblank


        ! Initialisation
        ! --------------

        ! standard output
        IF (linfos(3)>=2) WRITE(*,*) ' ... forward(picard)'

        ! initial values for variables
        ijk = i0*j0*k0
        mode = 0
        iter_nl = 0

        ! initialize nonlinear solver parameters: flow
        loc_nltolf = nltolf
#ifdef pres_base
!       convert [MPa] into [Pa]
        loc_nltolf = nltolf*pa_conv
#endif
        difmaxf = 1.0D9*loc_nltolf
        if (linfos(3)>=2) then
           WRITE(*,*) 'difmaxf = ', difmaxf
        end if
        relaxf = nlrelaxf
        difrmsf = big

        ! initialize nonlinear solver parameters: temp
        loc_nltolt = nltolt
        difmaxt = 1.0D9*loc_nltolt
        difrmst = big
        relaxt = nlrelaxt

        ! initialize nonlinear solver parameters: conc
        loc_nltolc = nltolc
        difmaxc = 1.0D9*loc_nltolc
        difrmsc = big
        relaxc = nlrelaxc

        ! initialize nonlinear solver parameters: satn
        loc_nltols = nltols
        difmaxs = 1.0D9*loc_nltols
        difrmss = big
        relaxs = nlrelaxs

        ! init history [1] for flow
        lcon(1,ismpl) = converged(difmaxf,loc_nltolf,-1,conv_hlen, &
          conv_hmax,conv_history(1,1,ismpl),conv_chlen(1,ismpl), &
          conv_ipos(1,ismpl))
        ! init history [2] for temperature
        lcon(2,ismpl) = converged(difmaxt,loc_nltolt,-2,conv_hlen, &
          conv_hmax,conv_history(1,1,ismpl),conv_chlen(1,ismpl), &
          conv_ipos(1,ismpl))
        ! init history for concentration
        DO i = 4, ntrans + 3
          lcon(i,ismpl) = converged(difmaxc,loc_nltolc,-i,conv_hlen, &
            conv_hmax,conv_history(1,1,ismpl),conv_chlen(1,ismpl), &
            conv_ipos(1,ismpl))
        END DO

        ! pre-set to enter the computation at least once
        lcon_sum = .FALSE.

! -------- begin nonlinear iteration
! ----------------------------------
!     LOOP = ITERATION lcon_sum - error !!!
!     LOOP = ITERATION difmaxf,difmaxt,difmaxc - wrong !!!
!$TAF LOOP = ITERATION head,temp,conc,difmaxf,difmaxt,difmaxc
        DO WHILE ( .NOT. lcon_sum .AND. iter_nl<maxiter_nl)
! --------

          ! Preprocessing
          ! -------------

          ! loop Counter for outer nonlinear iteration
          iter_nl = iter_nl + 1

          ! save old variable arrays for checking difference later
          CALL old_save(cgen_fw,ismpl)

          ! static relaxation (only flow and temperature)
          IF (transient .AND. tr_switch(ismpl)) THEN
            call static_relaxation(ijk,ismpl)
          END IF

          ! user directory functions
          CALL calc_user(ismpl)

#ifdef head_base
          ! update pressure for current head
          CALL head2pres(1,ismpl)
#endif
#ifdef pres_base
          ! update head from current pressure
          CALL pres2head(1,ismpl)
#endif

          ! Flow equation
          ! -------------

          ! Call solver
          IF (head_active .OR. pres_active) THEN
#ifdef head_base
            ! head computation
             CALL calc_head(ismpl)
#endif
#ifdef pres_base
             ! pressure computation
            CALL calc_pres(ismpl)
#endif
          ELSE
            lcon(1,ismpl) = .TRUE.
          END IF

          ! evaluate flow
          IF (head_active .OR. pres_active) THEN

            ! Neumann numbers
            IF (transient .AND. tr_switch(ismpl) .AND. linfos(3)>=1) THEN
#ifdef head_base
              CALL neumann_head(neu_max,ismpl)
#endif
#ifdef pres_base
              CALL neumann_pres(neu_max,ismpl)
#endif
!?            CALL Courant(Cou_max,ismpl)
            END IF

            ! compute difference
            difmaxfold = difmaxf
#ifdef head_base
            if (linfos(3)>=2) then
              WRITE(*,*) 'Check change with difmaxf = ', difmaxf
            end if
            CALL check_change(mode,pv_head,loc_nltolf,difrmsf,difmaxf,i0,j0,k0, &
                 head(1,1,1,ismpl),headold(1,cgen_fw,ismpl),ismpl)
#endif
#ifdef pres_base
            CALL check_change(mode,pv_pres,loc_nltolf,difrmsf,difmaxf,i0,j0,k0, &
              pres(1,1,1,ismpl),presold(1,cgen_fw,ismpl),ismpl)
!           Include non-linear residual in error
            difmaxf=  ABS(difmaxf)
#endif

            ! check for convergence head/pressure
            lcon(1,ismpl) = converged(difmaxf,loc_nltolf,1,conv_hlen, &
              conv_hmax,conv_history(1,1,ismpl),conv_chlen(1,ismpl), &
              conv_ipos(1,ismpl))

            ! adaptive relaxation
            IF (iter_nl>=2) THEN

              ! relaxation for better convergence
              IF (nladapt==1) THEN
                relaxold = relaxf
                ! DANGER not working for saturation !
                CALL nl_relax(iter_nl,difmaxf,difmaxfold,nlmaxf,relaxf, &
                  relaxold,ismpl)
              END IF
#ifdef head_base
              ! relaxing head
              CALL dscal(ijk,relaxf,head(1,1,1,ismpl),1)
              CALL daxpy(ijk,1.0D0-relaxf,headold(1,cgen_fw,ismpl),1, &
                head(1,1,1,ismpl),1)
#endif
#ifdef pres_base
              ! relaxing pressure
              CALL dscal(ijk,relaxf,pres(1,1,1,ismpl),1)
              CALL daxpy(ijk,1.0D0-relaxf,presold(1,cgen_fw,ismpl),1, &
                pres(1,1,1,ismpl),1)
#endif
            END IF

          END IF

          ! Heat equation
          ! -------------

          ! Call solver
          IF (temp_active) THEN
            CALL calc_temp(ismpl)
          ELSE
            lcon(2,ismpl) = .TRUE.
          END IF

          ! Evaluate temperature
          IF (temp_active) THEN

            IF (transient .AND. tr_switch(ismpl) .AND. linfos(3)>=1) &
                THEN
              CALL courant(cou_max,ismpl)
              CALL neumann_temp(neu_max,ismpl)
              CALL peclet_temp(pec_max,ismpl)
            END IF

            ! compute difference
            difmaxtold = difmaxt
            CALL check_change(mode,pv_temp,loc_nltolt,difrmst,difmaxt,i0,j0,k0, &
              temp(1,1,1,ismpl),tempold(1,cgen_fw,ismpl),ismpl)

            ! check for convergence
            lcon(2,ismpl) = converged(difmaxt,loc_nltolt,2,conv_hlen, &
              conv_hmax,conv_history(1,1,ismpl),conv_chlen(1,ismpl), &
              conv_ipos(1,ismpl))

            ! adaptive relaxation
            IF (iter_nl>=2) THEN
              ! relaxation for better convergence
              IF (nladapt==1) THEN
                relaxold = relaxt
                CALL nl_relax(iter_nl,difmaxt,difmaxtold,nlmaxt,relaxt, &
                  relaxold,ismpl)
              END IF
              CALL dscal(ijk,relaxt,temp(1,1,1,ismpl),1)
              CALL daxpy(ijk,1.0D0-relaxt,tempold(1,cgen_fw,ismpl),1, &
                temp(1,1,1,ismpl),1)
            END IF

          END IF

          ! Open spot for future implementation of another equation at
          ! lcon(3,ismpl)
          lcon(3,ismpl) = .TRUE.

          ! Transport equation
          ! -------------

          IF (trac_active) THEN
            IF (linfos(3)>=2) WRITE(*,*) ' ... calc_conc (tracer)'
          END IF
          DO species = 1, ntrac

            ! Call solver
            IF (trac_active) THEN
              CALL calc_conc(species,ismpl)
            ELSE
              lcon(3+species,ismpl) = .TRUE.
            END IF

            ! evaluate concentration
            IF (trac_active) THEN

              IF (transient .AND. tr_switch(ismpl) .AND. &
                  linfos(3)>=1 .AND. species==1) THEN
                CALL courant(cou_max,ismpl)
                CALL neumann_conc(neu_max,ismpl)
                CALL peclet_conc(pec_max,ismpl)
              END IF

              ! compute difference
              difmaxcold = difmaxc

              ! Check change
              CALL check_change(mode,pv_conc,loc_nltolc,difrmsc,difmaxc,i0,j0,k0, &
                conc(1,1,1,species,ismpl),concold(1,species,cgen_fw, &
                ismpl),ismpl)

              ! check for convergence
              lcon(3+species,ismpl) = converged(difmaxc,loc_nltolc, &
                3+species,conv_hlen,conv_hmax,conv_history(1,1,ismpl), &
                conv_chlen(1,ismpl),conv_ipos(1,ismpl))

              ! adaptive relaxation
              IF (iter_nl>=2) THEN
                ! relaxation for better convergence
                IF (nladapt==1) THEN
                  relaxold = relaxc
                  CALL nl_relax(iter_nl,difmaxc,difmaxcold,nlmaxc, &
                    relaxc,relaxold,ismpl)
                END IF
                CALL dscal(ijk,relaxc,conc(1,1,1,species,ismpl),1)
                CALL daxpy(ijk,1.0D0-relaxc,concold(1,species,cgen_fw, &
                  ismpl),1,conc(1,1,1,species,ismpl),1)
              END IF

              conc_conv(species,ismpl) = difmaxc

            END IF

          END DO

          ! Postprocessing
          ! -------------

!         check whether pres/temp/(conc) in domain of props validity
          CALL check_domain(ismpl)

!         summarise all convergency criteria
          lcon_sum = .TRUE.
          DO i = 1, conv_hmax
            IF ( .NOT. lcon(i,ismpl)) lcon_sum = .FALSE.
          END DO

          ! generate convergency output
          IF (linfos(3)>=1) THEN
#ifdef head_base
            WRITE(*,'(1A,1I6,2(1A,1e16.8))') '  [I] : iter_nl =', iter_nl, &
              ', difmaxh =', difmaxf, ', difmaxt =', difmaxt
#endif
#ifdef pres_base
            WRITE(*,'(1A,1I6,2(1A,1e16.8))') '  [I] : iter_nl =', iter_nl, &
              ', difmaxp =', difmaxf*Pa_conv1, ', difmaxt =', difmaxt, &
              ', difmaxs =', difmaxs
#endif
            IF (ntrac>=1) WRITE(*,'(1A,99e16.8)') '        difmaxc[*] =', (conc_conv(i,ismpl),i=1,ntrac)
          END IF
          IF ((runmode>=-1) .AND. ( .NOT. write_iter_disable)) THEN
            IF (linfos(3)>=2) WRITE(*,'(3A)') '  [W] : "', status_log(1:lblank(status_log)), '"'
            OPEN(76,file=status_log,status='unknown',position='append')
#ifdef head_base
            WRITE(76,'(I11,99e16.8)') iter_nl, difmaxf, difmaxt,(conc_conv(i,ismpl),i=1,ntrac)
#endif
#ifdef pres_base
            WRITE(76,'(I11,99e16.8)') iter_nl, difmaxf*Pa_conv1, difmaxt, difmaxs, (conc_conv(i,ismpl),i=1,ntrac)
#endif
            CLOSE(76)
          END IF

          ! calculate total salinity (used only in property module
          ! basc)
          CALL set_tsal(ismpl)

#ifdef DEBUG
          DO i = 1, n_debugout
            IF (debugout(1,i)==iter_time .AND. debugout(2,i)==iter_nl) THEN
              WRITE(project_sfx(ismpl),'(1A1,1I5.5,1A1,1I3.3)') '_', iter_time, 'x', iter_nl
              CALL write_tecdiff(-2,ismpl)
              project_sfx(ismpl) = ' '
            END IF
          END DO
#endif

! --------
        END DO
! -------- end nonlinear iteration
! --------------------------------

        ! Variable time step
        call set_var_deltat(iter_nl, ismpl)

        ! Standard output
        IF (transient .AND. tr_switch(ismpl)) THEN
          IF (linfos(3)>=2) WRITE(*,'(1A,1I10,1A)') 'timestep ', iter_time, ', leaving nonlinear iteration'
        END IF

        RETURN

      END SUBROUTINE forward_picard
