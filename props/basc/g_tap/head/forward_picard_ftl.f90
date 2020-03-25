!        Generated by TAPENADE     (INRIA, Ecuador team)
!  tapenade 3.x
!
!  Differentiation of forward_picard in forward (tangent) mode:
!   variations   of useful results: *d *e *f *concold *g *temp
!                *w *headold *x *head *dbc_data *tempold *tsal
!                *presold *conc *pres *a *b *c
!   with respect to varying inputs: *d *e *f *concold *g *temp
!                *w *headold *x *head *dbc_data *bcperiod *tempold
!                *propunit *tsal *presold *conc *pres *a *b *c
!   Plus diff mem management of: d:in e:in f:in concold:in g:in
!                temp:in w:in headold:in x:in head:in dbc_data:in
!                bcperiod:in tempold:in propunit:in tsal:in presold:in
!                conc:in pres:in simtime:in a:in b:in c:in
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
SUBROUTINE g_FORWARD_PICARD(iter_time, ismpl)
  USE ARRAYS

  USE g_ARRAYS

  USE MOD_GENRL
  USE MOD_GENRLC
  USE MOD_FLOW
  USE MOD_TEMP
  USE MOD_CONC
  USE MOD_TIME

  USE g_MOD_TIME

  USE MOD_DATA
  USE MOD_LINFOS
  IMPLICIT NONE
  INTEGER :: ismpl
  INTEGER :: i
  INTEGER :: iter_nl
! Maximal difference and RMS-difference between old and new
! variable array
  DOUBLE PRECISION :: difmaxf, difrmsf
  DOUBLE PRECISION :: g_difmaxf
  DOUBLE PRECISION :: difmaxt, difrmst
  DOUBLE PRECISION :: g_difmaxt
  DOUBLE PRECISION :: difmaxc, difrmsc
  DOUBLE PRECISION :: g_difmaxc
  DOUBLE PRECISION :: difmaxs, difrmss
  EXTERNAL CONVERGED, LBLANK
!
  INTEGER :: ijk, mode, iter_time, LBLANK, species
  DOUBLE PRECISION :: neu_max, pec_max, cou_max
  DOUBLE PRECISION :: difmaxfold, difmaxtold, difmaxcold
  DOUBLE PRECISION :: g_difmaxfold, g_difmaxtold, g_difmaxcold
  DOUBLE PRECISION :: relaxold, relaxf, relaxt, relaxc, relaxs
  DOUBLE PRECISION :: g_relaxold, g_relaxf, g_relaxt, g_relaxc
! Local non linear tolerance variables for convergence check
  DOUBLE PRECISION :: loc_nltolf, loc_nltolt, loc_nltolc, loc_nltols
!
  LOGICAL :: CONVERGED, lcon_sum
  EXTERNAL DUMMY
  DOUBLE PRECISION :: arg1
  DOUBLE PRECISION :: g_arg1
! Initialisation
! --------------
! standard output
  IF (linfos(3) .GE. 2) WRITE(*, *) ' ... forward(picard)'
! initial values for variables
  ijk = i0*j0*k0
  mode = 0
  iter_nl = 0
! initialize nonlinear solver parameters: flow
  loc_nltolf = nltolf
  difmaxf = 1.0d9*loc_nltolf
  IF (linfos(3) .GE. 2) WRITE(*, *) 'difmaxf = ', difmaxf
  relaxf = nlrelaxf
  difrmsf = big
! initialize nonlinear solver parameters: temp
  loc_nltolt = nltolt
  difmaxt = 1.0d9*loc_nltolt
  difrmst = big
  relaxt = nlrelaxt
! initialize nonlinear solver parameters: conc
  loc_nltolc = nltolc
  difmaxc = 1.0d9*loc_nltolc
  difrmsc = big
  relaxc = nlrelaxc
! initialize nonlinear solver parameters: satn
  loc_nltols = nltols
  difmaxs = 1.0d9*loc_nltols
  difrmss = big
  relaxs = nlrelaxs
! init history [1] for flow
  lcon(1, ismpl) = CONVERGED(difmaxf, loc_nltolf, -1, conv_hlen, &
&   conv_hmax, conv_history(1, 1, ismpl), conv_chlen(1, ismpl), &
&   conv_ipos(1, ismpl))
! init history [2] for temperature
  lcon(2, ismpl) = CONVERGED(difmaxt, loc_nltolt, -2, conv_hlen, &
&   conv_hmax, conv_history(1, 1, ismpl), conv_chlen(1, ismpl), &
&   conv_ipos(1, ismpl))
! init history for concentration
  DO i=4,ntrans+3
    lcon(i, ismpl) = CONVERGED(difmaxc, loc_nltolc, -i, conv_hlen, &
&     conv_hmax, conv_history(1, 1, ismpl), conv_chlen(1, ismpl), &
&     conv_ipos(1, ismpl))
  END DO
! pre-set to enter the computation at least once
  lcon_sum = .false.
  g_relaxc = 0.D0
  g_relaxf = 0.D0
  g_relaxt = 0.D0
  g_difmaxc = 0.D0
  g_difmaxf = 0.D0
  g_difmaxt = 0.D0
! -------- begin nonlinear iteration
! ----------------------------------
!     LOOP = ITERATION lcon_sum - error !!!
!     LOOP = ITERATION difmaxf,difmaxt,difmaxc - wrong !!!
!$TAF LOOP = ITERATION head,temp,conc,difmaxf,difmaxt,difmaxc
  DO WHILE (.NOT.lcon_sum .AND. iter_nl .LT. maxiter_nl)
! --------
! Preprocessing
! -------------
! loop Counter for outer nonlinear iteration
    iter_nl = iter_nl + 1
! save old variable arrays for checking difference later
    CALL g_OLD_SAVE(cgen_fw, ismpl)
! static relaxation (only flow and temperature)
    IF (transient .AND. tr_switch(ismpl)) CALL g_STATIC_RELAXATION(ijk&
&                                                              , ismpl)
! user directory functions
    CALL g_CALC_USER(ismpl)
#ifdef head_base
! update pressure for current head
    CALL g_HEAD2PRES(1, ismpl)
#endif
! Flow equation
! -------------
! Call solver
    IF (head_active .OR. pres_active) THEN
#ifdef head_base
! head computation
      CALL g_CALC_HEAD(ismpl)
#endif
    ELSE
      lcon(1, ismpl) = .true.
    END IF
! evaluate flow
    IF (head_active .OR. pres_active) THEN
! Neumann numbers
      IF (transient .AND. tr_switch(ismpl) .AND. linfos(3) .GE. 1) THEN
#ifdef head_base
        CALL NEUMANN_HEAD(neu_max, ismpl)
#endif
!?            CALL Courant(Cou_max,ismpl)
      END IF
! compute difference
      g_difmaxfold = g_difmaxf
      difmaxfold = difmaxf
#ifdef head_base
      IF (linfos(3) .GE. 2) WRITE(*, *) 'Check change with difmaxf = ', &
&                           difmaxf
      CALL g_CHECK_CHANGE(mode, pv_head, loc_nltolf, difrmsf, difmaxf&
&                     , g_difmaxf, i0, j0, k0, head(1, 1, 1, ismpl), &
&                     g_head(1, 1, 1, ismpl), headold(1, cgen_fw, &
&                     ismpl), g_headold(1, cgen_fw, ismpl), ismpl)
#endif
! check for convergence head/pressure
      lcon(1, ismpl) = CONVERGED(difmaxf, loc_nltolf, 1, conv_hlen, &
&       conv_hmax, conv_history(1, 1, ismpl), conv_chlen(1, ismpl), &
&       conv_ipos(1, ismpl))
! adaptive relaxation
      IF (iter_nl .GE. 2) THEN
! relaxation for better convergence
        IF (nladapt .EQ. 1) THEN
          g_relaxold = g_relaxf
          relaxold = relaxf
! DANGER not working for saturation !
          CALL g_NL_RELAX(iter_nl, difmaxf, g_difmaxf, difmaxfold, &
&                     g_difmaxfold, nlmaxf, relaxf, g_relaxf, &
&                     relaxold, g_relaxold, ismpl)
        END IF
#ifdef head_base
! relaxing head
        CALL g_DSCAL(ijk, relaxf, g_relaxf, head(1, 1, 1, ismpl), &
&                g_head(1, 1, 1, ismpl), 1)
        g_arg1 = -g_relaxf
        arg1 = 1.0d0 - relaxf
        CALL g_DAXPY(ijk, arg1, g_arg1, headold(1, cgen_fw, ismpl), &
&                g_headold(1, cgen_fw, ismpl), 1, head(1, 1, 1, ismpl)&
&                , g_head(1, 1, 1, ismpl), 1)
#endif
      END IF
    END IF
! Heat equation
! -------------
! Call solver
    IF (temp_active) THEN
      CALL g_CALC_TEMP(ismpl)
    ELSE
      lcon(2, ismpl) = .true.
    END IF
! Evaluate temperature
    IF (temp_active) THEN
      IF (transient .AND. tr_switch(ismpl) .AND. linfos(3) .GE. 1) THEN
        CALL COURANT(cou_max, ismpl)
        CALL NEUMANN_TEMP(neu_max, ismpl)
        CALL PECLET_TEMP(pec_max, ismpl)
      END IF
! compute difference
      g_difmaxtold = g_difmaxt
      difmaxtold = difmaxt
      CALL g_CHECK_CHANGE(mode, pv_temp, loc_nltolt, difrmst, difmaxt&
&                     , g_difmaxt, i0, j0, k0, temp(1, 1, 1, ismpl), &
&                     g_temp(1, 1, 1, ismpl), tempold(1, cgen_fw, &
&                     ismpl), g_tempold(1, cgen_fw, ismpl), ismpl)
! check for convergence
      lcon(2, ismpl) = CONVERGED(difmaxt, loc_nltolt, 2, conv_hlen, &
&       conv_hmax, conv_history(1, 1, ismpl), conv_chlen(1, ismpl), &
&       conv_ipos(1, ismpl))
! adaptive relaxation
      IF (iter_nl .GE. 2) THEN
! relaxation for better convergence
        IF (nladapt .EQ. 1) THEN
          g_relaxold = g_relaxt
          relaxold = relaxt
          CALL g_NL_RELAX(iter_nl, difmaxt, g_difmaxt, difmaxtold, &
&                     g_difmaxtold, nlmaxt, relaxt, g_relaxt, &
&                     relaxold, g_relaxold, ismpl)
        END IF
        CALL g_DSCAL(ijk, relaxt, g_relaxt, temp(1, 1, 1, ismpl), &
&                g_temp(1, 1, 1, ismpl), 1)
        g_arg1 = -g_relaxt
        arg1 = 1.0d0 - relaxt
        CALL g_DAXPY(ijk, arg1, g_arg1, tempold(1, cgen_fw, ismpl), &
&                g_tempold(1, cgen_fw, ismpl), 1, temp(1, 1, 1, ismpl)&
&                , g_temp(1, 1, 1, ismpl), 1)
      END IF
    END IF
! Open spot for future implementation of another equation at
! lcon(3,ismpl)
    lcon(3, ismpl) = .true.
! Transport equation
! -------------
    IF (trac_active) THEN
      IF (linfos(3) .GE. 2) WRITE(*, *) ' ... calc_conc (tracer)'
    END IF
    DO species=1,ntrac
! Call solver
      IF (trac_active) THEN
        CALL g_CALC_CONC(species, ismpl)
      ELSE
        lcon(3+species, ismpl) = .true.
      END IF
! evaluate concentration
      IF (trac_active) THEN
        IF (transient .AND. tr_switch(ismpl) .AND. linfos(3) .GE. 1 &
&           .AND. species .EQ. 1) THEN
          CALL COURANT(cou_max, ismpl)
          CALL NEUMANN_CONC(neu_max, ismpl)
          CALL PECLET_CONC(pec_max, ismpl)
        END IF
! compute difference
        g_difmaxcold = g_difmaxc
        difmaxcold = difmaxc
! Check change
        CALL g_CHECK_CHANGE(mode, pv_conc, loc_nltolc, difrmsc, &
&                       difmaxc, g_difmaxc, i0, j0, k0, conc(1, 1, 1, &
&                       species, ismpl), g_conc(1, 1, 1, species, &
&                       ismpl), concold(1, species, cgen_fw, ismpl), &
&                       g_concold(1, species, cgen_fw, ismpl), ismpl)
! check for convergence
        lcon(3+species, ismpl) = CONVERGED(difmaxc, loc_nltolc, 3 + &
&         species, conv_hlen, conv_hmax, conv_history(1, 1, ismpl), &
&         conv_chlen(1, ismpl), conv_ipos(1, ismpl))
! adaptive relaxation
        IF (iter_nl .GE. 2) THEN
! relaxation for better convergence
          IF (nladapt .EQ. 1) THEN
            g_relaxold = g_relaxc
            relaxold = relaxc
            CALL g_NL_RELAX(iter_nl, difmaxc, g_difmaxc, difmaxcold&
&                       , g_difmaxcold, nlmaxc, relaxc, g_relaxc, &
&                       relaxold, g_relaxold, ismpl)
          END IF
          CALL g_DSCAL(ijk, relaxc, g_relaxc, conc(1, 1, 1, species&
&                  , ismpl), g_conc(1, 1, 1, species, ismpl), 1)
          g_arg1 = -g_relaxc
          arg1 = 1.0d0 - relaxc
          CALL g_DAXPY(ijk, arg1, g_arg1, concold(1, species, &
&                  cgen_fw, ismpl), g_concold(1, species, cgen_fw, &
&                  ismpl), 1, conc(1, 1, 1, species, ismpl), g_conc(1&
&                  , 1, 1, species, ismpl), 1)
        END IF
        conc_conv(species, ismpl) = difmaxc
      END IF
    END DO
! Postprocessing
! -------------
!         check whether pres/temp/(conc) in domain of props validity
    CALL g_CHECK_DOMAIN(ismpl)
!         summarise all convergency criteria
    lcon_sum = .true.
    DO i=1,conv_hmax
      IF (.NOT.lcon(i, ismpl)) lcon_sum = .false.
    END DO
! generate convergency output
    IF (linfos(3) .GE. 1) THEN
#ifdef head_base
      WRITE(*, '(1A,1I6,2(1A,1e16.8))') '  [I] : iter_nl =', iter_nl, &
&     ', difmaxh =', difmaxf, ', difmaxt =', difmaxt
#endif
      IF (ntrac .GE. 1) WRITE(*, '(1A,99e16.8)') '        difmaxc[*] ='&
&                       , (conc_conv(i, ismpl), i=1,ntrac)
    END IF
    IF (runmode .GE. -1 .AND. (.NOT.write_iter_disable)) THEN
      IF (linfos(3) .GE. 2) WRITE(*, '(3A)') '  [W] : "', status_log(1:&
&                           LBLANK(status_log)), '"'
      OPEN(76, file=status_log, status='unknown', position='append') 
#ifdef head_base
      WRITE(76, '(I11,99e16.8)') iter_nl, difmaxf, difmaxt, (conc_conv(i&
&     , ismpl), i=1,ntrac)
#endif
      CLOSE(76) 
    END IF
! calculate total salinity (used only in property module
! basc)
    CALL g_SET_TSAL(ismpl)
#ifdef DEBUG
    DO i=1,n_debugout
      IF (debugout(1, i) .EQ. iter_time .AND. debugout(2, i) .EQ. &
&         iter_nl) THEN
        WRITE(project_sfx(ismpl), '(1A1,1I5.5,1A1,1I3.3)') '_', &
&       iter_time, 'x', iter_nl
!              CALL write_tecdiff(-2,ismpl)
        project_sfx(ismpl) = ' '
      END IF
    END DO
#endif
  END DO
! --------
! -------- end nonlinear iteration
! --------------------------------
! Variable time step
  CALL SET_VAR_DELTAT(iter_nl, ismpl)
! Standard output
  IF (transient .AND. tr_switch(ismpl)) THEN
    IF (linfos(3) .GE. 2) WRITE(*, '(1A,1I10,1A)') 'timestep ', &
&                         iter_time, ', leaving nonlinear iteration'
  END IF
  RETURN
END SUBROUTINE g_FORWARD_PICARD

