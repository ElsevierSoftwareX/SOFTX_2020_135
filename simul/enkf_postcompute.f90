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

!> @brief Postcomputation of samples or mean
!> @details
!> The simulation is restarted using the updated parameters.
subroutine enkf_postcompute()

  use arrays, only: &
       mpara,&
       head,&
       temp,&
       conc,&
       pres,&
       dbc_data,&
       dbc_dataold,&
       simtime,&
       main_input,&
       main_output,&
       nbc_data,&
       project_sfx

  use mod_genrl, only: &
       cgen_opti,&
       head_active,&
       temp_active,&
       chem_active,&
       pres_active,&
       nsmpl,&
       i0,&
       j0,&
       k0

  use mod_genrlc, only: &
       project

  use mod_conc, only: &
       ntrac

  use mod_time, only: &
       simtime_0,&
       max_simtime,&
       itimestep_0

  use mod_simul, only: &
       maxiter_enkf,&
       enkf_post,&
       ssample_outdir

  use mod_enkf, only: &
       nrens,&
       iter_out,&
       itimestep_ismpl,&
       project_org,&
       ttstart,&
       ttend
  
  implicit none

  double precision :: get_optip
  external get_optip

  integer :: rc_von, rc_bis, irens, i, j, l, l1

  rc_von = 0
  rc_bis = 0
  ! case switch 1: only samples
  IF (enkf_post==1) THEN
     rc_von = 1
     rc_bis = nrens
  END IF
  ! case switch 2: only mean
  IF (enkf_post==2) THEN
     rc_von = nsmpl
     rc_bis = nsmpl
  END IF
  ! case switch 2: samples and mean
  IF (enkf_post==3) THEN
     rc_von = 1
     rc_bis = nsmpl
  END IF
  ! sanity check
  IF (rc_bis<=0) THEN
     WRITE(*,'()') 'error: internal bug - ENKF postcompute switch not set!'
     STOP
  END IF

  !           switch to sample-subdir
  project_org = project
  project = trim(ssample_outdir)//trim(project_org)

  do irens = 1, nsmpl
     project_sfx(irens) = '_postcomp'
  end do

  !
  !           - INIT to recompute ensembles -
  !           !!! parallelisation needs to be inserted here !!!
  DO irens = rc_von, rc_bis
     ! init state-variables state (before)
     CALL old_restore(cgen_opti,irens)
     ! init "DBC" state (before)
     DO i = 1, nbc_data
        dbc_data(i,1,irens) = dbc_dataold(i)
     END DO
  END DO

  call enkf_log(6,0)

  !           set initial start time for simulation of all ensemble members
  ttstart = simtime_0
  ttend = max_simtime
  DO irens = rc_von, rc_bis
     simtime(irens) = simtime_0
  END DO
  !           restore start timestep
  itimestep_0 = itimestep_ismpl
  !
  IF (rc_von/=nsmpl) THEN
     !             update changes from the parameter space to the main-parameter input vector
     DO j = 1, rc_bis
        DO i = 1, mpara
           main_input(i,j) = get_optip(i,j)
        END DO
     END DO
     !             ---- compute all ensembles
     CALL forward_multi_compute(ttstart, ttend, iter_out, 'enkf+run+mean', rc_bis)
     !             -------------------
     !
     !             -> compute mean of state-variables
     !             - HEAD -
     IF (head_active) THEN
        CALL set_dval(i0*j0*k0,0.0d0,head(1,1,1,nsmpl))
        DO l = 1, nrens
           CALL daxpy(i0*j0*k0,1.0d0/dble(nrens),head(1,1,1,l),1,head(1,1,1,nsmpl),1)
        ENDDO
     ELSE
        !               when not active copy any value to the mean
        CALL dcopy(i0*j0*k0,head(1,1,1,1),1,head(1,1,1,nsmpl),1)
     END IF
     !             - TEMP -
     IF (temp_active) THEN
        CALL set_dval(i0*j0*k0,0.0d0,temp(1,1,1,nsmpl))
        DO l = 1, nrens
           CALL daxpy(i0*j0*k0,1.0d0/dble(nrens),temp(1,1,1,l),1,temp(1,1,1,nsmpl),1)
        ENDDO
     ELSE
        !               when not active copy any value to the mean
        CALL dcopy(i0*j0*k0,temp(1,1,1,1),1,temp(1,1,1,nsmpl),1)
     END IF
     !             - CONC -
     DO l1 = 1, ntrac
        IF (chem_active) THEN
           CALL set_dval(i0*j0*k0,0.0d0,conc(1,1,1,l1,nsmpl))
           DO l = 1, nrens
              CALL daxpy(i0*j0*k0,1.0d0/dble(nrens),conc(1,1,1,l1,l),1,conc(1,1,1,l1,nsmpl),1)
           ENDDO
        ELSE
           !                 when not active copy any value to the mean
           CALL dcopy(i0*j0*k0,conc(1,1,1,l1,1),1,conc(1,1,1,l1,nsmpl),1)
        END IF
     ENDDO
     !             - PRES -
     IF (pres_active) THEN
        CALL set_dval(i0*j0*k0,0.0d0,pres(1,1,1,nsmpl))
        DO l = 1, nrens
           CALL daxpy(i0*j0*k0,1.0d0/dble(nrens),pres(1,1,1,l),1,pres(1,1,1,nsmpl),1)
        ENDDO
     ELSE
        !               when not active copy any value to the mean
        CALL dcopy(i0*j0*k0,pres(1,1,1,1),1,pres(1,1,1,nsmpl),1)
     END IF
     !             write the average model to file, index number -3 for "mean"-prefix
     CALL forward_write(-5,nsmpl)
     !             -------------------
  ELSE
     !             update changes from the parameter space to the main-parameter input vector
     DO i = 1, mpara
        main_input(i,nsmpl) = get_optip(i,nsmpl)
     END DO
     !             compute only the mean of all realisations
     CALL omp_forward_multi_compute(main_input(1,nsmpl), main_output(1,nsmpl), ttstart, ttend, &
          iter_out, 'enkf+run+mean', nsmpl, nrens, nsmpl)
  END IF
  
  project = project_org

end subroutine enkf_postcompute
