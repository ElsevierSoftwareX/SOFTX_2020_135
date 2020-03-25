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

!> @brief Computation and output of ensemble means.
!> @details
!> The means of variables and/or parameters are calculated, saved
!> in propunit or the variable array under ensemble member nsmpl
!> and finally the forward_write wrapper for typical output is
!> called.
subroutine enkf_compute_mean()

  use arrays, only: &
       nprop,&
       nunits,&
       propunit,&
       mpara,&
       seed_para,&
       firstidx,&
       lastidx,&
       bc_firstidx,&
       bc_lastidx,&
       head,&
       temp,&
       conc,&
       pres

  use mod_genrl, only: &
       i0,&
       j0,&
       k0,&
       nsmpl,&
       head_active,&
       temp_active,&
       chem_active,&
       pres_active

  use mod_genrlc, only: &
       project

  use mod_conc, only: &
       ntrac

  use mod_linfos, only: &
       linfos

  use mod_simul, only: &
       ssample_outdir

  use mod_enkf, only: &
       nrens,&
       project_org

  implicit none

  integer :: s_k, s_u
  
  integer :: param_enabled
  external param_enabled

  integer :: i, j, l, l1
  
  WRITE(*,'(1A)') '  [I] : compute (lin/log) mean values'
  !         Determine average of ensemble for all parameters and
  !         variables and store in sm_max=nsmpl for Shemat-output
  !         -> compute mean of rock-props and BC-base values, !!
  !         -> currently no mean for time depnt BC's ([bcperiod])
  !         - default linear mean for all parameters -
  DO j = 1, nprop
     DO i = 1, nunits
        propunit(i,j,nsmpl) = 0.D0
        DO l = 1, nrens
           propunit(i,j,nsmpl) = propunit(i,j,nsmpl) + &
                propunit(i,j,l)/dble(nrens)
        END DO
     END DO
  END DO
  !         - logarithmic mean for specific parameters (overwrites the
  !         - default linear mean) -
  DO i = 1, mpara
     s_k = seed_para(1,i)
     s_u = seed_para(2,i)
     IF ((s_k>=firstidx .AND. s_k<=lastidx) .OR. (s_k>=bc_firstidx .AND. s_k<=bc_lastidx)) THEN
        IF (param_enabled(i) == 2) THEN
           IF (linfos(2)>=1) WRITE(*,'(1A,1I2,1A,1I6,1A)') &
                ' [I] : compute logarithmic mean for [component=',s_k,', unit=',s_u,']'
           !       only for active logarithmic parameters
           propunit(s_u,s_k,nsmpl) = 0.D0
           DO l = 1, nrens
              propunit(s_u,s_k,nsmpl) = propunit(s_u,s_k,nsmpl) +log(propunit(s_u,s_k,l))
           END DO
           propunit(s_u,s_k,nsmpl) = exp(propunit(s_u,s_k,nsmpl) /dble(nrens))
        END IF
     ELSE
        WRITE(*,'(1A,1I2,1A)') 'error: active parameter (code: ',s_k,') for ENKF not supported yet!'
        STOP
     END IF
  END DO
  !
  !         -> compute mean of state-variables
  !         - HEAD -
  IF (head_active) THEN
     CALL set_dval(i0*j0*k0,0.0d0,head(1,1,1,nsmpl))
     DO l = 1, nrens
        CALL daxpy(i0*j0*k0,1.0d0/dble(nrens),head(1,1,1,l),1,head(1,1,1,nsmpl),1)
     ENDDO
  ELSE
     !           when not active copy any value to the mean
     CALL dcopy(i0*j0*k0,head(1,1,1,1),1,head(1,1,1,nsmpl),1)
  END IF
  !         - TEMP -
  IF (temp_active) THEN
     CALL set_dval(i0*j0*k0,0.0d0,temp(1,1,1,nsmpl))
     DO l = 1, nrens
        CALL daxpy(i0*j0*k0,1.0d0/dble(nrens),temp(1,1,1,l),1,temp(1,1,1,nsmpl),1)
     ENDDO
  ELSE
     !           when not active copy any value to the mean
     CALL dcopy(i0*j0*k0,temp(1,1,1,1),1,temp(1,1,1,nsmpl),1)
  END IF
  !         - CONC -
  DO l1 = 1, ntrac
     IF (chem_active) THEN
        CALL set_dval(i0*j0*k0,0.0d0,conc(1,1,1,l1,nsmpl))
        DO l = 1, nrens
           CALL daxpy(i0*j0*k0,1.0d0/dble(nrens),conc(1,1,1,l1,l),1,conc(1,1,1,l1,nsmpl),1)
        ENDDO
     ELSE
        !             when not active copy any value to the mean
        CALL dcopy(i0*j0*k0,conc(1,1,1,l1,1),1,conc(1,1,1,l1,nsmpl),1)
     END IF
  ENDDO
  !         - PRES -
  IF (pres_active) THEN
     CALL set_dval(i0*j0*k0,0.0d0,pres(1,1,1,nsmpl))
     DO l = 1, nrens
        CALL daxpy(i0*j0*k0,1.0d0/dble(nrens),pres(1,1,1,l),1,pres(1,1,1,nsmpl),1)
     ENDDO
  ELSE
     !           when not active copy any value to the mean
     CALL dcopy(i0*j0*k0,pres(1,1,1,1),1,pres(1,1,1,nsmpl),1)
  END IF
  !         write the average model to file, index number -3 for "mean"-prefix


  !--------------------------------------------------
  project_org = project
  project = trim(ssample_outdir)//trim(project_org)

  CALL forward_write(-3,nsmpl)

  project = project_org
  !


end subroutine enkf_compute_mean
