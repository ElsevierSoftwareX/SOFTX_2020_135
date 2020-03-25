! MIT License
!
! Copyright (c) 2019 Geir Evensen
!
! Permission is hereby granted, free of charge, to any person
! obtaining a copy of this software and associated documentation
! files (the "Software"), to deal in the Software without
! restriction, including without limitation the rights to use, copy,
! modify, merge, publish, distribute, sublicense, and/or sell copies
! of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be
! included in all copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
! EXPRESS OR IMPLIED, INCLUDING BUT NOT NO EVENT SHALL THE AUTHORS OR
! COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
! ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE
! OR OTHER DEALINGS IN THE SOFTWARE.

      SUBROUTINE analysis(a,r,e,s,d,innov,ndim,nrens,nrobs, &
          truncation,mode,update_randrot, assidamp)
!-----------
! EnKF code from Geir Evensen   http://enkf.nersc.no/
!-----------
! Computes the analysed ensemble for A using the EnKF or square root schemes.

        USE mod_anafunc
        USE m_multa
        use mod_enkf, only:&
             cov_loc_switch,&
             cov_loc_assim_id,&
             rhos,&
             hybrid_switch,&
             pb_r,&
             pb_reps,&
             hybrid_alpha,&
             assimstp_switch,&
             ana_mat_out
             

        IMPLICIT NONE
        integer, INTENT (IN) :: & ! dimension of model 
          ndim
        integer, INTENT (IN) :: & ! number of ensemble 
          nrens
        integer, INTENT (IN) :: & 
          nrobs
! number of observati
        double precision, INTENT (INOUT), dimension (ndim,nrens) :: a ! ensemble matrix  
        double precision, INTENT (INOUT), dimension (nrobs,nrobs) :: r ! matrix holding R 
        double precision, INTENT (IN), dimension (nrobs,nrens) :: d ! matrix holding pe
        double precision, INTENT (IN), dimension (nrobs,nrens) :: e ! matrix holding pe
        double precision, INTENT (IN), dimension (nrobs,nrens) :: s ! matrix holding HA
        double precision, INTENT (IN), dimension (nrobs)  :: innov
! vector holding d-
        double precision, INTENT (IN) :: & 
          truncation
! The ratio of vari
        integer, INTENT (IN) :: &                                            ! Second integer is pseudo in
          mode
                                           !  1=eigen value pseudo inver
                                           !  2=SVD subspace pseudo inve
                                           !  3=SVD subspace pseudo inve

! first integer means
        LOGICAL, INTENT (IN) :: &                                            ! updates when using local an
          update_randrot
                                           ! points need to use the same

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Normally true; fals
        double precision, dimension (nrens,nrens) :: x5
        double precision :: assidamp
        integer ::  i, nrmin, iblkmax
        LOGICAL lreps

        integer :: n, m, j


        double precision, ALLOCATABLE, dimension (:) :: eig
        double precision, ALLOCATABLE, dimension (:,:) :: w
        double precision, ALLOCATABLE, dimension (:,:) :: x2
        double precision, ALLOCATABLE, dimension (:,:) :: x3
        double precision, ALLOCATABLE, dimension (:,:) :: reps


1001    FORMAT (I4,2X,10(D13.4,2X))

        lreps = .FALSE.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Pseudo inversion of C=SS' +(N-1)*R
        if(assimstp_switch) then
           WRITE(12,*)
           WRITE(12,*) &
                'Analysis: Compute Pseudo inversion of C=SS +(N-1)*R '
        end if
        IF (nrobs==1) THEN
          nrmin = 1
          ALLOCATE(w(1,1))
          ALLOCATE(eig(1))
          eig(1) = dot_product(s(1,:),s(1,:)) + &
            dble(real(nrens-1))*r(1,1)
          eig(1) = 1.0D0/eig(1)
          w(1,1) = 1.0D0
        ELSE
          SELECT CASE (mode)

          CASE (11,21)
            nrmin = nrobs
!        Evaluate R= S*S` + (nrens-1)*R
            if(assimstp_switch) then
               WRITE(12,*)
               WRITE(12,*) 'Analysis: Case 11,21'
               WRITE(12,*) '           R= S*S + (nrens-1)*R'
               IF (nrobs<12 .AND. nrens<6) THEN
                  WRITE(12,*) '*S*'
                  DO m = 1, nrobs
                     WRITE(12,1001) m, (s(m,j),j=1,nrens)
                  END DO
                  WRITE(12,*) '*R*'
                  DO m = 1, nrobs
                     WRITE(12,1001) m, (r(m,j),j=1,nrobs)
                  END DO
               END IF
               WRITE(12,*)
            end if
            if (hybrid_switch) then
               ! pb_r = (1-alpha)*ss' + alpha*pb_r
               CALL dgemm('n','t',nrobs,nrobs,nrens,1.0D0-hybrid_alpha,s,nrobs,s, &
                    nrobs,hybrid_alpha,pb_r,nrobs)
               ! r = (1-alpha)*ss' + alpha*pb_r + (N-1)*r
               r = pb_r + dble(nrens-1)*r
            else
               ! r = ss' + (N-1)*r
               CALL dgemm('n','t',nrobs,nrobs,nrens,1.0D0,s,nrobs,s, &
                    nrobs,dble(nrens-1),r,nrobs)
            end if
            IF (nrobs<12 .and. assimstp_switch) THEN
              WRITE(12,*) 'R=S*S + (nrens -1)*R'
              DO m = 1, nrobs
                WRITE(12,1001) m, (r(m,j),j=1,nrobs)
              END DO
            END IF
!        Compute eigenvalue decomposition of R -> W*eig*W`
            ALLOCATE(w(nrobs,nrobs))
            ALLOCATE(eig(nrobs))
            CALL eigc(r,nrobs,w,eig)
            CALL eigsign(eig,nrobs,truncation)
            if(assimstp_switch) WRITE(12,*) 'Analysis: finished pseudo inversion'

          CASE (12,22)
            nrmin = min(nrobs,nrens)
            ALLOCATE(w(nrobs,nrmin))
            ALLOCATE(eig(nrmin))
            if(assimstp_switch) then
               WRITE(12,*)
               WRITE(12,*) 'Analysis: Case 12,22'
               WRITE(12,*) '           R= S*S + (nrens-1)*R  lowrank'
               IF (nrobs<12 .AND. nrens<6) THEN
                  WRITE(12,*) '*S*'
                  DO m = 1, nrobs
                     WRITE(12,1001) m, (s(m,j),j=1,nrens)
                  END DO
                  WRITE(12,*) '*R*'
                  DO m = 1, nrobs
                     WRITE(12,1001) m, (r(m,j),j=1,nrobs)
                  END DO
               END IF
               WRITE(12,*)
            end if
            CALL lowrankcinv(s,r,nrobs,nrens,nrmin,w,eig,truncation)

          CASE (13,23)
            nrmin = min(nrobs,nrens)
            ALLOCATE(w(nrobs,nrmin))
            ALLOCATE(eig(nrmin))
            if(assimstp_switch) then
               WRITE(12,*)
               WRITE(12,*) 'Analysis: Case 13,23'
               WRITE(12,*) '           R= S*S + (nrens-1)*R  lowrank'
               IF (nrobs<12 .AND. nrens<6) THEN
                  WRITE(12,*) '*S*'
                  DO m = 1, nrobs
                     WRITE(12,1001) m, (s(m,j),j=1,nrens)
                  END DO
                  WRITE(12,*) '*R*'
                  DO m = 1, nrobs
                     WRITE(12,1001) m, (r(m,j),j=1,nrobs)
                  END DO
               END IF
               WRITE(12,*)
            end if
            CALL lowranke(s,e,nrobs,nrens,nrmin,w,eig,truncation)

          CASE DEFAULT
            PRINT *, 'Analysis: Unknown mode: ', mode
            STOP
          END SELECT

          if(assimstp_switch) then
             WRITE(12,*)
             WRITE(12,*) '*W*'
             DO m = 1, nrobs
                WRITE(12,1001) m, (w(m,j),j=1,nrobs)
             END DO
             WRITE(12,*)
             WRITE(12,*) ' Eigenwerte'
             WRITE(12,'(5(d13.4,2x))') (eig(m),m=1,nrobs)
             WRITE(12,*)
          end if
          
        END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Generation of X5 (or representers in EnKF case with few measurements)
        if(assimstp_switch) then
           WRITE(12,*)
           WRITE(12,*) 'Analysis: Generation of X5:'
        end if
        SELECT CASE (mode)
        CASE (11,12,13)
          if(assimstp_switch) WRITE(12,*) 'Case 11,12,13'
          ALLOCATE(x3(nrobs,nrens))
          IF (nrobs>1) THEN
            CALL genx3(nrens,nrobs,nrmin,eig,w,d,x3)
          ELSE
            x3 = d*eig(1)
          END IF
           ! In the case of covariance localisation, representer method
           ! has to be used
           if(cov_loc_switch) then
              if(assimstp_switch) then
                 WRITE(12,*) 'Analysis: Representer approach is used'
                 write(unit = 12, fmt = *) 'WITH: LOCALISATION'
              end if
              lreps = .TRUE.
              ALLOCATE(reps(ndim,nrobs))
              !        Reps=matmul(A,transpose(S))
              CALL dgemm('n','t',ndim,nrobs,nrens,1.0D0,a,ndim,s,nrobs, &
                   0.0D0,reps,ndim)              
              !SCHUR product with rhos
              do m = 1, nrobs
                 do n = 1, ndim
                    reps(n,m) = rhos(n,m,cov_loc_assim_id)*reps(n,m)
                 end do
              end do
           ! Hybrid EnKF: Add stationary covariance to representers
           else if (hybrid_switch) THEN
              if(assimstp_switch) then
                 WRITE(12,*) 'Analysis: Representer approach is used'
                 write(unit = 12, fmt = *) 'WITH: HYBRID EnKF'
              end if
              lreps = .TRUE.
              ALLOCATE(reps(ndim,nrobs))
              ! (1-alpha)*reps + alpha*pb_reps = reps_hybrid
              CALL dgemm('n','t',ndim,nrobs,nrens,1.0D0-hybrid_alpha,a,ndim,s,nrobs, &
                   hybrid_alpha,pb_reps,ndim)
              reps = pb_reps
           ! No localisation or hybrid: proceed normally according to how many
           ! observations exits with repect to the number of states
           else
              IF (2_8*ndim*nrobs<1_8*nrens*(nrobs+ndim)) THEN
                 !        Code for few observations ( m<nN/(2n-N) )
                 if(assimstp_switch) WRITE(12,*) 'Analysis: Representer approach is used'
                 lreps = .TRUE.
                 ALLOCATE(reps(ndim,nrobs))
                 !        Reps=matmul(A,transpose(S))
                 CALL dgemm('n','t',ndim,nrobs,nrens,1.0D0,a,ndim,s,nrobs, &
                      0.0D0,reps,ndim)
              ELSE
                 if(assimstp_switch) WRITE(12,*) 'Analysis: X5 approach is used'
                 !        X5=matmul(transpose(S),X3)
                 CALL dgemm('t','n',nrens,nrens,nrobs,1.0D0*assidamp,s, &
                      nrobs,x3,nrobs,0.0D0,x5,nrens)
                 DO i = 1, nrens
                    x5(i,i) = x5(i,i) + 1.0D0
                 END DO
              END IF
           end if


        CASE (21,22,23)
          if(assimstp_switch) WRITE(12,*) 'Case 21,22,23'
! Mean part of X5
          CALL meanx5(nrens,nrobs,nrmin,s,w,eig,innov,x5)

! Generating X2
          ALLOCATE(x2(nrmin,nrens))
          CALL genx2(nrens,nrobs,nrmin,s,w,eig,x2)

! Generating X5 matrix
          CALL x5sqrt(x2,nrobs,nrens,nrmin,x5,update_randrot,mode)

        CASE DEFAULT
          if(assimstp_switch) WRITE(12,*) 'Analysis: Unknown flag for mode: ', mode
          STOP
        END SELECT


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Generation of inflation
!   call inflationTEST(X5,nrens)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Final ensemble update
        if(assimstp_switch) then
           WRITE(12,*)
           WRITE(12,*) '      Analysis: ensemble update:'
        end if
        IF (lreps) THEN
!     A=A+matmul(Reps,X3)
          if(assimstp_switch) WRITE(12,*) 'using A=A+matmul(Reps,X3)'
          CALL dgemm('n','n',ndim,nrens,nrobs,1.0D0*assidamp,reps, &
            ndim,x3,nrobs,1.0D0,a,ndim)
          if(ana_mat_out) CALL dumpx3(x3,s,nrobs,nrens)
        ELSE
          if(assimstp_switch) WRITE(12,*) 'using sub multa with A=A*X5'
          iblkmax = min(ndim,200)
          CALL multa(a,x5,ndim,nrens,iblkmax)
          if(ana_mat_out) CALL dumpx5(x5,nrens)
        END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF (allocated(x2)) DEALLOCATE(x2)
        IF (allocated(x3)) DEALLOCATE(x3)
        IF (allocated(eig)) DEALLOCATE(eig)
        IF (allocated(w)) DEALLOCATE(w)
        IF (allocated(reps)) DEALLOCATE(reps)
      END SUBROUTINE
