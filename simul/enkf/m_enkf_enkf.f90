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

      MODULE m_enkf_enkf

      CONTAINS
        SUBROUTINE enkf(mem,s,nx,nrens,obs,obsvar,obspos,nrobs,irobs, &
            fixsamp,mode_analysis,truncation,covmodel,clh,clv,rexact, &
            assidamp,lstate)
!-----------
! after EnKF code from Geir Evensen   http://enkf.nersc.no/
!-----------
          USE m_obs_pert_enkf
          use mod_genrlc
          use mod_enkf, only:&
               assimstp_switch
          IMPLICIT NONE
          integer, INTENT (IN) :: nx
          integer, INTENT (IN) :: nrens
          integer, INTENT (IN) :: nrobs, irobs
          integer, INTENT (IN) :: mode_analysis
          integer, INTENT (IN) :: lstate

          double precision, INTENT (INOUT) :: mem(nx,nrens)
          double precision, INTENT (INOUT) :: s(nrobs,nrens)
          double precision, INTENT (IN) :: obs(nrobs)
          double precision, INTENT (IN) :: obsvar(nrobs)
          integer, INTENT (IN) :: obspos(nrobs)
          LOGICAL, INTENT (IN) :: fixsamp
          LOGICAL, INTENT (IN) :: rexact
          character (len=100), INTENT (IN) :: covmodel
          double precision, INTENT (IN) :: truncation

          double precision, ALLOCATABLE :: r(:,:)
          double precision, ALLOCATABLE :: e(:,:)
          double precision, ALLOCATABLE :: d(:,:)
!       double precision, allocatable :: S(:,:)
          double precision, ALLOCATABLE :: means(:)
          double precision, ALLOCATABLE :: innovation(:)

          double precision lh, lv, control
          double precision clh, clv, assidamp

          integer :: iens, m, i, j
          LOGICAL :: update_randrot = .TRUE.

!-------------------------------------
!          CHARACTER*3 afile
!          integer  ifile
!
!          ifile = irobs + 100
!1000      FORMAT (I3)
!          WRITE(afile,1000) ifile
!          afile = trim(afile)
!          OPEN(unit=12,file=senkf_outdir//'assimstp'//afile)
!---------------------------------------

          ALLOCATE(e(nrobs,nrens))
          ALLOCATE(d(nrobs,nrens))
!       allocate(S(nrobs,nrens))
          ALLOCATE(means(nrobs))

!======================================
! horizontal and vertical correlation length (should be given in input)
!        clh=500.d0
!        clv=10.d0

1001      FORMAT (5(E11.4,2X))
1002      FORMAT (I4,2X,6(E11.4,2X))
          IF (nrens<6 .and. assimstp_switch) THEN
            WRITE(12,*) &
              '   e1          e2          e3          e4          e5'
            DO j = 1, nx
              WRITE(12,1001) (mem(j,i),i=1,nrens)
            END DO
            WRITE(12,*)
            WRITE(12,*) 'obsvar = ', obsvar
            WRITE(12,*) ' nrobs  obs      obspos'
            WRITE(12,'(i4,2x,e11.4,2x,i4)') (j,obs(j),obspos(j),j=1, &
              nrobs)
            WRITE(12,*)
            WRITE(12,*) 'modeana = ', mode_analysis, 'truncation = ', &
              truncation, 'cov = ', trim(covmodel), ' Rexact = ', &
              rexact
            WRITE(12,*)
          END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Observe ensemble to construct the matrix S=HA
!        do iens =1, nrens
!          do m =1, nrobs
!            S(m,iens)      =  mem(obspos(m),iens)
!          enddo
!        enddo

!======================================
          if(assimstp_switch) then
             WRITE(12,*) 'construct S=HA'
             IF (nrens<6) THEN
                DO m = 1, nrobs
                   WRITE(12,1002) m, obs(m), (s(m,iens),iens=1,nrens)
                END DO
             END IF
             WRITE(12,*)
          end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Construct observation perturbations E
!   here only covmodel='diagonal' is possible
!  (for covmodel='gaussian' a field generator would be needed)
          CALL obs_pert(e,nrens,nrobs,fixsamp)
! Introduce correct variances
          DO iens = 1, nrens
            DO m = 1, nrobs
              e(m,iens) = dsqrt(obsvar(m))*e(m,iens)
            END DO
          END DO

!==========================================
          if(assimstp_switch) then
             WRITE(12,*) 'construct E'
             IF (nrens<6) THEN
                DO m = 1, nrobs
                   WRITE(12,1002) m, obs(m), (e(m,iens),iens=1,nrens)
                END DO
             END IF
             WRITE(12,*)
          end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Construct ensemble of measurements D=d+E
          DO iens = 1, nrens
            DO m = 1, nrobs
              d(m,iens) = obs(m) + e(m,iens)
            END DO
          END DO

!==============================================
          if(assimstp_switch) then
             WRITE(12,*) 'construct D'
             IF (nrens<6) THEN
                DO m = 1, nrobs
                   WRITE(12,1002) m, obs(m), (d(m,iens),iens=1,nrens)
                END DO
             END IF
             WRITE(12,*)
          end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute innovation D'=D-HA
          d = d - s

!=================================================
          if(assimstp_switch) then
             WRITE(12,*) 'construct innovation D=D-HA'
             IF (nrens<6) THEN
                DO m = 1, nrobs
                   WRITE(12,1002) m, obs(m), (d(m,iens),iens=1,nrens)
                END DO
             END IF
             WRITE(12,*)
          end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute mean(HA)
          means = 0.0D0
          DO iens = 1, nrens
            DO m = 1, nrobs
              means(m) = means(m) + s(m,iens)
            END DO
          END DO
          means = (1.0D0/dble(float(nrens)))*means
          if(assimstp_switch) then
             WRITE(12,*) 'meanS'
             WRITE(12,*) means
          end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute HA'=HA-mean(HA)
          DO iens = 1, nrens
            s(:,iens) = s(:,iens) - means(:)
          END DO

!=====================================================
          if(assimstp_switch) then
             WRITE(12,*) ' Compute HA =HA-mean(HA) and save on S'
             IF (nrens<6) THEN
                DO m = 1, nrobs
                   WRITE(12,1002) m, obs(m), (s(m,iens),iens=1,nrens)
                END DO
             END IF
             WRITE(12,*)
          end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute R
          ALLOCATE(r(nrobs,nrobs))
          r = 0.0D0
          IF (rexact) THEN
             if(assimstp_switch) WRITE(12,*) ' Exact R using covariance model: ', &
                  trim(covmodel)
            SELECT CASE (trim(covmodel))
            CASE ('diagonal')
              DO m = 1, nrobs
                r(m,m) = obsvar(m)
              END DO
            CASE ('gaussian')
              DO i = 1, nrobs
                DO j = 1, nrobs
                  CALL length(lh,lv,obspos(i),obspos(j),control,lstate)
                  if(assimstp_switch) WRITE(12,*) 'obspos and distance ', obspos(i), &
                       obspos(j), lh, lv
                  r(i,j) = obsvar(j)*control*dexp(-lh**2/clh**2)* &
                    dexp(-lv**2/clv**2)
!original     &                      dexp(-dble(real(i-j))**2/20.0d0**2)
                END DO
              END DO
            CASE DEFAULT
              if(assimstp_switch) WRITE(12,*) 'Covmodel is invalid : ', trim(covmodel)
            END SELECT
          ELSE
            if(assimstp_switch) WRITE(12,*) ' enkf: Lowrank R using covariance model: ', &
                 trim(covmodel)
            r = matmul(e,transpose(e))/dble(float(nrens))
          END IF
!=====================================================
          if(assimstp_switch) then
             WRITE(12,*) ' Compute R ', trim(covmodel)
             IF (nrobs<10) THEN
                DO m = 1, nrobs
                   WRITE(12,'(i3,2x,10(e11.4,2x))') m, (r(m,j),j=1,nrobs)
                END DO
             END IF
             WRITE(12,*)
          end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute innovation
          ALLOCATE(innovation(nrobs))
          DO m = 1, nrobs
            innovation(m) = obs(m) - means(m)
          END DO
!=====================================================
          if(assimstp_switch) then
             WRITE(12,*) ' Compute innovation'
             WRITE(12,'(10(e11.4,2x))') (innovation(j),j=1,nrobs)
             WRITE(12,*)
          end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          if(assimstp_switch) WRITE(12,*) '  enkf: calling analysis with mode: ', &
               mode_analysis

          CALL analysis(mem,r,e,s,d,innovation,nx,nrens,nrobs, &
            truncation,mode_analysis,update_randrot, assidamp)

          DEALLOCATE(innovation)
          DEALLOCATE(r)

!          CLOSE(12)
        END SUBROUTINE


        SUBROUTINE length(lh,lv,l1,l2,control,ijk)
!     sub to determine horizontal and vertical distances between two points given in grid indices

          use arrays
        use mod_genrl
        use mod_enkf, only:&
             enkf_log_out
          IMPLICIT NONE

          integer :: ijk
          integer :: i1, j1, k1, i2, j2, k2, l1, l2, mm1, mm2, m1, m2, &
            h1, h2
          double precision lh, lv, x1, y1, z1, x2, y2, z2, control

          control = 1.D0
!          ijk = i0*j0*k0
!    Test if variables are of the same type
          mm1 = mod(l1,ijk)
          h1 = int(l1/ijk)
          IF (mm1/=0) h1 = h1 + 1
          mm2 = mod(l2,ijk)
          h2 = int(l2/ijk)
          IF (mm2/=0) h2 = h2 + 1

          IF (h1==h2) THEN
!     Determine indices of the two observations of same type (h1=h2)
            IF (mm1==0) mm1 = ijk
            m1 = mod(mm1,i0*j0)
            m2 = mod(m1,i0)
            k1 = int(mm1/(i0*j0))
            IF (m1/=0) k1 = k1 + 1
            j1 = int((mm1-(k1-1)*(i0*j0))/i0)
            IF (m2/=0) j1 = j1 + 1
            i1 = mm1 - (k1-1)*i0*j0 - (j1-1)*i0
            if (enkf_log_out) WRITE(37,*) l1, mm1, i1, j1, k1

            IF (mm2==0) mm2 = ijk
            m1 = mod(mm2,i0*j0)
            m2 = mod(m1,i0)
            k2 = int(mm2/(i0*j0))
            IF (m1/=0) k2 = k2 + 1
            j2 = int((mm2-(k2-1)*(i0*j0))/i0)
            IF (m2/=0) j2 = j2 + 1
            i2 = mm2 - (k2-1)*i0*j0 - (j2-1)*i0
            if (enkf_log_out) WRITE(37,*) l1, mm2, i2, j2, k2
            if (enkf_log_out) WRITE(37,*)
            x1 = delxa(i1)
            y1 = delya(j1)
            z1 = delza(k1)
            x2 = delxa(i2)
            y2 = delya(j2)
            z2 = delza(k2)

!     Determine horizontal and vertical distances
            lh = dsqrt((x1-x2)**2+(y1-y2)**2)
            lv = dabs(z1-z2)
          ELSE
!     Set distances to zero for observations of different type (h1.ne.h2)
            lh = 0.D0
            lv = 0.D0
            control = 0.D0
          END IF
          RETURN
        END SUBROUTINE


      END MODULE m_enkf_enkf
