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

      MODULE mod_anafunc

      CONTAINS

        SUBROUTINE lowranke(s,e,nrobs,nrens,nrmin,w,eig,truncation)
!-----------
! EnKF code after original by Geir Evensen   http://enkf.nersc.no/
!-----------
          use mod_enkf, only:&
               assimstp_switch
          
          IMPLICIT NONE
          integer, INTENT (IN) :: nrobs
          integer, INTENT (IN) :: nrens
          integer, INTENT (IN) :: nrmin
          double precision, INTENT (IN), dimension (nrobs,nrens) :: s
          double precision, INTENT (IN), dimension (nrobs,nrens) :: e
          double precision, INTENT (OUT), dimension (nrobs,nrmin) :: w
          double precision, INTENT (OUT), dimension (nrmin) :: eig
          double precision, INTENT (IN) :: truncation

          double precision, dimension (nrobs,nrmin) :: u0
          double precision, dimension (nrmin) :: sig0
          double precision, dimension (nrmin,nrens) ::  x0
          integer :: i, j

          double precision, dimension (nrmin,nrmin) ::  u1
          double precision, dimension (1,1) :: vt1
          double precision, ALLOCATABLE, dimension (:) :: work
          integer :: lwork
          integer :: ierr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute SVD of S=HA`  ->  U0, sig0
          CALL svds(s,nrobs,nrens,nrmin,u0,sig0,truncation)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute X0=sig0^{*T} U0^T E

! X0= U0^T R
          CALL dgemm('t','n',nrmin,nrens,nrobs,1.0,u0,nrobs,e,nrobs, &
            0.0,x0,nrmin)

          DO j = 1, nrens
            DO i = 1, nrmin
              x0(i,j) = sig0(i)*x0(i,j)
            END DO
          END DO


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute singular value decomposition  of X0(nrmin,nrens)
          lwork = 2*max(3*nrens+nrobs,5*nrens)
          ALLOCATE(work(lwork))
          eig = 0.0D0

          CALL dgesvd('S','N',nrmin,nrens,x0,nrmin,eig,u1,nrmin,vt1,1, &
            work,lwork,ierr)
          DEALLOCATE(work)
          IF (ierr/=0) THEN
            if(assimstp_switch) WRITE(12,*) &
              'mod_anafunc (lowrankE): ierr from call dgesvd 1= ', &
              ierr
            STOP
          END IF

          DO i = 1, nrmin
            eig(i) = 1.0D0/(1.0D0+eig(i)**2)
          END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! W = U0 * sig0^{-1} * U1
          DO j = 1, nrmin
            DO i = 1, nrmin
              u1(i,j) = sig0(i)*u1(i,j)
            END DO
          END DO

          CALL dgemm('n','n',nrobs,nrmin,nrmin,1.0,u0,nrobs,u1,nrmin, &
            0.0,w,nrobs)
        END SUBROUTINE


        SUBROUTINE eigc(r,nrobs,z,eig)
! Compute eigenvalue decomposition of R -> Z*eig*Z`
          integer, INTENT (IN) :: nrobs
          double precision, INTENT (IN), dimension (nrobs,nrobs) :: r
          double precision, INTENT (OUT), dimension (nrobs,nrobs) :: z
          double precision, INTENT (OUT), dimension (nrobs) :: eig

#ifdef IBM
          double precision, ALLOCATABLE, dimension (:) :: ap
          integer ::  k
#endif
          double precision, dimension (nrobs,nrobs) :: rr

          double precision, dimension (8*nrobs) :: fwork
          integer, dimension (5*nrobs) :: iwork
          integer, dimension (nrobs) ::  ifail
          double precision abstol, ddum
          integer ::  idum, neig, ierr
          double precision, EXTERNAL :: dlamch

          idum = 1

#ifdef IBM
! Upper packed storage as in ESSL manual
          ALLOCATE(ap(nrobs*(nrobs+1)/2))
          k = 0
          DO j = 1, nrobs
            DO i = 1, j
              k = k + 1
              ap(k) = r(i,j)
            END DO
          END DO
          CALL dspev(21,ap,eig,z,nrobs,nrobs,fwork,2*nrobs)
          DEALLOCATE(ap)
#else
          abstol = 2.0D0*dlamch('S')
          rr = r
          CALL dsyevx('V','A','U',nrobs,rr,nrobs,ddum,ddum,idum,idum, &
            abstol,neig,eig,z,nrobs,fwork,8*nrobs,iwork,ifail,ierr)
#endif


        END SUBROUTINE



        SUBROUTINE eigsign(eig,nrobs,truncation)
! Returns the inverse of the truncated eigenvalue spectrum
        use mod_simul
        use mod_enkf, only:&
             assimstp_switch,&
             ana_mat_out
          IMPLICIT NONE
          integer, INTENT (IN) :: nrobs
          double precision, INTENT (INOUT), dimension (nrobs) :: eig
          double precision, INTENT (IN) :: truncation

          integer :: i, nrsigma
          double precision sigsum, sigsum1
          LOGICAL ex

          INQUIRE (file=senkf_outdir//'eigenvalues.dat',exist=ex)
          if(ana_mat_out) then
             IF (ex) THEN
                OPEN(10,file=senkf_outdir//'eigenvalues.dat', &
                     position='append')
                WRITE(10,'(a,i5,a)') ' ZONE  F=POINT, I=', nrobs, &
                     ' J=1 K=1'
                DO i = 1, nrobs
                   WRITE(10,'(i3,g13.5)') i, eig(nrobs-i+1)
                END DO
                CLOSE(10)
             ELSE
                OPEN(10,file=senkf_outdir//'eigenvalues.dat')
                WRITE(10,*) 'TITLE = "Eigenvalues of C"'
                WRITE(10,*) 'VARIABLES = "obs" "eigenvalues"'
                WRITE(10,'(a,i5,a)') ' ZONE  F=POINT, I=', nrobs, &
                     ' J=1 K=1'
                DO i = 1, nrobs
                   WRITE(10,'(i3,g13.5)') i, eig(nrobs-i+1)
                END DO
                CLOSE(10)
             END IF
          end if

! Significant eigenvalues
          sigsum = sum(eig(1:nrobs))
          sigsum1 = 0.0D0
          nrsigma = 0
          DO i = nrobs, 1, -1
!      print '(a,i5,g13.5)','Eigen values: ',i,eig(i)
            IF (sigsum1/sigsum<truncation) THEN
              nrsigma = nrsigma + 1
              sigsum1 = sigsum1 + eig(i)
              eig(i) = 1.0D0/eig(i)
            ELSE
              eig(1:i) = 0.0D0
              EXIT
            END IF
          END DO
          if(assimstp_switch) WRITE(12,'(2(a,i5))') &
               ' analysis: Number of dominant eigenvalues: ', nrsigma, &
               ' of ', nrobs
          if(assimstp_switch) WRITE(12,'(2(a,g13.4),a)') &
               ' analysis: Share (and truncation)  : ', sigsum1/sigsum, &
               ' (', truncation, ')'


        END SUBROUTINE



        SUBROUTINE genx2(nrens,nrobs,idim,s,w,eig,x2)
!! Generate X2= (I+eig)^{-0.5} * W^T * S
          IMPLICIT NONE
          integer, INTENT (IN) :: nrens
          integer, INTENT (IN) :: nrobs
          integer, INTENT (IN) :: idim
!! idim=nrobs for A4 and nrmin for A5
          double precision, INTENT (IN), dimension (idim,nrens) :: w
          double precision, INTENT (IN), dimension (nrobs,nrens) :: s
          double precision, INTENT (IN), dimension (idim) :: eig
          double precision, INTENT (OUT), dimension (idim,nrens) :: x2
          ! double precision, intent(in), dimension (nrobs,nrobs) :: w
          ! double precision, intent(in), dimension (nrobs,nrens) :: s
          ! double precision, intent(in), dimension (nrobs) :: eig
          ! double precision, intent(out), dimension (nrobs,nrens) :: x2
          integer :: i, j

!! Johannes: w seems to get wrong dimensions (idim,nrens) -
!! (nrobs,nrobs) expected
          CALL dgemm('t','n',idim,nrens,nrobs,1.0d0,w,nrobs,s,nrobs,0.0d0, &
            x2,idim)

          DO j = 1, nrens
            DO i = 1, idim
              x2(i,j) = dsqrt(eig(i))*x2(i,j)
            END DO
          END DO

        END SUBROUTINE



        SUBROUTINE genx3(nrens,nrobs,nrmin,eig,w,d,x3)
          IMPLICIT NONE
          integer, INTENT (IN) :: nrens
          integer, INTENT (IN) :: nrobs
          integer, INTENT (IN) :: nrmin
          double precision, INTENT (IN), dimension (nrmin) :: eig
          double precision, INTENT (IN), dimension (nrobs,nrmin) :: w
          double precision, INTENT (IN), dimension (nrobs,nrens) :: d
          double precision, INTENT (OUT), dimension (nrobs,nrmin) :: x3

          double precision, dimension (nrmin,nrobs) :: x1
          double precision, dimension (nrmin,nrens) :: x2
          integer :: i, j

          DO i = 1, nrmin
            DO j = 1, nrobs
              x1(i,j) = eig(i)*w(j,i)
            END DO
          END DO

!       X2=matmul(X1,D)
          CALL dgemm('n','n',nrmin,nrens,nrobs,1.D0,x1,nrmin,d,nrobs, &
            0.D0,x2,nrmin)

!     X3=matmul(W,X2)
          CALL dgemm('n','n',nrobs,nrens,nrmin,1.D0,w,nrobs,x2,nrmin, &
            0.D0,x3,nrobs)

        END SUBROUTINE



        SUBROUTINE meanx5(nrens,nrobs,nrmin,s,w,eig,innov,x5)
          IMPLICIT NONE
          integer, INTENT (IN) :: nrens
          integer, INTENT (IN) :: nrobs
          integer, INTENT (IN) :: nrmin
          double precision, INTENT (IN), dimension (nrmin,nrmin) :: w
          double precision, INTENT (IN), dimension (nrobs,nrens) :: s
          double precision, INTENT (IN), dimension (nrmin) :: eig
          double precision, INTENT (IN), dimension (nrobs) :: innov
          double precision, INTENT (OUT), dimension (nrens,nrens) :: x5

          double precision, dimension (nrmin) :: y1
          double precision, dimension (nrmin) :: y2
          double precision, dimension (nrobs) :: y3
          double precision, dimension (nrens) :: y4
          integer :: i

          IF (nrobs==1) THEN
            y1(1) = w(1,1)*innov(1)
            y2(1) = eig(1)*y1(1)
            y3(1) = w(1,1)*y2(1)
            y4(:) = y3(1)*s(1,:)
          ELSE
            CALL dgemv('t',nrobs,nrmin,1.D0,w,nrobs,innov,1,0.D0,y1,1)
            y2 = eig*y1
            CALL dgemv('n',nrobs,nrmin,1.D0,w,nrobs,y2,1,0.D0,y3,1)
            CALL dgemv('t',nrobs,nrens,1.D0,s,nrobs,y3,1,0.D0,y4,1)
          END IF

          DO i = 1, nrens
            x5(:,i) = y4(:)
          END DO

! X5=enN + (I - enN) X5  = enN + X5
          x5 = 1.0D0/dble(real(nrens)) + x5

        END SUBROUTINE



        SUBROUTINE x5sqrt(x2,nrobs,nrens,nrmin,x5,update_randrot,mode)
          USE m_randrot
          USE m_mean_preserving_rotation
          use mod_enkf, only:&
               assimstp_switch
          IMPLICIT NONE
          integer, INTENT (IN) :: nrobs
          integer, INTENT (IN) :: nrens
          integer, INTENT (INOUT) :: & ! note that nrmin=nrobs in a4 
            nrmin
          ! double precision, INTENT (IN), dimension (nrmin,nrens) :: x2
          ! double precision, INTENT (INOUT), dimension (nrens,nrens) :: x5
          double precision, intent(in), dimension (nrmin,nrens) :: x2
          double precision, intent(inout), dimension (nrens,nrens) :: x5
          LOGICAL, INTENT (IN) :: update_randrot
          integer, INTENT (IN) :: mode

          double precision, dimension (nrens,nrens) :: x3
          double precision, dimension (nrens,nrens) :: x33
          double precision, dimension (nrens,nrens) :: x4
          double precision, dimension (nrens,nrens) ::  ienn
          double precision, SAVE , ALLOCATABLE, dimension (:,:) :: rot

          double precision, dimension (nrmin,1) :: u
          double precision, dimension (nrmin) :: sig
          double precision, dimension (nrens,nrens) :: vt
          double precision, ALLOCATABLE, DIMENSION (:) :: work
          double precision, ALLOCATABLE, DIMENSION (:) :: isigma
          integer :: i, j, lwork, ierr

          
          if(assimstp_switch) WRITE(12,*) '      analysis (X5sqrt): update_randrot= ', &
               update_randrot
          IF (update_randrot) THEN
            IF (allocated(rot)) DEALLOCATE(rot)
            ALLOCATE(rot(nrens,nrens))
            CALL mean_preserving_rotation(rot,nrens)
          END IF

! SVD of X2
          lwork = 2*max(3*nrens+nrens,5*nrens)
          ALLOCATE(work(lwork))
          sig = 0.0D0
          CALL dgesvd('N','A',nrmin,nrens,x2,nrmin,sig,u,nrmin,vt, &
            nrens,work,lwork,ierr)
          DEALLOCATE(work)
          IF (ierr/=0) THEN
            if(assimstp_switch) WRITE(12,*) 'X5sqrt: ierr from call dgesvd = ', ierr
            STOP
          END IF

          IF (mode==21) nrmin = min(nrens,nrobs)
          ALLOCATE(isigma(nrmin))
          isigma = 1.0D0
          DO i = 1, nrmin
             IF (sig(i)>1.0 .and. assimstp_switch) WRITE(12,*) &
                  'X5sqrt: WARNING (m_X5sqrt): sig > 1', i, sig(i)
            isigma(i) = dsqrt(max(1.0D0-sig(i)**2,0.0D0))
          END DO

          DO j = 1, nrens
            x3(:,j) = vt(j,:)
          END DO

          DO j = 1, nrmin
            x3(:,j) = x3(:,j)*isigma(j)
          END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Multiply  X3* V' = (V*sqrt(I-sigma*sigma) * V' to ensure symmetric sqrt and 
! mean preserving rotation.   Sakov paper eq 13
          CALL dgemm('n','n',nrens,nrens,nrens,1.0d0,x3,nrens,vt,nrens, &
            0.0d0,x33,nrens)
! Check mean preservation X33*1_N = a* 1_N (Sakov paper eq 15)
!   do i=1,nrens
!      print *,'sum(X33)= ',i,sum(X33(i,:))
!   enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!   print '(a)','X5sqrt: sig: '
!   print '(5g11.3)',sig(1:min(nrmin,nrens))

          CALL dgemm('n','n',nrens,nrens,nrens,1.0d0,x33,nrens,rot, &
            nrens,0.0d0,x4,nrens)

          ienn = -1.0D0/dble(real(nrens))
          DO i = 1, nrens
            ienn(i,i) = ienn(i,i) + 1.0D0
          END DO

          CALL dgemm('n','n',nrens,nrens,nrens,1.D0,ienn,nrens,x4, &
            nrens,1.D0,x5,nrens)

          DEALLOCATE(isigma)

        END SUBROUTINE



        SUBROUTINE dumpx3(x3,s,nrobs,nrens)
          IMPLICIT NONE
          integer, INTENT (IN) :: nrens
          integer, INTENT (IN) :: nrobs
!  orig        double precision, INTENT (IN), dimension (nrens,nrens) :: x3
          double precision, INTENT (IN), dimension (nrobs,nrens) :: x3   ! corrected gm
          double precision, INTENT (IN), dimension (nrobs,nrens) :: s
          character (len=2) :: tag2

          tag2(1:2) = 'X3'
          OPEN(10,file='enkf_output/X3.uf',form='unformatted')
          WRITE(10) tag2, nrens, nrobs, x3, s
          CLOSE(10)

        END SUBROUTINE



        SUBROUTINE dumpx5(x5,nrens)
          IMPLICIT NONE
          integer, INTENT (IN) :: nrens
          double precision, INTENT (IN), dimension (nrens,nrens) :: x5
          integer :: j
          character (len=2) :: tag2

          tag2(1:2) = 'X5'
          OPEN(10,file='enkf_output/X5.uf',form='unformatted')
          WRITE(10) tag2, nrens, x5
          CLOSE(10)

          OPEN(10,file='X5col.dat')
          DO j = 1, nrens
            WRITE(10,'(i5,f10.4)') j, sum(x5(:,j))
          END DO
          CLOSE(10)

          OPEN(10,file='X5row.dat')
          DO j = 1, nrens
            WRITE(10,'(i5,f10.4)') j, sum(x5(j,:))/real(nrens)
          END DO
          CLOSE(10)
        END SUBROUTINE


        SUBROUTINE lowrankcinv(s,r,nrobs,nrens,nrmin,w,eig,truncation)
          IMPLICIT NONE
          integer, INTENT (IN) :: nrobs
          integer, INTENT (IN) :: nrens
          integer, INTENT (IN) :: nrmin
          double precision, INTENT (IN), dimension (nrobs,nrens) :: s
          double precision, INTENT (IN), dimension (nrobs,nrobs) :: r
          double precision, INTENT (OUT), dimension (nrobs,nrmin) :: w
          double precision, INTENT (OUT), dimension (nrmin) :: eig
          double precision, INTENT (IN) :: truncation

          double precision, dimension (nrobs,nrmin) ::  u0
          double precision, dimension (nrmin) :: sig0
          double precision, dimension (nrmin,nrmin) :: b
          double precision, dimension (nrmin,nrmin) :: z
          integer :: i, j

! Compute SVD of S=HA`  ->  U0, sig0
          CALL svds(s,nrobs,nrens,nrmin,u0,sig0,truncation)

! Compute B=sig0^{-1} U0^T R U0 sig0^{-1}
          CALL lowrankcee(b,nrmin,nrobs,nrens,r,u0,sig0)

! Compute eigenvalue decomposition  of B(nrmin,nrmin)
          CALL eigc(b,nrmin,z,eig)

!   print *,'eig:',nrmin
!   print '(6g11.3)',eig

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute inverse diagonal of (I+Lamda)
          DO i = 1, nrmin
            eig(i) = 1.0D0/(1.0D0+eig(i))
          END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! W = U0 * sig0^{-1} * Z
          DO j = 1, nrmin
            DO i = 1, nrmin
              z(i,j) = sig0(i)*z(i,j)
            END DO
          END DO

          CALL dgemm('n','n',nrobs,nrmin,nrmin,1.D0,u0,nrobs,z,nrmin, &
            0.D0,w,nrobs)

        END SUBROUTINE



        SUBROUTINE lowrankcee(b,nrmin,nrobs,nrens,r,u0,sig0)
          IMPLICIT NONE
          integer, INTENT (IN) :: nrmin
          integer, INTENT (IN) :: nrobs
          integer, INTENT (IN) :: nrens
          double precision, INTENT (INOUT), dimension (nrmin,nrmin) :: b
          double precision, INTENT (IN), dimension (nrobs,nrobs) :: r
          double precision, INTENT (IN), dimension (nrobs,nrmin) :: u0
          double precision, INTENT (IN), dimension (nrmin) :: sig0
          double precision, dimension (nrmin,nrobs) ::  x0
          integer :: i, j

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute B=sig0^{-1} U0^T R U0 sig0^{-1}

! X0= U0^T R
          CALL dgemm('t','n',nrmin,nrobs,nrobs,1.D0,u0,nrobs,r,nrobs, &
            0.D0,x0,nrmin)

! B= X0 U0
          CALL dgemm('n','n',nrmin,nrmin,nrobs,1.D0,x0,nrmin,u0,nrobs, &
            0.D0,b,nrmin)

          DO j = 1, nrmin
            DO i = 1, nrmin
              b(i,j) = sig0(i)*b(i,j)
            END DO
          END DO

          DO j = 1, nrmin
            DO i = 1, nrmin
              b(i,j) = sig0(j)*b(i,j)
            END DO
          END DO

          b = dble(real(nrens-1))*b

        END SUBROUTINE


        SUBROUTINE svds(s,nrobs,nrens,nrmin,u0,sig0,truncation)
          use mod_enkf, only:&
               assimstp_switch
          integer, INTENT (IN) :: nrobs
          integer, INTENT (IN) :: nrens
          integer, INTENT (IN) :: nrmin
          double precision, INTENT (IN), dimension (nrobs,nrens) :: s
          double precision, INTENT (OUT), dimension (nrmin) :: sig0
          double precision, INTENT (IN), dimension (nrobs,nrmin) :: u0
          double precision, INTENT (IN) :: truncation

          double precision, dimension (nrobs,nrens) :: s0
          double precision, dimension (1,1) :: vt0
          integer :: ierr
          integer :: lwork
          double precision, ALLOCATABLE, DIMENSION (:) :: work
          integer :: nrsigma, i

          double precision sigsum, sigsum1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute SVD of S=HA`  ->  U0, sig0
          lwork = 2*max(3*nrens+nrobs,5*nrens)
          ALLOCATE(work(lwork))

          s0 = s
          sig0 = 0.0D0
          CALL dgesvd('S','N',nrobs,nrens,s0,nrobs,sig0,u0,nrobs,vt0, &
            nrens,work,lwork,ierr)
          DEALLOCATE(work)
          IF (ierr/=0) THEN
            PRINT *, 'svdS: ierr from call dgesvd 0= ', ierr
            STOP
          END IF

          sigsum = 0.0D0
          DO i = 1, nrmin
            sigsum = sigsum + sig0(i)**2
          END DO

          sigsum1 = 0.0D0
! Significant eigenvalues.
          nrsigma = 0
          DO i = 1, nrmin
            IF (sigsum1/sigsum<truncation) THEN
              nrsigma = nrsigma + 1
              sigsum1 = sigsum1 + sig0(i)**2
            ELSE
              sig0(i:nrmin) = 0.0D0
              EXIT
            END IF
          END DO

          if(assimstp_switch) WRITE(12,'(a,i5,g13.5)') &
               '      analysis svdS: dominant sing. values and share ', &
               nrsigma, sigsum1/sigsum
!   write(*,'(5g11.3)')sig0

          DO i = 1, nrsigma
            sig0(i) = 1.0D0/sig0(i)
          END DO

        END SUBROUTINE

        SUBROUTINE inflationtest(x5,nrens)
          USE m_multa
          USE m_random
          use mod_enkf, only:&
               ana_mat_out
          IMPLICIT NONE
          integer, INTENT (IN) :: nrens
          double precision, INTENT (INOUT), dimension (nrens,nrens) :: x5

          integer, PARAMETER :: ndim = 100
          integer, PARAMETER :: nrnn = 250
          double precision aveverens
          double precision stdverens
          double precision, SAVE, dimension (ndim,nrnn) :: verens
          double precision, dimension (ndim,nrens) :: tmp
          double precision, dimension (nrens,nrens) :: verx
          double precision, dimension (nrens,nrens) :: xx
          double precision, dimension (ndim) :: std
          integer :: i, j
          LOGICAL, SAVE :: lfirst = .TRUE.
          integer, SAVE :: it = 0

          IF (nrnn/=nrens) THEN
            PRINT *, 'nrens not equal to nrnn', nrens, nrnn
            STOP
          END IF

          IF (lfirst) THEN
            lfirst = .FALSE.
            CALL random(verens,ndim*nrens)

! subtract mean to get ensemble of mean=0.0
            DO i = 1, ndim
              aveverens = sum(verens(i,1:nrens))/dble(real(nrens))
              DO j = 1, nrens
                verens(i,j) = verens(i,j) - aveverens
              END DO
            END DO

! compute std dev and scale ensemble so that it has variance=1.0
            DO i = 1, ndim
              stdverens = 0.0D0
              DO j = 1, nrens
                stdverens = stdverens + verens(i,j)**2
              END DO
              stdverens = dsqrt(stdverens/dble(real(nrens)))
              DO j = 1, nrens
                verens(i,j) = verens(i,j)/stdverens
              END DO
            END DO

            if(ana_mat_out) OPEN(10,file='stddev.dat')
            std(:) = 0.0D0
            DO j = 1, nrens
              DO i = 1, ndim
                std(i) = std(i) + verens(i,j)**2
              END DO
            END DO
            DO i = 1, ndim
              std(i) = dsqrt(std(i)/dble(real(nrens)))
            END DO
            stdverens = sum(std(1:ndim))/dble(real(ndim))
            if(ana_mat_out) WRITE(10,'(i5,g13.5)') it, stdverens
            if(ana_mat_out) CLOSE(10)
          END IF

          it = it + 1

          tmp = verens
          CALL multa(verens,x5,ndim,nrens,ndim)

! subtract mean from verens
          DO i = 1, ndim
            aveverens = sum(verens(i,1:nrens))/dble(real(nrens))
            DO j = 1, nrens
              verens(i,j) = verens(i,j) - aveverens
            END DO
          END DO

! compute average variance over all ndim states
          std(:) = 0.0D0
          DO j = 1, nrens
            DO i = 1, ndim
              std(i) = std(i) + verens(i,j)**2
            END DO
          END DO
          DO i = 1, ndim
            std(i) = dsqrt(std(i)/dble(real(nrens)))
          END DO
          stdverens = sum(std(1:ndim))/dble(real(ndim))

! Compute correction to ensemble perturbations
          verx = -1.0D0/dble(real(nrens))
          DO j = 1, nrens
            verx(j,j) = verx(j,j) + 1.0D0
          END DO
          DO j = 1, nrens
            DO i = 1, nrens
              verx(i,j) = verx(i,j)/stdverens
            END DO
          END DO
          DO j = 1, nrens
            DO i = 1, nrens
              verx(i,j) = verx(i,j) + 1.0D0/dble(real(nrens))
            END DO
          END DO

          CALL dgemm('n','n',nrens,nrens,nrens,1.0D0,x5,nrens,verx, &
            nrens,0.0D0,xx,nrens)
          PRINT *, '      analysis: inflation factor= ', 1.0/stdverens

          if(ana_mat_out) then
             OPEN(10,file='stddev.dat',position='append')
             WRITE(10,'(i5,g13.5)') it, stdverens
             CLOSE(10)
          end if
          RETURN


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Updating verens with inflated X5 (just for testing)
          x5 = xx
          verens = tmp
          CALL multa(verens,x5,ndim,nrens,ndim)

! subtract mean from verens
          DO i = 1, ndim
            aveverens = sum(verens(i,1:nrens))/dble(real(nrens))
            DO j = 1, nrens
              verens(i,j) = verens(i,j) - aveverens
            END DO
          END DO

! compute average variance over all ndim states
          std(:) = 0.0D0
          DO j = 1, nrens
            DO i = 1, ndim
              std(i) = std(i) + verens(i,j)**2
            END DO
          END DO
          DO i = 1, ndim
            std(i) = dsqrt(std(i)/dble(real(nrens)))
          END DO
          stdverens = sum(std(1:ndim))/dble(real(ndim))

          if(ana_mat_out) then
             OPEN(unit=10,file='stddev.dat',position='append')
             WRITE(10,'(i5,g13.5)') it, stdverens
             CLOSE(10)
          end if

        END SUBROUTINE



        SUBROUTINE inflation(x5,nrens)
          USE m_multa
          USE m_random
          use mod_enkf, only:&
               assimstp_switch
          IMPLICIT NONE
          integer, INTENT (IN) :: nrens
          double precision, INTENT (INOUT), dimension (nrens,nrens) :: x5

          integer, PARAMETER :: ndim = 100
          double precision aveverens
          double precision stdverens
          double precision, dimension (ndim,nrens) :: verens
          double precision, dimension (nrens,nrens) :: verx
          double precision, dimension (nrens,nrens) :: xx
          double precision, dimension (ndim) :: std
          integer :: i, j

          CALL random(verens,ndim*nrens)

! subtract mean to get ensemble of mean=0.0
          DO i = 1, ndim
            aveverens = sum(verens(i,1:nrens))/dble(real(nrens))
            DO j = 1, nrens
              verens(i,j) = verens(i,j) - aveverens
            END DO
          END DO

! compute std dev and scale ensemble so that it has variance=1.0
          DO i = 1, ndim
            stdverens = 0.0D0
            DO j = 1, nrens
              stdverens = stdverens + verens(i,j)**2
            END DO
            stdverens = dsqrt(stdverens/dble(real(nrens)))
            DO j = 1, nrens
              verens(i,j) = verens(i,j)/stdverens
            END DO
          END DO

          CALL multa(verens,x5,ndim,nrens,ndim)

! subtract mean from verens
          DO i = 1, ndim
            aveverens = sum(verens(i,1:nrens))/dble(real(nrens))
            DO j = 1, nrens
              verens(i,j) = verens(i,j) - aveverens
            END DO
          END DO

! compute average variance over all ndim states
          std(:) = 0.0D0
          DO j = 1, nrens
            DO i = 1, ndim
              std(i) = std(i) + verens(i,j)**2
            END DO
          END DO
          DO i = 1, ndim
            std(i) = dsqrt(std(i)/dble(real(nrens)))
          END DO
          stdverens = sum(std(1:ndim))/dble(real(ndim))

! Compute correction to ensemble perturbations
          verx = -1.0D0/dble(real(nrens))
          DO j = 1, nrens
            verx(j,j) = verx(j,j) + 1.0D0
          END DO
          DO j = 1, nrens
            DO i = 1, nrens
              verx(i,j) = verx(i,j)/stdverens
            END DO
          END DO
          DO j = 1, nrens
            DO i = 1, nrens
              verx(i,j) = verx(i,j) + 1.0D0/dble(real(nrens))
            END DO
          END DO

          CALL dgemm('n','n',nrens,nrens,nrens,1.0D0,x5,nrens,verx, &
            nrens,0.0,xx,nrens)
          x5 = xx

          if(assimstp_switch) WRITE(12,*) '      analysis: inflation factor= ', &
               1.0/stdverens

        END SUBROUTINE

      END MODULE mod_anafunc






