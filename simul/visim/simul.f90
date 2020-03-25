!>Copyright (C) 1996, The Board of Trustees of the Leland Stanford     %\n
!>Junior University.  All rights reserved.                             %\n
!>
!>The programs in GSLIB are distributed in the hope that they will be  %\n
!>useful, but WITHOUT ANY WARRANTY.  No author or distributor accepts  %\n
!>responsibility to anyone for the consequences of using them or for   %\n
!>whether they serve any particular purpose or work at all, unless he  %\n
!>says so in writing.  Everyone is granted permission to copy, modify  %\n
!>and redistribute the programs in GSLIB, but only under the condition %\n
!>that this notice and the above copyright notice remain intact.       %\n

!> @brief SIMUL wrapper for VISIM call
!> @param[in] realz realisation (no use at this time : AW)
!> @param[in] ismpl local sample index
!> @param[in] imodel parameter file index
!> @details
!>
!>               Volume Integration SIMulation\n
!>               ******************************c\n
!>The program is executed either by specifying the \n
!>parameter file as a commandline parameter  \n
!>  visim visim.par\n
!>or with no command line arguments\n
!>  visim\n
!>In the latter case the user\n
!>will be prompted for the name of a parameter file.  The parameter\n
!>file is described in the documentation\n
!>The output file will be a GEOEAS file containing the simulated values\n
!>The file is ordered by x,y,z, and then simulation (i.e., x cycles\n
!>fastest, then y, then z, then simulation number).  \n
!>can be transformed to reproduce a target histogram.\n
      SUBROUTINE simul(ismpl,realz,imodel)
      USE simul_arrays
      use mod_genrl
      use mod_simul
      use mod_OMP_TOOLS
      IMPLICIT NONE
      include 'visim.inc'
      include 'gslib.inc'
      include 'OMP_TOOLS.inc'
!
      integer :: ismpl
      integer :: i, j
      character (len=80) :: tmpfl
      DOUBLE PRECISION float,temp
      LOGICAL testfl
!     realisation (no use at this time : AW)
      integer :: realz,imodel,llout
      INTRINSIC int,sin,cos,dble
!
!
!AW CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!           (re-)init RNG - makes the start of all realisations identical
!aw   Call set_ival(MAXOP1,0,ixv)
!aw   Call set_ival(MAXOP1,0,ixv2)
      DO i = 1,MAXOP1
        ixv(i) = int(65536.0d0*sin(dble(i)))
        ixv2(i) = int(65536.0d0*cos(dble(i)))
      END DO
!AW CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
! Input/Output units used:
!
      CALL omp_new_file_handler(lin,8)
      CALL omp_new_file_handler(lout,9)
      CALL omp_new_file_handler(ldbg,10)
      CALL omp_new_file_handler(llvm,11)
      CALL omp_new_file_handler(lkv,12)
      CALL omp_new_file_handler(lout_mean,13)
      CALL omp_new_file_handler(llout,1)

!
! Read the parameters and data (transform as required):
!
!$OMP critical
      CALL readparm(realz,imodel)
!$OMP end critical
!
! Create Conditional Prob Lookup Table
!
!      if ((idrawopt.eq.2).or.(idrawopt.eq.4)) then
      IF (idrawopt.EQ.1) THEN
        CALL create_condtab
      END IF


!
!  setup the krige estimation variance matrix for honoring the local data 
!  (Conditional simulation) if it is to be CalCulated in stead of read
!  from a seperate file.
!
      WRITE(ldbg,*) 'icond=',icond
      IF (icond.EQ.1 .AND. ivar.EQ.0) THEN
!        Call setup_krgvar
!     write(*,*) 'Kriging varianCe CalCulated'
      END IF

!
!  open the output file
!
      OPEN(lout,file=outfl,status='UNKNOWN')

      IF (doestimation.EQ.1) THEN
        WRITE(tmpfl,871) 'visim_estimation',outfl
        OPEN(lout_mean,file=tmpfl,status='UNKNOWN')
  871   FORMAT (A,'_',A)
        WRITE(lout_mean,110)
      END IF

      WRITE(lout,108)

  108 FORMAT ('VISIM Realizations',/,'1',/,'value')
  110 FORMAT ('VISIM ESTIMATION',/,'2',/,'mean',/,'std')


!     INITIALIZE Data2Volume lookup table
      CALL cov_data2vol(0,0,0,0,0,temp)
!      write(*,*) 'start  : CalCulating Data2Vol CovarianCe'
!      write(*,*) 'stop  : CalCulating Data2Vol CovarianCe'


!     INITIALIZE Volume2Volume lookup table
      CALL cov_vol2vol(0,0,temp)


! next line should be set via parameter file.
!      write(*,*) 'read_Covtable',read_Covtable
!      stop
      IF (read_covtable.EQ.1 .AND. nvol.GT.0) THEN
!     ALSO CHECK IF FILES EXIST
        WRITE(tmpfl,1871) 'cv2v',outfl
        INQUIRE (file=tmpfl,exist=testfl)
        IF (.NOT.testfl) THEN
          WRITE(*,*) 'Could not read Covariance lookup table : ',tmpfl, &
               ' - Continuing without reading...'
        ELSE
          IF (idbg.GT.0) THEN
            WRITE(*,*) 'Reading cv2v=',tmpfl
          END IF
!$OMP critical
          OPEN(llout,file=tmpfl,status='unknown',form='unformatted')
          DO i = 1,nvol
            READ (llout) (cv2v(i,j),j=1,nvol)
          END DO
          CLOSE(llout)
!$OMP end critical
        END IF
        WRITE(tmpfl,1871) 'cd2v',outfl
        INQUIRE (file=tmpfl,exist=testfl)
        IF (.NOT.testfl) THEN
          WRITE(*,*) 'Could not read Covariance lookup table : ',tmpfl, &
               ' - Continuing without reading...'
        ELSE
          IF (idbg.GT.0) THEN
            WRITE(*,*) 'Reading cd2v=',tmpfl
          END IF
!$OMP critical
          OPEN(llout,file=tmpfl,status='unknown',form='unformatted')
          DO i = 1,nxyz
            READ (llout) (cd2v(i,j),j=1,nvol)
          END DO
          CLOSE(llout)
!$OMP end critical
        END IF
      END IF
!
!      write(*,*) 'start : CalCulating Vol2Vol CovarianCe'
!      do i =1, nvol
!         do j =i, nvol
!            Call Cov_vol2vol(i,j,temp)
!            Cv2v(i,j) = dble(temp)
!         enddo
!      enddo
!      write(*,*) 'stop  : CalCulating Vol2Vol CovarianCe'


!        write(tmpfl,1871) 'cv2v_init',outfl
!        open(#9, file=tmpfl, status = 'unknown')

!        do i=1,nvol
!           do j=1,nvol
!              write(#9,*) i, j, cv2v(i,j)
!           end do
!        end do
!        close(#9)


!
!AW-c     Open File Handle for lambdas #99
!AW-c
!AW-      write(tmpfl,1871) 'lambda',outfl
!AW-      open(#99, file=tmpfl, status = 'unknown',form='unformatted')
!
!  begin the aCtual simulation
!

      DO isim = 1,nsim
!
! Call visim for the simulation(s):
        CALL visim(realz)
      END DO

      CLOSE(lout)
      IF (doestimation.EQ.1) CLOSE(lout_mean)

      IF (idbg.GT.-2) THEN

!         write(tmpfl,1871) 'txtCv2v',outfl
!         open(#9, file=tmpfl, status = 'unknown')
!         do i=1,nvol
!            do j=1,nvol
!               write(#9,*) i, j, Cv2v(i,j)
!            end do
!         end do
!         Close(#9)
!
        IF (nvol.GT.0) THEN
!
          WRITE(tmpfl,1871) 'cv2v',outfl
          OPEN(llout,file=tmpfl,status='unknown',form='unformatted')
          DO i = 1,nvol
            WRITE(llout) ((cv2v(i,j)),j=1,nvol)
          END DO
          CLOSE(llout)
!
          WRITE(tmpfl,1871) 'cd2v',outfl
          OPEN(llout,file=tmpfl,status='unknown',form='unformatted')
          DO i = 1,nxyz
            WRITE(llout) ((cd2v(i,j)),j=1,nvol)
          END DO
          CLOSE(llout)
        END IF
!
 1871   FORMAT (A,'_',A)
!
!     Close file handle #99 (lambdas)
!AW-not really needed         close(#99)
!
      END IF
!
! Finished:
!
      IF (idbg.GT.-1) THEN
        WRITE(*,9998) VERSION
 9998   FORMAT (/' VISIM Version: ',f6.4,' Finished'/)
      END IF
      RETURN
      END
