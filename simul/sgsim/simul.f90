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

!> @brief SIMUL wrapper for SGSIM call
!> @param[in] realz realisation (no use at this time : AW)
!> @param[in] ismpl local sample index
!> @param[in] imodel parameter file index
!> @details
!>
!>               Sequential Gaussian Simulation\n
!>               ******************************\n
!>The program is executed with no command line arguments.  The user\n
!>will be prompted for the name of a parameter file.  The parameter\n
!>file is described in the documentation (see the example sgsim.par)\n
!>The output file will be a GEOEAS file containing the simulated values\n
!>The file is ordered by x,y,z, and then simulation (i.e., x cycles\n
!>fastest, then y, then z, then simulation number).  The values will be\n
!>backtransformed to the original data values if a normal scores\n
!>transform was performed.\n
      SUBROUTINE simul(ismpl,realz,imodel)
        use mod_genrl
        use mod_simul
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        IMPLICIT integer (I-N)
        INCLUDE 'sgsim.inc'
!     realisation (no use at this time : AW)
        integer :: realz, imodel
        INTRINSIC int, sin, dble


!AW !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     (re-)init RNG - makes the start of all realisations identical
!aw   Call set_ival(MAXOP1,0,ixv)
        DO i = 1, maxop1
          ixv(i) = int(65536.0D0*sin(dble(i)))
        END DO
!AW !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Input/Output units used:

        CALL omp_new_file_handler(lin,1)
        CALL omp_new_file_handler(lout,2)
        CALL omp_new_file_handler(ldbg,3)
        CALL omp_new_file_handler(llvm,4)

! Read the parameters and data (transform as required):

!$OMP critical
        CALL readparm(realz,imodel)
!$OMP end critical

! Call sgsim for the simulation(s):

        CALL sgsim(realz)

! Finished:

!        WRITE(*,9998) version
9998    FORMAT (/' SGSIM Version: ',F5.3,' Finished'/)
        RETURN
      END
