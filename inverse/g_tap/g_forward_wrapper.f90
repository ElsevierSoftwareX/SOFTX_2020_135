!>    @brief wrapper for the forward simulation call
!>    @param[in] iseed seeding index
!>    @param[in] iter time iteration counter
!>    @param[in] ismpl local sample index
!>    @details
!> using delayed propagation of derivatives,\n
!> code to use reading/writing data files instead of compute values\n
      SUBROUTINE g_forward_wrapper(iter,iseed,ismpl)
        use arrays
        use mod_genrl
        use mod_genrlc
        use mod_time
        IMPLICIT NONE
        integer :: ismpl
        INTEGER iter, iseed

        IF (lread_joutt) THEN
          WRITE(project_sfx(ismpl),'(1A,1I7.7)') '_all', iter
!         read data of each time step
          CALL read_joutt_hdf(iseed,ismpl)
          project_sfx(ismpl) = ' '
        ELSE
! ######### Forward Iteration ######
!       - normal Picard method -
!         !!! AD runtime optimisation !!!
          CALL forward_picard(iter,ismpl)
!         call gradient computation
          CALL g_forward_picard(iter,ismpl)
! ##################################
        END IF
!
        IF (lwrite_joutt) THEN
          WRITE(project_sfx(ismpl),'(1A,1I7.7)') '_all', iter
!          write data for each time step
          CALL write_joutt_hdf(iseed,ismpl)
          project_sfx(ismpl) = ' '
        END IF
!       write data for specific time steps
        CALL write_joutt(iseed,ismpl)
!
        RETURN
      END
