!>    @brief (reads the state variables) and the Jacoby matrix
!>    @param[in] iseed seeding index number
!>    @param[in] ismpl local sample index
!>    @details
!> Read state data from hdf5 output file (and Jacoby matrix).\n
!> Only the OpenMP-Master reads the data, all other threads waiting for his reading\n
!> -> special sorted reading -> special dependency\n
!> -> needs OpenMP "special ordered" usage (no adaption from standard).\n
      SUBROUTINE read_joutt_hdf(iseed,ismpl)
#ifndef noHDF
        USE hdf5
        use mod_hdf5_vars, only: error
#endif
        use arrays
#ifdef AD
        use g_arrays
#endif
#ifdef AD_RM
        use arrays_ad
#endif
        use mod_genrl
        use mod_genrlc
        use mod_flow
        use mod_conc
        use mod_linfos
        use mod_inverse
        use mod_data
        use mod_OMP_TOOLS
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j
        INCLUDE 'OMP_TOOLS.inc'

        character (len=256) :: filename
        character (len=10) :: snumber
        DOUBLE PRECISION, ALLOCATABLE :: ltmp(:,:)
        INTEGER iseed, i1, i2, i1s, i2s, i3, i4, imaster

#ifndef noHDF
!       compute the master index
        imaster = ((iseed -1) /Tlevel_0) *Tlevel_0 +1
!       each thread needs to place here its seeding index
        seed_index(ismpl) = iseed
!
        CALL omp_ordered_dep_begin(iseed, Tlevel_0)
!$OMP   flush(seed_index)
        IF (iseed == imaster) THEN
!         read the master copy of all state variables
          CALL read_outt_hdf(ismpl)
!         mark the master sample for all other active threads
          seed_index(nsmpl+1) = ismpl
        ELSE
!         copy from the master
          IF (head_active) CALL dcopy(i0*j0*k0,head(1,1,1,seed_index(nsmpl+1)),1,head(1,1,1,ismpl),1)
          IF (temp_active) CALL dcopy(i0*j0*k0,temp(1,1,1,seed_index(nsmpl+1)),1,temp(1,1,1,ismpl),1)
          IF (trans_active) CALL dcopy(i0*j0*k0*ntrans,conc(1,1,1,1,seed_index(nsmpl+1)),1,conc(1,1,1,1,ismpl),1)
          IF (pres_active) CALL dcopy(i0*j0*k0,pres(1,1,1,seed_index(nsmpl+1)),1,pres(1,1,1,ismpl),1)
        END IF
!$OMP   flush(seed_index)
        CALL omp_ordered_end(iseed)
!
        CALL omp_ordered_dep_begin(iseed, Tlevel_0)
        IF (iseed == imaster) THEN
          CALL chln(project,i1,i2)
          CALL chln(project_sfx(ismpl),i1s,i2s)
          WRITE(snumber,'(1I7)') iter_inv
          CALL chln(snumber,i3,i4)
          filename = project(i1:i2) // project_sfx(ismpl) (i1s:i2s) // &
            '_J_' // snumber(i3:i4) // '.h5'
          CALL h5open_f(error)
          ALLOCATE(ltmp(i0*j0*k0,mpara))
!$OMP     flush(seed_index)
!
#ifndef AD_RM
          IF (head_active) THEN
            CALL read4_hdf5(i0,j0,k0,mpara,ltmp,'j_head',filename)
            DO i = 1, Tlevel_0
              CALL dcopy(i0*j0*k0,ltmp(1,seed_index(i)),1,g_head(1,1,1,i),1)
            END DO
!$OMP       flush(g_head)
          END IF

          IF (temp_active) THEN
            CALL read4_hdf5(i0,j0,k0,mpara,ltmp,'j_temp',filename)
            DO i = 1, Tlevel_0
              CALL dcopy(i0*j0*k0,ltmp(1,seed_index(i)),1,g_temp(1,1,1,i),1)
            END DO
!$OMP       flush(g_temp)
          END IF

          IF (trans_active) THEN
            DO j = 1, ntrans
              WRITE(snumber,'(1A6,1I4.4)') 'j_conc', j
              CALL read4_hdf5(i0,j0,k0,mpara,ltmp,snumber,filename)
              DO i = 1, Tlevel_0
                CALL dcopy(i0*j0*k0,ltmp(1,seed_index(i)),1,g_conc(1,1,1,j,i),1)
              END DO
            END DO
!$OMP       flush(g_conc)
          END IF

          IF (pres_active) THEN
            CALL read4_hdf5(i0,j0,k0,mpara,ltmp,'j_pres',filename)
            DO i = 1, Tlevel_0
              CALL dcopy(i0*j0*k0,ltmp(1,seed_index(i)),1,g_pres(1,1,1,i),1)
            END DO
!$OMP       flush(g_pres)
          END IF

          DEALLOCATE(ltmp)
          CALL h5close_f(error)
        ELSE
!$OMP      flush(g_head)
!$OMP      flush(g_temp)
!$OMP      flush(g_conc)
!$OMP      flush(g_pres)
        END IF
#else
          IF (head_active) THEN
            CALL read4_hdf5(i0,j0,k0,mpara,ltmp,'j_head',filename)
            DO i = 1, Tlevel_0
              CALL dcopy(i0*j0*k0,ltmp(1,seed_index(i)),1,head_ad(1,1,1,i),1)
            END DO
!$OMP       flush(head_ad)
          END IF

          IF (temp_active) THEN
            CALL read4_hdf5(i0,j0,k0,mpara,ltmp,'j_temp',filename)
            DO i = 1, Tlevel_0
              CALL dcopy(i0*j0*k0,ltmp(1,seed_index(i)),1,temp_ad(1,1,1,i),1)
            END DO
!$OMP       flush(temp_ad)
          END IF

          IF (trans_active) THEN
            DO j = 1, ntrans
              WRITE(snumber,'(1A6,1I4.4)') 'j_conc', j
              CALL read4_hdf5(i0,j0,k0,mpara,ltmp,snumber,filename)
              DO i = 1, Tlevel_0
                CALL dcopy(i0*j0*k0,ltmp(1,seed_index(i)),1,conc_ad(1,1,1,j,i),1)
              END DO
            END DO
!$OMP       flush(conc_ad)
          END IF

          IF (pres_active) THEN
            CALL read4_hdf5(i0,j0,k0,mpara,ltmp,'j_pres',filename)
            DO i = 1, Tlevel_0
              CALL dcopy(i0*j0*k0,ltmp(1,seed_index(i)),1,pres_ad(1,1,1,i),1)
            END DO
!$OMP       flush(pres_ad)
          END IF


          DEALLOCATE(ltmp)
          CALL h5close_f(error)
        ELSE
!$OMP      flush(head_ad)
!$OMP      flush(temp_ad)
!$OMP      flush(conc_ad)
!$OMP      flush(pres_ad)
        END IF
#endif
        CALL omp_ordered_end(iseed)
#else
        WRITE(*,*) &
          'error: need HDF5 support for reading time step data!'
        STOP
#endif
        RETURN
      END
