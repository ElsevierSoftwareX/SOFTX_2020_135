!>    @brief reads the state variables
!>    @param[in] ismpl local sample index
!>    @details
!>    read hdf5 output file (from jacobian collection file)\n
      SUBROUTINE read_outt_hdf(ismpl)
#ifndef noHDF
        USE hdf5
        use mod_hdf5_vars, only: error
#endif
        use arrays
        use mod_genrl
        use mod_genrlc
        use mod_flow
        use mod_conc
        use mod_data
        use mod_inverse
        use mod_linfos
        IMPLICIT NONE
        character (len=256) :: filename
        character (len=10) :: snumber
        ! DOUBLE PRECISION, ALLOCATABLE :: ltmp(:,:)
        INTEGER i1, i2, i1s, i2s, i3, i4, i, ismpl

#ifndef noHDF
        CALL chln(project,i1,i2)
        CALL chln(project_sfx(ismpl),i1s,i2s)
        WRITE(snumber,'(1I7)') iter_inv
        CALL chln(snumber,i3,i4)
        filename = project(i1:i2) // project_sfx(ismpl) (i1s:i2s) // &
          '_J_' // snumber(i3:i4) // '.h5'

        CALL h5open_f(error)

        IF (head_active) THEN
          CALL read3_hdf5(i0,j0,k0,head(1,1,1,ismpl),'head',filename)
        END IF

        IF (temp_active) THEN
          CALL read3_hdf5(i0,j0,k0,temp(1,1,1,ismpl),'temp',filename)
        END IF

        IF (trans_active) THEN
          DO i = 1, ntrans
            WRITE(snumber,'(1A4,1I4.4)') 'conc', i
            CALL read3_hdf5(i0,j0,k0,conc(1,1,1,i,ismpl),snumber, &
              filename)
          END DO
        END IF

        IF (pres_active) THEN
          CALL read3_hdf5(i0,j0,k0,pres(1,1,1,ismpl),'pres',filename)
!         convert [MPa] into [Pa]
          CALL dscal(i0*j0*k0,pa_conv,pres(1,1,1,ismpl),1)
        END IF

        CALL h5close_f(error)
#else
        WRITE(*,*) &
          'error: need HDF5 support for reading time step data!'
        STOP
#endif
        RETURN
      END
