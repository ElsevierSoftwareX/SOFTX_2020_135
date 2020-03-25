!>    @brief writes out some state variables and Jacoby matrix
!>    @param[in] iseed seeding index number
!>    @param[in] ismpl local sample index
!>    @details
!> Create an hdf5 output file for eacht time step,\n
!> with additional Jacoby information -> Jacoby matrix needs to be sorted -> needs OpenMP "ordered"!\n
      SUBROUTINE write_joutt_hdf(iseed,ismpl)
#ifndef noHDF
        USE hdf5
        use mod_hdf5_vars, only: error
#endif
        use arrays
#ifndef AD_RM
        use g_arrays
#else
        use arrays_ad
#endif
        use mod_genrl
        use mod_genrlc
        use mod_flow
        use mod_temp
        use mod_conc
        use mod_inverse
        use mod_data
        use mod_linfos
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j
        DOUBLE PRECISION dx, dy, dz
        INTEGER i1, i2, i3, i4, i1s, i2s, iseed

        character (len=80) :: filename
        character (len=10) :: snumber
        character (len=16) :: pname

#ifndef noHDF
!      File identifiers
        INTEGER (hid_t) file_id

!      Dataset identifier
        INTEGER (hid_t) dset_id

!      Data space identifier
        INTEGER (hid_t) dataspace
#endif

#ifdef HDF6432
!     need for 32Bit HDF5 library
        INTEGER (kind=4) :: rank
        INTEGER (kind=4) :: gzlevel
        INTEGER (kind=4), allocatable :: inttmp(:,:,:)
#else
#ifdef HDF64
        INTEGER (kind=8) :: rank
        INTEGER (kind=8) :: gzlevel
        INTEGER (kind=8), allocatable :: inttmp(:,:,:)
#else
        INTEGER rank
        INTEGER gzlevel
        INTEGER, ALLOCATABLE :: inttmp(:,:,:)
#endif
#endif

        PARAMETER (rank=3)

!     gzip compression level
        PARAMETER (gzlevel=9)

#ifndef noHDF
!      Data buffers
        INTEGER (hsize_t) dims(rank)
        INTEGER (hsize_t) data_dims(7)

!     for chunk size and compression
        INTEGER (hid_t) plist_id
        INTEGER (hsize_t) chunk_dims(rank)
#endif

!     temp. needed for 'uindex' & 'idata'
        DOUBLE PRECISION, ALLOCATABLE :: dp4tmp(:,:,:,:)
        character (len=16), dimension (:), allocatable :: ctmp

        INTEGER lblank
        EXTERNAL lblank

#ifndef noHDF

#ifdef NOJAC
        RETURN
#endif

        CALL chln(project,i1,i2)
        CALL chln(project_sfx(ismpl),i1s,i2s)
        WRITE(snumber,'(1I7)') iter_inv
        CALL chln(snumber,i3,i4)
        IF (i1s==0) THEN
          filename = project(i1:i2) // '_J_' // snumber(i3:i4) // &
            '.h5'
        ELSE
          filename = project(i1:i2) // project_sfx(ismpl) (i1s:i2s) // &
            '_J_' // snumber(i3:i4) // '.h5'
        END IF

        IF (linfos(2)>=1) THEN
          WRITE(*,'(3A,1I4,1A,1I4)') '  [W] : HDF5 Jacobian to "', &
            filename(1:lblank(filename)), '", seeding=', iseed, &
            ', sample=', ismpl
        END IF


#ifdef fOMP
! !!! to avoid compiler bugs
!--- C$OMP ordered
        CALL omp_ordered_begin(iseed)
!$OMP critical
#endif
!      Initialize FORTRAN interface.
        CALL h5open_f(error)


! ############
        IF (iseed==1) THEN
! ############

!      Create a new file, late only open it for read and writes
          CALL h5fcreate_f(filename,h5f_acc_trunc_f,file_id,error, &
            h5p_default_f,h5p_default_f)

          dims(1) = i0
          dims(2) = j0
          dims(3) = k0
          data_dims(1) = i0
          data_dims(2) = j0
          data_dims(3) = k0
          chunk_dims(1) = min(20,i0)
          chunk_dims(2) = min(20,j0)
          chunk_dims(3) = min(20,k0)

!     for compression only
          CALL h5pcreate_f(h5p_dataset_create_f,plist_id,error)
          CALL h5pset_chunk_f(plist_id,rank,chunk_dims,error)
          CALL h5pset_deflate_f(plist_id,gzlevel,error)

          ALLOCATE(dp4tmp(i0,j0,k0,3))
!     not very performant !!!
          DO i3 = 1, k0
            dz = delza(i3)
            DO i2 = 1, j0
              dy = delya(i2)
              DO i1 = 1, i0
                dx = delxa(i1)
                dp4tmp(i1,i2,i3,1) = dx
                dp4tmp(i1,i2,i3,2) = dy
                dp4tmp(i1,i2,i3,3) = dz
              END DO
            END DO
          END DO

          IF (out_ijk(cout_i)) THEN
!      Create data space for the dataset.
            CALL h5screate_simple_f(rank,dims,dataspace,error)
!      Create dataset "A" inside file "f_name".
            CALL h5dcreate_f(file_id,'x',h5t_native_double,dataspace, &
              dset_id,error,plist_id)
!      Write 'A' to the dataset
            CALL h5dwrite_f(dset_id,h5t_native_double,dp4tmp(1,1,1,1), &
              data_dims,error)
!      Close file, dataset and dataspace identifiers.
            CALL h5sclose_f(dataspace,error)
            CALL h5dclose_f(dset_id,error)
          END IF

          IF (out_ijk(cout_j)) THEN
            CALL h5screate_simple_f(rank,dims,dataspace,error)
            CALL h5dcreate_f(file_id,'y',h5t_native_double,dataspace, &
              dset_id,error,plist_id)
            CALL h5dwrite_f(dset_id,h5t_native_double,dp4tmp(1,1,1,2), &
              data_dims,error)
            CALL h5sclose_f(dataspace,error)
            CALL h5dclose_f(dset_id,error)
          END IF

          IF (out_ijk(cout_k)) THEN
            CALL h5screate_simple_f(rank,dims,dataspace,error)
            CALL h5dcreate_f(file_id,'z',h5t_native_double,dataspace, &
              dset_id,error,plist_id)
            CALL h5dwrite_f(dset_id,h5t_native_double,dp4tmp(1,1,1,3), &
              data_dims,error)
            CALL h5sclose_f(dataspace,error)
            CALL h5dclose_f(dset_id,error)
          END IF
          DEALLOCATE(dp4tmp)

          ALLOCATE(inttmp(i0,j0,k0))
          DO i3 = 1, k0
            DO i2 = 1, j0
              DO i1 = 1, i0
                inttmp(i1,i2,i3) = uindex(i1,i2,i3)
              END DO
            END DO
          END DO
          IF (out_ijk(cout_uindex)) THEN
            CALL h5screate_simple_f(rank,dims,dataspace,error)
            CALL h5dcreate_f(file_id,'uindex',h5t_native_integer, &
              dataspace,dset_id,error,plist_id)
            CALL h5dwrite_f(dset_id,h5t_native_integer,inttmp, &
              data_dims,error)
            CALL h5sclose_f(dataspace,error)
            CALL h5dclose_f(dset_id,error)
          END IF
          DEALLOCATE(inttmp)

! - main arrays -
          CALL h5screate_simple_f(rank,dims,dataspace,error)
          CALL h5dcreate_f(file_id,'head',h5t_native_double,dataspace, &
            dset_id,error,plist_id)
          CALL h5dwrite_f(dset_id,h5t_native_double,head(1,1,1,ismpl), &
            data_dims,error)
          CALL h5sclose_f(dataspace,error)
          CALL h5dclose_f(dset_id,error)

          CALL h5screate_simple_f(rank,dims,dataspace,error)
          CALL h5dcreate_f(file_id,'temp',h5t_native_double,dataspace, &
            dset_id,error,plist_id)
          CALL h5dwrite_f(dset_id,h5t_native_double,temp(1,1,1,ismpl), &
            data_dims,error)
          CALL h5sclose_f(dataspace,error)
          CALL h5dclose_f(dset_id,error)

!     convert [Pa] into [MPa]
          DO i3 = 1, k0
            DO i2 = 1, j0
              DO i1 = 1, i0
                x(i1,i2,i3,ismpl) = pres(i1,i2,i3,ismpl)*pa_conv1
              END DO
            END DO
          END DO
          CALL h5screate_simple_f(rank,dims,dataspace,error)
          CALL h5dcreate_f(file_id,'pres',h5t_native_double,dataspace, &
            dset_id,error,plist_id)
          CALL h5dwrite_f(dset_id,h5t_native_double,x(1,1,1,ismpl), &
            data_dims,error)
          CALL h5sclose_f(dataspace,error)
          CALL h5dclose_f(dset_id,error)

          IF (trans_active) THEN
            DO i = 1, ntrans
              WRITE(snumber,'(1A4,1I4.4)') 'conc', i
              CALL h5screate_simple_f(rank,dims,dataspace,error)
              CALL h5dcreate_f(file_id,snumber,h5t_native_double, &
                dataspace,dset_id,error,plist_id)
              CALL h5dwrite_f(dset_id,h5t_native_double, &
                conc(1,1,1,i,ismpl),data_dims,error)
              CALL h5sclose_f(dataspace,error)
              CALL h5dclose_f(dset_id,error)
            END DO
          END IF

! ------------------------------
          CALL h5pclose_f(plist_id,error)
          CALL h5fclose_f(file_id,error)

! ############
        END IF
! ############


! ------------------------------
!     jacoby matrix
#ifndef AD_RM
        CALL add_cube(i0,j0,k0,iseed,g_head(1,1,1,ismpl),'j_head', filename)
        CALL add_cube(i0,j0,k0,iseed,g_temp(1,1,1,ismpl),'j_temp', filename)
        IF (trans_active) THEN
          DO i = 1, ntrans
            WRITE(snumber,'(1A6,1I4.4)') 'j_conc', i
            CALL add_cube(i0,j0,k0,iseed,g_conc(1,1,1,i,ismpl), snumber,filename)
          END DO
        END IF
        CALL add_cube(i0,j0,k0,iseed,g_pres(1,1,1,ismpl),'j_pres', filename)
#else
        CALL add_cube(i0,j0,k0,iseed,head_ad(1,1,1,ismpl),'j_head', filename)
        CALL add_cube(i0,j0,k0,iseed,temp_ad(1,1,1,ismpl),'j_temp', filename)
        IF (trans_active) THEN
          DO i = 1, ntrans
            WRITE(snumber,'(1A6,1I4.4)') 'j_conc', i
            CALL add_cube(i0,j0,k0,iseed,conc_ad(1,1,1,i,ismpl), snumber,filename)
          END DO
        END IF
        CALL add_cube(i0,j0,k0,iseed,pres_ad(1,1,1,ismpl),'j_pres', filename)

#endif
! ------------------------------
        IF (iseed==1) THEN
          ALLOCATE(ctmp(mpara))
!        collect names
          DO i = 1, mpara
            CALL param_name(i,pname,ismpl)
            WRITE(ctmp(i),'(1A16)') pname
          END DO
          CALL write2_hdf5_char(16,mpara,ctmp,'para_name',filename)
          DEALLOCATE(ctmp)
        END IF

! ------------------------------


!     close interface
        CALL h5close_f(error)
#ifdef fOMP
!$OMP end critical
! !!! to avoid compiler bugs
!--- C$OMP end ordered
        CALL omp_ordered_end(iseed)
#endif

#endif
        RETURN
      END
