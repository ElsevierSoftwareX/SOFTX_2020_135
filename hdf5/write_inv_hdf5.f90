! ******************************************************
!   WARNING: need 32 Bit version of the HDF5 library !
! ******************************************************

!>    @brief writes out all jacoby matrices into a HDF5 file
!>    @param[in] ismpl local sample index
!>    @details
!>    create an hdf5 output file (from jacobian)\n
      SUBROUTINE write_inv_hdf(ismpl)
#ifndef noHDF
        USE hdf5
        use mod_hdf5_vars, only: error
#endif
        use arrays
!        use g_arrays
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
        integer :: ia, ib, ic
        DOUBLE PRECISION dx, dy, dz
        INTEGER i1, i2, i3, i4

        DOUBLE PRECISION vxc, vyc, vzc, kx, ky, kz, lx, ly, lz, por, &
          rhof, visf
        EXTERNAL vxc, vyc, vzc, kx, ky, kz, lx, ly, lz, por, rhof, &
          visf

        character (len=256) :: filename
        character (len=8) :: snumber
        character (len=16) :: pname

#ifndef noHDF
!      File identifiers
        INTEGER (hid_t) file_id

!      Dataset identifier
        INTEGER (hid_t) dset_id

!      Data space identifier
        INTEGER (hid_t) dataspace
#endif

!     error: flag to check operation success
#ifdef HDF6432
        INTEGER (kind=4) :: rank
        INTEGER (kind=4) :: gzlevel
        INTEGER (kind=4), allocatable :: inttmp(:)
#else
#ifdef HDF64
        INTEGER (kind=8) :: rank
        INTEGER (kind=8) :: gzlevel
        INTEGER (kind=8), allocatable :: inttmp(:)
#else
        INTEGER rank
        INTEGER gzlevel
        INTEGER, ALLOCATABLE :: inttmp(:)
#endif
#endif

!     data buffers
        PARAMETER (rank=3)

!     gzip compression level
        PARAMETER (gzlevel=9)

#ifndef noHDF
        INTEGER (hsize_t) dims(rank)
        INTEGER (hsize_t) data_dims(7)

!     for chunk size and compression
        INTEGER (hid_t) plist_id
        INTEGER (hsize_t) chunk_dims(rank)
#endif

!     temp. array, need for 'uindex' & 'idata'
        DOUBLE PRECISION, ALLOCATABLE :: dp4tmp(:,:,:,:)
        CHARACTER (len=16), dimension (:), allocatable :: ctmp

        INTEGER lblank
        EXTERNAL lblank

#ifndef noHDF

        IF (write_disable) RETURN
#ifdef NOJAC
        RETURN
#endif

        CALL chln(project,i1,i2)
        WRITE(snumber,'(1I7)') iter_inv
        CALL chln(snumber,i3,i4)
        filename = project(i1:i2) // '_J_' // snumber(i3:i4) // '.h5'

        IF (linfos(2)>=1) THEN
          WRITE(*,'(3A)') '  [W] : HDF5-Inverse to "', &
            filename(1:lblank(filename)), '"'
        END IF

! ------------------------------
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


#ifdef fOMP
!$OMP critical
#endif
!      Initialize FORTRAN interface.
        CALL h5open_f(error)

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

        CALL h5screate_simple_f(rank,dims,dataspace,error)
        CALL h5dcreate_f(file_id,'y',h5t_native_double,dataspace, &
          dset_id,error,plist_id)
        CALL h5dwrite_f(dset_id,h5t_native_double,dp4tmp(1,1,1,2), &
          data_dims,error)
        CALL h5sclose_f(dataspace,error)
        CALL h5dclose_f(dset_id,error)

        CALL h5screate_simple_f(rank,dims,dataspace,error)
        CALL h5dcreate_f(file_id,'z',h5t_native_double,dataspace, &
          dset_id,error,plist_id)
        CALL h5dwrite_f(dset_id,h5t_native_double,dp4tmp(1,1,1,3), &
          data_dims,error)
        CALL h5sclose_f(dataspace,error)
        CALL h5dclose_f(dset_id,error)
        DEALLOCATE(dp4tmp)
! ------------------------------

        ALLOCATE(inttmp(i0*j0*k0))
        DO i = 1, i0*j0*k0
          CALL ijk_m(i,ia,ib,ic)
          inttmp(i) = uindex(ia,ib,ic)
        END DO
        CALL h5screate_simple_f(rank,dims,dataspace,error)
        CALL h5dcreate_f(file_id,'uindex',h5t_native_integer, &
          dataspace,dset_id,error,plist_id)
        CALL h5dwrite_f(dset_id,h5t_native_integer,inttmp,data_dims, &
          error)
        CALL h5sclose_f(dataspace,error)
        CALL h5dclose_f(dset_id,error)
        DEALLOCATE(inttmp)

! ------------------------------------------------
        IF (ndata>=1) THEN

          dims(1) = ndata
          dims(2) = mpara
          dims(3) = 1
          data_dims(1) = ndata
          data_dims(2) = mpara
          data_dims(3) = 1
          chunk_dims(1) = min(8000,ndata)
          chunk_dims(2) = 1
          chunk_dims(3) = 1

          IF (mpara>0) THEN
!     for compression only
            CALL h5pclose_f(plist_id,error)
            CALL h5pcreate_f(h5p_dataset_create_f,plist_id,error)
            CALL h5pset_chunk_f(plist_id,rank,chunk_dims,error)
            CALL h5pset_deflate_f(plist_id,gzlevel,error)

            CALL h5screate_simple_f(rank,dims,dataspace,error)
            CALL h5dcreate_f(file_id,'jacoby',h5t_native_double, &
              dataspace,dset_id,error,plist_id)
            CALL h5dwrite_f(dset_id,h5t_native_double,jac,data_dims, &
              error)
            CALL h5sclose_f(dataspace,error)
            CALL h5dclose_f(dset_id,error)
          END IF

! ------------------------------------------------
          dims(1) = ndata
          dims(2) = n_idata
          dims(3) = 1
          data_dims(1) = ndata
          data_dims(2) = n_idata
          data_dims(3) = 1
          chunk_dims(1) = min(8000,ndata)
          chunk_dims(2) = 1
          chunk_dims(3) = 1

!     for compression only
          CALL h5pclose_f(plist_id,error)
          CALL h5pcreate_f(h5p_dataset_create_f,plist_id,error)
          CALL h5pset_chunk_f(plist_id,rank,chunk_dims,error)
          CALL h5pset_deflate_f(plist_id,gzlevel,error)


          ALLOCATE(inttmp(ndata*n_idata))
          DO j = 1, n_idata
            DO i = 1, ndata
              inttmp(i+(j-1)*ndata) = idata(i,j)
            END DO
          END DO
          CALL h5screate_simple_f(rank,dims,dataspace,error)
          CALL h5dcreate_f(file_id,'idata',h5t_native_integer, &
            dataspace,dset_id,error,plist_id)
          CALL h5dwrite_f(dset_id,h5t_native_integer,inttmp,data_dims, &
            error)
          CALL h5sclose_f(dataspace,error)
          CALL h5dclose_f(dset_id,error)
          DEALLOCATE(inttmp)

          dims(1) = ndata
          dims(2) = n_ddata
          dims(3) = 1
          data_dims(1) = ndata
          data_dims(2) = n_ddata
          data_dims(3) = 1
          chunk_dims(1) = min(8000,ndata)
          chunk_dims(2) = 1
          chunk_dims(3) = 1

!     for compression only
          CALL h5pclose_f(plist_id,error)
          CALL h5pcreate_f(h5p_dataset_create_f,plist_id,error)
          CALL h5pset_chunk_f(plist_id,rank,chunk_dims,error)
          CALL h5pset_deflate_f(plist_id,gzlevel,error)

          CALL h5screate_simple_f(rank,dims,dataspace,error)
          CALL h5dcreate_f(file_id,'data',h5t_native_double,dataspace, &
            dset_id,error,plist_id)
          CALL h5dwrite_f(dset_id,h5t_native_double,ddata,data_dims, &
            error)
          CALL h5sclose_f(dataspace,error)
          CALL h5dclose_f(dset_id,error)

        END IF
! ------------------------------------------------

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
        CALL h5pclose_f(plist_id,error)
        CALL h5pcreate_f(h5p_dataset_create_f,plist_id,error)
        CALL h5pset_chunk_f(plist_id,rank,chunk_dims,error)
        CALL h5pset_deflate_f(plist_id,gzlevel,error)


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


        CALL h5pclose_f(plist_id,error)
        CALL h5fclose_f(file_id,error)

        IF (covar/=0 .OR. resmat/=0) THEN
          CALL write2_hdf5(mpara,mpara,covar_p,'covariance_p', &
            filename)
          CALL write2_hdf5(mpara,mpara,resmat_p,'resolution_p', &
            filename)
        END IF

        ALLOCATE(ctmp(mpara))
!     collect names
        DO i = 1, mpara
          CALL param_name(i,pname,ismpl)
          WRITE(ctmp(i),'(1A16)') pname
        END DO
        CALL write2_hdf5_char(16,mpara,ctmp,'para_name',filename)
        DEALLOCATE(ctmp)

!      Close FORTRAN interface.
        CALL h5close_f(error)
#ifdef fOMP
!$OMP end critical
#endif

#endif
        RETURN
      END
