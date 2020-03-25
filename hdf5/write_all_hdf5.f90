! MIT License
!
! Copyright (c) 2020 SHEMAT-Suite
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.

! ******************************************************
!   WARNING: need 32 Bit version of the HDF5 library !
! ******************************************************

!>    @brief write out the physical state and all properties
!>    @param[in] ident file/iteration index number
!>    @param[in] ismpl local sample index
!>    @details
!>    create an hdf5 output file\n
      SUBROUTINE write_hdf(ident,ismpl)
#ifndef noHDF
        USE hdf5
        use mod_hdf5_vars, only: error
#endif
        use arrays
        use mod_genrl
        use mod_genrlc
        use mod_flow
        use mod_temp
        use mod_conc
        use mod_time
        use mod_linfos
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k
        double precision :: dx, dy, dz
        integer :: i1s, i2s, i1, i2, i3, i4, ident

        double precision, external :: vxc, vyc, vzc, kx, ky, kz, lx, ly, lz, por, &
             rhof, visf

        character (len=80) :: filename
        character (len=8) :: snumber

        logical, dimension (3) :: out_bc

#ifndef noHDF
!      File identifiers
        integer (kind=hid_t) :: file_id

!      Dataset identifier
        integer (kind=hid_t) :: dset_id

!      Data space identifier
        integer (kind=hid_t) :: dataspace
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
        INTEGER (kind=8),  allocatable :: inttmp(:,:,:)
#else
        integer :: rank
        integer :: gzlevel
        integer, allocatable, dimension (:,:,:) :: inttmp
#endif
#endif
        PARAMETER (rank=3)

!     gzip compression level
        PARAMETER (gzlevel=9)

#ifndef noHDF
!      Data buffers
        integer (kind=hsize_t), dimension (rank) :: dims
        integer (kind=hsize_t), dimension (7) :: data_dims

!     for chunk size and compression
        integer (kind=hid_t) :: plist_id
        integer (kind=hsize_t), dimension (rank) :: chunk_dims
#endif

        DOUBLE PRECISION, ALLOCATABLE :: dp3tmp(:,:,:), &
          dp4tmp(:,:,:,:)

!     bc-naming
        character (len=3), dimension (2) :: c_bt
        DATA c_bt/'bcd', 'bcn'/
        character (len=8) :: c_name

        integer, external :: lblank

#ifndef noHDF

#ifdef NOOUT
        RETURN
#endif

        CALL chln(project,i1,i2)
        CALL chln(project_sfx(ismpl),i1s,i2s)
        IF (ident>=0) THEN
          WRITE(snumber,'(1I7)') ident
        ELSE IF (ident==-1) THEN
          WRITE(snumber,'(A8)') 'final'
        ELSE IF (ident==-2) THEN
          WRITE(snumber,'(A8)') 'debug'
        ELSE IF (ident==-3) THEN
          WRITE(snumber,'(A8)') 'ens_mean'
        ELSE IF (ident==-4) THEN
          WRITE(snumber,'(A8)') 'mean'
        ELSE IF (ident==-5) THEN
          WRITE(snumber,'(A8)') 'ens_mean'
        END IF
        CALL chln(snumber,i3,i4)
        IF (i1s==0) THEN
          filename = project(i1:i2) // '_' // snumber(i3:i4) // '.h5'
        ELSE
          filename = project(i1:i2) // project_sfx(ismpl) (i1s:i2s) // &
            '_' // snumber(i3:i4) // '.h5'
        END IF

        IF (linfos(3)>=1) THEN
          WRITE(*,'(3A)') '  [W] : HDF5 to "', &
            filename(1:lblank(filename)), '"'
        END IF


#ifdef fOMP
!$OMP critical
#endif
!     Initialize FORTRAN interface.
        CALL h5open_f(error)

!     Create a new file, later only open it for read and writes
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

        ALLOCATE(inttmp(i0,j0,k0))
        ALLOCATE(dp4tmp(i0,j0,k0,7))
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
              dp4tmp(i1,i2,i3,4) = vxc(i1,i2,i3,ismpl)
              dp4tmp(i1,i2,i3,5) = vyc(i1,i2,i3,ismpl)
              dp4tmp(i1,i2,i3,6) = vzc(i1,i2,i3,ismpl)
            END DO
          END DO
        END DO

        IF (ident<=-1 .OR. out_ijk(cout_i)) THEN
!        Create data space for the dataset.
          CALL h5screate_simple_f(rank,dims,dataspace,error)
!        Create dataset "A" inside file "f_name".
          CALL h5dcreate_f(file_id,'x',h5t_native_double,dataspace, &
            dset_id,error,plist_id)
!        Write 'A' to the dataset
          CALL h5dwrite_f(dset_id,h5t_native_double,dp4tmp(1,1,1,1), &
            data_dims,error)
!        Close file, dataset and dataspace identifiers.
          CALL h5sclose_f(dataspace,error)
          CALL h5dclose_f(dset_id,error)
        END IF

        IF (ident<=-1 .OR. out_ijk(cout_j)) THEN
          CALL h5screate_simple_f(rank,dims,dataspace,error)
          CALL h5dcreate_f(file_id,'y',h5t_native_double,dataspace, &
            dset_id,error,plist_id)
          CALL h5dwrite_f(dset_id,h5t_native_double,dp4tmp(1,1,1,2), &
            data_dims,error)
          CALL h5sclose_f(dataspace,error)
          CALL h5dclose_f(dset_id,error)
        END IF

        IF (ident<=-1 .OR. out_ijk(cout_k)) THEN
          CALL h5screate_simple_f(rank,dims,dataspace,error)
          CALL h5dcreate_f(file_id,'z',h5t_native_double,dataspace, &
            dset_id,error,plist_id)
          CALL h5dwrite_f(dset_id,h5t_native_double,dp4tmp(1,1,1,3), &
            data_dims,error)
          CALL h5sclose_f(dataspace,error)
          CALL h5dclose_f(dset_id,error)
        END IF

        DO i3 = 1, k0
          DO i2 = 1, j0
            DO i1 = 1, i0
              inttmp(i1,i2,i3) = uindex(i1,i2,i3)
            END DO
          END DO
        END DO
        IF (ident<=-1 .OR. out_ijk(cout_uindex)) THEN
          CALL h5screate_simple_f(rank,dims,dataspace,error)
          CALL h5dcreate_f(file_id,'uindex',h5t_native_integer, &
            dataspace,dset_id,error,plist_id)
          CALL h5dwrite_f(dset_id,h5t_native_integer,inttmp,data_dims, &
            error)
          CALL h5sclose_f(dataspace,error)
          CALL h5dclose_f(dset_id,error)
        END IF

! - main arrays -
        IF (ident<=-1 .OR. out_pv(pv_head)) THEN
          CALL h5screate_simple_f(rank,dims,dataspace,error)
          CALL h5dcreate_f(file_id,'head',h5t_native_double,dataspace, &
            dset_id,error,plist_id)
          CALL h5dwrite_f(dset_id,h5t_native_double,head(1,1,1,ismpl), &
            data_dims,error)
          CALL h5sclose_f(dataspace,error)
          CALL h5dclose_f(dset_id,error)
        END IF

        IF (ident<=-1 .OR. out_pv(pv_temp)) THEN
          CALL h5screate_simple_f(rank,dims,dataspace,error)
          CALL h5dcreate_f(file_id,'temp',h5t_native_double,dataspace, &
            dset_id,error,plist_id)
          CALL h5dwrite_f(dset_id,h5t_native_double,temp(1,1,1,ismpl), &
            data_dims,error)
          CALL h5sclose_f(dataspace,error)
          CALL h5dclose_f(dset_id,error)
        END IF

        IF (ident<=-1 .OR. out_pv(pv_pres)) THEN
!        convert [Pa] into [MPa]
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
        END IF

        IF ((ident<=-1 .OR. out_pv(pv_conc)) .AND. trans_active) THEN
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

! ---------------

        IF (ident<=-1 .OR. out_ijk(cout_vx)) THEN
          CALL h5screate_simple_f(rank,dims,dataspace,error)
          CALL h5dcreate_f(file_id,'vx',h5t_native_double,dataspace, &
            dset_id,error,plist_id)
          CALL h5dwrite_f(dset_id,h5t_native_double,dp4tmp(1,1,1,4), &
            data_dims,error)
          CALL h5sclose_f(dataspace,error)
          CALL h5dclose_f(dset_id,error)
        END IF

        IF (ident<=-1 .OR. out_ijk(cout_vy)) THEN
          CALL h5screate_simple_f(rank,dims,dataspace,error)
          CALL h5dcreate_f(file_id,'vy',h5t_native_double,dataspace, &
            dset_id,error,plist_id)
          CALL h5dwrite_f(dset_id,h5t_native_double,dp4tmp(1,1,1,5), &
            data_dims,error)
          CALL h5sclose_f(dataspace,error)
          CALL h5dclose_f(dset_id,error)
        END IF

        IF (ident<=-1 .OR. out_ijk(cout_vz)) THEN
          CALL h5screate_simple_f(rank,dims,dataspace,error)
          CALL h5dcreate_f(file_id,'vz',h5t_native_double,dataspace, &
            dset_id,error,plist_id)
          CALL h5dwrite_f(dset_id,h5t_native_double,dp4tmp(1,1,1,6), &
            data_dims,error)
          CALL h5sclose_f(dataspace,error)
          CALL h5dclose_f(dset_id,error)
        END IF


        DO i3 = 1, k0
          DO i2 = 1, j0
            DO i1 = 1, i0
              dp4tmp(i1,i2,i3,7) = por(i1,i2,i3,ismpl)
            END DO
          END DO
        END DO
        IF (ident<=-1 .OR. out_prop(idx_por)) THEN
          CALL h5screate_simple_f(rank,dims,dataspace,error)
          CALL h5dcreate_f(file_id,'por',h5t_native_double,dataspace, &
            dset_id,error,plist_id)
          CALL h5dwrite_f(dset_id,h5t_native_double,dp4tmp(1,1,1,7), &
            data_dims,error)
          CALL h5sclose_f(dataspace,error)
          CALL h5dclose_f(dset_id,error)
        END IF

        DO i3 = 1, k0
          DO i2 = 1, j0
            DO i1 = 1, i0
              dp4tmp(i1,i2,i3,1) = kx(i1,i2,i3,ismpl)
              dp4tmp(i1,i2,i3,2) = ky(i1,i2,i3,ismpl)
              dp4tmp(i1,i2,i3,3) = kz(i1,i2,i3,ismpl)
              dp4tmp(i1,i2,i3,4) = lx(i1,i2,i3,ismpl)
              dp4tmp(i1,i2,i3,5) = ly(i1,i2,i3,ismpl)
              dp4tmp(i1,i2,i3,6) = lz(i1,i2,i3,ismpl)
            END DO
          END DO
        END DO
        IF (ident<=-1 .OR. out_prop(idx_an_kx)) THEN
          CALL h5screate_simple_f(rank,dims,dataspace,error)
          CALL h5dcreate_f(file_id,'kx',h5t_native_double,dataspace, &
            dset_id,error,plist_id)
          CALL h5dwrite_f(dset_id,h5t_native_double,dp4tmp(1,1,1,1), &
            data_dims,error)
          CALL h5sclose_f(dataspace,error)
          CALL h5dclose_f(dset_id,error)
        END IF

        IF (ident<=-1 .OR. out_prop(idx_an_ky)) THEN
          CALL h5screate_simple_f(rank,dims,dataspace,error)
          CALL h5dcreate_f(file_id,'ky',h5t_native_double,dataspace, &
            dset_id,error,plist_id)
          CALL h5dwrite_f(dset_id,h5t_native_double,dp4tmp(1,1,1,2), &
            data_dims,error)
          CALL h5sclose_f(dataspace,error)
          CALL h5dclose_f(dset_id,error)
        END IF

        IF (ident<=-1 .OR. out_prop(idx_kz)) THEN
          CALL h5screate_simple_f(rank,dims,dataspace,error)
          CALL h5dcreate_f(file_id,'kz',h5t_native_double,dataspace, &
            dset_id,error,plist_id)
          CALL h5dwrite_f(dset_id,h5t_native_double,dp4tmp(1,1,1,3), &
            data_dims,error)
          CALL h5sclose_f(dataspace,error)
          CALL h5dclose_f(dset_id,error)
        END IF

        IF (ident<=-1 .OR. out_prop(idx_an_lx)) THEN
          CALL h5screate_simple_f(rank,dims,dataspace,error)
          CALL h5dcreate_f(file_id,'lx',h5t_native_double,dataspace, &
            dset_id,error,plist_id)
          CALL h5dwrite_f(dset_id,h5t_native_double,dp4tmp(1,1,1,4), &
            data_dims,error)
          CALL h5sclose_f(dataspace,error)
          CALL h5dclose_f(dset_id,error)
        END IF

        IF (ident<=-1 .OR. out_prop(idx_an_ly)) THEN
          CALL h5screate_simple_f(rank,dims,dataspace,error)
          CALL h5dcreate_f(file_id,'ly',h5t_native_double,dataspace, &
            dset_id,error,plist_id)
          CALL h5dwrite_f(dset_id,h5t_native_double,dp4tmp(1,1,1,5), &
            data_dims,error)
          CALL h5sclose_f(dataspace,error)
          CALL h5dclose_f(dset_id,error)
        END IF

        IF (ident<=-1 .OR. out_prop(idx_lz)) THEN
          CALL h5screate_simple_f(rank,dims,dataspace,error)
          CALL h5dcreate_f(file_id,'lz',h5t_native_double,dataspace, &
            dset_id,error,plist_id)
          CALL h5dwrite_f(dset_id,h5t_native_double,dp4tmp(1,1,1,6), &
            data_dims,error)
          CALL h5sclose_f(dataspace,error)
          CALL h5dclose_f(dset_id,error)
        END IF

        DO i3 = 1, k0
          DO i2 = 1, j0
            DO i1 = 1, i0
              i = uindex(i1,i2,i3)
              dp4tmp(i1,i2,i3,1) = propunit(i,idx_comp,ismpl)
              dp4tmp(i1,i2,i3,2) = propunit(i,idx_q,ismpl)
              dp4tmp(i1,i2,i3,3) = propunit(i,idx_rc,ismpl)
              dp4tmp(i1,i2,i3,4) = propunit(i,idx_df,ismpl)
              dp4tmp(i1,i2,i3,5) = propunit(i,idx_ec,ismpl)
              dp4tmp(i1,i2,i3,6) = propunit(i,idx_lc,ismpl)
            END DO
          END DO
        END DO

        IF (ident<=-1 .OR. out_prop(idx_comp)) THEN
          CALL h5screate_simple_f(rank,dims,dataspace,error)
          CALL h5dcreate_f(file_id,'comp',h5t_native_double,dataspace, &
            dset_id,error,plist_id)
          CALL h5dwrite_f(dset_id,h5t_native_double,dp4tmp(1,1,1,1), &
            data_dims,error)
          CALL h5sclose_f(dataspace,error)
          CALL h5dclose_f(dset_id,error)
        END IF

        IF (ident<=-1 .OR. out_prop(idx_q)) THEN
          CALL h5screate_simple_f(rank,dims,dataspace,error)
          CALL h5dcreate_f(file_id,'q',h5t_native_double,dataspace, &
            dset_id,error,plist_id)
          CALL h5dwrite_f(dset_id,h5t_native_double,dp4tmp(1,1,1,2), &
            data_dims,error)
          CALL h5sclose_f(dataspace,error)
          CALL h5dclose_f(dset_id,error)
        END IF

        IF (ident<=-1 .OR. out_prop(idx_rc)) THEN
          CALL h5screate_simple_f(rank,dims,dataspace,error)
          CALL h5dcreate_f(file_id,'rc',h5t_native_double,dataspace, &
            dset_id,error,plist_id)
          CALL h5dwrite_f(dset_id,h5t_native_double,dp4tmp(1,1,1,3), &
            data_dims,error)
          CALL h5sclose_f(dataspace,error)
          CALL h5dclose_f(dset_id,error)
        END IF

        IF (ident<=-1 .OR. out_prop(idx_df)) THEN
          CALL h5screate_simple_f(rank,dims,dataspace,error)
          CALL h5dcreate_f(file_id,'df',h5t_native_double,dataspace, &
            dset_id,error,plist_id)
          CALL h5dwrite_f(dset_id,h5t_native_double,dp4tmp(1,1,1,4), &
            data_dims,error)
          CALL h5sclose_f(dataspace,error)
          CALL h5dclose_f(dset_id,error)
        END IF

        IF (ident<=-1 .OR. out_prop(idx_ec)) THEN
          CALL h5screate_simple_f(rank,dims,dataspace,error)
          CALL h5dcreate_f(file_id,'ec',h5t_native_double,dataspace, &
            dset_id,error,plist_id)
          CALL h5dwrite_f(dset_id,h5t_native_double,dp4tmp(1,1,1,5), &
            data_dims,error)
          CALL h5sclose_f(dataspace,error)
          CALL h5dclose_f(dset_id,error)
        END IF

        IF (ident<=-1 .OR. out_prop(idx_lc)) THEN
          CALL h5screate_simple_f(rank,dims,dataspace,error)
          CALL h5dcreate_f(file_id,'lc',h5t_native_double,dataspace, &
            dset_id,error,plist_id)
          CALL h5dwrite_f(dset_id,h5t_native_double,dp4tmp(1,1,1,6), &
            data_dims,error)
          CALL h5sclose_f(dataspace,error)
          CALL h5dclose_f(dset_id,error)
        END IF


        DO i3 = 1, k0
          DO i2 = 1, j0
            DO i1 = 1, i0
              dp4tmp(i1,i2,i3,1) = rhof(i1,i2,i3,ismpl)
              dp4tmp(i1,i2,i3,2) = visf(i1,i2,i3,ismpl)
            END DO
          END DO
        END DO
        IF (ident<=-1 .OR. out_ijk(cout_rhof)) THEN
          CALL h5screate_simple_f(rank,dims,dataspace,error)
          CALL h5dcreate_f(file_id,'rhof',h5t_native_double,dataspace, &
            dset_id,error,plist_id)
          CALL h5dwrite_f(dset_id,h5t_native_double,dp4tmp(1,1,1,1), &
            data_dims,error)
          CALL h5sclose_f(dataspace,error)
          CALL h5dclose_f(dset_id,error)
        END IF

        IF (ident<=-1 .OR. out_ijk(cout_visf)) THEN
          CALL h5screate_simple_f(rank,dims,dataspace,error)
          CALL h5dcreate_f(file_id,'visf',h5t_native_double,dataspace, &
            dset_id,error,plist_id)
          CALL h5dwrite_f(dset_id,h5t_native_double,dp4tmp(1,1,1,2), &
            data_dims,error)
          CALL h5sclose_f(dataspace,error)
          CALL h5dclose_f(dset_id,error)
        END IF

        dims(1) = i0
        dims(2) = 1
        dims(3) = 1
        data_dims(1) = i0
        data_dims(2) = 1
        data_dims(3) = 1
        chunk_dims(1) = min(8000,i0)
        chunk_dims(2) = 1
        chunk_dims(3) = 1
!     for compression only
        CALL h5pclose_f(plist_id,error)
        CALL h5pcreate_f(h5p_dataset_create_f,plist_id,error)
        CALL h5pset_chunk_f(plist_id,rank,chunk_dims,error)
        CALL h5pset_deflate_f(plist_id,gzlevel,error)

        IF (out_ijk(cout_i)) THEN
!     Create data space for the dataset.
          CALL h5screate_simple_f(rank,dims,dataspace,error)
!     Create dataset "A" inside file "f_name".
          CALL h5dcreate_f(file_id,'delx',h5t_native_double,dataspace, &
               dset_id,error,plist_id)
!     Write 'A' to the dataset
          CALL h5dwrite_f(dset_id,h5t_native_double,delx,data_dims, &
               error)
!     Close file, dataset and dataspace identifiers.
          CALL h5sclose_f(dataspace,error)
          CALL h5dclose_f(dset_id,error)
        END IF

        dims(1) = 1
        dims(2) = j0
        dims(3) = 1
        data_dims(1) = 1
        data_dims(2) = j0
        data_dims(3) = 1
        chunk_dims(1) = 1
        chunk_dims(2) = min(8000,j0)
        chunk_dims(3) = 1
!     for compression only
        CALL h5pclose_f(plist_id,error)
        CALL h5pcreate_f(h5p_dataset_create_f,plist_id,error)
        CALL h5pset_chunk_f(plist_id,rank,chunk_dims,error)
        CALL h5pset_deflate_f(plist_id,gzlevel,error)

        IF (out_ijk(cout_j)) THEN
          CALL h5screate_simple_f(rank,dims,dataspace,error)
          CALL h5dcreate_f(file_id,'dely',h5t_native_double,dataspace, &
               dset_id,error,plist_id)
          CALL h5dwrite_f(dset_id,h5t_native_double,dely,data_dims, &
               error)
          CALL h5sclose_f(dataspace,error)
          CALL h5dclose_f(dset_id,error)
        END IF

        dims(1) = 1
        dims(2) = 1
        dims(3) = k0
        data_dims(1) = 1
        data_dims(2) = 1
        data_dims(3) = k0
        chunk_dims(1) = 1
        chunk_dims(2) = 1
        chunk_dims(3) = min(8000,k0)
!     for compression only
        CALL h5pclose_f(plist_id,error)
        CALL h5pcreate_f(h5p_dataset_create_f,plist_id,error)
        CALL h5pset_chunk_f(plist_id,rank,chunk_dims,error)
        CALL h5pset_deflate_f(plist_id,gzlevel,error)

        IF (out_ijk(cout_k)) THEN
          CALL h5screate_simple_f(rank,dims,dataspace,error)
          CALL h5dcreate_f(file_id,'delz',h5t_native_double,dataspace, &
               dset_id,error,plist_id)
          CALL h5dwrite_f(dset_id,h5t_native_double,delz,data_dims, &
               error)
          CALL h5sclose_f(dataspace,error)
          CALL h5dclose_f(dset_id,error)
        END IF

        DEALLOCATE(inttmp)
        DEALLOCATE(dp4tmp)
!   ------------------------------------------------------------------
!     full tables of boundary conditions

        IF (nbc_data>0) THEN

          out_bc(1) = out_prop(idx_hbc)
          out_bc(2) = out_prop(idx_tbc)
          out_bc(3) = out_prop(idx_cbc)

          DO j = 1, 3
            if (out_bc(j)) then
            DO i = 1, 2
              c_name = pv_name(j) // '_' // c_bt(i)
!              count the number of bc with this type (i,j)
              k = 0
              DO i1 = 1, nbc_data
                IF (ibc_data(i1,cbc_pv)==j .AND. &
                  ibc_data(i1,cbc_bt)==i) k = k + 1
              END DO
              ALLOCATE(dp3tmp(k,ndbc,1))
              ALLOCATE(inttmp(k,nibc,1))
!              copy the counted bc
              k = 0
              DO i1 = 1, nbc_data
                IF (ibc_data(i1,cbc_pv)==j .AND. &
                    ibc_data(i1,cbc_bt)==i) THEN
                  k = k + 1
                  DO i2 = 1, nibc
                    inttmp(k,i2,1) = ibc_data(i1,i2)
                  END DO
                  DO i2 = 1, ndbc
                    dp3tmp(k,i2,1) = dbc_data(i1,i2,ismpl)
                  END DO
                END IF
              END DO

              IF (k>0) THEN
                dims(1) = k
                dims(2) = nibc
                dims(3) = 1
                data_dims(1) = k
                data_dims(2) = nibc
                data_dims(3) = 1
!                 for compression only
                chunk_dims(1) = min(2000,k)
                chunk_dims(2) = 1
                chunk_dims(3) = 1
                CALL h5pclose_f(plist_id,error)
                CALL h5pcreate_f(h5p_dataset_create_f,plist_id,error)
                CALL h5pset_chunk_f(plist_id,rank,chunk_dims,error)
                CALL h5pset_deflate_f(plist_id,gzlevel,error)

                CALL h5screate_simple_f(rank,dims,dataspace,error)
                CALL h5dcreate_f(file_id,'i'//c_name, &
                  h5t_native_integer,dataspace,dset_id,error,plist_id)
                CALL h5dwrite_f(dset_id,h5t_native_integer,inttmp, &
                  data_dims,error)
                CALL h5sclose_f(dataspace,error)
                CALL h5dclose_f(dset_id,error)

                dims(2) = ndbc
                data_dims(2) = ndbc
!                 for compression only
                chunk_dims(1) = min(4000,k)
                CALL h5pclose_f(plist_id,error)
                CALL h5pcreate_f(h5p_dataset_create_f,plist_id,error)
                CALL h5pset_chunk_f(plist_id,rank,chunk_dims,error)
                CALL h5pset_deflate_f(plist_id,gzlevel,error)

                CALL h5screate_simple_f(rank,dims,dataspace,error)
                CALL h5dcreate_f(file_id,c_name,h5t_native_double, &
                  dataspace,dset_id,error,plist_id)
                CALL h5dwrite_f(dset_id,h5t_native_double,dp3tmp, &
                  data_dims,error)
                CALL h5sclose_f(dataspace,error)
                CALL h5dclose_f(dset_id,error)
              END IF
              DEALLOCATE(inttmp)
              DEALLOCATE(dp3tmp)
            END DO
          end if
          END DO
        END IF

!   ------------------------------------------------------------------

!     close interface
        CALL h5pclose_f(plist_id,error)
        CALL h5fclose_f(file_id,error)

!       write additional property array
        IF (index(project_sfx(ismpl),'_postcomp')>0) THEN
          CALL write2_hdf5(nunits,nprop,propunit(1,1,ismpl),'props_full',filename)
        END IF

!     close FORTRAN interface
        CALL h5close_f(error)

#ifdef fOMP
!$OMP end critical
#endif

#endif
        RETURN
      END
