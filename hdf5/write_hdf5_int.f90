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

!>    @brief writes an 3-dimensional array to a existing HDF5 file
!>    @param[in] NI 1.dimension
!>    @param[in] NJ 2.dimension
!>    @param[in] NK 3.dimension
!>    @param[in] A integer array
!>    @param[in] A_name name of the array
!>    @param[in] f_name HDF5 file name
!>    @details
!>    Used only for restart output.
      SUBROUTINE write_hdf5_int(ni,nj,nk,a,a_name,f_name)
#ifndef noHDF
        USE hdf5
        use mod_hdf5_vars, only: error
#endif
        IMPLICIT NONE

!      arrayname and filename
        character (len=*) :: a_name, f_name

        INTEGER ni, nj, nk, i, j, k
        INTEGER a(ni,nj,nk)

#ifndef noHDF
!      File identifiers
        INTEGER (hid_t) file_id

!      Dataset identifier
        INTEGER (hid_t) dset_id

!      Data space identifier
        INTEGER (hid_t) dataspace
#endif

!      Data buffers
#ifdef HDF6432
        INTEGER (kind=4) :: rank
#else
#ifdef HDF64
        INTEGER (kind=8) :: rank
#else
        INTEGER rank
#endif
#endif
        PARAMETER (rank=3)

#ifndef noHDF
        INTEGER (hsize_t) dims(rank)
        INTEGER (hsize_t) data_dims(7)

!     for chunk size and compression
        INTEGER (hid_t) plist_id
        INTEGER (hsize_t) chunk_dims(rank)
#endif

!     gzip compression level
#ifdef HDF6432
        INTEGER (kind=4) :: gzlevel
#else
#ifdef HDF64
        INTEGER (kind=8) :: gzlevel
#else

        INTEGER gzlevel
#endif
#endif
        PARAMETER (gzlevel=9)

!     neeed for 32Bit HDF5 library
#ifdef HDF6432
        INTEGER (kind=4), allocatable :: inttmp(:,:,:)
#else
#ifdef HDF64
        INTEGER (kind=8), allocatable :: inttmp(:,:,:)
#else

        INTEGER, ALLOCATABLE :: inttmp(:,:,:)
#endif
#endif

#ifndef noHDF

        ALLOCATE(inttmp(ni,nj,nk))

!      Initialize FORTRAN interface.
!aw      Call h5open_f(error)

        dims(1) = ni
        dims(2) = nj
        dims(3) = nk
        data_dims(1) = ni
        data_dims(2) = nj
        data_dims(3) = nk
        chunk_dims(1) = min(ni,20)
        chunk_dims(2) = min(nj,20)
        chunk_dims(3) = min(nk,20)

!      Reopen hdf5 file.
        CALL h5fopen_f(f_name,h5f_acc_rdwr_f,file_id,error)


!      Create data space for the dataset.
        CALL h5screate_simple_f(rank,dims,dataspace,error)

!      for compression only
        CALL h5pcreate_f(h5p_dataset_create_f,plist_id,error)
        CALL h5pset_chunk_f(plist_id,rank,chunk_dims,error)
        CALL h5pset_deflate_f(plist_id,gzlevel,error)

!      Create dataset "A" inside file "f_name".
        CALL h5dcreate_f(file_id,a_name,h5t_native_integer,dataspace, &
          dset_id,error,plist_id)

        DO k = 1, nk
          DO j = 1, nj
            DO i = 1, ni
              inttmp(i,j,k) = a(i,j,k)
            END DO
          END DO
        END DO
!      Write 'A' to the dataset
        CALL h5dwrite_f(dset_id,h5t_native_integer,inttmp,data_dims, &
          error)

!      Close file, dataset, propertylist and dataspace identifiers.
        CALL h5pclose_f(plist_id,error)
        CALL h5sclose_f(dataspace,error)
        CALL h5dclose_f(dset_id,error)
        CALL h5fclose_f(file_id,error)

!      Close FORTRAN interface.
!aw      Call h5Close_f(error)

        DEALLOCATE(inttmp)

#else
        WRITE(*,*) 'error: HDF5 support was not compiled in'
        STOP
#endif
        RETURN
      END
