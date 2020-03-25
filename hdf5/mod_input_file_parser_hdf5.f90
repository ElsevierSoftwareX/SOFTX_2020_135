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

!>    @brief hdf5 input file parsing
!>    @details
!> Note: To be able to use input file parsing with hdf5, the
!> hdf5-input-files have to be generated using the script:
!> `convert_to_hdf5.py`. This script can be found in the repository
!> `SHEMAT-Suite_Scripts` under
!> `python/preprocessing/convert_to_hdf5.py`.
module mod_input_file_parser_hdf5
    ! HDF5 input data format parser
    ! Developed as part of the EoCoE project 2016
    ! author: Sebastian Luehrs, JSC, Forschungszentrum Juelich
#ifndef noHDF
    use hdf5
#ifdef fMPI
    use mpi
#endif
    implicit none
    ! globally available information, if HDF5 input file is used
    logical :: h5parse_use_hdf5_datafile = .false.
    logical :: h5parse_hdf5_environment = .false.
    integer(kind=HID_T) :: file_id
    private
    public :: h5parse_read_dimension_size_for_dataset, h5parse_read_2d_double_dataset, &
        & h5parse_read_1d_double_dataset, h5parse_use_hdf5_datafile, &
        & h5parse_open_datafile, h5parse_close_datafile, h5parse_check_attr_exist, &
        & h5parse_read_double_attribute, h5parse_read_integer_attribute, &
        & h5parse_check_dataset_exist, h5parse_read_2d_integer_dataset, &
        & h5parse_read_3d_double_dataset, h5parse_read_3d_integer_dataset, &
        & h5parse_hdf5_environment
contains

    subroutine h5parse_init()
        ! Initalize HDF5 Fortran interface
        integer(kind=4) :: hdferr
        call H5open_f(hdferr)
        h5parse_hdf5_environment = .true.
    end subroutine h5parse_init

    subroutine h5parse_open_file(filename)
        ! Open HDF5-datafile and store the filehandle
        ! @param[in] filename: Input filename
        character(len=80), intent(in) :: filename
        integer(kind=HID_T) :: plist_id
        integer(kind=4) :: hdferr
#ifdef fMPI
        !  set up file access propterty for parallel HDF5
        call H5Pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdferr)
        call H5Pset_fapl_mpio_f(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL,hdferr)
#else
        plist_id = H5P_DEFAULT_F
#endif
        ! create new file collectively and release property list identifier
        call H5Fopen_f(filename,H5F_ACC_RDONLY_F,file_id,hdferr,plist_id)
#ifdef fMPI
        call H5Pclose_f(plist_id,hdferr)
#endif
    end subroutine h5parse_open_file

    function h5parse_read_dimension_size_for_dataset(dataset_name,dimension) result(dimension_size)
        ! Return the dimension size for a specific dataset (given by name).
        ! If the rank of the choosen dataset is > 1, the specific dimension can be selected by number.
        ! @param[in] dataset_name: The name of the dataset
        ! @param[in] dimension: The dimension number (default: 1)
        ! @param[out] dimension_size: Size of the selected dimension in dataset
        character(len=*), intent(in) :: dataset_name
        integer, intent(in), optional :: dimension
        integer :: dimension_size
        integer(kind=HID_T) :: dataset_id
        integer(kind=HID_T) :: dataspace_id
        integer(kind=HSIZE_T), dimension(:), allocatable :: dimsr, maxdimsr
        integer(kind=4) :: rank
        integer(kind=4) :: hdferr
        call H5Dopen_f(file_id,dataset_name,dataset_id,hdferr)
        ! Open dataset
        call H5Dget_space_f(dataset_id, dataspace_id, hdferr)
        ! Get rank
        call H5Sget_simple_extent_ndims_f(dataspace_id, rank, hdferr)
        allocate(dimsr(rank))
        allocate(maxdimsr(rank))
        ! Get dimension size
        call H5Sget_simple_extent_dims_f(dataspace_id, dimsr, maxdimsr, hdferr)
        if (present(dimension)) then
            dimension_size = int(dimsr(dimension))
        else
            dimension_size = int(dimsr(1))
        end if
        deallocate(dimsr)
        deallocate(maxdimsr)
        call H5Sclose_f(dataspace_id, hdferr)
        call H5Dclose_f(dataset_id, hdferr)
    end function h5parse_read_dimension_size_for_dataset

    subroutine h5parse_read_3d_double_dataset(dataset_name,storage,data_shape_ext)
        ! Read a 3D double dataset
        ! @param[in] dataset_name: The name of the dataset
        ! @param[inout] storage: Storage array which will hold the the contents of the dataset
        ! @param[in] data_shape_ext: Selectable file data_shape (default: shape of storage)
        character(len=*), intent(in) :: dataset_name
        double precision, dimension(:,:,:), intent(inout) :: storage
        integer, optional, dimension(3), intent(in) :: data_shape_ext
        integer(kind=HSIZE_T), dimension(3) :: data_shape
        if (present(data_shape_ext)) then
            data_shape = data_shape_ext
        else
            data_shape = shape(storage)
        end if
        call h5parse_read_nd_double_dataset(dataset_name,storage,data_shape)
    end subroutine h5parse_read_3d_double_dataset

    subroutine h5parse_read_3d_integer_dataset(dataset_name,storage,data_shape_ext)
        ! Read a 3D integer dataset
        ! @param[in] dataset_name: The name of the dataset
        ! @param[inout] storage: Storage array which will hold the the contents of the dataset
        ! @param[in] data_shape_ext: Selectable file data_shape (default: shape of storage)
        character(len=*), intent(in) :: dataset_name
        integer, dimension(:,:,:), intent(inout) :: storage
        integer, optional, dimension(3), intent(in) :: data_shape_ext
        integer(kind=HSIZE_T), dimension(3) :: data_shape
        if (present(data_shape_ext)) then
            data_shape = data_shape_ext
        else
            data_shape = shape(storage)
        end if
        call h5parse_read_nd_integer_dataset(dataset_name,storage,data_shape)
    end subroutine h5parse_read_3d_integer_dataset

    subroutine h5parse_read_2d_double_dataset(dataset_name,storage,data_shape_ext)
        ! Read a 2D double dataset
        ! @param[in] dataset_name: The name of the dataset
        ! @param[inout] storage: Storage array which will hold the the contents of the dataset
        ! @param[in] data_shape_ext: Selectable file data_shape (default: shape of storage)
        character(len=*), intent(in) :: dataset_name
        double precision, dimension(:,:), intent(inout) :: storage
        integer, optional, dimension(2), intent(in) :: data_shape_ext
        integer(kind=HSIZE_T), dimension(2) :: data_shape
        if (present(data_shape_ext)) then
            data_shape = data_shape_ext
        else
            data_shape = shape(storage)
        end if
        call h5parse_read_nd_double_dataset(dataset_name,storage,data_shape)
    end subroutine h5parse_read_2d_double_dataset

    subroutine h5parse_read_2d_integer_dataset(dataset_name,storage,data_shape_ext)
        ! Read a 2D integer dataset
        ! @param[in] dataset_name: The name of the dataset
        ! @param[inout] storage: Storage array which will hold the the contents of the dataset
        ! @param[in] data_shape_ext: Selectable file data_shape (default: shape of storage)
        character(len=*), intent(in) :: dataset_name
        integer, dimension(:,:), intent(inout) :: storage
        integer, optional, dimension(2), intent(in) :: data_shape_ext
        integer(kind=HSIZE_T), dimension(2) :: data_shape
        if (present(data_shape_ext)) then
            data_shape = data_shape_ext
        else
            data_shape = shape(storage)
        end if
        call h5parse_read_nd_integer_dataset(dataset_name,storage,data_shape)
    end subroutine h5parse_read_2d_integer_dataset

    subroutine h5parse_read_1d_double_dataset(dataset_name,storage,data_shape_ext)
        ! Read a 1D double dataset
        ! @param[in] dataset_name: The name of the dataset
        ! @param[inout] storage: Storage array which will hold the the contents of the dataset
        ! @param[in] data_shape_ext: Selectable file data_shape (default: shape of storage)
        character(len=*), intent(in) :: dataset_name
        double precision, dimension(:), intent(inout) :: storage
        integer, optional, intent(in) :: data_shape_ext
        integer(kind=HSIZE_T), dimension(1) :: data_shape
        if(present(data_shape_ext))then
            data_shape = data_shape_ext
        else
            data_shape = shape(storage)
        end if
        call h5parse_read_nd_double_dataset(dataset_name,storage,data_shape)
    end subroutine h5parse_read_1d_double_dataset

    subroutine h5parse_read_nd_double_dataset(dataset_name,storage,data_shape)
        ! Read a nD double dataset
        ! @param[in] dataset_name: The name of the dataset
        ! @param[inout] storage: Storage array which will hold the the contents of the dataset
        ! @param[in] data_shape_ext: Selectable file data_shape (default: shape of storage)
        character(len=*), intent(in) :: dataset_name
        double precision, dimension(*), intent(inout) :: storage
        integer(kind=HSIZE_T), dimension(:), intent(in) :: data_shape
        integer(kind=HSIZE_T), dimension(:), allocatable :: offset
        integer(kind=HID_T) :: dataset_id
        integer(kind=HID_T) :: dataspace_id, memspace_id
        integer(kind=4) :: hdferr
        allocate(offset(size(data_shape)))
        offset = 0
        call H5Dopen_f(file_id, dataset_name, dataset_id, hdferr)
        ! Open dataset
        call H5Dget_space_f(dataset_id, dataspace_id, hdferr)
        ! Select data_shape in dataset
        call H5Sselect_hyperslab_f(dataspace_id, H5S_SELECT_SET_F, offset, data_shape, hdferr)
        ! Create memory dataspace
        call H5Screate_simple_f(size(data_shape), data_shape, memspace_id, hdferr)
        ! Read the data
        call H5Dread_f(dataset_id, H5T_NATIVE_DOUBLE, storage, data_shape, hdferr, &
            & memspace_id, dataspace_id)
        call H5Sclose_f(memspace_id, hdferr)
        call H5Sclose_f(dataspace_id, hdferr)
        call H5Dclose_f(dataset_id, hdferr)
        deallocate(offset)
    end subroutine h5parse_read_nd_double_dataset

    subroutine h5parse_read_nd_integer_dataset(dataset_name,storage,data_shape)
        ! Read a nD integer dataset
        ! @param[in] dataset_name: The name of the dataset
        ! @param[inout] storage: Storage array which will hold the the contents of the dataset
        ! @param[in] data_shape_ext: Selectable file data_shape (default: shape of storage)
        character(len=*), intent(in) :: dataset_name
        integer, dimension(*), intent(inout) :: storage
        integer(kind=HSIZE_T), dimension(:), intent(in) :: data_shape
        integer(kind=HSIZE_T), dimension(:), allocatable :: offset
        integer(kind=HID_T) :: dataset_id
        integer(kind=HID_T) :: dataspace_id, memspace_id
        integer(kind=4) :: hdferr
        allocate(offset(size(data_shape)))
        offset = 0
        call H5Dopen_f(file_id,dataset_name,dataset_id,hdferr)
        ! Open dataset
        call H5Dget_space_f(dataset_id, dataspace_id, hdferr)
        ! Select data_shape in dataset
        call H5Sselect_hyperslab_f(dataspace_id, H5S_SELECT_SET_F, offset, data_shape, hdferr)
        ! Create memory dataspace
        call H5Screate_simple_f(size(data_shape), data_shape, memspace_id, hdferr)
        ! Read the data
        call H5Dread_f(dataset_id, H5T_NATIVE_INTEGER, storage, data_shape, hdferr, &
            & memspace_id, dataspace_id)
        call H5Sclose_f(memspace_id, hdferr)
        call H5Sclose_f(dataspace_id, hdferr)
        call H5Dclose_f(dataset_id, hdferr)
        deallocate(offset)
    end subroutine h5parse_read_nd_integer_dataset

    function h5parse_check_attr_exist(attr_name,obj_name_ext) result(attr_exists)
        ! Check if a specifc attribute, given by name, exists
        ! @param[in] attr_name: The name of the attribute
        ! @param[in] obj_name_ext: The relative context the attribute belongs to (default: root)
        ! @param[out] attr_exists: Attribute exists?
        character(len=*), intent(in) :: attr_name
        character(len=*), optional, intent(in) :: obj_name_ext
        character(len=255) :: obj_name
        logical(kind=4) :: attr_exists
        integer(kind=4) :: hdferr
        if(present(obj_name_ext)) then
            obj_name = obj_name_ext
        else
            obj_name = "/"
        end if
        call h5aexists_by_name_f(file_id, obj_name, attr_name, attr_exists, hdferr)
    end function h5parse_check_attr_exist

    function h5parse_check_dataset_exist(dataset_name) result(dataset_exists)
        ! Check if a specifc dataset, given by name, exists
        ! @param[in] dataset_name: The name of the dataset
        ! @param[out] dataset_exists: Dataset exists?
        character(len=*), intent(in) :: dataset_name
        logical(kind=4) :: dataset_exists
        integer(kind=4) :: hdferr
        call h5lexists_f(file_id, dataset_name, dataset_exists, hdferr)
    end function h5parse_check_dataset_exist

    subroutine h5parse_read_double_attribute(attr_name,storage,obj_name_ext)
        ! Read the content of a double attribute
        ! @param[in] attr_name: The name of the attribute
        ! @param[in] obj_name_ext: The relative context the attribute belongs to (default: root)
        ! @param[inout] storage: Storage variable which will hold the the content of the attribute
        character(len=*), intent(in) :: attr_name
        double precision, intent(inout) :: storage
        character(len=*), optional, intent(in) :: obj_name_ext
        character(len=255) :: obj_name
        integer(kind=4) :: hdferr
        integer(kind=HID_T) :: attr_id
        integer(kind=HSIZE_T), dimension(1) :: dims = 0
        if(present(obj_name_ext)) then
            obj_name = obj_name_ext
        else
            obj_name = "/"
        end if
        call h5aopen_by_name_f(file_id, obj_name, attr_name, attr_id, hdferr)
        call h5aread_f(attr_id, H5T_NATIVE_DOUBLE, storage, dims, hdferr)
        call h5aclose_f(attr_id, hdferr)
    end subroutine h5parse_read_double_attribute

    subroutine h5parse_read_integer_attribute(attr_name,storage,obj_name_ext)
        ! Read the content of a integer attribute
        ! @param[in] attr_name: The name of the attribute
        ! @param[in] obj_name_ext: The relative context the attribute belongs to (default: root)
        ! @param[inout] storage: Storage variable which will hold the the content of the attribute
        character(len=*), intent(in) :: attr_name
        integer, intent(inout) :: storage
        character(len=*), optional, intent(in) :: obj_name_ext
        character(len=255) :: obj_name
        integer(kind=4) :: hdferr
        integer(kind=HID_T) :: attr_id
        integer(kind=HSIZE_T), dimension(1) :: dims = 0
        if(present(obj_name_ext)) then
            obj_name = obj_name_ext
        else
            obj_name = "/"
        end if
        call h5aopen_by_name_f(file_id, obj_name, attr_name, attr_id, hdferr)
        call h5aread_f(attr_id, H5T_NATIVE_INTEGER, storage, dims, hdferr)
        call h5aclose_f(attr_id, hdferr)
    end subroutine h5parse_read_integer_attribute

    subroutine h5parse_close_file()
        ! Close the HDF5 internal file handle
        integer(kind=4) :: hdferr
        call H5Fclose_f(file_id, hdferr)
    end subroutine h5parse_close_file

    subroutine h5parse_finalize()
        ! Close the HDF5 Fortran interface
        integer(kind=4) :: hdferr
        h5parse_hdf5_environment = .false.
        call H5close_f(hdferr)
    end subroutine h5parse_finalize

    subroutine h5parse_open_datafile(filename)
        ! Find the HDF5 datafile filenam in the standard ASCII input file
        ! @param[in] filename: Name of the standard ASCII input file
        use mod_genrl
        character(len=*), intent(in) :: filename
        character(len=255) :: hdf5_data_filename
        character(len=80) :: line
        logical :: found
        integer :: lblank
        open(79,file=filename,status='old')
        call h5parse_init()
        ! Search for the HDF5 datafile filename
        if (found(79,key_char//' h5parse data file',line,.false.)) then
            read(79,'(1A)',err=200,end=200) hdf5_data_filename
            write(*,*) ' '
            write(*,*) '  reading model input data:'
            write(*,*) '    from file "', hdf5_data_filename(&
                &:lblank(hdf5_data_filename)),'"'
            write(*,*) ' '
            ! Use HDF5 input file parser
            h5parse_use_hdf5_datafile = .true.
            ! Open the datafile
            call h5parse_open_file(hdf5_data_filename)
        else
            ! Do not use HDF5 input file parser (use old input format)
            h5parse_use_hdf5_datafile = .false.
        end if
        close(79)
        return
200     write(*,'(1A)') 'error: can not read "h5parse data file"!'
        stop
    end subroutine h5parse_open_datafile

    subroutine h5parse_close_datafile()
        ! Close the HDF5 internal file handle and the HDF5 interface
        if (h5parse_use_hdf5_datafile) then
            call h5parse_close_file()
            h5parse_use_hdf5_datafile = .false.
        end if
        call h5parse_finalize()
    end subroutine h5parse_close_datafile
#endif

end module mod_input_file_parser_hdf5
