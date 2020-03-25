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

!>    @brief shared hdf5 variables
module mod_hdf5_vars

#ifndef noHDF

  use hdf5, only: hid_t
  
  !> @brief Name of opened hdf file
  !> @details
  !> Name of opened hdf file. \n
  !> Saves the name of the hdf file opened in open_hdf.
  character (len=256) :: default_hdf_file

  !> @brief HDF5 file id
  !> @details
  !> HDF5 file id. \n
  !> HDF5 file id.
  integer (kind=hid_t) :: file_id


  !> @brief HDF5 Error flag
  !> @details
  !> HDF5 Error flag. \n
  !> HDF5 Error flag to check operation success.
#ifdef HDF6432
  integer (kind=4) :: error
#else
#ifdef HDF64
  integer (kind=8) :: error
#else
  integer :: error
#endif
#endif

#endif

end module mod_hdf5_vars
