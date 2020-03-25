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

!>    @brief creates a HDF5 file
!>    @param[in] f_name : hdf5 file name
!>    @details
!>    create a new hdf5-file, used only for restart output.\n
      SUBROUTINE create_hdf5(f_name)
#ifndef noHDF
        USE hdf5
        use mod_hdf5_vars, only: error
#endif
        IMPLICIT NONE

!      arrayname and filename
        character (len=*) :: f_name

#ifndef noHDF
!      File identifiers
        INTEGER (hid_t) file_id
#endif

#ifndef noHDF

!      Initialize FORTRAN interface.
!aw      Call h5open_f(error)

!      Create a new file, later only open it for read and writes
!org      call h5fcreate_f(f_name, H5F_ACC_TRUNC_F, file_id, error)
        CALL h5fcreate_f(f_name,h5f_acc_trunc_f,file_id,error, &
          h5p_default_f,h5p_default_f)
        CALL h5fclose_f(file_id,error)

!      Close FORTRAN interface.
!aw      Call h5Close_f(error)

#else
        WRITE(*,*) 'error: HDF5 support was not compiled in'
        STOP
#endif
        RETURN
      END

