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

!>    @brief open hdf5-file
!>    @param[in] f_name hdf5 file name
!>    @details
!>    Opening hdf5 files, used in read-routines.
!> Note: To be able to use input file parsing with hdf5, the
!> hdf5-input-files have to be generated using the script:
!> `convert_to_hdf5.py`. This script can be found in the repository
!> `SHEMAT-Suite_Scripts` under
!> `python/preprocessing/convert_to_hdf5.py`.
      SUBROUTINE open_hdf5(f_name)
#ifndef noHDF
        USE hdf5
        use mod_input_file_parser_hdf5
        use mod_hdf5_vars, only: default_hdf_file, file_id, error
#endif
        IMPLICIT NONE

!      arrayname and filename
        character (len=*) :: f_name


#ifndef noHDF

!     Initialize FORTRAN interface.
        if (.not. h5parse_hdf5_environment) then
           CALL h5open_f(error)
        end if

!     Reopen hdf5 file.
        default_hdf_file = f_name
        IF (f_name/=' ') THEN
          CALL h5fopen_f(f_name,h5f_acc_rdwr_f,file_id,error)
        END IF

#endif
        RETURN
      END

