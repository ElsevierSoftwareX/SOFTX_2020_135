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

!>    @brief extract extern file name and read the double precision array 'A'
!>    @param[in] NI i-direction array dimension
!>    @param[in] NJ j-direction array dimension
!>    @param[in] NK k-direction array dimension
!>    @param[in] A_name array name (ascii)
!>    @param[in] name_link line for extracting
!>    @param[out] A array with all data
!>    @return "true" it was not an external file, FALSE: data from external file readed
      LOGICAL FUNCTION no_ext_link(ni,nj,nk,a,a_name,name_link)
        use mod_genrlc
        IMPLICIT NONE
        INTEGER i1, i2
        character (len=80) :: defaultname
!     arrayname and filename
        character (len=*) :: a_name
!
        INTEGER ni, nj, nk
        INTEGER j, k
        DOUBLE PRECISION a(ni,nj,nk)
!
        character (len=*) :: name_link
!
        INTEGER sfirst, end_sec
        EXTERNAL sfirst
!
        LOGICAL is_txt

!     default HDF5 input file name
        CALL chln(project,i1,i2)
!
        is_txt = .FALSE.
        CALL get_arg('HDF5',name_link,j,k)
        defaultname = project(i1:i2) // '.h5'
        IF (j<1 .OR. k<j) THEN
          is_txt = .TRUE.
          CALL get_arg('TXT',name_link,j,k)
          defaultname = project(i1:i2) // '.txt'
          IF (j<1 .OR. k<j) THEN
            no_ext_link = .TRUE.
            RETURN
          END IF
        END IF
!     end of the section name
        end_sec = sfirst(name_link)
!
!     found extern file declaration
        IF (name_link(j:k)=='default') THEN
          IF (is_txt) THEN
            CALL read_open_txt(ni,nj,nk,a,name_link(1:end_sec), &
              defaultname)
          ELSE
            CALL read_hdf5(ni,nj,nk,a,a_name,defaultname)
          END IF
        ELSE
          IF (is_txt) THEN
            CALL read_open_txt(ni,nj,nk,a,name_link(1:end_sec), &
              name_link(j:k))
          ELSE
            CALL read_hdf5(ni,nj,nk,a,a_name,name_link(j:k))
          END IF
        END IF
        no_ext_link = .FALSE.
!
        RETURN
      END

!>    @brief extract extern file name and read the integer array 'A'
!>    @param[in] NI i-direction array dimension
!>    @param[in] NJ j-direction array dimension
!>    @param[in] NK k-direction array dimension
!>    @param[in] A_name array name (ascii)
!>    @param[in] name_link line for extracting
!>    @param[out] A array with all data
!>    @return "true" it was not an external file, FALSE: data from external file readed
      LOGICAL FUNCTION no_ext_link_int(ni,nj,nk,a,a_name,name_link)
        use mod_genrlc
        IMPLICIT NONE
        INTEGER i1, i2
        character (len=80) :: defaultname
!     arrayname and filename
        character (len=*) :: a_name
!
        INTEGER ni, nj, nk
        INTEGER j, k
        INTEGER a(ni,nj,nk)
!
        character (len=*) :: name_link
!
        INTEGER sfirst, end_sec
        EXTERNAL sfirst
!
        LOGICAL is_txt

!     default HDF5 input file name
        CALL chln(project,i1,i2)
        defaultname = project(i1:i2) // '.h5'
!
        is_txt = .FALSE.
        CALL get_arg('HDF5',name_link,j,k)
        IF (j<1 .OR. k<j) THEN
          is_txt = .TRUE.
          CALL get_arg('TXT',name_link,j,k)
          IF (j<1 .OR. k<j) THEN
            no_ext_link_int = .TRUE.
            RETURN
          END IF
        END IF
!     end of the section name
        end_sec = sfirst(name_link)
!
!     found extern file declaration
        IF (name_link(j:k)=='default') THEN
          IF (is_txt) THEN
            CALL read_open_txt_int(ni,nj,nk,a,name_link(1:end_sec), &
              defaultname)
          ELSE
            CALL read_hdf5_int(ni,nj,nk,a,a_name,defaultname)
          END IF
        ELSE
          IF (is_txt) THEN
            CALL read_open_txt_int(ni,nj,nk,a,name_link(1:end_sec), &
              name_link(j:k))
          ELSE
            CALL read_hdf5_int(ni,nj,nk,a,a_name,name_link(j:k))
          END IF
        END IF
        no_ext_link_int = .FALSE.
!
        RETURN
      END

!>    @brief read double precision data from file
!>    @param[in] NI i-direction array dimension
!>    @param[in] NJ j-direction array dimension
!>    @param[in] NK k-direction array dimension
!>    @param[in] A_name array name (ascii)
!>    @param[in] f_name file name
!>    @param[out] A array with all data
      SUBROUTINE read_open_txt(ni,nj,nk,a,a_name,f_name)
        IMPLICIT NONE
        character (len=*) :: a_name, f_name
        character (len=80) :: line
        INTEGER ni, nj, nk, i, j, k
        DOUBLE PRECISION a(ni,nj,nk)
        LOGICAL found
        EXTERNAL found

        WRITE(*,*) '  open TXT file: ', f_name
        OPEN(89,file=f_name,status='old',err=99)
        IF (found(89,a_name,line,.TRUE.)) THEN
          READ(89,*) (((a(i,j,k),i=1,ni),j=1,nj),k=1,nk)
        END IF
        CLOSE(89)
        RETURN
!
99      WRITE(*,*) 'error: can not open file "', f_name, '" !'
        STOP
      END

!>    @brief read integer data from file
!>    @param[in] NI i-direction array dimension
!>    @param[in] NJ j-direction array dimension
!>    @param[in] NK k-direction array dimension
!>    @param[in] A_name array name (ascii)
!>    @param[in] f_name file name
!>    @param[out] A array with all data
      SUBROUTINE read_open_txt_int(ni,nj,nk,a,a_name,f_name)
        IMPLICIT NONE
        character (len=*) :: a_name, f_name
        character (len=80) :: line
        INTEGER ni, nj, nk, i, j, k
        INTEGER a(ni,nj,nk)
        LOGICAL found
        EXTERNAL found

        WRITE(*,*) '  open TXT file: ', f_name
        OPEN(89,file=f_name,status='old',err=99)
        IF (found(89,a_name,line,.TRUE.)) THEN
          READ(89,*) (((a(i,j,k),i=1,ni),j=1,nj),k=1,nk)
        END IF
        CLOSE(89)
        RETURN
!
99      WRITE(*,*) 'error: can not open file "', f_name, '" !'
        STOP
      END
