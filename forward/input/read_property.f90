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

!>    @brief read and overwrites a specific rock property (only full dimension is supported)
!>    @param[in] filename "model" (or other) file name
!>    @param[in] prop_idx property index number
!>    @param[in] ismpl local sample index
      SUBROUTINE read_property(filename,prop_idx,ismpl)
        use arrays
        use mod_genrl
        IMPLICIT NONE
!
        INTEGER ismpl, i, j, k
        character (len=80) ::  filename
        character (len=80) :: line
        character (len=10) :: pname
        INTEGER prop_idx
!
        LOGICAL found, no_ext_link
        EXTERNAL found, no_ext_link
        INTRINSIC trim, adjustl

!
        IF (prop_idx<firstidx .OR. prop_idx>lastidx) THEN
          WRITE(*,'(1A,3(1I2,1A))') 'error: property index number ',prop_idx,' out of range (',firstidx,',',lastidx,')!'
          STOP
        END IF
!
        pname = trim(adjustl(properties(prop_idx)))
        WRITE(*,*) ' [R] : '//trim(doc_properties(prop_idx))//' ('//trim(pname)//') from file "', trim(filename),'"'
!     read file
        OPEN(79,file=filename,status='old')
!     init HDF5 support, when available
        CALL open_hdf5(' ')
!
        IF (found(79,key_char//' '//trim(pname),line,.TRUE.)) THEN
          IF (no_ext_link(i0,j0,k0,x(1,1,1,ismpl),trim(pname),line)) &
            READ(79,*,err=200,end=200) (((x(i,j,k,ismpl),i=1,i0),j=1,J0),k=1,K0)
          DO k = 1, K0
            DO j = 1, J0
              DO i = 1, I0
                propunit(uindex(i,j,k),prop_idx,ismpl) = x(i,j,k,ismpl)
              END DO
            END DO
          END DO
        END IF
!
!     finish HDF5 support, when available
        CALL close_hdf5()
!     close project config file
        CLOSE(79)
        RETURN

!     error handler
200     WRITE(*,'(3A)') 'error: to few values in section "',trim(pname),'"!'
        STOP
      END
