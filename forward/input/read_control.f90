
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

!>    @brief read control sections for the external file-infrastructure
!>    @param[in] filename
!>    @param[in] ismpl local sample index
      SUBROUTINE read_control(filename,ismpl)
        use mod_genrl
        use mod_genrlc
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j

        character (len=*) :: filename
        character (len=5000) :: line
        LOGICAL found
        EXTERNAL found
        INTRINSIC trim

!       External input file switch not set
        if (.not. read_external_input) then
          filename_data = filename
          filename_simul = filename
          filename_enkf = filename
          filename_inverse = filename
          return
        end if

!       open file
        OPEN(79,file=filename,status='old')
        WRITE(*,*) ' '
        WRITE(*,*) '  reading external file control'
        WRITE(*,*) ' '
!
        filename_data = filename
        IF (runmode>=1) THEN
          IF (found(79,key_char//' data: external file',line,.FALSE.)) THEN
            CALL get_arg('file',line,i,j)
            IF (i<1 .OR. j<i) THEN
              READ(79,*,err=100,end=100) filename_data
            ELSE
              READ(line(i:j),*) filename_data
            END IF
            WRITE(*,*) ' [R] : external file "'//trim(filename_data)//'" for DATA'
!          ELSE
!            WRITE(*,*) ' <D> : no external file for DATA'
          END IF
        END IF
!
        filename_simul = filename
        IF (def_binary=='simul') THEN
          IF (found(79,key_char//' simul: external file',line,.FALSE.)) THEN
            CALL get_arg('file',line,i,j)
            IF (i<1 .OR. j<i) THEN
              READ(79,*,err=101,end=101) filename_simul
            ELSE
              READ(line(i:j),*) filename_simul
            END IF
            WRITE(*,*) ' [R] : external file "'//trim(filename_simul)//'" for SIMUL'
!          ELSE
!            WRITE(*,*) ' <D> : no external file for SIMUL'
          END IF
        END IF
!
        filename_enkf = filename
        IF (runmode>=2.AND.def_binary=='simul') THEN
          IF (found(79,key_char//' enkf: external file',line,.FALSE.)) THEN
            CALL get_arg('file',line,i,j)
            IF (i<1 .OR. j<i) THEN
              READ(79,*,err=102,end=102) filename_enkf
            ELSE
              READ(line(i:j),*) filename_enkf
            END IF
            WRITE(*,*) ' [R] : external file "'//trim(filename_enkf)//'" for ENKF'
!          ELSE
!            WRITE(*,*) ' <D> : no external file for ENKF'
          END IF
        END IF
!
        filename_inverse = filename
        IF (runmode>=2.AND.def_binary=='inverse') THEN
          IF (found(79,key_char//' inverse: external file',line,.FALSE.)) THEN
            CALL get_arg('file',line,i,j)
            IF (i<1 .OR. j<i) THEN
              READ(79,*,err=103,end=103) filename_inverse
            ELSE
              READ(line(i:j),*) filename_inverse
            END IF
            WRITE(*,*) ' [R] : external file "'//trim(filename_inverse)//'" for INVERSE'
!          ELSE
!            WRITE(*,*) ' <D> : no external file for INVERSE'
          END IF
        END IF
!
        CLOSE(79)
        RETURN
!
!       error handling
100     WRITE(*,'(1A)') 'error: no external file name found in section "data: external file"!'
        STOP
101     WRITE(*,'(1A)') 'error: no external file name found in section "simul: external file"!'
        STOP
102     WRITE(*,'(1A)') 'error: no external file name found in section "enkf: external file"!'
        STOP
103     WRITE(*,'(1A)') 'error: no external file name found in section "inverse: external file"!'
        STOP
      END
