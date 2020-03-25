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

!> @brief read additional ENKF variables and create subdirectories
!> @param[in] filename ENKF parameter file name
!> @param[in] ismpl local sample index
      SUBROUTINE read_enkf(filename,ismpl)

        use mod_genrl, only: &
             key_char

        use mod_simul, only: &
             enkf_post,&
             maxiter_enkf

        IMPLICIT NONE
        integer :: ismpl
        integer :: i
        character (len=*) :: filename
        character (len=5000) :: line
        integer :: locstr
        LOGICAL found
        EXTERNAL locstr, found

!     open file
        OPEN(79,file=filename,status='old')
        WRITE(*,*)
        WRITE(*,*) '  reading ENKF parameter'
        WRITE(*,*)
!
!     init HDF5 support, when available
        CALL open_hdf5(' ')

        ! # enkf postcompute
        enkf_post = 0
        IF (found(79,key_char//' enkf postcompute',line,.FALSE.)) THEN
          enkf_post = 0
          i = locstr(line,'samples')
          IF (i>=1) THEN
            enkf_post = 1
            WRITE(*,*) ' [R] : postcomputing ENKF ensembles'
          END IF
          i = locstr(line,'mean')
          IF (i>=1) THEN
            enkf_post = enkf_post +2
            WRITE(*,*) ' [R] : postcomputing ENKF mean (currently, recomputed from all ensembles)'
          END IF
          i = locstr(line,'none')
          IF (i>=1) enkf_post = 0
          IF (enkf_post == 0) WRITE(*,*) ' [R] : no ENKF postcomputing'
        ELSE
          WRITE(*,*) ' <D> : no ENKF postcomputing'
        END IF

        ! # enkf iter
        maxiter_enkf = 1
        IF (found(79,key_char//' enkf iter',line,.FALSE.)) THEN
          READ(79,*) maxiter_enkf
          WRITE(*,'(1A,1I3)') '  [R] : number of ENKF global assimilation iterations = ',maxiter_enkf
        ELSE
          WRITE(*,'(1A)') '  <D> : number of ENKF global assimilation iterations = 1'
        END IF
!
!     finish HDF5 support, when available
        CALL close_hdf5()
        CLOSE(79)
        !
        call enkf_make_output_dirs(.false.)
        
        RETURN
      END
