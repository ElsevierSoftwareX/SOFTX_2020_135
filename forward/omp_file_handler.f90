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

!>    @brief create a new file handler (number)
!>    @param[out] fh file handler (number)
!>    @param[in] i 0: reset file handler table, >0 : handler index, specific for each thread
      SUBROUTINE omp_new_file_handler(fh,i)
        use arrays
        use mod_OMP_TOOLS
        IMPLICIT NONE
        INCLUDE 'OMP_TOOLS.inc'
        INTEGER fh, i, o, l2

!     reset the file-handler table
        IF (i==0) THEN
          l2 = tlevel_0
          DO fh = 1, c_fhandler
            DO o = 1, l2
              fh_table(fh,o) = 0
            END DO
          END DO
          fh = -1
          RETURN
        END IF
!
        o = omp_get_his_thread_num()
!     looking for an already used file-handler (for this thread)
        DO l2 = 1, c_fhandler
          IF (fh_table(l2,o+1)==i) THEN
            fh = l2 + o*c_fhandler + c_foffset
            RETURN
          END IF
        END DO
!     looking for a free file-handler (not used from any thread)
        DO l2 = 1, c_fhandler
          IF (fh_table(l2,o+1)==0) THEN
            fh_table(l2,o+1) = i
            fh = l2 + o*c_fhandler + c_foffset
            RETURN
          END IF
        END DO
!     error handler
        WRITE(*,'(1A,1I4,1A)') 'error: max number (', c_fhandler, &
          ') of open files per thread, in "omp_file_handler.f" !'
        STOP
      END

!>    @brief remove used file handler (number)
!>    @param[in,out] fh file handler (number), gets 0 when successful
      SUBROUTINE omp_del_file_handler(fh)
        use arrays
        use mod_OMP_TOOLS
        IMPLICIT NONE
        INCLUDE 'OMP_TOOLS.inc'
        INTEGER fh, o, l2

        o = omp_get_his_thread_num()
        l2 = fh - o*c_fhandler - c_foffset
!     sanity check
        IF (l2<=0 .OR. l2>c_fhandler) THEN
          WRITE(*,'(1A)') &
            'error: wrong file-handler, in "omp_file_handler.f" !'
          STOP
        END IF
!     free the already used file-handler
        IF (fh_table(l2,o+1)/=0) THEN
          fh_table(l2,o+1) = 0
          fh = 0
          RETURN
        END IF
!     error handler
        WRITE(*,'(1A,1I4,1A)') 'error: file-handler ', fh, &
          ' already closed, in "omp_file_handler.f" !'
        STOP
      END
