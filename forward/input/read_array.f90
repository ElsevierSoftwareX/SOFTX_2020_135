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

!>    @brief read MM*NN double precision values from file
!>    @param[in] fid file id
!>    @param[out] AA array to fill
!>    @param[in] MM number of values expected in each line
!>    @param[in] NN number of lines
!>    @param[in] ustr section name
!>    @param[in] ismpl local sample index (ignored here)
      SUBROUTINE read_array(fid,mm,nn,aa,ustr,ismpl)
        use arrays
        use mod_genrl
        IMPLICIT NONE
!
        INTEGER mm, nn, i, j, ismpl
!       warning flag
        LOGICAL w_flag
!       file id
        INTEGER fid
!       array to fill
        DOUBLE PRECISION aa(mm,nn)
!       unit string, for error message
        character (len=*) :: ustr
!       for reading each line separated
        character (len=1024) :: line

        w_flag = .false.
        IF (ustr==key_char//' units' .OR. &
            ustr==key_char//' errors' .OR. &
            ustr==key_char//' apriori') THEN
          DO j = 1, nn
            READ(fid,'(A)',err=200,end=200) line
            READ(line,*,err=100,end=100) (aa(i,j),i=1,mm)
            GOTO 10
100         IF (i-1<mm) w_flag = .true.
10          CONTINUE
          END DO
          IF (w_flag) GOTO 300
        ELSE
          DO j = 1, nn
            READ(fid,*,err=200,end=200) (aa(i,j),i=1,mm)
          END DO
        ENDIF
        RETURN
!
!       error handler
200     WRITE(*,'(3A,1I5,1A,1I5,1A)') 'error: while reading section "', ustr, &
          '", to few values (need ', mm, ' in each of the ', nn, ' lines) !'
        STOP
300     WRITE(*,'(3A,1I5,1A,1I5,1A)') '  <D> : WARNING while reading section "', ustr, &
          '", to few values (need ', mm, ' in each of the ', nn, ' lines) !'
        IF (mm==nprop_load) THEN
          IF (ustr==key_char//' units') WRITE(*,'(1I3,2A,'//c_npropunit//'(",",1A),1A)') mm, &
            ' properties each line : [', (properties(i),i=1,nprop_load), ']'
          IF (ustr==key_char//' errors') WRITE(*,'(1I3,2A,'//c_npropunit//'(",",1A),1A)') mm, &
            ' property errors each line : [', (properties(i),i=1,nprop_load), ']'
          IF (ustr==key_char//' apriori') WRITE(*,'(1I3,2A,'//c_npropunit//'(",",1A),1A)') &
            mm, ' apriori values each line : [', (properties(i),i=1,nprop_load), ']'
          WRITE(*,*) ' '
          WRITE(*,'(1A)') ' properties:'
          DO i = 1, nprop_load
            WRITE(*,'(4A)') '  ', properties(i), ': ', doc_properties(i)
          END DO
          WRITE(*,*) ' '
        END IF
        RETURN
      END
