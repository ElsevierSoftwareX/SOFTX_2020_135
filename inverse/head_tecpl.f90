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

!>    @brief write a Tecplot header to file "fileid"
!>    @param[in] fileid file handler
!>    @param[in] ni grid size in I0 direction
!>    @param[in] nj grid size in J0 direction
!>    @param[in] nk grid size in K0 direction
!>    @param[in] feld_name array/component name
!>    @param[in] ismpl local sample index
      SUBROUTINE head_tecpl(fileid,ni,nj,nk,feld_name,ismpl)
        use arrays
        use mod_inverse
        use mod_genrl
        use mod_genrlc
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j
        CHARACTER*25 zeilen_form, zz_form
        INTEGER, PARAMETER :: zbuff_size=16768
        CHARACTER*(zbuff_size) zeilen_buff
        CHARACTER*4 feld_name
        CHARACTER*16 pname
        INTEGER l10, maxp, sollp, fileid, ni, nj, nk


        CALL chln(project,i,j)
        WRITE(fileid,'(3A)') 'title = "', project(i:j), '"'

        maxp = 1
        sollp = 4*15 + 9 + 19*mpara
        IF (zbuff_size<sollp) THEN
          WRITE(*,*) &
            'error, "zeilen_buff" to small in "head_tecpl" !!!'
          STOP
        END IF

        l10 = int(log(dble(sollp))/log(10.0D0)) + 1
        zz_form = '(A13,I' // achar(l10+48) // ',A4)'
        WRITE(zeilen_form,zz_form) '(A', sollp, ')'

        zeilen_buff(maxp:maxp+15) = 'variables = "x"'
        maxp = maxp + 15
        zeilen_buff(maxp:maxp+15) = ',           "y"'
        maxp = maxp + 15
        zeilen_buff(maxp:maxp+15) = ',           "z"'
        maxp = maxp + 15
        zeilen_buff(maxp:maxp+9) = ',"uindex"'
        maxp = maxp + 9
        zeilen_buff(maxp:maxp+15) = ',        "' // feld_name // '"'
        maxp = maxp + 15

100     FORMAT (1A2,1A16,1A1)
        DO i = 1, mpara
          CALL param_name(i,pname,ismpl)
          WRITE(zeilen_buff(maxp:maxp+19),100) ',"', pname, '"'
          maxp = maxp + 19
        END DO

        IF ((maxp-1)/=sollp) THEN
          WRITE(*,*) 'error, "maxP"<>"sollP" in "head_tecpl" !!!'
          STOP
        END IF

        WRITE(fileid,'(A)') zeilen_buff(1:maxp-1)

        WRITE(fileid,'(3(a,i5),a)') 'zone i=', ni, ', j=', nj, &
          ', k=', nk, ', f=point'

        RETURN
      END
