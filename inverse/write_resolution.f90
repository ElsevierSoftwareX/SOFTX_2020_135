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

!>    @brief writes resolution matrix
!>    @param[in] ident inverse iteration counter
!>    @param[in] ismpl local sample index
      SUBROUTINE write_resolution(ident,ismpl)
        use arrays
        use mod_genrl
        use mod_genrlc
        use mod_inverse
        use mod_data
        use mod_linfos
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j

        INTEGER ident, lout
        INTEGER i1_, i2_, i3_, i4_

        CHARACTER filename*80, snumber*8
        INTEGER lblank
        EXTERNAL lblank


        IF ((resmat<=0) .OR. write_disable) RETURN
#ifdef NOOUT
        RETURN
#endif

        IF (mpara<1) THEN
          WRITE(*,*) 'error : in "write_resolution" need "mpara">0 !'
          STOP
        END IF

        CALL omp_new_file_handler(lout,16)

        CALL chln(project,i1_,i2_)
        IF (ident>=0) THEN
          WRITE(snumber,'(1I7)') ident
        ELSE IF (ident==-1) THEN
          WRITE(snumber,'(A8)') 'final'
        ELSE IF (ident==-2) THEN
          WRITE(snumber,'(A8)') 'debug'
        ELSE IF (ident==-3) THEN
          WRITE(snumber,'(A8)') 'ens_mean'
        ELSE IF (ident==-4) THEN
          WRITE(snumber,'(A8)') 'mean'
        ELSE IF (ident==-5) THEN
          WRITE(snumber,'(A8)') 'ens_mean'
        END IF
        CALL chln(snumber,i3_,i4_)
        filename = project(i1_:i2_) // '_res-p_' // &
          snumber(i3_:i4_) // '.dat'

        OPEN(lout,file=filename,status='unknown',blank='null')

        IF (linfos(2)>=0) THEN
          WRITE(*,'(3A)') '  [W] : Resolution-Parameter to "', &
            filename(1:lblank(filename)), '"'
        END IF

        IF (mpara>200) WRITE(*,*) &
          'warning: line to short in "write_resolution"'
        DO i = 1, mpara
          WRITE(lout,'(200(e15.7,1x))') (resmat_p(i,j),j=1,mpara)
        END DO

        CLOSE(lout)
        CALL omp_del_file_handler(lout)

        RETURN
      END
