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

!>    @brief writess data in text-format
!>    @param[in] ident index/iteration number
!>    @param[in] ismpl local sample index
      SUBROUTINE write_text(ident,ismpl)
        use arrays
        use mod_genrl
        use mod_genrlc
        use mod_data
        use mod_time
        use mod_flow
        use mod_temp
        use mod_conc
        use mod_linfos
        use mod_OMP_TOOLS
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k





        INCLUDE 'OMP_TOOLS.inc'

        INTEGER i1, i2, i3, i4, i1s, i2s, ident, lblank, lout, tracer
        EXTERNAL lblank
        character (len=80) :: filename
        character (len=20) :: snumber
        character (len=32) :: strng

        DOUBLE PRECISION vxc, vyc, vzc
        EXTERNAL vxc, vyc, vzc


        IF ( .NOT. txt_out) RETURN
#ifdef NOTXT
        RETURN
#endif

!     get his own file discriptor index
        CALL omp_new_file_handler(lout,16)

        CALL chln(project,i1,i2)
        CALL chln(project_sfx(ismpl),i1s,i2s)
        IF (ident>=0) THEN
          WRITE(snumber,'(1I7)') ident
        ELSE IF (ident==-1) THEN
          WRITE(snumber,'(A20)') 'final'
        ELSE IF (ident==-2) THEN
          WRITE(snumber,'(A20)') 'debug'
        ELSE IF (ident==-3) THEN
          WRITE(snumber,'(A20)') 'ens_mean'
        ELSE IF (ident==-4) THEN
          WRITE(snumber,'(A20)') 'mean'
        ELSE IF (ident==-5) THEN
          WRITE(snumber,'(A20)') 'ens_mean'
        END IF
        CALL chln(snumber,i3,i4)

        IF (i1s==0) THEN
          filename = project(i1:i2) // '_' // snumber(i3:i4) // '.txt'
          OPEN(lout,file=filename,status='unknown',blank='null')
          WRITE(lout,'(1A)') key_char//' title'
          WRITE(lout,'(1A)') project(i1:i2) // '_' // snumber(i3:i4)
        ELSE
          filename = project(i1:i2) // project_sfx(ismpl) (i1s:i2s) // &
            '_' // snumber(i3:i4) // '.txt'
          OPEN(lout,file=filename,status='unknown',blank='null')
          WRITE(lout,'(1A)') key_char//' title'
          WRITE(lout,'(1A)') project(i1:i2) // &
            project_sfx(ismpl) (i1s:i2s) // '_' // snumber(i3:i4)
        END IF

        IF (linfos(3)>=1) THEN
          WRITE(*,'(3A)') '  [W] : Text to "', &
            filename(1:lblank(filename)), '"'
        END IF

! --------
        WRITE(lout,'(1A)') key_char//' delx'
        CALL write_dense3d(delx,i0,1,1,lout,ismpl)

        WRITE(lout,'(1A)') key_char//' dely'
        CALL write_dense3d(dely,1,j0,1,lout,ismpl)

        WRITE(lout,'(1A)') key_char//' delz'
        CALL write_dense3d(delz,1,1,k0,lout,ismpl)

! --------
        WRITE(lout,'(1A)') key_char//' head init'
        CALL write_dense3d(head(1,1,1,ismpl),i0,j0,k0,lout,ismpl)

        WRITE(lout,'(1A)') key_char//' temp init'
        CALL write_dense3d(temp(1,1,1,ismpl),i0,j0,k0,lout,ismpl)

        DO i = 1, ntrans
          WRITE(lout,'(1A,1I4.4,1A)') key_char//' tracer', i, ' init'
          CALL write_dense3d(conc(1,1,1,i,ismpl),i0,j0,k0,lout,ismpl)
        END DO

        WRITE(lout,'(1A)') key_char//' pres init'
        DO k = 1, k0
          DO j = 1, j0
            DO i = 1, i0
              x(i,j,k,ismpl) = pres(i,j,k,ismpl)*pa_conv1
            END DO
          END DO
        END DO
        CALL write_dense3d(x(1,1,1,ismpl),i0,j0,k0,lout,ismpl)

        IF (trac_active) THEN
          DO tracer = 1, ntrac
            WRITE(strng,'(1A,1I4.4)') 'tracer', tracer
            CALL chln(strng,i1,i2)
            WRITE(lout,'(3A)') key_char//' ', strng(i1:i2), ' init'
            CALL write_dense3d(conc(1,1,1,tracer,ismpl),i0,j0,k0,lout, &
              ismpl)
          END DO
        END IF

! --------
        WRITE(lout,'(1A)') key_char//' vx init'
        DO k = 1, k0
          DO j = 1, j0
            DO i = 1, i0
              x(i,j,k,ismpl) = vxc(i,j,k,ismpl)
            END DO
          END DO
        END DO
        CALL write_dense3d(x(1,1,1,ismpl),i0,j0,k0,lout,ismpl)

        WRITE(lout,'(1A)') key_char//' vy init'
        DO k = 1, k0
          DO j = 1, j0
            DO i = 1, i0
              x(i,j,k,ismpl) = vyc(i,j,k,ismpl)
            END DO
          END DO
        END DO
        CALL write_dense3d(x(1,1,1,ismpl),i0,j0,k0,lout,ismpl)

        WRITE(lout,'(1A)') key_char//' vz init'
        DO k = 1, k0
          DO j = 1, j0
            DO i = 1, i0
              x(i,j,k,ismpl) = vzc(i,j,k,ismpl)
            END DO
          END DO
        END DO
        CALL write_dense3d(x(1,1,1,ismpl),i0,j0,k0,lout,ismpl)

! --------
        CLOSE(lout)
        CALL omp_del_file_handler(lout)
        CALL compress_file(compress_out,filename)

        RETURN
      END
