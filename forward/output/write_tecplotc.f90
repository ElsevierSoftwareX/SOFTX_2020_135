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

!>    @brief writess conc-data in tecplot-format
!>    @param[in] ident index number for file name
!>    @param[in] ismpl local sample index
      SUBROUTINE write_tecplotc(ident,ismpl)
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
        integer :: i, j, k, l
        INCLUDE 'OMP_TOOLS.inc'

        DOUBLE PRECISION v_abs
        INTEGER i1, i2, i3, i4, i1s, i2s, ident, lblank, lout, locstr, &
          clast
        EXTERNAL lblank, locstr, clast

        DOUBLE PRECISION vxc, vyc, vzc, kx, ky, kz, lx, ly, lz, por, &
          rhof, visf
        EXTERNAL vxc, vyc, vzc, kx, ky, kz, lx, ly, lz, por, rhof, &
          visf
        DOUBLE PRECISION qxc, qyc, qzc
        EXTERNAL qxc, qyc, qzc

        character (len=80) :: filename
        character (len=20) :: snumber
        character (len=1024) :: sline
        INTRINSIC dsqrt


        IF ( .NOT. tec_out .OR. .NOT. trans_active) RETURN
#ifdef NOPLT
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
          filename = project(i1:i2) // '_' // snumber(i3:i4) // &
            '_chem.plt'
          OPEN(lout,file=filename,status='unknown',blank='null')
          WRITE(lout,'(3A)') 'title = "', project(i1:i2) // '_' // &
            snumber(i3:i4), '_chem"'
        ELSE
          filename = project(i1:i2) // project_sfx(ismpl) (i1s:i2s) // &
            '_' // snumber(i3:i4) // '_chem.plt'
          OPEN(lout,file=filename,status='unknown',blank='null')
          WRITE(lout,'(3A)') 'title = "', project(i1:i2) // &
            project_sfx(ismpl) (i1s:i2s) // '_' // snumber(i3:i4), '"'
        END IF

        IF (linfos(3)>=1) THEN
          WRITE(*,'(3A)') '  [W] : Tecplot to "', &
            filename(1:lblank(filename)), '"'
        END IF


        WRITE(lout,'(7A,9999(1A12,1I4.4,1A1))') &
          'variables = "x",            "y",            "z","uindex",', &
          '  "i",  "j",  "k", ', '           "head",', &
          '           "temp",', '           "pres",', &
          '           "v_abs",', &
          '           "tsal"', (',      "conc',i,'"',i=1,ntrans)

!     use sline for ZONE string
        i1 = locstr(project_sfx(ismpl),'_time_')
        i2 = locstr(project_sfx(ismpl),'_out')
        IF (i1>=1 .AND. i2>=1) THEN
          WRITE(sline,'(1A,1I7,3A,1I7)') 'zone T="', ident, &
            '.time step", SolutionTime=', project_sfx(ismpl) (i1+6:i2- &
            1), ', StrandID=', ident
          WRITE(lout,'(1A,3(A,I5),1A)') sline(1:clast(sline)), &
            ', i=', i0, ', j=', j0, ', k=', k0, ', f=point'
        ELSE
          WRITE(lout,'(3(A,I5),1A)') 'zone i=', i0, ', j=', j0, &
            ', k=', k0, ', f=point'
        END IF

!  output of box-Centred data:
        DO k = 1, k0
          DO j = 1, j0
            DO i = 1, i0
              v_abs = dsqrt(vxc(i,j,k,ismpl)*vxc(i,j,k,ismpl)+ &
                vyc(i,j,k,ismpl)*vyc(i,j,k,ismpl)+ &
                vzc(i,j,k,ismpl)*vzc(i,j,k,ismpl))
              WRITE(lout,'(3(g15.6,1x),1I9,3I6,1X,5(e17.8,1x),9999(e17.8,1A))') & 
                delxa(i), delya(j), delza(k), &
                uindex(i,j,k), i, j, k, &
                head(i,j,k,ismpl), &
                temp(i,j,k,ismpl), pres(i,j,k,ismpl)*pa_conv1, &
                v_abs, tsal(i,j,k,ismpl), &
                (conc(i,j,k,l,ismpl),' ',l=1,ntrans)
            END DO
          END DO
        END DO

        CLOSE(lout)
        CALL omp_del_file_handler(lout)
        CALL compress_file(compress_out,filename)

        RETURN
      END
