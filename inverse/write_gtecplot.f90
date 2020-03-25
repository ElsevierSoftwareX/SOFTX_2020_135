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

!>    @brief writes gradients in tecplot-format
!>    @param[in] ident index number for file name
!>    @param[in] ismpl local sample index
      SUBROUTINE write_gtecplot(ident,ismpl)
        use arrays
#ifdef AD
        use g_arrays
#endif
#ifdef AD_RM
        use arrays_ad
#endif
        use mod_genrl
        use mod_genrlc
        use mod_data
        use mod_time
        use mod_flow
        use mod_temp
        use mod_inverse
        use mod_linfos
        use mod_OMP_TOOLS
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k
        INCLUDE 'OMP_TOOLS.inc'

        DOUBLE PRECISION tu
        INTEGER i1, i2, i3, i4, i1s, i2s, ident, lblank, lout, locstr, &
          clast
        EXTERNAL lblank, locstr, clast

        CHARACTER filename*256, snumber*8, sline*1024


        IF ( .NOT. tec_out) RETURN
#ifdef NOPLT
        RETURN
#endif

!     get his own file discriptor index
        CALL omp_new_file_handler(lout,16)

        tu = tunit

        CALL chln(project,i1,i2)
        CALL chln(project_sfx(ismpl),i1s,i2s)
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
        CALL chln(snumber,i3,i4)

        IF (i1s==0) THEN
          filename = project(i1:i2) // '_G_p' // snumber(i3:i4) // '.plt'
          OPEN(lout,file=filename,status='unknown',blank='null')
          WRITE(lout,'(3A)') 'title = "', project(i1:i2) // '_' // &
            snumber(i3:i4), '"'
        ELSE
          filename = project(i1:i2) // project_sfx(ismpl) (i1s:i2s) // &
            '_G_p' // snumber(i3:i4) // '.plt'
          OPEN(lout,file=filename,status='unknown',blank='null')
          WRITE(lout,'(3A)') 'title = "', project(i1:i2) // &
            project_sfx(ismpl) (i1s:i2s) // '_' // snumber(i3:i4), '"'
        END IF

        IF (linfos(3)>=1) THEN
          WRITE(*,'(3A)') '  [W] : gradients Tecplot to "', &
            filename(1:lblank(filename)), '"'
        END IF

        WRITE(lout,'(13A)') 'variables = "x",            "y",            "z", "uindex",', &
          '  "i",  "j",  "k",', &
#ifndef AD_RM
          '         "g_head",         "g_temp",','         "g_conc"','         "g_pres"'
#else
          '         "head_ad",         "temp_ad",','         "conc_ad"','         "pres_ad"'

#endif
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

!       output of box-centred data:
        DO k = 1, k0
          DO j = 1, j0
            DO i = 1, i0
              WRITE(lout, '(3(g15.6,1x),i10,3i6,4(e17.8,1x))') &
                delxa(i), delya(j), delza(k), uindex(i,j,k), i, j, k, &
!#ifndef AD_RM
                g_head(i,j,k,ismpl), g_temp(i,j,k,ismpl), g_conc(i,j,k,1,ismpl), g_pres(i,j,k,ismpl)
!#else
!                head_ad(i,j,k,ismpl), temp_ad(i,j,k,ismpl), conc_ad(i,j,k,1,ismpl), pres_ad(i,j,k,ismpl)
!#endif

            END DO
          END DO
        END DO

        CLOSE(lout)
        CALL omp_del_file_handler(lout)
        CALL compress_file(compress_out,filename)

        RETURN
      END
