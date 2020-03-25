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

!>    @brief writes data in vtk-format
!>    @param[in] ident output file index number
!>    @param[in] ismpl local sample index
!>    @details
!>writes data in vtk-format  \n
!>see:         http://mayavi.sourceforge.net\n
!>             http://www.vtk.org/pdf/file-formats.pdf\n
      SUBROUTINE write_vtk(ident,ismpl)
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

        DOUBLE PRECISION, ALLOCATABLE :: val(:)
        ! DOUBLE PRECISION dx, dy, dz
        INTEGER i1s, i2s, i1, i2, i3, i4, species, ident, lblank, lout
        EXTERNAL lblank

        DOUBLE PRECISION vxc, vyc, vzc, kx, ky, kz, lx, ly, lz, por, & 
          qt, rhof,visf, visn
        EXTERNAL vxc, vyc, vzc, kx, ky, kz, lx, ly, lz, por, qt, rhof, & 
          visf, visn

        DOUBLE PRECISION qxc, qyc, qzc
        EXTERNAL qxc, qyc, qzc

        character (len=20) :: snumber
        character (len=256) :: filename, prname, strng
        

        IF ( .NOT. vtk_out) RETURN
#ifdef NOVTK
        RETURN
#endif

!     get his own file discriptor index
        CALL omp_new_file_handler(lout,16)

        ALLOCATE(val(max(i0,j0,k0)))

        CALL chln(project,i1,i2)
        CALL chln(project_sfx(ismpl),i1s,i2s)
        IF (ident>=0) THEN
          WRITE(snumber,'(I7)') ident
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
          prname = project(i1:i2) // '_' // snumber(i3:i4)
          filename = project(i1:i2) // '_' // snumber(i3:i4) // '.vtk'
        ELSE
          prname = project(i1:i2) // project_sfx(ismpl) (i1s:i2s) // & 
            '_' // snumber(i3:i4)
          filename = project(i1:i2) // project_sfx(ismpl) (i1s:i2s) // & 
            '_' // snumber(i3:i4) // '.vtk'
        END IF

        OPEN(lout,file=filename,status='unknown',blank='null')

        IF (linfos(3)>=1) THEN
          WRITE(*,'(3A)') '  [W] : VTK to "', & 
            filename(1:lblank(filename)), '"'
        END IF

! HEADER
        WRITE(lout,'(a/a/a/a/a,3I7)') '# vtk DataFile Version 2.0', & 
          prname, 'ASCII', 'DATASET RECTILINEAR_GRID', 'DIMENSIONS  ', & 
          i0, j0, k0
          !i0+1, j0+1, k0


! X_COORDINATES
        WRITE(lout,'(a,I7,a)') 'X_COORDINATES  ', i0, ' float'
        val(1) = 0.5D0*delx(1)
        DO i = 2, i0
          val(i) = val(i-1) + 0.5D0*(delx(i-1)+delx(i))
        END DO
        WRITE(lout,'(100e16.6)') (val(i),i=1,i0)

        !WRITE(lout,'(a,I7,a)') 'X_COORDINATES  ', i0+1, ' float'
        !WRITE(lout,'(100e16.6)') (delxa(i)-0.5d0*delx(i),i=1,i0)
        !write(unit = lout, fmt = '(e16.6)') delxa(i0) + 0.5d0*delx(i0) 

! Y_COORDINATES
        WRITE(lout,'(a,I7,a)') 'Y_COORDINATES  ', j0, ' float'
        val(1) = 0.5D0*dely(1)
        DO i = 2, j0
          val(i) = val(i-1) + 0.5D0*(dely(i-1)+dely(i))
        END DO
        WRITE(lout,'(100e16.6)') (val(i),i=1,j0)

        !WRITE(lout,'(a,I7,a)') 'Y_COORDINATES  ', j0+1, ' float'
        !WRITE(lout,'(100e16.6)') (delya(i)-0.5d0*dely(i),i=1,j0)
        !write(unit = lout, fmt = '(e16.6)') delya(j0) + 0.5d0*dely(j0) 


! Z_COORDINATES
        WRITE(lout,'(a,I7,a)') 'Z_COORDINATES  ', k0, ' float'
        val(1) = 0.5D0*delz(1)
        DO i = 2, k0
          val(i) = val(i-1) + 0.5D0*(delz(i-1)+delz(i))
        END DO
        WRITE(lout,'(100e16.6)') (val(i),i=1,k0)

        !WRITE(lout,'(a,I7,a)') 'Z_COORDINATES  ', k0+1, ' float'
        !WRITE(lout,'(100e16.6)') (delza(i)-0.5d0*delz(i),i=1,k0)
        !write(unit = lout, fmt = '(e16.6)') delza(k0) + 0.5d0*delz(k0)

        !WRITE(lout,'(a,I7,a)') 'Z_COORDINATES  ', k0, ' float'
        !write(unit = lout, fmt = '(100e16.6)') (delza(i),i=1,k0)


        WRITE(lout,'(/a,i8)') 'POINT_DATA', i0*j0*k0
        !WRITE(lout,'(/a,i8)') 'CELL_DATA', i0*j0*k0


        if(out_ijk(cout_uindex)) then
           WRITE(lout,'(a)') 'SCALARS uindex int 1'
           WRITE(lout,'(a)') 'LOOKUP_TABLE default'
           WRITE(lout,'(25I10)') (((uindex(i,j,k),i=1,i0),j=1,j0),k=1,k0)
        end if

        if(out_prop(idx_por)) then
           WRITE(lout,'(a)') 'SCALARS por float 1'
           WRITE(lout,'(a)') 'LOOKUP_TABLE default'
           WRITE(lout,'(100e16.6)') (((por(i,j,k,ismpl),i=1,i0),j=1,j0),k=1,k0)
        end if

        IF (temp_active) THEN
           if(out_pv(pv_temp)) then
              WRITE(lout,'(a)') 'SCALARS temp float 1'
              WRITE(lout,'(a)') 'LOOKUP_TABLE default'
              WRITE(lout,'(100e16.6)') (((temp(i,j,k,ismpl),i=1,i0),j=1,j0),k=1,k0)
           end if
           if(out_prop(idx_an_lx)) then
              WRITE(lout,'(a)') 'SCALARS lx float 1'
              WRITE(lout,'(a)') 'LOOKUP_TABLE default'
              WRITE(lout,'(100e16.6)') (((lx(i,j,k,ismpl),i=1,i0),j=1,j0),k=1,k0)
           end if
           if(out_prop(idx_an_ly)) then
              WRITE(lout,'(a)') 'SCALARS ly float 1'
              WRITE(lout,'(a)') 'LOOKUP_TABLE default'
              WRITE(lout,'(100e16.6)') (((ly(i,j,k,ismpl),i=1,i0),j=1,j0),k=1,k0)
           end if
           if(out_prop(idx_lz)) then
              WRITE(lout,'(a)') 'SCALARS lz float 1'
              WRITE(lout,'(a)') 'LOOKUP_TABLE default'
              WRITE(lout,'(100e16.6)') (((lz(i,j,k,ismpl),i=1,i0),j=1,j0),k=1,k0)
           end if
           if(out_prop(idx_q)) then
              WRITE(lout,'(a)') 'VECTORS q float'
              WRITE(lout,'(3e16.6)') (((qxc(i,j,k,ismpl),qyc(i,j,k,ismpl),qzc(i,j,k,ismpl),i=1,i0),j=1,j0),k=1,k0)
              WRITE(lout,'(a)') 'SCALARS h float 1'
              WRITE(lout,'(a)') 'LOOKUP_TABLE default'
              WRITE(lout,'(100e16.6)') (((qt(i,j,k,ismpl),i=1,i0),j=1,j0),k=1,k0)
           end if
        END IF

        IF (head_active) THEN
           if(out_pv(pv_head)) then
              WRITE(lout,'(a)') 'SCALARS head float 1'
              WRITE(lout,'(a)') 'LOOKUP_TABLE default'
              WRITE(lout,'(100e16.6)') (((head(i,j,k,ismpl),i=1,i0),j=1,j0),k=1,k0)
           end if
        END IF
        IF (pres_active) THEN
           if(out_pv(pv_pres)) then
              WRITE(lout,'(a)') 'SCALARS pres float 1'
              WRITE(lout,'(a)') 'LOOKUP_TABLE default'
              WRITE(lout,'(100e16.6)') (((pres(i,j,k,ismpl)*pa_conv1,i=1,i0),j=1,j0),k=1,k0)
#ifdef pres_base
              CALL pres2head(1,ismpl)
              WRITE(lout,'(a)') 'SCALARS head float 1'
              WRITE(lout,'(a)') 'LOOKUP_TABLE default'
              WRITE(lout,'(100e16.6)') (((head(i,j,k,ismpl),i=1,i0),j=1,j0),k=1,k0)
#endif
           end if
        END IF

        IF (head_active .OR. pres_active) THEN
          IF (klogflag) THEN
             if(out_prop(idx_an_kx)) then
                WRITE(lout,'(a)') 'SCALARS kx float 1'
                WRITE(lout,'(a)') 'LOOKUP_TABLE default'
                WRITE(lout,'(100e16.6)') (((log10(kx(i,j,k,ismpl)),i=1, &
                     i0),j=1,j0),k=1,k0)
             end if
             if(out_prop(idx_an_ky)) then
                WRITE(lout,'(a)') 'SCALARS ky float 1'
                WRITE(lout,'(a)') 'LOOKUP_TABLE default'
                WRITE(lout,'(100e16.6)') (((log10(ky(i,j,k,ismpl)),i=1, &
                     i0),j=1,j0),k=1,k0)
             end if
             if(out_prop(idx_kz)) then
                WRITE(lout,'(a)') 'SCALARS kz float 1'
                WRITE(lout,'(a)') 'LOOKUP_TABLE default'
                WRITE(lout,'(100e16.6)') (((log10(kz(i,j,k,ismpl)),i=1, &
                     i0),j=1,j0),k=1,k0)
             end if
             if(out_ijk(cout_rhof)) then
                WRITE(lout,'(a)') 'SCALARS rhof float 1'
                WRITE(lout,'(a)') 'LOOKUP_TABLE default'
                WRITE(lout,'(100e16.6)') (((rhof(i,j,k,ismpl),i=1, &
                     i0),j=1,j0),k=1,k0)
             end if
          ELSE
             if(out_prop(idx_an_kx)) then
                WRITE(lout,'(a)') 'SCALARS kx float 1'
                WRITE(lout,'(a)') 'LOOKUP_TABLE default'
                WRITE(lout,'(100e16.6)') (((kx(i,j,k,ismpl),i=1, &
                     i0),j=1,j0),k=1,k0)
             end if
             if(out_prop(idx_an_ky)) then
                WRITE(lout,'(a)') 'SCALARS ky float 1'
                WRITE(lout,'(a)') 'LOOKUP_TABLE default'
                WRITE(lout,'(100e16.6)') (((ky(i,j,k,ismpl),i=1, &
                     i0),j=1,j0),k=1,k0)
             end if
             if(out_prop(idx_kz)) then
                WRITE(lout,'(a)') 'SCALARS kz float 1'
                WRITE(lout,'(a)') 'LOOKUP_TABLE default'
                WRITE(lout,'(100e16.6)') (((kz(i,j,k,ismpl),i=1, &
                     i0),j=1,j0),k=1,k0)
             end if
          END IF
          ! WRITE(lout,'(a)') 'VECTORS v float'
          ! WRITE(lout,'(3e16.6)') (((vxc(i,j,k,ismpl), & 
          !   vyc(i,j,k,ismpl),vzc(i,j,k,ismpl),i=1, & 
          !   i0),j=1,j0),k=1,k0)
        END IF

        IF (trac_active) THEN
           if(out_pv(pv_conc)) then
              DO species = 1, ntrans
                 WRITE(strng,'(I7)') species
                 CALL chln(strng,i1,i2)
                 WRITE(lout,'(a)') 'SCALARS tracer' // strng(i1:i2) // ' float 1'
                 WRITE(lout,'(a)') 'LOOKUP_TABLE default'
                 WRITE(lout,'(100e16.6)') (((conc(i,j,k,species,ismpl),i=1, &
                      i0),j=1,j0),k=1,k0)
              END DO
           end if
        END IF

        !Supposed to be written behind tracer output
        IF (head_active .OR. pres_active) THEN
          WRITE(lout,'(a)') 'VECTORS v float'
          WRITE(lout,'(3e16.6)') (((vxc(i,j,k,ismpl), &
            vyc(i,j,k,ismpl),vzc(i,j,k,ismpl),i=1, &
            i0),j=1,j0),k=1,k0)
        END IF

        DEALLOCATE(val)

        CLOSE(lout)
        CALL omp_del_file_handler(lout)
        CALL compress_file(compress_out,filename)

        RETURN
      END
