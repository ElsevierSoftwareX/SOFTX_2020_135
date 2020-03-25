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

!>    @brief writes data in tecplot-format
!>    @param[in] ident index number for file name
!>    @param[in] ismpl local sample index
      SUBROUTINE write_tecplot(ident,ismpl)
        use arrays
        use mod_genrl
        use mod_genrlc
        use mod_data
        use mod_time
        use mod_flow
        use mod_temp
        use mod_linfos
        use mod_OMP_TOOLS
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k
        INCLUDE 'OMP_TOOLS.inc'

        INTEGER i1, i2, i3, i4, i1s, i2s, ident, lblank, lout, locstr, &
          clast
        EXTERNAL lblank, locstr, clast

        DOUBLE PRECISION vxc, vyc, vzc, kx, ky, kz, lx, ly, lz, por, &
          rhof, visf
        EXTERNAL vxc, vyc, vzc, kx, ky, kz, lx, ly, lz, por, rhof, &
          visf
        DOUBLE PRECISION qxc, qyc, qzc
        EXTERNAL qxc, qyc, qzc

        DOUBLE PRECISION px, py, pz
        character (len=256) :: filename
        character (len=20) :: snumber
        character (len=1024) :: sline

        double precision, parameter :: zero = 0.0d0

        logical :: cell_centered
        integer :: line_break_cc

        IF ( .NOT. tec_out) RETURN
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
          filename = project(i1:i2) // '_' // snumber(i3:i4) // '.plt'
          OPEN(lout,file=filename,status='unknown',blank='null')
          WRITE(lout,'(3A)') 'title = "', project(i1:i2) // '_' // &
            snumber(i3:i4), '"'
        ELSE
          filename = project(i1:i2) // project_sfx(ismpl) (i1s:i2s) // &
            '_' // snumber(i3:i4) // '.plt'
          OPEN(lout,file=filename,status='unknown',blank='null')
          WRITE(lout,'(3A)') 'title = "', project(i1:i2) // &
            project_sfx(ismpl) (i1s:i2s) // '_' // snumber(i3:i4), '"'
        END IF

        IF (linfos(3)>=1) THEN
          WRITE(*,'(3A)') '  [W] : Tecplot to "', &
            filename(1:lblank(filename)), '"'
        END IF


        cell_centered = .false.
           
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if(cell_centered) then
           
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           WRITE(lout,'(10A)') 'variables = "xnode","ynode","znode","i",  "j",  "k", "x",           &
                & "y",            "z", "uindex",', &
                '   "obsz","obst","cindex",', &
                '           "head",           "temp",', &
                '           "pres",', &
                '           "vx",           "vy",           "vz",', &
                '          "por",           "kx",           "ky",', &
                '           "kz",           "lx",           "ly",', &
                '           "lz",         "rhof",         "visf",', &
                '          "qxc",          "qyc",          "qzc",', &
                '        "sigma",           "lc",          "src",'


           ! DIFFERENCE FOR CELLCENTERED
           i1 = locstr(project_sfx(ismpl),'_time_')
           i2 = locstr(project_sfx(ismpl),'_out')
           IF (i1>=1 .AND. i2>=1) THEN
              WRITE(sline,'(1A,1I7,3A,1I7)') 'zone T="', ident, &
                   '.time step", SolutionTime=', project_sfx(ismpl) (i1+6:i2- &
                   1), ', StrandID=', ident
              WRITE(lout,'(1A,3(A,I5),1A)') sline(1:clast(sline)), &
                   ', i=', i0+1, ', j=', j0+1, ', k=', k0+1, ', f=block'
           ELSE
              WRITE(lout,'(3(A,I5),1A)') 'zone i=', i0+1, ', j=', j0+1, &
                   ', k=', k0+1, ', f=block'  
           END IF

           ! DIFFERENCE FOR CELLCENTERED
           write(lout,*) 'VARLOCATION=([1,2,3]=NODAL,&
                & [4,5,6,7,8,9,10,11,12,13,14,15,16,17&
                &,18,19,20,21,22,23,24,25,26,27,28,29,&
                &30,31,32,33,34,35,36,37,38,39]=CELLCENTERED)'
           line_break_cc = 0
           do k = 1, k0+1
              do j = 1, j0+1
                 do i = 1, i0
                    write(unit = lout, fmt = '(g15.6,1x)', advance='no') delxa(i)-0.5d0*delx(i) 
                    if(i==i0) write(unit = lout, fmt = '(g15.6,1x)', advance='no') delxa(i)+0.5d0*delx(i) 
                    ! Line break after 100 entries
                    line_break_cc = line_break_cc + 1
                    if(line_break_cc >= 100) then
                       write(unit = lout, fmt = *) 
                       line_break_cc = 0
                    end if
                 end do
              end do
           end do
           line_break_cc = 0
           write(unit = lout, fmt = *) 
           do k = 1, k0+1
              do j = 1, j0
                 do i = 1, i0+1
                    write(unit = lout, fmt = '(g15.6,1x)', advance='no') delya(j)- 0.5d0*dely(j) 
                    ! Line break after 100 entries
                    line_break_cc = line_break_cc + 1
                    if(line_break_cc >= 100) then
                       write(unit = lout, fmt = *) 
                       line_break_cc = 0
                    end if
                 end do
                 if(j==j0) then
                    do i = 1, i0+1
                       write(unit = lout, fmt = '(g15.6,1x)', advance='no') delya(j)+0.5d0*dely(j) 
                       ! Line break after 100 entries
                       line_break_cc = line_break_cc + 1
                       if(line_break_cc >= 100) then
                          write(unit = lout, fmt = *) 
                          line_break_cc = 0
                       end if
                    end do
                 end if
              end do
           end do
           line_break_cc = 0
           write(unit = lout, fmt = *) 
           do k = 1, k0
              do j = 1, j0+1
                 do i = 1, i0+1
                    write(unit = lout, fmt = '(g15.6,1x)', advance='no') delza(k)-0.5d0*delz(k) 
                    ! Line break after 100 entries
                    line_break_cc = line_break_cc + 1
                    if(line_break_cc >= 100) then
                       write(unit = lout, fmt = *) 
                       line_break_cc = 0
                    end if
                 end do
              end do
              if(k==k0) then
                 do j = 1, j0+1
                    do i = 1, i0+1
                       write(unit = lout, fmt = '(g15.6,1x)', advance='no') delza(k)+0.5d0*delz(k) 
                       ! Line break after 100 entries
                       line_break_cc = line_break_cc + 1
                       if(line_break_cc >= 100) then
                          write(unit = lout, fmt = *) 
                          line_break_cc = 0
                       end if
                    end do
                 end do
              end if
           end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           write(unit = lout, fmt = *) 
           write(unit = lout, fmt = *) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           line_break_cc = 0
           do k = 1, k0
              do j = 1, j0
                 do i = 1, i0
                    write(unit = lout, fmt = '(i6)', advance='no') i
                    ! Line break after 100 entries
                    line_break_cc = line_break_cc + 1
                    if(line_break_cc >= 100) then
                       write(unit = lout, fmt = *) 
                       line_break_cc = 0
                    end if
                 end do
              end do
           end do

           line_break_cc = 0
           write(unit = lout, fmt = *)


           do k = 1, k0
              do j = 1, j0
                 do i = 1, i0
                    write(unit = lout, fmt = '(i6)', advance='no') j
                    ! Line break after 100 entries
                    line_break_cc = line_break_cc + 1
                    if(line_break_cc >= 100) then
                       write(unit = lout, fmt = *) 
                       line_break_cc = 0
                    end if
                 end do
              end do
           end do
           
           line_break_cc = 0
           write(unit = lout, fmt = *) 

           do k = 1, k0
              do j = 1, j0
                 do i = 1, i0
                    write(unit = lout, fmt = '(i6)', advance='no') k
                    ! Line break after 100 entries
                    line_break_cc = line_break_cc + 1
                    if(line_break_cc >= 100) then
                       write(unit = lout, fmt = *) 
                       line_break_cc = 0
                    end if
                 end do
              end do
           end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           write(unit = lout, fmt = *) 
           write(unit = lout, fmt = *) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           line_break_cc = 0
           do k = 1, k0
              do j = 1, j0
                 do i = 1, i0
                    write(unit = lout, fmt = '(g15.6,1x)', advance='no') delxa(i) 
                    ! Line break after 100 entries
                    line_break_cc = line_break_cc + 1
                    if(line_break_cc >= 100) then
                       write(unit = lout, fmt = *) 
                       line_break_cc = 0
                    end if
                 end do
              end do
           end do

           line_break_cc = 0
           write(unit = lout, fmt = *)


           do k = 1, k0
              do j = 1, j0
                 do i = 1, i0
                    write(unit = lout, fmt = '(g15.6,1x)', advance='no') delya(j)
                    ! Line break after 100 entries
                    line_break_cc = line_break_cc + 1
                    if(line_break_cc >= 100) then
                       write(unit = lout, fmt = *) 
                       line_break_cc = 0
                    end if
                 end do
              end do
           end do

           line_break_cc = 0
           write(unit = lout, fmt = *) 

           do k = 1, k0
              do j = 1, j0
                 do i = 1, i0
                    write(unit = lout, fmt = '(g15.6,1x)', advance='no') delza(k)
                    ! Line break after 100 entries
                    line_break_cc = line_break_cc + 1
                    if(line_break_cc >= 100) then
                       write(unit = lout, fmt = *) 
                       line_break_cc = 0
                    end if
                 end do
              end do
           end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           write(unit = lout, fmt = *) 
           write(unit = lout, fmt = *) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           line_break_cc = 0
           do k = 1, k0
              do j = 1, j0
                 do i = 1, i0
                    write(unit = lout, fmt = '(i10)', advance='no') uindex(i,j,k)
                    ! Line break after 100 entries
                    line_break_cc = line_break_cc + 1
                    if(line_break_cc >= 100) then
                       write(unit = lout, fmt = *) 
                       line_break_cc = 0
                    end if
                 end do
              end do
           end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           write(unit = lout, fmt = *) 
           write(unit = lout, fmt = *) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           line_break_cc = 0
           do k = 1, k0
              do j = 1, j0
                 do i = 1, i0
                    write(unit = lout, fmt = '(es17.8,1x)', advance='no') head(i,j,k,ismpl)
                    ! Line break after 100 entries
                    line_break_cc = line_break_cc + 1
                    if(line_break_cc >= 100) then
                       write(unit = lout, fmt = *) 
                       line_break_cc = 0
                    end if
                 end do
              end do
           end do
           line_break_cc = 0
           write(unit = lout, fmt = *) 
           do k = 1, k0
              do j = 1, j0
                 do i = 1, i0
                    write(unit = lout, fmt = '(es17.8,1x)', advance='no') temp(i,j,k,ismpl)
                    ! Line break after 100 entries
                    line_break_cc = line_break_cc + 1
                    if(line_break_cc >= 100) then
                       write(unit = lout, fmt = *) 
                       line_break_cc = 0
                    end if
                 end do
              end do
           end do
           line_break_cc = 0
           write(unit = lout, fmt = *) 
           do k = 1, k0
              do j = 1, j0
                 do i = 1, i0
                    write(unit = lout, fmt = '(es17.8,1x)', advance='no') pres(i,j,k,ismpl)*pa_conv1
                    ! Line break after 100 entries
                    line_break_cc = line_break_cc + 1
                    if(line_break_cc >= 100) then
                       write(unit = lout, fmt = *) 
                       line_break_cc = 0
                    end if
                 end do
              end do
           end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           write(unit = lout, fmt = *) 
           write(unit = lout, fmt = *) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           line_break_cc = 0
           do k = 1, k0
              do j = 1, j0
                 do i = 1, i0
                    write(unit = lout, fmt = '(es15.6,1x)', advance='no') vxc(i,j,k,ismpl)
                    ! Line break after 100 entries
                    line_break_cc = line_break_cc + 1
                    if(line_break_cc >= 100) then
                       write(unit = lout, fmt = *) 
                       line_break_cc = 0
                    end if
                 end do
              end do
           end do
           line_break_cc = 0
           write(unit = lout, fmt = *) 
           do k = 1, k0
              do j = 1, j0
                 do i = 1, i0
                    write(unit = lout, fmt = '(es15.6,1x)', advance='no') vyc(i,j,k,ismpl)
                    ! Line break after 100 entries
                    line_break_cc = line_break_cc + 1
                    if(line_break_cc >= 100) then
                       write(unit = lout, fmt = *) 
                       line_break_cc = 0
                    end if
                 end do
              end do
           end do
           line_break_cc = 0
           write(unit = lout, fmt = *) 
           do k = 1, k0
              do j = 1, j0
                 do i = 1, i0
                    write(unit = lout, fmt = '(es15.6,1x)', advance='no') vzc(i,j,k,ismpl)
                    ! Line break after 100 entries
                    line_break_cc = line_break_cc + 1
                    if(line_break_cc >= 100) then
                       write(unit = lout, fmt = *) 
                       line_break_cc = 0
                    end if
                 end do
              end do
           end do
           line_break_cc = 0
           write(unit = lout, fmt = *) 
           do k = 1, k0
              do j = 1, j0
                 do i = 1, i0
                    write(unit = lout, fmt = '(es15.6,1x)', advance='no') por(i,j,k,ismpl)
                    ! Line break after 100 entries
                    line_break_cc = line_break_cc + 1
                    if(line_break_cc >= 100) then
                       write(unit = lout, fmt = *) 
                       line_break_cc = 0
                    end if
                 end do
              end do
           end do
           line_break_cc = 0
           write(unit = lout, fmt = *) 
           do k = 1, k0
              do j = 1, j0
                 do i = 1, i0
                    px = kx(i,j,k,ismpl)
                    if(klogflag) then
                       px = log10(kx(i,j,k,ismpl))
                    end if
                    write(unit = lout, fmt = '(es15.6,1x)', advance='no') px
                    ! Line break after 100 entries
                    line_break_cc = line_break_cc + 1
                    if(line_break_cc >= 100) then
                       write(unit = lout, fmt = *) 
                       line_break_cc = 0
                    end if
                 end do
              end do
           end do
           line_break_cc = 0
           write(unit = *, fmt = *) 
           do k = 1, k0
              do j = 1, j0
                 do i = 1, i0
                    py = ky(i,j,k,ismpl)
                    if(klogflag) then
                       py = log10(ky(i,j,k,ismpl))
                    end if
                    write(unit = lout, fmt = '(es15.6,1x)', advance='no') py
                    ! Line break after 100 entries
                    line_break_cc = line_break_cc + 1
                    if(line_break_cc >= 100) then
                       write(unit = lout, fmt = *) 
                       line_break_cc = 0
                    end if
                 end do
              end do
           end do
           line_break_cc = 0
           write(unit = *, fmt = *) 
           do k = 1, k0
              do j = 1, j0
                 do i = 1, i0
                    pz = kz(i,j,k,ismpl)
                    if(klogflag) then
                       pz = log10(kz(i,j,k,ismpl))
                    end if
                    write(unit = lout, fmt = '(es15.6,1x)', advance='no') pz
                    ! Line break after 100 entries
                    line_break_cc = line_break_cc + 1
                    if(line_break_cc >= 100) then
                       write(unit = lout, fmt = *) 
                       line_break_cc = 0
                    end if
                 end do
              end do
           end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           write(unit = lout, fmt = *) 
           write(unit = lout, fmt = *) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           line_break_cc = 0
           do k = 1, k0
              do j = 1, j0
                 do i = 1, i0
                    write(unit = lout, fmt = '(es15.6,1x)', advance='no') lx(i,j,k,ismpl)
                    ! Line break after 100 entries
                    line_break_cc = line_break_cc + 1
                    if(line_break_cc >= 100) then
                       write(unit = lout, fmt = *) 
                       line_break_cc = 0
                    end if
                 end do
              end do
           end do
           line_break_cc = 0
           write(unit = lout, fmt = *) 
           do k = 1, k0
              do j = 1, j0
                 do i = 1, i0
                    write(unit = lout, fmt = '(es15.6,1x)', advance='no') ly(i,j,k,ismpl)
                    ! Line break after 100 entries
                    line_break_cc = line_break_cc + 1
                    if(line_break_cc >= 100) then
                       write(unit = lout, fmt = *) 
                       line_break_cc = 0
                    end if
                 end do
              end do
           end do
           line_break_cc = 0
           write(unit = lout, fmt = *) 
           do k = 1, k0
              do j = 1, j0
                 do i = 1, i0
                    write(unit = lout, fmt = '(es15.6,1x)', advance='no') lz(i,j,k,ismpl)
                    ! Line break after 100 entries
                    line_break_cc = line_break_cc + 1
                    if(line_break_cc >= 100) then
                       write(unit = lout, fmt = *) 
                       line_break_cc = 0
                    end if
                 end do
              end do
           end do
           line_break_cc = 0
           write(unit = lout, fmt = *) 
           do k = 1, k0
              do j = 1, j0
                 do i = 1, i0
                    write(unit = lout, fmt = '(es15.6,1x)', advance='no') rhof(i,j,k,ismpl)
                    ! Line break after 100 entries
                    line_break_cc = line_break_cc + 1
                    if(line_break_cc >= 100) then
                       write(unit = lout, fmt = *) 
                       line_break_cc = 0
                    end if
                 end do
              end do
           end do
           line_break_cc = 0
           write(unit = lout, fmt = *) 
           do k = 1, k0
              do j = 1, j0
                 do i = 1, i0
                    write(unit = lout, fmt = '(es15.6,1x)', advance='no') visf(i,j,k,ismpl)
                    ! Line break after 100 entries
                    line_break_cc = line_break_cc + 1
                    if(line_break_cc >= 100) then
                       write(unit = lout, fmt = *) 
                       line_break_cc = 0
                    end if
                 end do
              end do
           end do
           write(unit = lout, fmt = *) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           write(unit = lout, fmt = *) 
           write(unit = lout, fmt = *) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           line_break_cc = 0
           do k = 1, k0
              do j = 1, j0
                 do i = 1, i0
                    write(unit = lout, fmt = '(es15.6,1x)', advance='no') qxc(i,j,k,ismpl)
                    ! Line break after 100 entries
                    line_break_cc = line_break_cc + 1
                    if(line_break_cc >= 100) then
                       write(unit = lout, fmt = *) 
                       line_break_cc = 0
                    end if
                 end do
              end do
           end do
           line_break_cc = 0
           write(unit = lout, fmt = *) 
           do k = 1, k0
              do j = 1, j0
                 do i = 1, i0
                    write(unit = lout, fmt = '(es15.6,1x)', advance='no') qyc(i,j,k,ismpl)
                    ! Line break after 100 entries
                    line_break_cc = line_break_cc + 1
                    if(line_break_cc >= 100) then
                       write(unit = lout, fmt = *) 
                       line_break_cc = 0
                    end if
                 end do
              end do
           end do
           line_break_cc = 0
           write(unit = lout, fmt = *) 
           do k = 1, k0
              do j = 1, j0
                 do i = 1, i0
                    write(unit = lout, fmt = '(es15.6,1x)', advance='no') qzc(i,j,k,ismpl)
                    ! Line break after 100 entries
                    line_break_cc = line_break_cc + 1
                    if(line_break_cc >= 100) then
                       write(unit = lout, fmt = *) 
                       line_break_cc = 0
                    end if
                 end do
              end do
           end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           write(unit = lout, fmt = *) 
           write(unit = lout, fmt = *) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           line_break_cc = 0
           do k = 1, k0
              do j = 1, j0
                 do i = 1, i0
                    write(unit = lout, fmt = '(es15.6,1x)', advance='no') w(i,j,k,ismpl)
                    ! Line break after 100 entries
                    line_break_cc = line_break_cc + 1
                    if(line_break_cc >= 100) then
                       write(unit = lout, fmt = *) 
                       line_break_cc = 0
                    end if
                 end do
              end do
           end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        else
           
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           
           
           WRITE(lout,'(9A)') 'variables = "i",  "j",  "k", "x",           &
                & "y",            "z", "uindex",', &
                '           "head",           "temp",', &
                '           "pres",', &
                '           "vx",           "vy",           "vz",', &
                '          "por",           "kx",           "ky",', &
                '           "kz",           "lx",           "ly",', &
                '           "lz",         "rhof",         "visf",', &
                '          "qxc",          "qyc",          "qzc",', &
                '        "sigma",           "lc",          "src",'

           !     use sline for ZONE string
           i1 = locstr(project_sfx(ismpl),'_time_')
           i2 = locstr(project_sfx(ismpl),'_out')
           IF (i1>=1 .AND. i2>=1) THEN
!!$              WRITE(sline,'(1A,1I7,3A,1I7)') 'zone T="', ident, &
!!$                   '.time step", SolutionTime=', project_sfx(ismpl) (i1+6:i2- &
!!$                   1), ', StrandID=', 1
              WRITE(sline,'(1A,1I7)') 'zone StrandID=', 1
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
                    px = kx(i,j,k,ismpl)
                    py = ky(i,j,k,ismpl)
                    pz = kz(i,j,k,ismpl)
                    IF (klogflag) THEN
                       px = log10(kx(i,j,k,ismpl))
                       py = log10(ky(i,j,k,ismpl))
                       pz = log10(kz(i,j,k,ismpl))
                    END IF
                    WRITE(lout, &
                         '(3i6,3(g15.6,1x),i10,3(e17.8,1x),19(e15.6,1x))') &
                         i, j, k, delxa(i), delya(j), delza(k), uindex(i,j,k), &
                         head(i,j,k,ismpl), temp(i,j,k,ismpl), &
                         pres(i,j,k,ismpl)*pa_conv1, &
                         vxc(i,j,k,ismpl), vyc(i,j,k,ismpl), &
                         vzc(i,j,k,ismpl), por(i,j,k,ismpl), px, py, pz, &
                         lx(i,j,k,ismpl), ly(i,j,k,ismpl), lz(i,j,k,ismpl), &
                         rhof(i,j,k,ismpl), visf(i,j,k,ismpl), &
                         qxc(i,j,k,ismpl), qyc(i,j,k,ismpl), qzc(i,j,k,ismpl), &
                         w(i,j,k,ismpl)
                 END DO
              END DO
           END DO
           
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        end if
           
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CLOSE(lout)
        CALL omp_del_file_handler(lout)
        CALL compress_file(compress_out,filename)

        RETURN
      END
