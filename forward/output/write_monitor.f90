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
!>    @param[in] otype 1: new file, 2 append
!>    @param[in] ismpl local sample index
      SUBROUTINE write_monitor(otype,ismpl)
        use arrays
        use mod_genrl
        use mod_genrlc
        use mod_data
        use mod_time
        use mod_flow
        use mod_temp
        use mod_conc
        use mod_linfos
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k, m
        integer :: ib
        DOUBLE PRECISION dx, dy, dz, val
!     otype:
!       1  new file
!       2  append
        INTEGER otype
!     mm_*: outer iteration loop (file for each monitor point)
!     mt_*: inner iteration loop (file for each time step)
        INTEGER mm_b, mm_e, mt_b, mt_e
        INTEGER i1, i2, i3, m2, lblank, lout, in1, in2
        EXTERNAL lblank
!     external functions
        DOUBLE PRECISION vxc, vyc, vzc, kx, ky, kz, lx, ly, lz, por, &
          rhof, visf
        EXTERNAL vxc, vyc, vzc, kx, ky, kz, lx, ly, lz, por, rhof, &
          visf
        DOUBLE PRECISION qxc, qyc, qzc, bhpr
        EXTERNAL qxc, qyc, qzc, bhpr
!     normal or transposed, monitor point files, time step files
        LOGICAL ltransposed, lmonip, ltimes

        character (len=256) :: filename
        character (len=65536) :: line
        character (len=8) :: snumber
        INTRINSIC trim


        IF ( .NOT. (transient .AND. tr_switch(ismpl))) RETURN
        IF ( .NOT. monitor) RETURN
        IF ( .NOT. write_smonitor) THEN
          IF (write_disable) RETURN
        END IF
#ifdef NOMON
        RETURN
#endif

        snumber = '0'
        IF (write_smonitor) THEN 
          WRITE(snumber,'(1I7)') smon_idx(ismpl)
          IF (smon_idx(ismpl)>=0) THEN
            WRITE(snumber,'(1I7)') smon_idx(ismpl)
          ELSE IF (smon_idx(ismpl)==-1) THEN
            WRITE(snumber,'(A8)') 'final'
          ELSE IF (smon_idx(ismpl)==-2) THEN
            WRITE(snumber,'(A8)') 'debug'
          ELSE IF (smon_idx(ismpl)==-3) THEN
            WRITE(snumber,'(A8)') 'ens_mean'
          ELSE IF (smon_idx(ismpl)==-4) THEN
            WRITE(snumber,'(A8)') 'mean'
          ELSE IF (smon_idx(ismpl)==-5) THEN
            WRITE(snumber,'(A8)') 'ens_mean'
          END IF
        END IF
        CALL chln(snumber,in1,in2)

        CALL omp_new_file_handler(lout,16)

!     out_orientation:
!       0  normal (pv column wise), each time step has its own block (mp row wise)
!       1  transposed (pv row wise), each time step has its own block (mp column wise)
!       2  normal (pv column wise), new file for each time step (mp row wise)
!       3  transposed (pv row wise), new file for each time step (mp column wise)
!       4  normal (pv column wise), new file for each monitor point (time row wise)
        IF (out_orientation==1 .OR. out_orientation==3) THEN
          ltransposed = .TRUE.
        ELSE
          ltransposed = .FALSE.
        END IF
        IF (out_orientation==4) THEN
!        setup outer iteration
          mm_b = 1
          mm_e = nmon
!        inner loop makes only one step, computed later
          mt_b = 0
          mt_e = 0
          lmonip = .TRUE.
        ELSE
!        outer loop makes only one dummy step
          mm_b = 1
          mm_e = 1
!        setup inner iteration
          mt_b = 1
          mt_e = nmon
          lmonip = .FALSE.
        END IF
        IF (out_orientation==2 .OR. out_orientation==3) THEN
          ltimes = .TRUE.
        ELSE
          ltimes = .FALSE.
        END IF

        IF (nmon>9999 .AND. ltransposed) THEN
          WRITE(*,'(1A)') &
            'error: number of monitor points limited to 9999!'
          STOP
        END IF

        CALL chln(project,i1,i2)
!     steady state - file name
        WRITE(filename,'(3A)') project(i1:i2),trim(project_sfx(ismpl)), '_monitor.dat'
        IF (write_smonitor) WRITE(filename,'(5A)') project(i1:i2),trim(project_sfx(ismpl)), &
          '_monitor_', snumber(in1:in2), '.dat'
!     transient - file name
        IF (ltimes) THEN
          WRITE(filename,'(3A,1e15.9,1A)') project(i1:i2),trim(project_sfx(ismpl)), '_', &
            (simtime(ismpl))/tunit, '_monitor.dat'
          IF (write_smonitor) WRITE(filename,'(3A,1e15.9,3A)') &
            project(i1:i2),trim(project_sfx(ismpl)), '_', (simtime(ismpl))/tunit, &
            '_monitor_', snumber(in1:in2), '.dat'
        END IF

        IF (linfos(3)>=2 .AND. .NOT. lmonip) THEN
          WRITE(*,'(3A)') '  [W] : Monitor points to "', &
            filename(1:lblank(filename)), '"'
        END IF

! ----------------------------------------------
        DO m = mm_b, mm_e
          IF (lmonip) THEN
!           compute inner loop, one step with the right index
            mt_b = m
            mt_e = m
!           monitor point - file name
            CALL chln(project,i1,i2)
            WRITE(filename,'(3A,1I6.6,1A)') project(i1:i2),trim(project_sfx(ismpl)), '_', m, &
              '_monitor.dat'
            IF (write_smonitor) WRITE(filename,'(3A,1I6.6,3A)') &
              project(i1:i2),trim(project_sfx(ismpl)), '_', m, '_monitor_', snumber(in1:in2), &
              '.dat'
            IF (linfos(3)>=2) THEN
              WRITE(*,'(3A)') '  [W] : Monitor points to "', &
                filename(1:lblank(filename)), '"'
            END IF
          END IF

!        Header
          IF ((otype==1 .OR. ltimes) .AND. .NOT. ltransposed) THEN
            WRITE(line,'(9999(1A11,1I4.4,1A2))') &
              ('      "conc', i, '",',i=1,ntrans)
            OPEN(lout,file=filename,status='replace')
            WRITE(lout,'(12A)') '%       "time"           "x" ', &
              '           "y"           "z"      "uindex" ', &
              '   "i"  "j"  "k"', '          "head" ', &
              '          "temp" ', '          "pres" ', &
              line(1:ntrans*17), &
              '            "vx" ', '            "vy" ', &
              '            "vz" ', '         "bhpr"', &
              '            "kz" '
            CLOSE(lout)
          END IF
          IF ((otype==1 .OR. ltimes) .AND. ltransposed) THEN
            OPEN(lout,file=filename,status='replace')
            WRITE(lout,'(1A)') '% "x"'
            WRITE(lout,'(1A)') '% "y"'
            WRITE(lout,'(1A)') '% "z"'
            WRITE(lout,'(1A)') '% "uindex"'
            WRITE(lout,'(1A)') '% "i"'
            WRITE(lout,'(1A)') '% "j"'
            WRITE(lout,'(1A)') '% "k"'
            DO m2 = mt_b, mt_e
              i1 = imon(m2,1)
              dx = 0.5D0*delx(1)
              DO i = 2, i1
                dx = dx + 0.5D0*(delx(i-1)+delx(i))
              END DO
              WRITE(lout,'(1e14.6,1X)',advance='NO') dx
            END DO
            WRITE(lout,*) ' '
            DO m2 = mt_b, mt_e
              i2 = imon(m2,2)
              dy = 0.5D0*dely(1)
              DO j = 2, i2
                dy = dy + 0.5D0*(dely(j-1)+dely(j))
              END DO
              WRITE(lout,'(1e14.6,1X)',advance='NO') dy
            END DO
            WRITE(lout,*) ' '
            DO m2 = mt_b, mt_e
              i3 = imon(m2,3)
              dz = 0.5D0*delz(1)
              DO k = 2, i3
                dz = dz + 0.5D0*(delz(k-1)+delz(k))
              END DO
              WRITE(lout,'(1e14.6,1X)',advance='NO') dz
            END DO
            WRITE(lout,*) ' '
            WRITE(lout,'(9999(I14,1X))') (uindex(imon(m2,1),imon(m2,2),imon(m2,3)),m2=mt_b,mt_e)
            WRITE(lout,'(9999(I14,1X))') (imon(m2,1),m2=mt_b,mt_e)
            WRITE(lout,'(9999(I14,1X))') (imon(m2,2),m2=mt_b,mt_e)
            WRITE(lout,'(9999(I14,1X))') (imon(m2,3),m2=mt_b,mt_e)
            WRITE(lout,'(1A)') '% "head"'
            WRITE(lout,'(1A)') '% "temp"'
            WRITE(lout,'(1A)') '% "pres"'
            DO i = 1, ntrans
              WRITE(lout,'(1A7,1I4.4,1A1)') '% "conc', i, '"'
            END DO
            WRITE(lout,'(1A)') '%   "vx"'
            WRITE(lout,'(1A)') '%   "vy"'
            WRITE(lout,'(1A)') '%   "vz"'
            CLOSE(lout)
          END IF

!        append data block
          OPEN(lout,file=filename,status='unknown',position='append')

!        Body
          IF ( .NOT. ltransposed) THEN
!          - normal orientation -
!          separate each time steps
            IF ( .NOT. ltimes .AND. .NOT. lmonip) WRITE(lout,'(A)') '%'
!          - inner loop -
            DO m2 = mt_b, mt_e
              i3 = imon(m2,3)
              i2 = imon(m2,2)
              i1 = imon(m2,1)
              dz = 0.5D0*delz(1)
!
              val = 0.0d0
              DO ib = first_flow, last_flow
                i = ibc_data(ib,cbc_i)
                j = ibc_data(ib,cbc_j)
                k = ibc_data(ib,cbc_k)
!              "neumann"?, skip otherwise
                IF (ibc_data(ib,cbc_bt)==bt_neuw.AND.i==i1.AND.j==i2.AND.k==i3) THEN
                  IF (ibc_data(ib,cbc_bcu)>0) THEN
                    WRITE(*,*) 'error: well function can not be defined in a BC-unit!'
                    STOP
                  END IF
!better-recompute-all-times instead of: val = dbc_data(ib,3,ismpl)
                  val = bhpr(i,j,k,ismpl)
                END IF
              END DO
!
              WRITE(lout, &
                '(1e17.9,3(e14.6,1X),1I9,3I6,9999(e16.8,1X))') &
                (simtime(ismpl))/tunit, delxa(i1), delya(i2), delza(i3), &
                uindex(i1,i2,i3), i1, i2, i3, head(i1,i2,i3,ismpl), &
                temp(i1,i2,i3,ismpl), pres(i1,i2,i3,ismpl)*pa_conv1, &
                (conc(i1,i2,i3,i,ismpl),i=1, &
                ntrans), vxc(i1,i2,i3,ismpl), &
                vyc(i1,i2,i3,ismpl), vzc(i1,i2,i3,ismpl),val*pa_conv1, &
                kz(i1,i2,i3,ismpl)
            END DO
          ELSE
!           - transposed orientation -
!           time comment
            WRITE(lout,'(1A,1e17.9)') '% "time: "', (simtime(ismpl))/tunit
            WRITE(lout,'(9999(e14.6,1X))') (head(imon(m2,1),imon(m2,2),imon(m2,3),ismpl),m2=mt_b,mt_e)
            WRITE(lout,'(9999(e14.6,1X))') (temp(imon(m2,1),imon(m2,2),imon(m2,3),ismpl),m2=mt_b,mt_e)
            WRITE(lout,'(9999(e14.6,1X))') (pres(imon(m2,1),imon(m2,2),imon(m2,3),ismpl)*pa_conv1,m2=mt_b,mt_e)
            DO i = 1, ntrans
              WRITE(lout,'(9999(e14.6,1X))') (conc(imon(m2,1),imon(m2,2),imon(m2,3),i,ismpl),m2=mt_b,mt_e)
            END DO
            WRITE(lout,'(9999(e14.6,1X))') (vxc(imon(m2,1),imon(m2,2),imon(m2,3),ismpl),m2=mt_b,mt_e)
            WRITE(lout,'(9999(e14.6,1X))') (vyc(imon(m2,1),imon(m2,2),imon(m2,3),ismpl),m2=mt_b,mt_e)
            WRITE(lout,'(9999(e14.6,1X))') (vzc(imon(m2,1),imon(m2,2),imon(m2,3),ismpl),m2=mt_b,mt_e)
          END IF

          CLOSE(lout)
        END DO
! ----------------------------------------------
        CALL omp_del_file_handler(lout)

        RETURN
      END
