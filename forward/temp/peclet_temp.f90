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

!>    @brief OpenMP wrapper for "omp_peclet_temp"
!>    @param[out] peclet_max peclet criteria
!>    @param[in] ismpl local sample index
      SUBROUTINE peclet_temp(peclet_max,ismpl)
        use mod_genrl
        use mod_OMP_TOOLS
        IMPLICIT NONE
        integer :: ismpl
        INCLUDE 'OMP_TOOLS.inc'
        DOUBLE PRECISION peclet_max

#ifdef fOMP
!$OMP parallel num_threads(Tlevel_1)
!$      call omp_binding(ismpl)
#endif
        CALL omp_peclet_temp(peclet_max,ismpl)
#ifdef fOMP
!$OMP end parallel
#endif
!
        RETURN
      END

!>    @brief calculate grid peclet numbers (temperature)
!>    @param[out] peclet_max peclet criteria
!>    @param[in] ismpl local sample index
      SUBROUTINE omp_peclet_temp(peclet_max,ismpl)
        use arrays
        use mod_genrl
        use mod_genrlc
        use mod_temp
        use mod_linfos
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k


        INTEGER c1, c2, c3
        DOUBLE PRECISION peclet_maxx, peclet_minx, peclet_avgx, &
          peclet_maxy, peclet_miny, peclet_avgy, peclet_maxz, &
          peclet_minz, peclet_avgz, val, davg, peclet_max
        INTEGER ipt, jpt, kpt
        DOUBLE PRECISION li, lj, lk, rhocf, vx, vy, vz, por
        EXTERNAL li, lj, lk, rhocf, vx, vy, vz, por

        IF (linfos(3)>=2) THEN
!$OMP master
          WRITE(*,*)
          WRITE(*,'(A,1e16.8)') '  ... peclet-temp'
          WRITE(*,*)
!$OMP end master
        END IF
!
! temperature-val in x
        peclet_maxx = small
        peclet_minx = big
        peclet_avgx = 0.0D0
        val = 0.0D0
        ipt = 0
        jpt = 0
        kpt = 0
        c1 = 0
!$OMP do schedule(static)
        DO k = 1, k0
          DO j = 1, j0
            DO i = 2, i0 - 1
              c1 = c1 + 1
              davg = 0.5D0*(delx(i)+delx(i+1))
              val = abs(vx(i,j,k,ismpl))*davg/li(i,j,k,ismpl)
              IF (val>peclet_maxx) THEN
                peclet_maxx = val
                ipt = i
                jpt = j
                kpt = k
              END IF
              IF (val<peclet_minx) peclet_minx = val
              peclet_avgx = peclet_avgx + val
            END DO
          END DO
        END DO
!$OMP end do nowait
!      if (linfos(3).ge.2)write(*,*)
!     &      "max. temp-val in x:  ", ipt,jpt,kpt
!
! temperature-val in y
        peclet_maxy = small
        peclet_miny = big
        peclet_avgy = 0.0D0
        val = 0.0D0
        ipt = 0
        jpt = 0
        kpt = 0
        c2 = 0
!$OMP do schedule(static)
        DO k = 1, k0
          DO j = 2, j0 - 1
            DO i = 1, i0
              c2 = c2 + 1
              davg = 0.5D0*(dely(j)+dely(j+1))
              val = abs(vy(i,j,k,ismpl))*davg/lj(i,j,k,ismpl)
              IF (val>peclet_maxy) THEN
                peclet_maxy = val
                ipt = i
                jpt = j
                kpt = k
              END IF
              IF (val<peclet_miny) peclet_miny = val
              peclet_avgy = peclet_avgy + val
            END DO
          END DO
        END DO
!$OMP end do nowait
!      if(linfos(3).ge.2)
!     &    write(*,*)"max. temp-val in y:  " ,ipt,jpt,kpt
!
!  temperature-val in z
        peclet_maxz = small
        peclet_minz = big
        peclet_avgz = 0.0D0
        val = 0.0D0
        ipt = 0
        jpt = 0
        kpt = 0
        c3 = 0
!$OMP do schedule(static)
        DO k = 2, k0 - 1
          DO j = 1, j0
            DO i = 1, i0
              c3 = c3 + 1
              davg = 0.5D0*(delz(k)+delz(k+1))
              val = abs(vz(i,j,k,ismpl))*davg/lk(i,j,k,ismpl)
              IF (val>peclet_maxz) THEN
                peclet_maxz = val
                ipt = i
                jpt = j
                kpt = k
              END IF
              IF (val<peclet_minz) peclet_minz = val
              peclet_avgz = peclet_avgz + val
            END DO
          END DO
        END DO
!$OMP end do nowait
!      if(linfos(3).ge.2)
!     &    write(*,*)"max. temp-val in z:  " ,ipt,jpt,kpt
!
!     compute global sum for all values
        CALL omp_summe(peclet_maxx,peclet_minx,peclet_avgx, &
          peclet_maxy,peclet_miny,peclet_avgy,peclet_maxz,peclet_minz, &
          peclet_avgz,c1,c2,c3,ismpl)
!
!$OMP master
        IF (i0>2) THEN
          peclet_avgx = peclet_avgx/dble(c1)
        ELSE
          peclet_maxx = 0.0D0
          peclet_minx = 0.0D0
          peclet_avgx = 0.0D0
        END IF
        IF (j0>2) THEN
          peclet_avgy = peclet_avgy/dble(c2)
        ELSE
          peclet_maxy = 0.0D0
          peclet_miny = 0.0D0
          peclet_avgy = 0.0D0
        END IF
        IF (k0>2) THEN
          peclet_avgz = peclet_avgz/dble(c3)
        ELSE
          peclet_maxz = 0.0D0
          peclet_minz = 0.0D0
          peclet_avgz = 0.0D0
        END IF
!
        peclet_max = max(peclet_avgx,peclet_avgy,peclet_avgz)
!
        IF (linfos(3)>=2) THEN
          WRITE(*,*) 'peclet number for temperature in x,y,z:'
          WRITE(*,'(a,1e12.3,a,1e10.3,a,1e10.3)') '  max. : ', &
            peclet_maxx, ',    ', peclet_maxy, ',    ', peclet_maxz
          WRITE(*,'(a,1e12.3,a,1e10.3,a,1e10.3)') '  min. : ', &
            peclet_minx, ',    ', peclet_miny, ',    ', peclet_minz
          WRITE(*,'(a,1e12.3,a,1e10.3,a,1e10.3)') '  avg. : ', &
            peclet_avgx, ',    ', peclet_avgy, ',    ', peclet_avgz
        END IF
!
        IF (peclet_max>2.0D0 .AND. linfos(3)>=1) THEN
          WRITE(*,'(1A)') '!!!: peclet number for temperature > 2 :'
          WRITE(*,'(a,1e12.3,a,1e10.3,a,1e10.3)') 'x: ', peclet_maxx, &
            'y: ', peclet_maxy, 'z: ', peclet_maxz
        END IF
!$OMP end master
!
        RETURN
      END
