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

!>    @brief OpenMP wrapper for "omp_courant"
!>    @param[out] courant_max global courant number
!>    @param[in] ismpl local sample index
      SUBROUTINE courant(courant_max,ismpl)
        use mod_genrl
        use mod_OMP_TOOLS
        IMPLICIT NONE
        INCLUDE 'OMP_TOOLS.inc'
        DOUBLE PRECISION courant_max
        integer :: ismpl

#ifdef fOMP
!$OMP parallel num_threads(Tlevel_1)
!$      call omp_binding(ismpl)
#endif
        CALL omp_courant(courant_max,ismpl)
#ifdef fOMP
!$OMP end parallel
#endif

        RETURN
      END

!>    @brief calculate grid courant number
!>    @param[out] courant_max global courant number
!>    @param[in] ismpl local sample index
      SUBROUTINE omp_courant(courant_max,ismpl)
        use arrays
        use mod_genrl
        use mod_genrlc
        use mod_temp
        use mod_time
        use mod_linfos
        IMPLICIT NONE
        integer :: i, j, k
        integer :: ismpl

        DOUBLE PRECISION courant_maxx, courant_minx, courant_avgx
        DOUBLE PRECISION courant_maxy, courant_miny, courant_avgy, courant_maxz
        DOUBLE PRECISION courant_minz, courant_avgz
        DOUBLE PRECISION courant_max, fac, delt, deltat, davg, val, min_veloc
        EXTERNAL deltat
!     min. value, ignore courant numbers for lower velocities
        PARAMETER (min_veloc=1.0D-10)
!
        INTEGER c1, c2, c3
        DOUBLE PRECISION vx, vy, vz, por
        EXTERNAL vx, vy, vz, por

        delt = deltat(simtime(ismpl),ismpl)
!
        IF ( .NOT. (transient .AND. tr_switch(ismpl))) THEN
!$OMP master
          WRITE(*,*) ' courant: not defined for steady state'
!$OMP end master
          RETURN
        END IF
!
        c1 = 0
        courant_maxx = small
        courant_minx = big
        courant_avgx = 0.0D0
!$OMP do schedule(static)
        DO k = 1, k0
          DO j = 1, j0
            DO i = 2, i0 - 1
              val = abs(vx(i,j,k,ismpl))
              IF (val>=min_veloc) THEN
                c1 = c1 + 1
                davg = 0.5D0*(delx(i)+delx(i+1))
                fac = delt/(davg*por(i,j,k,ismpl))
                val = val*fac
                IF (val>courant_maxx) courant_maxx = val
                IF (val<courant_minx) courant_minx = val
                courant_avgx = courant_avgx + val
              END IF
            END DO
          END DO
        END DO
!$OMP end do nowait
!
        c2 = 0
        courant_maxy = small
        courant_miny = big
        courant_avgy = 0.0D0
!$OMP do schedule(static)
        DO k = 1, k0
          DO j = 2, j0 - 1
            DO i = 1, i0
              val = abs(vy(i,j,k,ismpl))
              IF (val>=min_veloc) THEN
                c2 = c2 + 1
                davg = 0.5D0*(dely(j)+dely(j+1))
                fac = delt/(davg*por(i,j,k,ismpl))
                val = val*fac
                IF (val>courant_maxy) courant_maxy = val
                IF (val<courant_miny) courant_miny = val
                courant_avgy = courant_avgy + val
              END IF
            END DO
          END DO
        END DO
!$OMP end do nowait
!
        c3 = 0
        courant_maxz = small
        courant_minz = big
        courant_avgz = 0.0D0
!$OMP do schedule(static)
        DO k = 2, k0 - 1
          DO j = 1, j0
            DO i = 1, i0
              val = abs(vz(i,j,k,ismpl))
              IF (val>=min_veloc) THEN
                c3 = c3 + 1
                davg = 0.5D0*(delz(k)+delz(k+1))
                fac = delt/(davg*por(i,j,k,ismpl))
                val = val*fac
                IF (val>courant_maxz) courant_maxz = val
                IF (val<courant_minz) courant_minz = val
                courant_avgz = courant_avgz + val
              END IF
            END DO
          END DO
        END DO
!$OMP end do nowait
!
!     compute global sum for all values
        CALL omp_summe(courant_maxx,courant_minx,courant_avgx, &
          courant_maxy,courant_miny,courant_avgy,courant_maxz, &
          courant_minz,courant_avgz,c1,c2,c3,ismpl)
!
!$OMP master
        IF (i0>2) THEN
          courant_avgx = courant_avgx/dble(c1)
        ELSE
          courant_maxx = 0.0D0
          courant_minx = 0.0D0
          courant_avgx = 0.0D0
        END IF
        IF (j0>2) THEN
          courant_avgy = courant_avgy/dble(c2)
        ELSE
          courant_maxy = 0.0D0
          courant_miny = 0.0D0
          courant_avgy = 0.0D0
        END IF
        IF (k0>2) THEN
          courant_avgz = courant_avgz/dble(c3)
        ELSE
          courant_maxz = 0.0D0
          courant_minz = 0.0D0
          courant_avgz = 0.0D0
        END IF
!
        courant_max = max(courant_maxx,courant_maxy,courant_maxz)
!
        IF (linfos(3)>=2) THEN
          WRITE(*,*) 'courant-number in x,y,z:'
          WRITE(*,'(a,1e10.3,a,1e10.3,a,1e10.3)') '  max. : ', &
            courant_maxx, ',    ', courant_maxy, ',    ', courant_maxz
          WRITE(*,'(a,1e10.3,a,1e10.3,a,1e10.3)') '  min. : ', &
            courant_minx, ',    ', courant_miny, ',    ', courant_minz
          WRITE(*,'(a,1e10.3,a,1e10.3,a,1e10.3,1a,3I8)') '  avg. : ', &
            courant_avgx, ',    ', courant_avgy, ',    ', &
            courant_avgz, ', #', c1, c2, c3
        END IF
!
        IF (linfos(3)>=1 .AND. courant_max>1.D0) THEN
          WRITE(*,'(a)') '!!!: courant number(s) greater than 1 :'
          WRITE(*,'(a,1e12.3,a,1e10.3,a,1e10.3)') 'x: ', &
            courant_maxx, 'y: ', courant_maxy, 'z: ', courant_maxz
          WRITE(*,*)
        END IF
!$OMP end master

        RETURN
      END
