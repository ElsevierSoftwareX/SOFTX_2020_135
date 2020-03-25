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

!>    @brief OpenMP wrapper for "omp_neumann_temp"
!>    @param[out] neumann_max neumann criteria
!>    @param[in] ismpl local sample index
      SUBROUTINE neumann_temp(neumann_max,ismpl)
        use mod_genrl
        use mod_OMP_TOOLS
        IMPLICIT NONE
        integer :: ismpl

        INCLUDE 'OMP_TOOLS.inc'
        DOUBLE PRECISION neumann_max

#ifdef fOMP
!$OMP parallel num_threads(Tlevel_1)
!$      call omp_binding(ismpl)
#endif
        CALL omp_neumann_temp(neumann_max,ismpl)
#ifdef fOMP
!$OMP end parallel
#endif

        RETURN
      END

!>    @brief calculate grid neuman numbers (temp)
!>    @param[out] neumann_max maximal neuman number
!>    @param[in] ismpl local sample index
      SUBROUTINE omp_neumann_temp(neumann_max,ismpl)
        use arrays
        use mod_genrl
        use mod_genrlc
        use mod_temp
        use mod_time
        use mod_linfos
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k



        INTEGER c1, c2, c3
        DOUBLE PRECISION neumann_maxx, neumann_minx, neumann_avgx
        DOUBLE PRECISION neumann_maxy, neumann_miny, neumann_avgy
        DOUBLE PRECISION neumann_maxz, neumann_minz, neumann_avgz
        DOUBLE PRECISION val, neumann_max, delt, fac, davg
        DOUBLE PRECISION deltat
        EXTERNAL deltat
        DOUBLE PRECISION li, lj, lk, rhocf, por, rhoceff
        EXTERNAL li, lj, lk, rhocf, por, rhoceff


        delt = deltat(simtime(ismpl),ismpl)

        IF ( .NOT. (transient .AND. tr_switch(ismpl))) THEN
!$OMP master
          WRITE(*,*) ' neumann-temp: not defined for steady state'
!$OMP end master
          RETURN
        ELSE IF (linfos(3)>=2) THEN
!$OMP master
          WRITE(*,*)
          WRITE(*,'(A,1e16.8)') '  ... neumann-temp:  delt/tunit = ', &
            delt/tunit
          WRITE(*,*)
!$OMP end master
        END IF

! val in x
        c1 = 0
        neumann_maxx = small
        neumann_minx = big
        neumann_avgx = 0.0D0
!$OMP do schedule(static)
        DO k = 1, k0
          DO j = 1, j0
            DO i = 2, i0 - 1
              c1 = c1 + 1
              davg = 0.5D0*(delx(i)+delx(i+1))
              fac = delt/rhoceff(i,j,k,ismpl)
              val = fac*li(i,j,k,ismpl)/(davg*davg)
              IF (val>neumann_maxx) neumann_maxx = val
              IF (val<neumann_minx) neumann_minx = val
              neumann_avgx = neumann_avgx + val
            END DO
          END DO
        END DO
!$OMP end do nowait

! val in y
        c2 = 0
        neumann_maxy = small
        neumann_miny = big
        neumann_avgy = 0.0D0
!$OMP do schedule(static)
        DO k = 1, k0
          DO j = 2, j0 - 1
            DO i = 1, i0
              c2 = c2 + 1
              davg = 0.5D0*(dely(j)+dely(j+1))
              fac = delt/rhoceff(i,j,k,ismpl)
              val = fac*lj(i,j,k,ismpl)/(davg*davg)
              IF (val>neumann_maxy) neumann_maxy = val
              IF (val<neumann_miny) neumann_miny = val
              neumann_avgy = neumann_avgy + val
            END DO
          END DO
        END DO
!$OMP end do nowait

! val in z
        c3 = 0
        neumann_maxz = small
        neumann_minz = big
        neumann_avgz = 0.0D0
!$OMP do schedule(static)
        DO k = 2, k0 - 1
          DO j = 1, j0
            DO i = 1, i0
              c3 = c3 + 1
              davg = 0.5D0*(delz(k)+delz(k+1))
              fac = delt/rhoceff(i,j,k,ismpl)
              val = fac*lk(i,j,k,ismpl)/(davg*davg)
              IF (val>neumann_maxz) neumann_maxz = val
              IF (val<neumann_minz) neumann_minz = val
              neumann_avgz = neumann_avgz + val
            END DO
          END DO
        END DO
!$OMP end do nowait

!     compute global sum for all values
        CALL omp_summe(neumann_maxx,neumann_minx,neumann_avgx, &
          neumann_maxy,neumann_miny,neumann_avgy,neumann_maxz, &
          neumann_minz,neumann_avgz,c1,c2,c3,ismpl)

!$OMP master
        IF (i0>2) THEN
          neumann_avgx = neumann_avgx/dble(c1)
        ELSE
          neumann_maxx = 0.0D0
          neumann_minx = 0.0D0
          neumann_avgx = 0.0D0
        END IF
        IF (j0>2) THEN
          neumann_avgy = neumann_avgy/dble(c2)
        ELSE
          neumann_maxy = 0.0D0
          neumann_miny = 0.0D0
          neumann_avgy = 0.0D0
        END IF
        IF (k0>2) THEN
          neumann_avgz = neumann_avgz/dble(c3)
        ELSE
          neumann_maxz = 0.0D0
          neumann_minz = 0.0D0
          neumann_avgz = 0.0D0
        END IF

        neumann_max = max(neumann_maxx,neumann_maxy,neumann_maxz)

        IF (linfos(3)>=2) THEN
          WRITE(*,*) 'neumann number for temperature in x,y,z:'
          WRITE(*,'(a,1e10.3,a,1e10.3,a,1e10.3)') '  max. : ', &
            neumann_maxx, ',    ', neumann_maxy, ',    ', neumann_maxz
          WRITE(*,'(a,1e10.3,a,1e10.3,a,1e10.3)') '  min. : ', &
            neumann_minx, ',    ', neumann_miny, ',    ', neumann_minz
          WRITE(*,'(a,1e10.3,a,1e10.3,a,1e10.3)') '  avg. : ', &
            neumann_avgx, ',    ', neumann_avgy, ',    ', neumann_avgz
        END IF

        IF (linfos(3)>=1 .AND. neumann_max>1.D0) THEN
          WRITE(*,'(a)') &
            '!!!: neumann temp number(s) greater than 1 :'
          WRITE(*,'(a,1e12.3,a,1e10.3,a,1e10.3)') 'x: ', &
            neumann_maxx, 'y: ', neumann_maxy, 'z: ', neumann_maxz
          WRITE(*,*)
        END IF
!$OMP end master

        RETURN
      END
