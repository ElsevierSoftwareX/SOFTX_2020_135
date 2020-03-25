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

!>    @brief domain of validity for module ice
!>    @param[in] ismpl local sample index
      SUBROUTINE check_domain(ismpl)
        use arrays
        use mod_genrl
        use mod_genrlc
        use mod_conc
        use mod_linfos
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k, l


!     counts the directions (thickness large enough)
        INTEGER icountp, icountt, icountc
!     [csmin] reasonable physical value when not [cmin]
        DOUBLE PRECISION pmin, pmax, tmin, tmax, cmin, csmin, cmax, &
          hmin, hmax, dpmax, dtmax, dcmax, dhmax, dpmin, dtmin, dcmin, &
          dhmin
        PARAMETER (pmin=0.01D6,pmax=150.D6,tmin=-60.0D0,tmax=350.0, &
          cmin=0.D0,csmin=1.D-30,cmax=1.D0)
        INTRINSIC trim


        icountp = 0
        icountt = 0
        icountc = 0
        dpmax = pmax
        dpmin = pmin
        dtmax = tmax
        dtmin = tmin
        dcmax = cmax
        dcmin = cmin

        DO k = 1, k0
          DO j = 1, j0
            DO i = 1, i0
              IF (pres(i,j,k,ismpl)<pmin) THEN
                icountp = icountp + 1
                dpmin = min(dpmin,pres(i,j,k,ismpl))
                pres(i,j,k,ismpl) = pmin
              END IF
              IF (pres(i,j,k,ismpl)>pmax) THEN
                icountp = icountp + 1
                dpmax = max(dpmax,pres(i,j,k,ismpl))
                pres(i,j,k,ismpl) = pmax
              END IF
              IF (temp(i,j,k,ismpl)<tmin) THEN
                icountt = icountt + 1
                dtmin = min(dtmin,temp(i,j,k,ismpl))
                temp(i,j,k,ismpl) = tmin
              END IF
              IF (temp(i,j,k,ismpl)>tmax) THEN
                icountt = icountt + 1
                dtmax = max(dtmax,temp(i,j,k,ismpl))
                temp(i,j,k,ismpl) = tmax
              END IF
              DO l = 1, ntrac
                IF (conc(i,j,k,l,ismpl)<cmin) THEN
                  icountc = icountc + 1
                  dcmin = min(dcmin,conc(i,j,k,l,ismpl))
                  conc(i,j,k,l,ismpl) = cmin
                END IF
!aw-later                  if (conc(i,j,k,l,ismpl).gt.cmax) then
!aw-later                     icountc = icountc +1
!aw-later                     dcmax = max(dcmax, conc(i,j,k,l,ismpl))
!aw-later                     conc(i,j,k,l,ismpl) = cmax
!aw-later                  endif
!                 use reasonable physical value [csmin]
                IF (conc(i,j,k,l,ismpl)<csmin) THEN
!                    to avoid numerically instabilities
                  conc(i,j,k,l,ismpl) = cmin
                END IF
              END DO
            END DO
          END DO
        END DO

        IF (icountp/=0) WRITE(*,'(3A,1I8,1A,1e16.7,1A,1e16.7,1A)') &
          'warning: pres not in domain of validity of module <', &
          trim(def_props), '> at ', icountp, ' points (min', dpmin, &
          ', max', dpmax, ')!'
        IF (icountt/=0) WRITE(*,'(3A,1I8,1A,1e16.7,1A,1e16.7,1A)') &
          'warning: temp not in domain of validity of module <', &
          trim(def_props), '> at ', icountt, ' points (min', dtmin, &
          ', max', dtmax, ')!'
        IF (icountc/=0) WRITE(*,'(3A,1I8,1A,1e16.7,1A,1e16.7,1A)') &
          'warning: conc not in domain of validity of module <', &
          trim(def_props), '> at ', icountc, ' points (min', dcmin, &
          ', max', dcmax, ')!'

        RETURN
      END
