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

!> @brief Save simulated variable values for later comparison with input data
!> @param[in] ismpl local sample index
!> @details
!> Set the values of sdata from computed (example: head) and old
!> variable arrays (example: headold). \n
!> sdata will be compared to the read in values in ddata. \n\n
!>
!> Two linear interpolations of the values are implemented:\n
!> 1. A linear interpolation of the position of the data (px,py,pz)
!>    inside the grid
!> 2. according to how the time specified for the data
!>    (ddata(l,cdd_time)) is located between the previous simulation
!>    time (simtime(ismpl)) and the current simulation time
!>    (simtime(ismpl)+deltt). \n\n
!>
!> collect and save the computed values for a comparison with
!> 'ddata(:,cid_pv)'\n
!> -> usage in 'write_data.f' and 'forward/j_*-array(inversion)'\n
      SUBROUTINE save_data(ismpl)

        use arrays
        use mod_genrl
        use mod_data
        use mod_time
        use mod_linfos

        IMPLICIT NONE

        integer :: ismpl
        integer :: i, j, k, l

        DOUBLE PRECISION deltt, deltat, numdiff, dalfa, dbeta, interpolatelin, bhpr
        DOUBLE PRECISION px, py, pz, vals, vals_old

        INTEGER i_type, i_si

        EXTERNAL deltat, interpolatelin, bhpr

        ! get current time step, set 1.0d0 for steady state
        deltt = deltat(simtime(ismpl),ismpl)
        IF ( .NOT. transient) deltt = 1.D0

!       allowed numerical difference
        numdiff = 1.0D2*const_dble(1)*simtime(ismpl)
!
        DO l = 1, ndata

!         correct time interval, or steady-state
          IF (ddata(l,cdd_time)>simtime(ismpl)+numdiff .AND. &
              ddata(l,cdd_time)<=simtime(ismpl)+deltt+numdiff .OR. &
              .NOT. transient) THEN

            i = idata(l,cid_i)
            j = idata(l,cid_j)
            k = idata(l,cid_k)
            px = ddata(l,cdd_i)
            py = ddata(l,cdd_j)
            pz = ddata(l,cdd_k)

            i_type = idata(l,cid_pv)
            i_si = idata(l,cid_si)

!           interpolation: m=(a*n+b*o)/(a+b)
            dalfa = ddata(l,cdd_time) - simtime(ismpl)
            dbeta = deltt - dalfa
            IF ( .NOT. transient) THEN
              dalfa = 1.D0
              dbeta = 0.D0
            END IF

            vals = 0.0D0
            vals_old = 0.0D0

!           choose physical value to save
            IF (i_type==pv_head) THEN

              vals = interpolatelin(i0,j0,k0, i,j,k, &
                head(1,1,1,ismpl), &
                px,py,pz, delx,dely,delz, delxa,delya,delza)

              vals_old = interpolatelin(i0,j0,k0, i,j,k, &
                headold(1,cgen_time,ismpl), &
                px,py,pz, delx,dely,delz, delxa,delya,delza)

            ELSE IF (i_type==pv_pres) THEN

              vals = interpolatelin(i0,j0,k0, i,j,k, &
                pres(1,1,1,ismpl), &
                px,py,pz, delx,dely,delz, delxa,delya,delza)

              vals_old = interpolatelin(i0,j0,k0, i,j,k, &
                presold(1,cgen_time,ismpl), &
                px,py,pz, delx,dely,delz, delxa,delya,delza)

            ELSE IF (i_type==pv_temp) THEN

              vals = interpolatelin(i0,j0,k0, i,j,k, &
                temp(1,1,1,ismpl), &
                px,py,pz, delx,dely,delz, delxa,delya,delza)

              vals_old = interpolatelin(i0,j0,k0, i,j,k, &
                tempold(1,cgen_time,ismpl), &
                px,py,pz, delx,dely,delz, delxa,delya,delza)

            ELSE IF (i_type==pv_conc) THEN

              vals = interpolatelin(i0,j0,k0, i,j,k, &
                conc(1,1,1,i_si,ismpl), &
                px,py,pz, delx,dely,delz, delxa,delya,delza)

              vals_old = interpolatelin(i0,j0,k0, i,j,k, &
                concold(1,i_si,cgen_time,ismpl), &
                px,py,pz, delx,dely,delz, delxa,delya,delza)

            ELSE IF (i_type==pv_bhpr) THEN

              vals = bhpr(i,j,k,ismpl)
              vals_old = vals

            END IF

            ! Write time-interpolated value to sdata
            sdata(l,ismpl) = (dalfa*vals+dbeta*vals_old)/deltt

          END IF

        END DO

        RETURN
      END
