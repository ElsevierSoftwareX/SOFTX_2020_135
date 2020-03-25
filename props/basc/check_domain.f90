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

!>    @brief domain of validity for module basc
!>    @param[in] ismpl local sample index
!>    @details
!>    Checking whether pres/temp/(conc) are in domain of props
!>    validity. Version for property module basc. \n
!>    \n
!>    For concentration, an error is thrown and the execution is
!>    stopped if the concentration is outside the physical values.
      subroutine check_domain(ismpl)
        use arrays, only: pres, temp, conc
        use mod_genrl, only: i0, j0, k0
        use mod_genrlc, only: def_props
        use mod_conc, only: ntrac
        use mod_linfos, only: linfos

        implicit none

        ! Sample index
        integer :: ismpl

        ! Iteration counters
        integer :: i, j, k, l

        ! counters for the values outside domain of validity
        ! pres
        integer :: icountp
        ! temp
        integer :: icountt
        ! conc
        integer :: icountc

        ! min/max boundaries of the domain of validity
        ! pres
        double precision, parameter ::  pmin = 0.01d6
        double precision, parameter ::  pmax = 110.0d6
        ! temp
        double precision, parameter ::  tmin = 0.0d0
        double precision, parameter ::  tmax = 350.0d0
        ! conc
        double precision, parameter ::  cmin = 0.0d0
        double precision, parameter ::  cmax = 1.0d5
        ! numerical boundary
        double precision, parameter ::  csmin = 1.0d-30

        ! records the overall min/max of values if they are outside
        ! domain of validity
        double precision :: dpmax, dtmax, dcmax, dhmax
        double precision :: dpmin, dtmin, dcmin, dhmin

        intrinsic trim


        ! Set counters to zero
        icountp = 0
        icountt = 0
        icountc = 0

        ! Set overall min/max to boundaries of the domain of validity
        dpmax = pmax
        dpmin = pmin
        dtmax = tmax
        dtmin = tmin
        dcmax = cmax
        dcmin = cmin

        ! Check pres
        do k = 1, k0
          do j = 1, j0
            do i = 1, i0
              if (pres(i,j,k,ismpl)<pmin) then
                ! Set min counter
                icountp = icountp + 1
                ! Set new overall minimum
                dpmin = min(dpmin,pres(i,j,k,ismpl))
                ! Change pres value to minimum of the domain of validity
                pres(i,j,k,ismpl) = pmin
              end if
              if (pres(i,j,k,ismpl)>pmax) then
                ! Set max counter
                icountp = icountp + 1
                ! Set new overall maximum
                dpmax = max(dpmax,pres(i,j,k,ismpl))
                ! Change pres value to maximum of the domain of validity
                pres(i,j,k,ismpl) = pmax
              end if
            end do
          end do
        end do

        ! Check temp
        do k = 1, k0
          do j = 1, j0
            do i = 1, i0
              if (temp(i,j,k,ismpl)<tmin) then
                icountt = icountt + 1
                dtmin = min(dtmin,temp(i,j,k,ismpl))
                temp(i,j,k,ismpl) = tmin
              end if
              if (temp(i,j,k,ismpl)>tmax) then
                icountt = icountt + 1
                dtmax = max(dtmax,temp(i,j,k,ismpl))
                temp(i,j,k,ismpl) = tmax
              end if
            end do
          end do
        end do

        ! Check conc
        do k = 1, k0
          do j = 1, j0
            do i = 1, i0
              do l = 1, ntrac
                if (conc(i,j,k,l,ismpl).gt.cmax) then
                  icountc = icountc +1
                  dcmax = max(dcmax, conc(i,j,k,l,ismpl))
                  conc(i,j,k,l,ismpl) = cmax
                end if
                if (conc(i,j,k,l,ismpl)<cmin .and. &
                  conc(i,j,k,l,ismpl)<-csmin) then
                  icountc = icountc + 1
                  dcmin = min(dcmin,conc(i,j,k,l,ismpl))
                  conc(i,j,k,l,ismpl) = cmin
                end if
                if (conc(i,j,k,l,ismpl)<csmin) then
                  ! very small conc values set to zero to avoid
                  ! numerically instabilities
                  conc(i,j,k,l,ismpl) = cmin
                end if
              end do
            end do
          end do
        end do

!       disable the warning output for linfos(3)==-1
        if (linfos(3)>=0) then
          if (icountp/=0) write(*,'(3A,1I8,1A,1e16.7,1A,1e16.7,1A)') &
            'warning: pres not in domain of validity of module <', &
            trim(def_props), '> at ', icountp, ' points (min', dpmin, &
            ', max', dpmax, ')!'
          if (icountt/=0) write(*,'(3A,1I8,1A,1e16.7,1A,1e16.7,1A)') &
            'warning: temp not in domain of validity of module <', &
            trim(def_props), '> at ', icountt, ' points (min', dtmin, &
            ', max', dtmax, ')!'
          if (icountc/=0) write(*,'(3A,1I8,1A,1e16.7,1A,1e16.7,1A)') &
            'warning: conc not in domain of validity of module <', &
            trim(def_props), '> at ', icountc, ' points (min', dcmin, &
            ', max', dcmax, ')!'

          ! error outputs for hard physical concentration boundaries
          if (dcmax > cmax) then
            write(unit = *, fmt = *) "[E1] Error in check_domain.f90:", &
                "  maximum concentration dcmax= ", dcmax, &
                " larger than allowed maximum value cmax=", cmax
            stop
          end if
          if (dcmin > cmin) then
            write(unit = *, fmt = *) "[E2] Error in check_domain.f90:", &
                "  minimum concentration dcmin= ", dcmin, &
                " smaller than allowed minimum value cmin=", cmin
            stop
          end if
        end if

        return

      end subroutine check_domain
