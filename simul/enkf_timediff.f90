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

!> @brief Set observation arrays for temperature difference observations
!> @param[in] nrobs Number of observation locations for current update
!> @param[inout] obsvec Observation values
!> @param[inout] obsvar_assim Observation variances
!> @param[inout] obspos Observation locations
!> @param[inout] s S-Matrix for observation locations
!> @details
!> Set observation arrays for temperature difference
!> observations. Dummy time-differences arrays are introduced, set
!> according to the original temperature observations and finally
!> output as observation data for the update.
subroutine enkf_timediff(nrobs,obsvec,obsvar_assim,obspos,s)

  use mod_enkf, only: &
       nrens,&
       tempdiff_inds,&
       num_tempdiff

  implicit none

  integer, intent(in) :: nrobs
  double precision, intent(inout), dimension(nrobs) :: obsvec
  double precision, intent(inout), dimension(nrobs) :: obsvar_assim
  double precision, intent(inout), dimension(nrobs) :: obspos
  double precision, intent(inout), dimension(nrobs,nrens) :: s

  double precision, allocatable :: td_obsvec(:)
  double precision, allocatable :: td_obsvar_assim(:)
  double precision, allocatable :: td_s(:,:)
  integer, allocatable :: td_obspos(:)

  integer :: itd, irens
  
  ALLOCATE(td_obsvec(nrobs))
  ALLOCATE(td_obsvar_assim(nrobs))
  ALLOCATE(td_obspos(nrobs))
  ALLOCATE(td_s(nrobs,nrens))
  td_obsvec(:) = 0.0d0
  td_obsvar_assim(:) = 0.0d0
  td_obspos(:) = 0
  td_s(:,:)=0.0d0

  ! Set time difference data
  do itd = 1, num_tempdiff
     td_obsvec(itd) = obsvec(tempdiff_inds(itd,1)) &
          - obsvec(tempdiff_inds(itd,2))
     td_obsvar_assim(itd) = obsvar_assim(tempdiff_inds(itd,1))
     td_obspos(itd) = obspos(tempdiff_inds(itd,1))
     do irens = 1, nrens
        td_s(itd,irens) = s(tempdiff_inds(itd,1),irens) &
             - s(tempdiff_inds(itd,2),irens)
     end do
  end do

  ! Check array sizes
  if(size(obsvec) .ne. size(td_obsvec)) then
     write(unit = *, fmt = *) '[E1] Error in enkf_timediff'
     stop 1
  end if
  if(size(obsvar_assim) .ne. size(td_obsvar_assim)) then
     write(unit = *, fmt = *) '[E2] Error in enkf_timediff'
     stop 1
  end if
  if(size(obspos) .ne. size(td_obspos)) then
     write(unit = *, fmt = *) '[E3] Error in enkf_timediff'
     stop 1
  end if
  if(size(s) .ne. size(td_s)) then
     write(unit = *, fmt = *) '[E4] Error in enkf_timediff'
     stop 1
  end if

  ! Set observation arrays
  obsvec = td_obsvec
  obsvar_assim = td_obsvar_assim
  obspos = td_obspos
  s = td_s

  DEALLOCATE(td_obsvec)
  DEALLOCATE(td_obsvar_assim)
  DEALLOCATE(td_obspos)
  DEALLOCATE(td_s)
  
end subroutine enkf_timediff
