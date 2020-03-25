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

!> @brief Call Evensen enkf() subroutine after final preparations
!> @param[in] assim_id Observation variable index
!> @details
!> Call Evensen `enkf()` subroutine after final preparations
subroutine enkf_assim_single(assim_id)

  use m_enkf_enkf
  use m_random

  use mod_enkf, only: &
       mem,&
       irobs,&
       nrobs_loc,&
       lstate,&
       nstate,&
       nrens,&
       i_obs,&
       j_obs,&
       k_obs,&
       var_obs,&
       fixsamp,&
       mode_analysis,&
       truncation,&
       covmodel,&
       clh,&
       clv,&
       rexact,&
       assidamp,&
       tracdamp,&
       obsvar,&
       ns_switch,&
       obsvar_ns,&
       hybrid_switch,&
       obs_ns,&
       tempdiff_switch

  implicit none

  integer, intent(in) :: assim_id

  double precision, allocatable :: obsvec(:)
  double precision, allocatable :: obsvar_assim(:)
  double precision, allocatable :: obswork(:)
  double precision, allocatable :: s(:,:)
  integer, allocatable :: obspos(:)
  
  double precision, external :: varpar

  integer, external :: index_obs_to_mem
  
  double precision :: kz, lz, por
  external kz, lz, por

  integer :: nrobs, lloc, nrobsex, iloc, irens
  double precision :: assisave

!-----------------------------------------
! Assimilate
!-----------------------------------------
     ! Set number of observation locations
     nrobs = nrobs_loc(irobs)

     ! Exclude measurements (NAN=1.D60)
     nrobsex = 0
     DO iloc = 1, nrobs_loc(irobs)
        if(var_obs(irobs,iloc,assim_id)>1.0d59) then
           nrobsex = nrobsex + 1
        end if
     END DO
     nrobs = nrobs - nrobsex

     ! Hybrid EnKF only
     if (hybrid_switch) then
        call enkf_hybrid_subcovs(assim_id,nrobs)
     end if
     
     IF (nrobs>0) THEN
        ALLOCATE(obsvec(nrobs))
        ALLOCATE(obsvar_assim(nrobs))
        ALLOCATE(obspos(nrobs))
        ALLOCATE(obswork(nrobs))
        ALLOCATE(s(nrobs,nrens))
        obsvec(:) = 0.0d0
        obsvar_assim(:) = 0.0d0
        obswork(:) = 0.0d0
        obspos(:) = 0
        s(:,:)=0.0d0
        lloc=0
        DO iloc = 1, nrobs_loc(irobs)
           if( var_obs(irobs,iloc,assim_id)<1.0d59) then
              lloc=lloc + 1
              obspos(lloc) = index_obs_to_mem(iloc,assim_id)

              ! Set observation values and variance
              if(ns_switch) then
                 obsvec(lloc) = var_obs(irobs,iloc,assim_id)
                 obsvar_assim(lloc) = obsvar_ns(iloc,assim_id)
                 do irens = 1, nrens
                    s(lloc,irens) = obs_ns(irens,iloc,assim_id)
                 end do
              else
                 obsvec(lloc) = var_obs(irobs,iloc,assim_id)
                 obsvar_assim(lloc) = obsvar(assim_id)
                 ! Determine the S-matrix for m_enkf_enkf.f90 from
                 ! original ensemble
                 do irens = 1, nrens
                    s(lloc,irens) = varpar(i_obs(irobs,iloc),j_obs(irobs,iloc),&
                         k_obs(irobs,iloc),irens,assim_id)
                 end do
              end if

           end if
        END DO

        ! Temperature differences as observations
        if(tempdiff_switch) then
           call enkf_timediff(nrobs,obsvec,obsvar_assim,obspos,s)
        end if

        ! Add observation noise
        CALL random(obswork,nrobs)
        DO iloc = 1, nrobs
           obsvec(iloc) = obsvec(iloc) + dsqrt(obsvar_assim(iloc))*obswork(iloc)
        END DO

        ! Special damping coefficient for concentration
        assisave = assidamp
        if(assim_id == 3) then
           assidamp = tracdamp
        end if

        call enkf_log(4,assim_id)

!------------------------------------------------------------------
!------------------------------------------------------------------
        CALL enkf(mem,s,nstate,nrens,obsvec,obsvar_assim,obspos,nrobs, &
             irobs,fixsamp,mode_analysis,truncation,covmodel,clh,clv, &
             rexact,assidamp,lstate)
        if(assim_id == 3) then
           assidamp = assisave
        end if
!------------------------------------------------------------------
!------------------------------------------------------------------

        if (hybrid_switch) then
           call enkf_hybrid_subcovs_dealloc()
        end if

        DEALLOCATE(obsvec)
        DEALLOCATE(obsvar_assim)
        DEALLOCATE(obspos)
        DEALLOCATE(obswork)
        DEALLOCATE(s)
     END IF

end subroutine enkf_assim_single
