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

!> @brief Wrapper for combined update
!> @details
!> Preparations for combined update and then calling subroutine
!> `enkf()` from Evensen's EnKF code.
subroutine enkf_assim_comb()

  use m_enkf_enkf
  use m_random

  use mod_enkf, only: &
       mem,&
       act_o,&
       normobs,&
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
       obsvar,&
       ns_switch,&
       obs_ns,&
       obsvar_ns,&
       num_enkf_vars
  
  implicit none

  double precision, allocatable :: obsvec(:)
  double precision, allocatable :: obsvar_assim(:)
  double precision, allocatable :: obswork(:)
  double precision, allocatable :: s(:,:)
  integer, allocatable :: obspos(:)

  double precision, external :: varpar

  integer, external :: index_obs_to_mem

  double precision :: kz, lz, por
  external kz, lz, por

  integer :: m, lloc, nrobs, iloc, irens
  integer :: ivar, nrobsex
  
  integer :: n_act_o

  call enkf_output_assimstp(0)
  
  ! Number of active observation variables
  n_act_o = sum(act_o)
  nrobs = nrobs_loc(irobs)*n_act_o

  ! Exclude measurements (NAN=1.D60)
  nrobsex = 0
  DO iloc = 1, nrobs_loc(irobs)
     do ivar = 1, num_enkf_vars
        if((act_o(irobs,ivar)==1) .and.&
             (var_obs(irobs,iloc,ivar)>1.0d59)) then
           nrobsex = nrobsex + 1
        end if
     end do
  END DO
  nrobs=nrobs-nrobsex

  ALLOCATE(obsvec(nrobs))
  ALLOCATE(obsvar_assim(nrobs))
  ALLOCATE(obspos(nrobs))
  ALLOCATE(obswork(nrobs))
  ALLOCATE(s(nrobs,nrens))
  obsvec(:) = 0.D0
  obsvar_assim(:) = 0.D0
  obswork(:) = 0.D0
  obspos(:) = 0
  s(:,:)=0.d0
  !-----------------------------------------
  lloc = 0
  DO iloc = 1, nrobs_loc(irobs)
     do ivar = 1, num_enkf_vars
        IF ((act_o(irobs,ivar)==1) .and. &
             (var_obs(irobs,iloc,ivar)<1.D59) )  THEN
           lloc = lloc + 1
           obspos(lloc) = index_obs_to_mem(iloc,ivar)

           if(ns_switch) then
              obsvec(lloc) = var_obs(irobs,iloc,ivar)
              obsvar_assim(lloc) = obsvar_ns(iloc,ivar)
              DO irens = 1, nrens
                 s(lloc,irens) = obs_ns(irens,iloc,ivar)
              END DO
           else
              obsvec(lloc) = var_obs(irobs,iloc,ivar)
              obsvar_assim(lloc) = obsvar(ivar)
              DO irens = 1, nrens
                 s(lloc,irens) = varpar(i_obs(irobs,iloc),j_obs(irobs,iloc),&
                      k_obs(irobs,iloc),irens,ivar)
              END DO
           end if

        END IF
     end do
  END DO

  ! add OBSERVATION NOISE (here better correlated noise may be added)
  CALL random(obswork,nrobs)
  DO iloc = 1, nrobs
     obsvec(iloc) = obsvec(iloc) + dsqrt(obsvar_assim(iloc))*obswork(iloc)
  END DO

  call enkf_log(5,0)

  IF (normobs) then
     !  normalize the  observation vector for different kind of observations
     DO irens =1, nrens
        DO m =1, nrobs
           S(m,irens)      =  S(m,irens) /dsqrt(obsvar_assim(m))
        END DO
     END DO
     DO m =1, nrobs
        obsvec(m)      =  obsvec(m) /dsqrt(obsvar_assim(m))
        obsvar_assim(m)=1.d0
     END DO
  END IF

  !------------------------------------------------------------------
  !------------------------------------------------------------------
  CALL enkf(mem,s,nstate,nrens,obsvec,obsvar_assim,obspos,nrobs, &
       irobs,fixsamp,mode_analysis,truncation,covmodel,clh,clv, &
       rexact,assidamp,lstate)
  !------------------------------------------------------------------
  !------------------------------------------------------------------
  DEALLOCATE(obsvec)
  DEALLOCATE(obsvar_assim)
  DEALLOCATE(obspos)
  DEALLOCATE(obswork)
  DEALLOCATE(s)

  call enkf_output_assimstp(1)

end subroutine enkf_assim_comb
