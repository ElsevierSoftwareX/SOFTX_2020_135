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

!> @brief Generate Normal Score arrays at observation locations
!> @details
!> Set `obs_ns`, `obs_ns_original`, `obs_ns_gauss`,
!> `obs_ns_orig_vars`. See `mod_enkf.f90` for more information.
subroutine enkf_normal_score_obsvars()

  use mod_enkf,only:&
       act_o,&
       irobs,&
       num_enkf_vars,&
       obs_ns,&
       obs_ns_original,&
       obs_ns_gauss,&
       obs_ns_orig_vars,&
       nrobs_loc,&
       nrens,&
       i_obs,&
       j_obs,&
       k_obs  

  implicit none
  
  integer :: ivar, i, iloc, irens, ii, jj, kk
  integer :: nrobs

  integer, dimension(nrens) :: obs_ind
  double precision :: avetmp, vartmp

  double precision, external :: varpar
  
  double precision :: kz, lz, por
  external kz, lz, por
  
  nrobs = nrobs_loc(irobs)

  allocate(obs_ns(nrens,nrobs,num_enkf_vars))
  allocate(obs_ns_original(nrens,nrobs,num_enkf_vars))
  allocate(obs_ns_gauss(nrens,nrobs,num_enkf_vars))
  allocate(obs_ns_orig_vars(nrobs,num_enkf_vars))
  obs_ns(:,:,:) = 0.0d0
  obs_ns_original(:,:,:) = 0.0d0
  obs_ns_gauss(:,:,:) = 0.0d0
  obs_ns_orig_vars(:,:) = 0.0d0

  obs_ind(:) = 0.0d0
  
  do ivar = 1, num_enkf_vars
     if (act_o(irobs,ivar) == 1) then

        do iloc = 1, nrobs

           ! Set obs_ns_original
           ii = i_obs(irobs,iloc)
           jj = j_obs(irobs,iloc)
           kk = k_obs(irobs,iloc)
           do irens = 1, nrens
              obs_ns_original(irens,iloc,ivar) = varpar(ii,jj,kk,irens,ivar)
           end do
           
           ! Set obs_ns_orig_vars
           call enkf_vector_average(obs_ns_original(:,iloc,ivar),nrens,avetmp)
           call enkf_vector_variance(obs_ns_original(:,iloc,ivar),nrens,avetmp,vartmp)
           obs_ns_orig_vars(iloc,ivar) = vartmp

           ! Sort obs_ns_original
           call enkf_sort_vector_ind(obs_ns_original(:,iloc,ivar),nrens,obs_ind)
           
           ! Set obs_ns, obs_ns_gauss
           call enkf_generate_gaussians('cumu',obs_ns(:,iloc,ivar),nrens)
           obs_ns_gauss(:,iloc,ivar) = obs_ns(:,iloc,ivar)

           ! Sort obs_ns
           call enkf_sort_vector_back(obs_ns(:,iloc,ivar),nrens,obs_ind)

        end do

     end if
  end do

end subroutine enkf_normal_score_obsvars
