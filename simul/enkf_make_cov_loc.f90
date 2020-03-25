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

!> @brief Sets the correlation matrix for covariance localisation
!> @details
!> The matrix `rhos` is set.
subroutine enkf_make_cov_loc()

  use mod_enkf, only:&
       cov_loc_lenx,&
       cov_loc_leny,&
       cov_loc_lenz,&
       rhos,&
       lstate,&
       irobs,&
       nrobs_loc,&
       num_enkf_vars,&
       cov_loc_assim_id,&
       act_s,&
       rank_mem,&
       i_obs,&
       j_obs,&
       k_obs

  implicit none


  integer :: imem, iloc, ivar, ivar_state_pos, istart, iend
!
  integer :: idu, iobs_mem, rmem
  integer :: io, jo, ko
  integer :: is, js, ks
  double precision :: d, dx, dy, dz
  double precision :: corr

  rhos(:,:,:) = 0.0d0
  
  do ivar = 1, num_enkf_vars

     if(act_s(ivar) == 1) then
        rmem = rank_mem(ivar)
        
        do iloc = 1, nrobs_loc(irobs)

           io = i_obs(irobs,iloc)
           jo = j_obs(irobs,iloc)
           ko = k_obs(irobs,iloc)
           
           istart = rmem*lstate + 1
           iend = rmem*lstate + lstate
           do imem = istart, iend
              
              ! Get the state vector position
              call enkf_index_mem_to_loc(imem,is,js,ks,idu) ! idu = ivar
              if(idu /= ivar) write(unit = *, fmt = *) '[E1] Error in enkf_make_cov_loc.f90'

              ! Euclidean distance
              call enkf_euclidean_dist(is,js,ks,io,jo,ko,d,dx,dy,dz)

              ! Correlation function
              call enkf_correlation_fcn(cov_loc_lenx,cov_loc_leny,cov_loc_lenz,dx,dy,dz,corr)
              rhos(imem,iloc,1) = corr
              cov_loc_assim_id = 1

           end do
        end do

     end if

  end do

end subroutine enkf_make_cov_loc
