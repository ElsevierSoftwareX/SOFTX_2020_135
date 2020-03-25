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

!> @brief Deallocate arrays used in EnKF simulation
!> @details
!> Deallocate all arrays used in the entire EnKF simulation.
subroutine enkf_dealloc()

  use mod_enkf, only: &
       stdvar,&
       resvar,&
       var_true,&
       obst,&
       nrobs_loc,&
       act_o,&
       i_obs,&
       j_obs,&
       k_obs,&
       var_obs,&
       sysindx,&
       mat_cov_ref,&
       mat_single_out,&
       pb,&
       pb_r,&
       pb_reps,&
       stoch_bc_nbc,&
       stoch_bc_ibc,&
       tempdiff_inds,&
       tcon_inds,&
       pp_inds,&
       gss

  implicit none

  if(allocated(stdvar)) deallocate(stdvar)

  if(allocated(resvar)) deallocate(resvar)

  if(allocated(var_true)) deallocate(var_true)

  if(allocated(obst)) deallocate(obst)
  if(allocated(nrobs_loc)) deallocate(nrobs_loc)

  if(allocated(act_o)) deallocate(act_o)
  
  if(allocated(i_obs)) deallocate(i_obs)
  if(allocated(j_obs)) deallocate(j_obs)
  if(allocated(k_obs)) deallocate(k_obs)

  if(allocated(var_obs)) deallocate(var_obs)
  
  if(allocated(sysindx)) deallocate(sysindx)

  if(allocated(mat_cov_ref)) deallocate(mat_cov_ref)
  if(allocated(mat_single_out)) deallocate(mat_single_out)

  if(allocated(pb)) deallocate(pb)
  if(allocated(pb_r)) deallocate(pb_r)
  if(allocated(pb_reps)) deallocate(pb_reps)

  if (allocated(stoch_bc_nbc)) deallocate(stoch_bc_nbc)
  if (allocated(stoch_bc_ibc)) deallocate(stoch_bc_ibc)

  if(allocated(tempdiff_inds)) deallocate(tempdiff_inds)

  if(allocated(tcon_inds)) deallocate(tcon_inds)
  if(allocated(pp_inds)) deallocate(pp_inds)

  if (allocated(gss)) deallocate(gss)

end subroutine enkf_dealloc
