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

!> @brief EnKF Input checks
!> @details
!> Variables read in from the EnKF input file in subroutine
!> `enkf_input()` are checked (mostly if they are inside acceptable
!> value ranges)
subroutine enkf_input_check()

  use arrays, only: &
       vdefaultswitch, first_flow, last_flow, ibc_data, cbc_bt
  
  use mod_enkf, only: &
       assidamp,&
       assitype,&
       clh,&
       clv,&
       covmodel,&
       ex_unit,&
       iassout,&
       iassout_single_start,&
       iassout_single,&
       iassout_cov_ref_start,&
       iassout_cov_ref,&
       mat_single_out,&
       num_single_out,&
       mode_analysis,&
       nrobs_int,&
       num_enkf_vars,&
       num_ref_dist,&
       ref_dist_var,&
       ref_dist_switch,&
       tracdamp,&
       truncation,&
       sysvar,&
       switch_cov_ref,&
       mat_cov_ref,&
       num_cov_ref,&
       vtk_out_covs,&
       cov_loc_switch,&
       cov_loc_lenx,&
       cov_loc_leny,&
       cov_loc_lenz,&
       hybrid_switch,&
       hybrid_alpha,&
       obsvar,&
       act_s,&
       enkf_variable_names,&
       dual_enkf_switch,&
       iterative_switch,&
       vtk_out_realisations,&
       num_out_realisations,&
       stoch_bc_switch,&
       pres_vel_switch,&
       num_pres_vel,&
       tempdiff_switch,&
       pp_switch,&
       ns_switch,&
       num_pp,&
       pp_inds,&
       pp_get_out,&
       pp_ivar

  use mod_genrl, only: &
       i0,&
       j0,&
       k0,&
       nsmpl,&
       head_active

  implicit none

  integer :: i, ivar, ib
  integer :: i_pp, j_pp, k_pp, i_pp_bef, j_pp_bef, k_pp_bef

  integer :: i_so, j_so, k_so, ivar_so
  character (len=3) :: single_out_lin_log

  integer :: i_cov_ref
  integer :: j_cov_ref
  integer :: k_cov_ref
  integer :: ivar_cov_ref

  write(unit = *, fmt = *) 
  write(unit = *, fmt = *) 
  write(unit = *, fmt = *) "---------- CHECK ENKF INPUT -------"
  write(unit = *, fmt = *) " [IE] means  Input Error"
  write(unit = *, fmt = *) " [OK] means OK"
  write(unit = *, fmt = *) 
  write(unit = *, fmt = *)

  !Number of observation times
  if(nrobs_int < 1) then
     write(unit = *, fmt = *) "[IE] nrobs_int less then 1"
     stop 1
  end if

  ! Truncation
  if(truncation > 1.0d0 .or. truncation < 0.0d0) then
     write(unit = *, fmt = *) "[IE] truncation out of range."
     stop 1
  end if

  ! Analysis mode switch
  if(mode_analysis /= 11 &
       .and. mode_analysis /= 12 &
       .and. mode_analysis /= 13 &
       .and. mode_analysis /= 21 &
       .and. mode_analysis /= 22 &
       .and. mode_analysis /= 23 ) then
     write(unit = *, fmt = *) "[IE] Wrong mode_analysis (11,12,13,21,22,23)"
     stop 1
  end if

  ! Covariance model
  if(trim(covmodel) /= 'gaussian' &
       .and. trim(covmodel) /= 'diagonal') then
     write(unit = *, fmt = *) "[IE] Wrong covmodel specification (gaussian/diagonal)"
     stop 1
  end if

  ! Correlation lengths
  if(clh < 0.0d0 .or. clv < 0.0d0) then
     write(unit = *, fmt = *) "[IE] Correlation length negative."
     stop 1
  end if

  ! Variances should be positive
  if(obsvar(1) < 0.0d0 &
       .or. obsvar(2) < 0.0d0 &
       .or. obsvar(3) < 0.0d0 &
       .or. obsvar(4) < 0.0d0 &
       .or. obsvar(5) < 0.0d0 &
       .or. obsvar(6) < 0.0d0 &
       .or. sysvar(1) < 0.0d0 &
       .or. sysvar(2) < 0.0d0 &
       .or. sysvar(3) < 0.0d0 &
       .or. sysvar(4) < 0.0d0 &
       .or. sysvar(5) < 0.0d0 &
       .or. sysvar(6) < 0.0d0) then
     write(unit = *, fmt = *) "[IE] Variances negative"
     stop 1
  end if

  ! Activity of parameters, 1:h, 2:t, 3:c, 4:kz, 5:lz, 6:por
  do ivar = 1, num_enkf_vars
     if(act_s(ivar) /= 0 .and. act_s(ivar) /= 1) then
        write(unit = *, fmt = *) "[IE] act_s out of range at variable ", ivar
     end if
  end do
  
  !Exclude unit
  if(ex_unit < 0) then
     write(unit = *, fmt = *) "[IE] ex_unit out of range"
     stop 1
  end if

  ! Damping
  if( assidamp < 0.0d0 .or. assidamp > 1.01d0) then
     write(unit = *, fmt = *) "[IE] assidamp out of range"
     stop 1
  end if

  if( tracdamp < 0.0d0 .or. tracdamp > 1.01d0) then
     write(unit = *, fmt = *) "[IE] tracdamp out of range"
     stop 1
  end if

  ! Measurement output
  if(iassout < 0) then
     write(unit = *, fmt = *) "[IE] iassout out of range"
     stop 1
  elseif( iassout == 0 ) then
     ! Zero as input should lead to no output at all
     iassout = nrobs_int + 1
  end if

  ! Single cell output times input
  if( iassout_single_start < 0 ) then
     ! Negative input leads to no output
     iassout_single_start = nrobs_int + 1
  elseif( iassout_single_start == 0 ) then
     ! "0" leads to 'ini' output, otherwise acts like "1"
     iassout_single_start = 1
  end if

  if( iassout_single < 0 ) then
     write(unit = *, fmt = *) "[IE] iassout_single out of range"
     stop 1
  elseif( iassout_single == 0) then
     ! Zero as input leads to maximally one output
     iassout_single = nrobs_int + 1
  end if
  
  ! Single cell output input
  do i = 1, num_single_out
     
     i_so = mat_single_out(i,1)
     j_so = mat_single_out(i,2)
     k_so = mat_single_out(i,3)
     ivar_so = mat_single_out(i,4)
     select case (mat_single_out(i,5))
     case(0)
        single_out_lin_log = 'lin'
     case(1)
        single_out_lin_log = 'log'
     case default
        write(unit = *, fmt = *) '[IE] single_out_lin_log in mat_single_out wrong'
     end select

     if( i_so < 0 .or. i_so > i0 ) then
        write(unit = *, fmt = *) "[IE]: i_so out of range"
        stop 1
     else if( j_so < 0 .or. j_so > j0 ) then
        write(unit = *, fmt = *) "[IE]: j_so out of range"
        stop 1
     else if( k_so < 0 .or. k_so > k0 ) then
        write(unit = *, fmt = *) "[IE]: k_so out of range"
        stop 1
     end if

     select case (ivar_so)
     case (1,2,3,4,5,6)
        if(act_s(ivar_so) /= 1) then
           write(unit = *, fmt = *) "[IE] ivar_so NOT in state_vector"
           stop 1
        end if
     case (11,12,13,14,15,16)
        if(act_s(ivar_so-10) == 1) then
           write(unit = *, fmt = *) "[IE] ivar_so IN state_vector"
           stop 1
        end if
     case default
        write(unit = *, fmt = *) "[IE]: ivar_so out of range"
        stop 1
     end select

     select case (single_out_lin_log)
     case ('lin')
        write(unit = *, fmt = *) "[OK] single_out_lin_log = lin"
        write(unit = *, fmt = *) 
     case('log')
        if(ivar_so /= 4 .and. ivar_so /= 3) then
           write(unit = *, fmt = *)&
                "[IE]: single_out_lin_log = log only for permeability/conc"
           stop 1
        end if
     case default
        write(unit = *, fmt = *)&
             "[IE]: single_out_lin_log has wrong character value"
        stop 1
     end select
  end do

  ! Reference distribution input
  if(ref_dist_switch) then
     do i = 1, num_ref_dist 
        if(ref_dist_var(i) < 1 &
             .or. ref_dist_var(i) > num_enkf_vars) then
           write(unit = *, fmt = *) "[IE]: ref_dist_var out of range"
           stop 1
        end if
     end do
  else
     write(unit = *, fmt = *) "[OK]: ref_dist_switch false"
  end if

  ! Reference correlation output times input
  if( iassout_cov_ref_start < 0 ) then
     write(unit = *, fmt = *) "[IE] iassout_cov_ref_start out of range"
     stop 1
  elseif( iassout_cov_ref_start == 0 ) then
     ! Zero as input should lead to no output at all
     iassout_cov_ref_start = nrobs_int + 1
  end if

  if( iassout_cov_ref < 0 ) then
     write(unit = *, fmt = *) "[IE] iassout_cov_ref out of range"
     stop 1
  elseif( iassout_cov_ref == 0) then
     ! Zero as input leads to maximally one output
     iassout_cov_ref = nrobs_int + 1
  end if

  ! Switches for covariance output not compatible
  if( switch_cov_ref .neqv. vtk_out_covs ) then
     write(unit = *, fmt = *) "[IE] switch_cov_ref and vtk_out_covs incompatible"
  end if

  !Reference correlation location
  do i = 1, num_cov_ref

     i_cov_ref = mat_cov_ref(i,1)
     j_cov_ref = mat_cov_ref(i,2)
     k_cov_ref = mat_cov_ref(i,3)
     ivar_cov_ref = mat_cov_ref(i,4)

     if(i_cov_ref < 1 .or. i_cov_ref > i0) then
        write(unit = *, fmt = *) "[IE] i_cov_ref out of range"
        stop 1
     else if(j_cov_ref < 1 .or. j_cov_ref > j0) then
        write(unit = *, fmt = *) "[IE] j_cov_ref out of range"
        stop 1
     else if(k_cov_ref < 1 .or. k_cov_ref > k0) then
        write(unit = *, fmt = *) "[IE] k_cov_ref out of range"
        stop 1
     end if
     
     ! Reference correlation variable
     select case (ivar_cov_ref)
     case (1,2,3,4,5,6)
        if(act_s(ivar_cov_ref) /= 1) then
           write(unit = *, fmt = *) "[IE] ivar_cov_ref not in state_vector"
           stop 1
        end if
     case (11,12,13)
        if(act_s(ivar_cov_ref-10) == 1) then
           write(unit = *, fmt = *) "[IE] ivar_cov_ref=", ivar_cov_ref ,&
           ", but ", enkf_variable_names(ivar_cov_ref-10) ," in state_vector"
           stop 1
        end if
     case default
        write(unit = *, fmt = *) "[IE]: ivar_cov_ref out of range"
        stop 1
     end select
  end do

  ! Covariance localisation length in reasonable length scale
  if ( (cov_loc_lenx < 0.0d0 ) ) then
     write(unit = *, fmt = *) "[IE] cov_loc_lenx negative"
     stop 1
  else if ( (cov_loc_leny < 0.0d0 ) ) then
     write(unit = *, fmt = *) "[IE] cov_loc_leny negative"
     stop 1
  else if ( (cov_loc_lenz < 0.0d0 ) ) then
     write(unit = *, fmt = *) "[IE] cov_loc_lenz negative"
     stop 1
  end if

  ! Hybrid EnKF
  ! Not supposed to be started together with covariance localisation
  ! This would have to be implemented in analysis.f90
  if(cov_loc_switch .and. hybrid_switch) then
     write(unit = *, fmt = *) "[IE]: cov_loc_switch and hybrid_switch on together"
     stop 1
  end if

  if(hybrid_alpha < -0.0001d0 .or. hybrid_alpha > 1.0001d0) then
     write(unit = *, fmt = *) "[IE]: hybrid_alpha ", hybrid_alpha, " outside interval [0,1]."
  end if
  
  ! Iterative and Dual EnKF together not possible
  if(dual_enkf_switch .and. iterative_switch) then
     write(unit = *, fmt = *) "[IE] Dual EnKF and Iterative EnKF both on"
     stop 1
  end if

  ! Assim variable output for less than 20 realisations
  if(vtk_out_realisations) then
     if (num_out_realisations < 1 .or. num_out_realisations > 20) then
        write(unit = *, fmt = *) "[IE] vtk output realisations: Number not in [1,20]" 
     end if
     if (num_out_realisations > nsmpl-1) then
        write(unit = *, fmt = *) "[IE] vtk_output realisations: More output than realisations"
        stop 1
     end if
  end if

  ! Prescribed velocity only possible if: No stochastic bc, no active
  ! head and prescribed velocity is on in forward
  if (pres_vel_switch) then
     if (stoch_bc_switch) then
        write(unit = *, fmt = *) "[IE] pres_vel_switch and stoch_bc_switch should not"
        write(unit = *, fmt = *) "     be turned on at the same time                 "
        stop 1
     end if
     if (head_active)  then
        do ib = first_flow, last_flow
           if (ibc_data(ib,cbc_bt) .ne. 1) then
              write(unit = *, fmt = *) "[IE] all head boundaries should be Dirichlet"
              stop 1
           end if
        end do
     end if
     if (.not. vdefaultswitch) then
        write(unit = *, fmt = *) "[IE] pres_vel_switch needs vdefaultswitch to be true"
     end if
     if (num_pres_vel>3 .or. num_pres_vel<0) then
        write(unit = *, fmt = *) "[IE] num_pres_vel with wrong value: ", num_pres_vel
        write(unit = *, fmt = *) "     It should be between 0 and 3."
     end if
  end if

  ! Temperature difference
  if ( (tempdiff_switch) .and. (.not. assitype) ) then
     write(unit = *, fmt = *) "[IE] Temperature differences only with"
     write(unit = *, fmt = *) "     assitype sequential"
     stop 1
  end if

  ! Pilot point should not interact with other methods
  if (pp_switch) then
     if (ns_switch) then
        write(unit = *, fmt = *) "[IE] Pilot Points and Normal score not together"
     end if
     if (cov_loc_switch) then
        write(unit = *, fmt = *) "[IE] Pilot Points and Covariance localisation not together"
     end if
     if (dual_enkf_switch) then
        write(unit = *, fmt = *) "[IE] Pilot Points and Dual EnKF not together"
     end if
     if (hybrid_switch) then
        write(unit = *, fmt = *) "[IE] Pilot Points and Hybrid EnKF not together"
     end if
     if (iterative_switch) then
        write(unit = *, fmt = *) "[IE] Pilot Points and Iterative EnKF not together"
     end if

     ! Check pp entries
     do i = 1, num_pp

        i_pp = pp_inds(i,1)
        j_pp = pp_inds(i,2)
        k_pp = pp_inds(i,3)

        ! Pilot point indices in range of model
        if(i_pp < 1 .or. i_pp > i0) then
           write(unit = *, fmt = *) "[IE] i_pp out of range"
           stop 1
        else if(j_pp < 1 .or. j_pp > j0) then
           write(unit = *, fmt = *) "[IE] j_pp out of range"
           stop 1
        else if(k_pp < 1 .or. k_pp > k0) then
           write(unit = *, fmt = *) "[IE] k_pp out of range"
           stop 1
        end if

        ! Pilot point variable in state vector
        select case (pp_ivar)
        case (1,2,3,4,5,6)
           if(act_s(pp_ivar) /= 1) then
              write(unit = *, fmt = *) "[IE] pp_ivar=", pp_ivar ,&
                   ", but ",enkf_variable_names(pp_ivar)," not in state_vector"
              stop 1
           end if
        case default
           write(unit = *, fmt = *) "[IE]: pp_ivar out of range"
           stop 1
        end select

        ! Pilot points in ascending order
        if(i > 1) then
           i_pp_bef = pp_inds(i-1,1)
           j_pp_bef = pp_inds(i-1,2)
           k_pp_bef = pp_inds(i-1,3)
           if(k_pp == k_pp_bef) then
              if(j_pp == j_pp_bef) then
                 if(i_pp == i_pp_bef) then
                    write(unit = *, fmt = *) "[IE]: pp_inds in ascending order",&
                         ", problem at i_pp_bef, i_pp=", i_pp_bef, i_pp
                    stop 1
                 else if(i_pp < i_pp_bef) then
                    write(unit = *, fmt = *) "[IE]: pp_inds in ascending order",&
                         ", problem at ipp_bef, i_pp=", i_pp_bef, i_pp
                    stop 1
                 end if
              else if (j_pp < j_pp_bef) then
                 write(unit = *, fmt = *) "[IE]: pp_inds in ascending order",&
                      ", problem at j_pp_bef, j_pp=", j_pp_bef, j_pp
                 stop 1
              end if
           else if(k_pp < k_pp_bef) then
              write(unit = *, fmt = *) "[IE]: pp_inds in ascending order",&
                   ", problem at k_pp_bef, k_pp=", k_pp_bef, k_pp
              stop 1
           end if
        end if
     end do

     select case (pp_get_out)
     case ('get','out','def')
        write (unit = *, fmt = *) "[OK]: pp_get_out"
     case default
        write (unit = *, fmt = *) "[IE]: pp_get_out should be 'get' or 'out' of 'def'."
     end select
  end if

  write(unit = *, fmt = *) " No input error detected"
  write(unit = *, fmt = *)

end subroutine enkf_input_check
