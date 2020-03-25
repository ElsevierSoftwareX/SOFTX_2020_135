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

!> @brief Read data from the observation file
!> @details
!> __USED__: Number of observation times, maximum number of
!> observation points at a given observation time and the name of the
!> observation file.
!>
!> __SET__: Observation times, number of observation points at each
!> observation time, activity of the variables, indices and values of
!> the variables at each observation point.
subroutine enkf_read_obs()

  use mod_enkf, only: &
       nrobs_int,&
       nrobs_int_pure,&
       iterative_nrobs_int,&
       max_obs_loc,&
       obs_name,&
       dual_enkf_switch,&
       iterative_switch,&
       act_o,&
       num_enkf_vars,&
       var_obs,&
       enkf_log_out,&
       iterative_doubleupdate,&
       obst,&
       nrobs_loc,&
       i_obs,&
       j_obs,&
       k_obs

  implicit none

  integer :: iobs, iloc, iobs_dual_enkf, remainder_dual_enkf, ienkfvar
  integer :: iiter, numiter, lenrep, lenbef, iobsend
  
  integer :: idu

  obs_name = trim(adjustl(obs_name))
  open(unit = 36, file = obs_name)
  read(unit = 36, fmt = *)

  !----------------------------------------------------------------------
  !----------------------------------------------------------------------
  !                Typical lines in an observational file 
  !  (First line: dummy; second line: iobs-loop; third line: iloc-loop)
  !
  !
  !         i                j             k               h               T              conc               kz     lz     por 
  !     1  0.0000000000E+00    16     0   0   1   0   0   0
  !             35             35              1  4.4000000000E+02  2.3000000000E+01  1.0000000000E-06  5.0000000000E-10  3.E+00  1.E-01
  !----------------------------------------------------------------------
  !----------------------------------------------------------------------

  ! read observation time obst and nr of observation location
  ! at that time and active parameters/variables 

  if(.not. dual_enkf_switch .and. .not. iterative_switch) then
     do iobs = 1, nrobs_int
        
        read(unit = 36, fmt = *) idu, obst(iobs), nrobs_loc(iobs), &
             act_o(iobs,1), act_o(iobs,2), act_o(iobs,3), &
             act_o(iobs,4), act_o(iobs,5), act_o(iobs,6)
        
        !----------------------------------------------------------------------
        if(enkf_log_out) then
           write(unit = 37, fmt = *) '*** read ', nrobs_loc(iobs), &
                ' data for interval ', iobs, ' of ', nrobs_int, &
                ' at time ', obst(iobs)
           write(unit = 37, fmt = *) 'Code for active in obs vector h, t, c, kz, &
                lz, por ', act_o(iobs,1), act_o(iobs,2), &
                act_o(iobs,3), act_o(iobs,4), act_o(iobs,5), &
                act_o(iobs,6)
           write(unit = 37, fmt = *)
        end if
        !---------------------------------------------------------------------

        ! read the observed values vor different variables/rock properties
        do iloc = 1, nrobs_loc(iobs)
           read(unit = 36, fmt = *) i_obs(iobs,iloc), j_obs(iobs,iloc), &
                k_obs(iobs,iloc), var_obs(iobs,iloc,1), &
                var_obs(iobs,iloc,2), var_obs(iobs,iloc,3), &
                var_obs(iobs,iloc,4), var_obs(iobs,iloc,5), &
                var_obs(iobs,iloc,6)
           if (var_obs(iobs,iloc,3)<=0.D0) var_obs(iobs,iloc,3) = 1.D-30
           if (var_obs(iobs,iloc,4)<=0.D0) var_obs(iobs,iloc,4) &
                = 1.D-30
        end do
     end do
     close(unit = 36)
  end if

  if(dual_enkf_switch) then
     do iobs = 1, nrobs_int
        
        remainder_dual_enkf = mod(iobs,2) !1: parameter update, 2:state update
        iobs_dual_enkf = (iobs+remainder_dual_enkf)/2 !Real observation, each is taken twice
        select case(remainder_dual_enkf)
           !Parameter update: Read in everything normally
        case (1)
           read(unit = 36, fmt = *) idu, obst(iobs), nrobs_loc(iobs), &
                act_o(iobs,1), act_o(iobs,2), act_o(iobs,3), &
                act_o(iobs,4), act_o(iobs,5), act_o(iobs,6)
        
           ! read the observed values for different variables/rock properties
           do iloc = 1, nrobs_loc(iobs_dual_enkf)
              read(unit = 36, fmt = *) i_obs(iobs,iloc), j_obs(iobs,iloc), &
                   k_obs(iobs,iloc), var_obs(iobs,iloc,1), &
                   var_obs(iobs,iloc,2), var_obs(iobs,iloc,3), &
                   var_obs(iobs,iloc,4), var_obs(iobs,iloc,5), &
                   var_obs(iobs,iloc,6)
              if (var_obs(iobs,iloc,3)<=0.D0) var_obs(iobs,iloc,3) = 1.D-30
              if (var_obs(iobs,iloc,4)<=0.D0) var_obs(iobs,iloc,4) &
                   = 1.D-30
           end do
           !State update: Take same values for observations as for parameter update
        case (0)
           obst(iobs) = obst(iobs-1)
           nrobs_loc(iobs) = nrobs_loc(iobs-1)
           act_o(iobs,:) = act_o(iobs-1,:)
           do iloc = 1, nrobs_loc(iobs_dual_enkf)
              i_obs(iobs,iloc) = i_obs(iobs-1,iloc)
              j_obs(iobs,iloc) = j_obs(iobs-1,iloc)
              k_obs(iobs,iloc) = k_obs(iobs-1,iloc)
              var_obs(iobs,iloc,:) = var_obs(iobs-1,iloc,:)
           end do
        end select
     end do
     close(unit = 36)
  end if


  if(iterative_switch) then
     numiter = (nrobs_int_pure-1)/iterative_nrobs_int
     
     do iiter = 1, numiter+1
        lenrep = (iiter-1)*iterative_nrobs_int
        lenbef = iterative_nrobs_int*(iiter-1)*iiter/2

        !Copy already read measurements
        if(lenbef > 0) then
           do iobs = lenbef+1, lenbef+lenrep 
              obst(iobs) = obst(iobs-lenrep)
              nrobs_loc(iobs) = nrobs_loc(iobs-lenrep)
              do ienkfvar = 1, num_enkf_vars
                 if(iterative_doubleupdate) then
                    act_o(iobs,ienkfvar) = act_o(iobs-lenrep,ienkfvar)
                 else
                    act_o(iobs,ienkfvar) = 0
                 end if
              end do
              do iloc = 1, nrobs_loc(iobs)
                 i_obs(iobs,iloc) = i_obs(iobs-lenrep,iloc)
                 j_obs(iobs,iloc) = j_obs(iobs-lenrep,iloc)
                 k_obs(iobs,iloc) = k_obs(iobs-lenrep,iloc)
                 var_obs(iobs,iloc,:) = var_obs(iobs-lenrep,iloc,:)
              end do
           end do
        end if

        ! Read new measurements
        if(iiter < numiter+1) then
           iobsend = lenbef + lenrep + iterative_nrobs_int
        else
           iobsend = nrobs_int
        end if
        do iobs = lenbef+lenrep+1, iobsend
           ! Read in info about obs
           read(unit = 36, fmt = *) idu, obst(iobs), nrobs_loc(iobs), &
                act_o(iobs,1), act_o(iobs,2), act_o(iobs,3), &
                act_o(iobs,4), act_o(iobs,5), act_o(iobs,6)
           ! Read in obs values
           do iloc = 1, nrobs_loc(iobs)
              read(unit = 36, fmt = *) i_obs(iobs,iloc), j_obs(iobs,iloc), &
                   k_obs(iobs,iloc), var_obs(iobs,iloc,1), &
                   var_obs(iobs,iloc,2), var_obs(iobs,iloc,3), &
                   var_obs(iobs,iloc,4), var_obs(iobs,iloc,5), &
                   var_obs(iobs,iloc,6)
              if (var_obs(iobs,iloc,3)<=0.D0) var_obs(iobs,iloc,3) = 1.D-30
              if (var_obs(iobs,iloc,4)<=0.D0) var_obs(iobs,iloc,4) &
                   = 1.D-30
           end do
        end do
        
     end do
     
     close(unit = 36)
  end if
  
  ! act_o_var entries in one array act_o
  do ienkfvar = 1, num_enkf_vars
     do iobs = 1, nrobs_int
        if( .not.( act_o(iobs,ienkfvar) == 1 .or. act_o(iobs,ienkfvar) == 0)) then
           write(unit = *, fmt = *) "[E1] Error in enkf_read_obs.f90 act_o(iobs,ienkfvar) not 0 or 1: "&
                , iobs, ienkfvar, act_o(iobs,ienkfvar)
        end if
     end do
  end do

  ! var_obs entries in one array
  do ienkfvar = 1, num_enkf_vars
     do iobs = 1, nrobs_int
        do iloc = 1, nrobs_loc(iobs)
           select case (ienkfvar)
           case(4)
              if(var_obs(iobs,iloc,ienkfvar) < 1.0d59) then
                 var_obs(iobs,iloc,ienkfvar) = dlog10(var_obs(iobs,iloc,ienkfvar) + 1.0d-30)
              else
                 var_obs(iobs,iloc,ienkfvar) = 1.0d60
              end if
           end select
        end do
     end do
  end do

end subroutine enkf_read_obs
