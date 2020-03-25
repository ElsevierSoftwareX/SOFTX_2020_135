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

!> @brief SIMUL/ENKF controller
!> @param[in] iter_enkf ENKF iteration counter
!> @details
!> Main loop over all ensemble members
!> - 1. Read-in/copy-in/initialise parameters for each member
!>   + 1.1 Read input file `MODEL.enkf`, read observations and possibly true.
!> - 2. generate parameter-initialisation/random-seed for each member
!> - 3. compute the normal forward simulation for all members until data is available
!>   + 3.1 If `runmode < 2`, just compute forward simulation
!> - 4. perform the enkf analysis based on Evensen-algorithm and routines
!>   + 4.1 Around EnKF-update: Outputs and EnKF-variant modifications
!> - 5. update the ensemble and continue forward simulation
!> - 6. collect/write the output/statistical values
      SUBROUTINE enkf_iter(iter_enkf)

        use arrays, only: &
             vdefaultswitch

        use mod_genrl, only: &
             runmode, tec_out, write_disable,&
             write_iter_disable,&
             write_smonitor, head_active

        use mod_simul, only: &
             enkf_post, maxiter_enkf

        use mod_enkf, only: &
             assitype,&
             checktrue,&
             err_cov,&
             irobs,&
             nrens,&
             nrobs_int,&
             ttstart,&
             ttend,&
             iter_out,&
             sflags,&
             iassout,&
             iassout_single_start,&
             iassout_single,&
             iassout_cov_ref_start,&
             iassout_cov_ref,&
             ns_switch,&
             switch_cov_ref,&
             cov_loc_switch,&
             dual_enkf_switch,&
             stoch_bc_switch,&
             hybrid_switch,&
             general_seed_switch,&
             iterative_switch,&
             pres_vel_switch,&
             tcon_switch,&
             pp_switch,&
             switch_so_ini,&
             general_seed_seed,&
             stoch_bc_seed_seed,&
             pres_vel_seed,&
             vtk_out_enkf,&
             vtk_out_stddev,&
             obs_standard_out,&
             compute_mean

        implicit none

        include 'gslib.inc'
        include 'OMP_TOOLS.inc'

        integer, intent(in) :: iter_enkf

        integer :: i
       
        ! Disable usual output
        write_disable = .true.
        write_iter_disable = .true.
        write_smonitor = .true.

!------------------------------------------------------
!------------- ENKF INPUT -----------------------------
!------------------------------------------------------
        IF (runmode>=2) THEN

           call enkf_log(0,0)
           
           ! EnKF-Input
           call enkf_input()

           ! Check EnKF-Input
           call enkf_input_check()

           ! Set nrobs_int for Dual EnKF/Iterative EnKF
           call enkf_nrobs_int()

           ! State vector length
           call enkf_state_len()

           ! TRUE-Input
           IF (checktrue) THEN
              call enkf_alloc_true()
              call enkf_read_true()
           END IF

           ! OBS-Input
           call enkf_read_obs_len()
           call enkf_alloc_obs()
           call enkf_read_obs()

        END IF

!------------------------------------------------------
!------------ ENKF INPUT END --------------------------
!------------------------------------------------------

!-------------------------------------------------------
!-------------- SGSIM / OLD RESTORE  -------------------
!-------------------------------------------------------

        if (vdefaultswitch .and. head_active) then
           call enkf_velocity_dbc(1)
           call enkf_velocity_dbcold()
        end if
           
        IF (iter_enkf == 1 .OR. runmode<2) THEN
           select case (runmode)
           case (0,1)
              ! Pure forward
              call enkf_only_forward(iter_enkf)
           case(2,3)
              ! Initialisation/simulation of stochastic ensembles
              call enkf_simulate(iter_enkf)
              if(switch_so_ini) then
                 call enkf_output_single_cell('ini')
              end if
              ! Output: enkf.log
              call enkf_log(1,0)
           case default
              write(unit = *, fmt = *) 'Wrong runmode specification (0,1,2,3).'
              stop
           end select
        ELSE
           ! Restoring from previous EnKF-iteration
           call enkf_old_restore()
        END IF

        IF (runmode>=2) THEN
           
           if(stoch_bc_switch) then
              call enkf_set_random_seed(stoch_bc_seed_seed(1),stoch_bc_seed_seed(2))
              call enkf_stoch_bc()
           end if

           if(tcon_switch) then
              call enkf_tcon()
           end if
           
           if(pres_vel_switch) then
              call enkf_set_random_seed(pres_vel_seed(1),pres_vel_seed(2))
              call enkf_pres_vel()
           end if
           
           if(general_seed_switch) then
              call enkf_set_random_seed(general_seed_seed(1),general_seed_seed(2))
           end if
           
        END IF
!-------------------------------------------------------
!--------END IF SGSIM / OLD RESTORE  -------------------
!-------------------------------------------------------

        ! save start timestep
        call enkf_simtime(1,iter_enkf)

!--------------------------------------------------------
!------------------ CHECK ENKF MODE ---------------------
!--------------------------------------------------------
        IF (runmode>=2) THEN
           
    !------------------------------------------------------
    !-------------------- MAIN ENKF LOOP ------------------
    !------------------------------------------------------
           DO irobs = 1, nrobs_int

              !-------------------------------------------------
              !---------- INITIALIZE ENKF VARIANTS -------------
              !-------------------------------------------------

              if(dual_enkf_switch) then
                 call enkf_dual_enkf()
              end if
              if(iterative_switch) then
                 call enkf_iterative_enkf()
              end if
              
              !-------------------------------------------------
              !---------- END INITIALIZE ENKF VARIANTS ---------
              !-------------------------------------------------

              !-------------------------------------------------
              !---------- FORWARD SIMULATION -------------------
              !-------------------------------------------------

              ! Forward until observation time
              call enkf_prepare_forward()
              call forward_multi_compute(ttstart, ttend, iter_out, sflags, nrens)
              call enkf_wrapup_forward()

              !-------------------------------------------------
              !-------- END FORWARD SIMULATION -----------------
              !-------------------------------------------------
              
              !-------------------------------------------------
              !---------- MAKE ENKF STATE VECTORS --------------
              !-------------------------------------------------
              ! Output: enkf.log
              call enkf_log(2,0)

              call enkf_alloc_state_vector()
              call enkf_make_state_vector()
              call enkf_add_system_noise()

              if(cov_loc_switch) then
                 call enkf_make_cov_loc()
              end if

              call enkf_calc_statistics()
              if(switch_cov_ref) then
                 call enkf_compute_covs()
              end if

              !-------------------------------------------------
              !-------- END MAKE ENKF STATE VECTORS ------------
              !-------------------------------------------------

              !-------------------------------------------------
              !------------ BEFORE ASSIMILATION ----------------
              !-------------------------------------------------

              ! Output every tenth irobs
              if( mod(irobs,10) == 0 .and. obs_standard_out) then
                 write(unit = *, fmt = *) 
                 write(unit = *, fmt = *) '----------------------------------------------------'
                 write(unit = *, fmt = *) '----- Before observation #  ',   irobs, '    -------'
                 write(unit = *, fmt = *) '----------------------------------------------------'
              end if

              ! SAVE: Stddev
              if(vtk_out_stddev) then
                 if( irobs == 1 ) then
                    call enkf_alloc_stddev()
                 end if
                 call enkf_stddev('bef')
              end if

              ! SAVE: Residual
              IF (checktrue) then
                 if( irobs == 1) then
                    call enkf_alloc_resid()
                 end if
                 call enkf_resid('bef')
              END IF

              ! OUTPUT: ensemble averages/deviations
              IF ( (irobs==1) .or. (mod(irobs,iassout)==0) ) then
                 !--- tec ---
                 if(tec_out) then
                    call enkf_output_tec('bef')
                 end if
                 !-- vtk ----
                 if(vtk_out_enkf) then
                    call enkf_output_vtk('bef')
                 end if
              END IF

              ! OUTPUT: Cell distributions
              if ( (irobs-iassout_single_start>=0) &
                   .and. (mod(irobs-iassout_single_start,iassout_single)==0) ) then
                 call enkf_output_single_cell('bef')
              end if

              ! OUTPUT: state covariance matrix
              IF ((irobs==nrobs_int) .and. (err_cov)) then 
                 CALL enkf_state_covari('bef')
              END IF

              ! OUTPUT: single covariances
              if ((irobs-iassout_cov_ref_start>=0) &
                   .and. (mod(irobs-iassout_cov_ref_start,iassout_cov_ref)==0) &
                   .and. switch_cov_ref) then
                 call enkf_output_covs('bef')
              end if

              if(tcon_switch) then
                 call enkf_output_tcon('bef')
              end if
              
              if(pres_vel_switch) then
                 call enkf_output_pres_vel('bef')
              end if
              
              !-------------------------------------------------
              !------------- END BEFORE ASSIMILATION -----------
              !-------------------------------------------------

              !-------------------------------------------------
              !------------- ENKF VARIANTS BEFORE --------------
              !-------------------------------------------------
              if(ns_switch .eqv. .true.) then
                 call enkf_normal_score()
                 call enkf_normal_score_obsvars()
                 call enkf_normal_score_observations_obsvars()
              end if

              if(irobs==1 .and. hybrid_switch) then
                 call enkf_hybrid_cov()
              end if

              if(pp_switch) then
                 call enkf_set_pp_inds()
                 if(irobs == 1) then
                    call enkf_cov_pp()
                 end if
                 call enkf_nstate_pp('bef')
                 call enkf_state_vector_pp()
              end if
              !-------------------------------------------------
              !------------- END ENKF VARIANTS BEFORE ----------
              !-------------------------------------------------

              !-------------------------------------------------
              !------------ ASSIMILATIONS ----------------------
              !-------------------------------------------------
              IF (assitype) THEN
                 ! Sequential assimilation
                 call enkf_assim_sequ()
              ELSE 
                 ! combined assimilation
                 call enkf_assim_comb()
              END IF
              !-------------------------------------------------
              !------------ END ASSIMILATIONS ------------------
              !-------------------------------------------------

              !-------------------------------------------------
              !------------- ENKF VARIANTS AFTER ---------------
              !-------------------------------------------------
              if(ns_switch .eqv. .true.) then
                 call enkf_normal_score_back()
                 call enkf_normal_score_dealloc()
              end if

              if(pp_switch) then
                 call enkf_nstate_pp('aft')
                 call enkf_proj_pp()
              end if

              !-------------------------------------------------
              !------------- END ENKF VARIANTS AFTER -----------
              !-------------------------------------------------

              !-------------------------------------------------
              !-------- UPDATE VARIABLES -----------------------
              !-------------------------------------------------
              call enkf_update_vars()

              call enkf_calc_statistics()
              if(switch_cov_ref) then
                 call enkf_compute_covs()
              end if
              !-------------------------------------------------
              !-------- END UPDATE VARIABLES -------------------
              !-------------------------------------------------

              !-------------------------------------------------
              !----------------- AFTER ASSIMILATION ------------
              !-------------------------------------------------

              ! SAVE: Stddev
              if(vtk_out_stddev) then
                 call enkf_stddev('aft')
              end if

              ! SAVE: Residual
              IF (checktrue) THEN
                 call enkf_resid('aft')
              END IF

              ! OUTPUT: ensemble averages/deviations
              IF ((mod(irobs,iassout)==0) .OR. (irobs==nrobs_int)) then
                 !--tec --
                 if(tec_out) then
                    call enkf_output_tec('aft')
                 end if
                 !-- vtk ---
                 if(vtk_out_enkf) then
                    call enkf_output_vtk('aft')
                 end if
              END IF

              ! OUTPUT: cell distributions
              if ( (irobs-iassout_single_start>=0) .and. &
                   (mod(irobs-iassout_single_start,iassout_single)==0) ) then
                 call enkf_output_single_cell('aft')
              end if

              ! OUTPUT: state covariance matrix
              IF ((irobs==nrobs_int) .and. (err_cov)) then 
                 CALL enkf_state_covari('aft')
              END IF

              ! OUTPUT: single covariances
              if ((irobs-iassout_cov_ref_start>=0) &
                   .and. (mod(irobs-iassout_cov_ref_start,iassout_cov_ref)==0) &
                   .and. switch_cov_ref) then
                 call enkf_output_covs('aft')
              end if

              if(tcon_switch) then
                 call enkf_output_tcon('aft')
              end if
              
              if(pres_vel_switch) then
                 call enkf_output_pres_vel('aft')
              end if
              
              !-------------------------------------------------
              !------------- END AFTER ASSIMILATION ------------
              !-------------------------------------------------

              ! deallocate system vector arrays for this assimilation step
              call enkf_dealloc_state_vector()

           END DO
    !-------------------------------------------------
    !----------- END MAIN ENKF LOOP -----------------
    !-------------------------------------------------

           ! OUTPUT: All over standard deviations
           if(vtk_out_stddev) then
              call enkf_output_stddev()
           end if

           ! OUTPUT: All over residuals
           IF (checktrue) THEN
              call enkf_output_resid()
           END IF

           ! FINAL ENSEMBLE AVERAGES
           if (compute_mean) then
              call enkf_compute_mean()
           end if

           ! RECOMPUTE full model (all ensembles) for the last ENKF iteration
           IF (iter_enkf==maxiter_enkf .AND. enkf_post>=1) THEN
              call enkf_postcompute()
           END IF

           ! Deallocate EnKF Arrays
           call enkf_dealloc()

        END IF

!--------------------------------------------------------
!---- END IF CHECK ENKF MODE ----------------------------
!--------------------------------------------------------
        
        ! output related
        write_disable = .false.
        write_iter_disable = .false.
        write_smonitor = .false.

        if(runmode >= 2) then
           ! Output: enkf.log
           call enkf_log(7,0)
        end if

        ! restore start timestep (itimestep_0)
        call enkf_simtime(5,iter_enkf)

      END SUBROUTINE enkf_iter
