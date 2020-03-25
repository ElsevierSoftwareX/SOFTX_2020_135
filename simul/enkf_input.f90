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

!> @brief read parameters from  ENKF input file
!> @details
!> Read EnKF input parameters. \n
!> `# nrobs_int`: Set nrobs_int. \n
!> `# mode_analysis`: Set mode_analysis. \n
!> `# truncation`: Set truncation. \n
!> `# covmodel`: Set and trim covmodel. \n
!> `# correlation lengths`: Set clh, clv. \n
!> `# rexact`: Set rexact. \n
!> `# fixsamp`: Set fixsamp. \n
!> `# checktrue`: Set checktrue. \n
!> `# observation variances`: Set obsvar array. \n
!> `# system variances`: Set sysvar array. \n
!> `# state vector activity`: Set act_s array. \n
!> `# ex_unit`: Set ex_unit. \n
!> `# assimod`: Set assitype by reading assimod. \n
!> `# normobs`: Set normobs. \n
!> `# damping factors`: Set assidamp and tracdamp. \n
!> `# iassout`: Set iassout. \n
!> `# err_cov`: Set err_cov. \n
!> `# true file names`: Set true_name and true_chem_name. \n
!> `# observation file names`: Set obs_name. \n
!> `# single cell output times`: Set iassout_single_start and
!> iassout_single. \n
!> `# single cell output`: Set mat_single_out and
!> num_single_out. \n
!> `# normal score`: Set ns_switch and ns_backfactor. \n
!> `# reference distribution`: Set num_ref_dist, ref_dist_switch,
!> ref_dist_var and ref_dist_file_name. \n
!> `# reference point covariance times`: `: Set
!> iassout_cov_ref_start and iassout_cov_ref. \n
!> `# reference point covariance`: Set mat_cov_ref and
!> switch_cov_ref. \n
!> `# covariance localisation`: Set cov_loc_switch, cov_loc_lenx,
!> cov_loc_leny and cov_loc_lenz. \n
!> `# cell centered`: Set cell_centered. \n
!> `# dualenkf`: Set dual_enkf_switch. \n
!> `# stochbc`: Set stoch_bc_switch, num_stoch_bc,
!> stoch_bc_seed_seed, stoch_bc_stddevs. \n
!> `# generalseed`: Set general_seed_switch and
!> general_seed_seed. \n
!> `# hybridenkf`: Set hybrid_switch and hybrid_alpha. \n
!> `# assimstp output`: Set assimstp_switch. \n
!> `# vtk output enkf`: Set vtk_out_enkf. \n
!> `# vtk output stddev`: Set vtk_out_stddev. \n
!> `# vtk output resid`: Set vtk_out_resid. \n
!> `# vtk output covs`: Set vtk_out_covs. \n
!> `# vtk output realisations`: Set vtk_out_realisations and
!> num_out_realisations. \n
!> `# obs standard out`: Set obs_standard_out. \n
!> `# enkf log output`: Set enkf_log_out. \n
!> `# analysis matrix output`: Set ana_mat_out. \n
!> `# compute realisation output`: Set comp_real_out. \n
!> `# iterativeenkf`: Set iterative_switch, iterative_nrobs_int,
!> iterative_doubleupdate. \n
!> `# prescribed velocity`: Set pres_vel_switch, num_pres_vel,
!> pres_vel_seed, vdefault_stddevs, vdefault_sysvarmem. \n
!> `# temperature difference`: Set tempdiff_switch,
!> num_tempdiff, tempdiff_inds. \n
!> `# tcon`: Set tcon_switch, tcon_stddev, tcon_sysvarmem,
!> tcon_inds. \n
!> `# pilot point`: Set pp_switch, num_pp, pp_get_out, pp_ivar
!> pp_inds. \n
!> `# compute mean`: Set compute_mean. \n
      SUBROUTINE enkf_input()

        use mod_genrl, only:&
             key_char

        use mod_genrlc, only:&
             project

        use mod_enkf, only:&
             covmodel,&
             rexact,&
             fixsamp,&
             mode_analysis, &
             truncation,&
             checktrue,&
             clh,&
             clv,&
             nrobs_int,  &
             assitype, &
             normobs,&
             assidamp,&
             tracdamp,&
             iassout,&
             err_cov,&
             ex_unit,&
             true_name,&
             true_chem_name,&
             obs_name,&
             iassout_single_start,&
             iassout_single,&
             num_single_out,&
             ns_switch,&
             ns_backfactor,&
             num_ref_dist,&
             ref_dist_file_name,&
             ref_dist_switch,&
             ref_dist_var,&
             switch_cov_ref,&
             num_cov_ref,&
             iassout_cov_ref_start,&
             iassout_cov_ref,&
             cov_loc_switch,&
             cov_loc_lenx,&
             cov_loc_leny,&
             cov_loc_lenz,&
             cell_centered,&
             dual_enkf_switch,&
             stoch_bc_switch,&
             num_stoch_bc,&
             stoch_bc_seed_seed,&
             stoch_bc_stddevs,&
             general_seed_seed,&
             general_seed_switch,&
             hybrid_switch,&
             hybrid_alpha,&
             assimstp_switch,&
             vtk_out_enkf,&
             vtk_out_stddev,&
             vtk_out_resid,&
             vtk_out_covs,&
             vtk_out_realisations,&
             num_out_realisations,&
             enkf_log_out,&
             ana_mat_out,&
             comp_real_out,&
             iterative_switch,&
             iterative_nrobs_int,&
             iterative_doubleupdate,&
             pres_vel_switch,&
             pres_vel_seed,&
             num_pres_vel,&
             tempdiff_switch,&
             num_tempdiff,&
             tcon_switch,&
             num_tcon,&
             tcon_stddev,&
             tcon_sysvarmem,&
             pp_switch,&
             num_pp,&
             pp_get_out,&
             pp_ivar,&
             mat_cov_ref,&
             mat_single_out,&
             num_enkf_vars,&
             act_s,&
             obsvar,&
             sysvar,&
             vdefault_stddevs,&
             vdefault_sysvarmem,&
             tempdiff_inds,&
             tcon_inds,&
             switch_so_ini,&
             pp_inds,&
             obs_standard_out,&
             compute_mean

        IMPLICIT NONE
        character (len=256) :: project_enkf
        character (len=100) :: assimod

        logical :: found
        external found
        character (len=80) :: line
        integer :: i_beg_arg, i_end_arg, i_arg, ienkfvar
        integer :: itp, itc, ipp
        integer :: i_pp, j_pp, k_pp
        integer :: i_dummy
        double precision :: r_dummy, s_dummy, t_dummy
        double precision :: j_dummy, k_dummy
        character (len=3) :: single_out_lin_log
        integer :: i_so, j_so, k_so, ivar_so
        integer :: i_cov_ref, j_cov_ref, k_cov_ref
        integer :: ivar_cov_ref

        write(project_enkf, fmt = *) trim(project), ".enkf"
        project_enkf = trim(adjustl(project_enkf))

        !Writing to enkf.log
        WRITE(37,*)
        WRITE(37,*) 'Enkf input:'
        WRITE(37,*) '----------'

        ! project_enkf = PROJECTNAME.enkf
        OPEN(unit=10,file=project_enkf)
        
        ! # nrobs_int
        if(found(10,key_char//' nrobs_int',line,.FALSE.)) then
           write(unit = *, fmt = *) '[R] : nrobs_int'
           read(unit = 10, fmt = *) nrobs_int
           write(unit = 37, fmt = *) 'nrobs_int', nrobs_int 
        else
           write(unit = *, fmt = *) '[E1] enkf_input.f90 (# nrobs_int) '
           stop
        end if

        ! # mode_analysis
        if(found(10,key_char//' mode_analysis',line,.FALSE.)) then
           write(unit = *, fmt = *) '[R] : mode_analysis'
           read(unit = 10, fmt = *) mode_analysis
           write(unit = 37, fmt = *) 'mode_analysis ', mode_analysis 
        else
           write(unit = *, fmt = *) '[E2] enkf_input.f90 (# mode_analysis) '
           stop
        end if

        ! # truncation
        if(found(10,key_char//' truncation',line,.FALSE.)) then
           write(unit = *, fmt = *) '[R] : truncation'
           read(unit = 10, fmt = *) truncation
           write(unit = 37, fmt = *) 'truncation ', truncation
        else
           write(unit = *, fmt = *) '[E3] enkf_input.f90 (# truncation) '
           stop
        end if

        ! # covmodel
        if(found(10,key_char//' covmodel',line,.FALSE.)) then
           write(unit = *, fmt = *) '[R] : covmodel'
           read(unit = 10, fmt = '(a100)') covmodel
           write(unit = 37, fmt = '(8Hcovmodel,12x,a20)') trim(covmodel) 
        else
           write(unit = *, fmt = *) '[E4] enkf_input.f90 (# covmodel) '
           stop
        end if

        ! # correlation lengths
        if(found(10,key_char//' correlation lengths',line,.FALSE.)) then
           write(unit = *, fmt = *) '[R] : clh, clv'
           read(unit = 10, fmt = *) clh, clv
           write(unit = 37, fmt = *) 'clh, clv ', clh, clv
        else
           write(unit = *, fmt = *) '[E5] enkf_input.f90 (# correlation lengths) '
           stop
        end if

        ! # rexact
        if(found(10,key_char//' rexact',line,.FALSE.)) then
           write(unit = *, fmt = *) '[R] : rexact'
           read(unit = 10, fmt = '(l1)') rexact
           write(unit = 37, fmt = *) 'rexact ', rexact
        else
           write(unit = *, fmt = *) '[E6] enkf_input.f90 (# rexact) '
           stop
        end if

        ! # fixsamp
        if(found(10,key_char//' fixsamp',line,.FALSE.)) then
           write(unit = *, fmt = *) '[R] : fixsamp'
           read(unit = 10, fmt = '(l1)') fixsamp
           write(unit = 37, fmt = *) 'fixsamp ', fixsamp
        else
           write(unit = *, fmt = *) '[E6] enkf_input.f90 (# fixsamp) '
           stop
        end if

        ! # checktrue
        if(found(10,key_char//' checktrue',line,.FALSE.)) then
           write(unit = *, fmt = *) '[R] : checktrue'
           read(unit = 10, fmt = '(l1)') checktrue
           write(unit = 37, fmt = *) 'checktrue ', checktrue
        else
           write(unit = *, fmt = *) '[E7] enkf_input.f90 (# checktrue) '
           stop
        end if

        ! # observation variances
        if(found(10,key_char//' observation variances',line,.FALSE.)) then
           write(unit = *, fmt = *) '[R] : obsvar(1), obsvar(2), obsvar(3), obsvar(4), obsvar(5), obsvar(6)'
           read(unit = 10, fmt = *) obsvar(1), obsvar(2), obsvar(3), obsvar(4), obsvar(5), obsvar(6)
           write(unit = 37, fmt = *) 'obsvar(1), obsvar(2), obsvar(3), obsvar(4), obsvar(5), obsvar(6)', &
                obsvar(1), obsvar(2), obsvar(3), obsvar(4), obsvar(5), obsvar(6)
        else
           write(unit = *, fmt = *) '[E8] enkf_input.f90 (# observation variances) '
           stop
        end if

        ! # system variances
        if(found(10,key_char//' system variances',line,.FALSE.)) then
           write(unit = *, fmt = *) '[R] : sysvar(1), sysvar(2), sysvar(3), sysvar(4), sysvar(5), sysvar(6)'
           read(unit = 10, fmt = *) sysvar(1), sysvar(2), sysvar(3), sysvar(4), sysvar(5), sysvar(6)
           write(unit = 37, fmt = *) 'sysvar(1), sysvar(2), sysvar(3), sysvar(4), sysvar(5), sysvar(6)', &
                sysvar(1), sysvar(2), sysvar(3), sysvar(4), sysvar(5), sysvar(6)
        else
           write(unit = *, fmt = *) '[E9] enkf_input.f90 (# system variances) '
           stop
        end if

        ! # state vector activity
        if(found(10,key_char//' state vector activity',line,.FALSE.)) then
           write(unit = *, fmt = *) '[R] : act_s'
           read(unit = 10, fmt = *) act_s(1), act_s(2), act_s(3), act_s(4), act_s(5), act_s(6)
           write(unit = 37, fmt = *) 'act_s ',&
                act_s(1), act_s(2), act_s(3), act_s(4), act_s(5), act_s(6)
        else
           write(unit = *, fmt = *) '[E10] enkf_input.f90 (# state vector activity) '
           stop
        end if

        ! # ex_unit
        if(found(10,key_char//' ex_unit',line,.FALSE.)) then
           write(unit = *, fmt = *) '[R] : ex_unit'
           read(unit = 10, fmt = *) ex_unit
           write(unit = 37, fmt = *) 'ex_unit', ex_unit
        else
           write(unit = *, fmt = *) '[E11] enkf_input.f90 (# ex_unit) '
           stop
        end if

        ! # assimod
        if(found(10,key_char//' assimod',line,.FALSE.)) then
           write(unit = *, fmt = *) '[R] : assimod'
           read(unit = 10, fmt = '(a100)') assimod
           write(unit = 37, fmt = '(7Hassimod,12x,a20)') assimod
        else
           write(unit = *, fmt = *) '[E12] enkf_input.f90 (# assimod) '
           stop
        end if
        
        if(trim(assimod) =='comb') then
           assitype= .FALSE.
        else if(trim(assimod) =='sequ') then
           assitype= .TRUE.
        else
           write(unit = *, fmt = *) "[E29] Error in enkf_input"
           stop
        end if

        ! # normobs
        if(found(10,key_char//' normobs',line,.FALSE.)) then
           write(unit = *, fmt = *) '[R] : normobs'
           read(unit = 10, fmt = '(l1)') normobs
           write(unit = 37, fmt = '(7Hnormobs,13x,l1)') normobs
        else
           write(unit = *, fmt = *) '[E13] enkf_input.f90 (# normobs) '
           stop
        end if

        ! # damping factors
        if(found(10,key_char//' damping factors',line,.FALSE.)) then
           write(unit = *, fmt = *) '[R] : assidamp, tracdamp'
           read(unit = 10, fmt = *) assidamp, tracdamp
           write(unit = 37, fmt = *) 'assidamp, tracdamp ', assidamp, tracdamp 
        else
           write(unit = *, fmt = *) '[E14] enkf_input.f90 (# damping factors) '
           stop
        end if

        ! # iassout
        if(found(10,key_char//' iassout',line,.FALSE.)) then
           write(unit = *, fmt = *) '[R] : iassout'
           read(unit = 10, fmt = *) iassout
           write(unit = 37, fmt = *) 'iassout ', iassout
        else
           write(unit = *, fmt = *) '[E15] enkf_input.f90 (# iassout) '
           stop
        end if

        ! err_cov
        if(found(10,key_char//' err_cov',line,.FALSE.)) then
           write(unit = *, fmt = *) '[R] : err_cov'
           read(unit = 10, fmt = '(l1)') err_cov
           write(unit = 37, fmt = '(7Herr_cov,13x,l1)') err_cov
        else
           write(unit = *, fmt = *) '[E16] enkf_input.f90 (# err_cov) '
           stop
        end if

        ! # true file names
        if(found(10,key_char//' true file names',line,.FALSE.)) then
           write(unit = *, fmt = *) '[R] : true_name'
           read(unit = 10, fmt = *) true_name
           write(unit = *, fmt = *) '[R] : true_chem_name'
           read(unit = 10, fmt = *) true_chem_name 
           write(unit = 37, fmt = *) true_name
           write(unit = 37, fmt = *) true_chem_name
        else
           write(unit = *, fmt = *) '[E17] enkf_input.f90 (# true file names) '
           stop
        end if

        ! # observation file names
        if(found(10,key_char//' observation file names',line,.FALSE.)) then
           write(unit = *, fmt = *) '[R] : obs_name'
           read(unit = 10, fmt = *) obs_name
           write(unit = 37, fmt = *) obs_name
        else
           write(unit = *, fmt = *) '[E18] enkf_input.f90 (# observation file names) '
           stop
        end if

        ! # single cell output times
        if(found(10,key_char//' single cell output times',line,.FALSE.)) then
           write(unit = *, fmt = *) '[R] : iassout_single_start, iassout_single'
           read(unit = 10, fmt = *) iassout_single_start, iassout_single
           write(unit = 37, fmt = *) 'iassout_single_start, iassout_single ',&
                iassout_single_start, iassout_single
           if(iassout_single_start==0) then
              switch_so_ini = .True.
           else
              switch_so_ini = .False.
           end if
        else
           write(unit = *, fmt = *) '[E18] enkf_input.f90 (# single cell output times) '
           stop
        end if

        ! # single cell output
        if(found(10,key_char//' single cell output',line,.FALSE.)) then

           ! Get the number of records
           call get_arg('records',line,i_beg_arg,i_end_arg)
           if( i_beg_arg < 1 .or. i_end_arg < i_beg_arg) then
              write(unit = *, fmt = *) '[E22] enkf_input.f90 (# single cell output, records=...'
              stop
           else
              read(line(i_beg_arg:i_end_arg),*) num_single_out
           end if           
           write(unit = 37, fmt = *) 'num_single_out ', num_single_out

           write(unit = *, fmt = *) '[R] : i_so, j_so, k_so'
           write(unit = *, fmt = *) '[R] : ivar_so, single_out_lin_log'
           allocate(mat_single_out(num_single_out,5))
           do i_arg = 1, num_single_out
              read(unit = 10, fmt = *) i_so, j_so, k_so,&
                   ivar_so, single_out_lin_log
              write(unit = 37, fmt = *) 'i_so, j_so, k_so ',&
                   i_so, j_so, k_so
              write(unit = 37, fmt = *) 'ivar_so, single_out_lin_log ',&
                   ivar_so, single_out_lin_log

              mat_single_out(i_arg,1) = i_so
              mat_single_out(i_arg,2) = j_so
              mat_single_out(i_arg,3) = k_so
              mat_single_out(i_arg,4) = ivar_so
              select case (single_out_lin_log)
              case('lin')
                 mat_single_out(i_arg,5) = 0
              case('log')
                 mat_single_out(i_arg,5) = 1
              case default
                 write(unit = *, fmt = *) '[E23] enkf_input.f90 (# single cell ouptut linlog)'
                 stop
              end select
              
           end do

        else
           write(unit = *, fmt = *) '[E19] enkf_input.f90 (# single cell output) '
           stop
        end if

        ! # normal score
        ns_switch = .false.
        if(found(10,key_char//' normal score',line,.FALSE.)) then
           write(unit = *, fmt = *) '[R] : ns_switch'
           read(unit = 10, fmt = *) ns_switch, ns_backfactor
           write(unit = 37, fmt = *) 'ns_switch', ns_switch
           write(unit = 37, fmt = *) 'ns_backfactor', ns_backfactor
        else
           write(unit = *, fmt = *) '[NR] : ns_switch'
           write(unit = 37, fmt = *) 'normal score false '
        end if

        ! # reference distribution
        ref_dist_switch = .false.
        if(found(10,key_char//' reference distribution',line,.FALSE.)) then

           ! Get the number of records
           call get_arg('records',line,i_beg_arg,i_end_arg)
           if( i_beg_arg < 1 .or. i_end_arg < i_beg_arg) then
              write(unit = *, fmt = *) '[E21] enkf_input.f90 (# reference distribution, records=...'
              stop
           else
              read(line(i_beg_arg:i_end_arg),*) num_ref_dist
           end if           
           write(unit = 37, fmt = *) 'num_ref_dist ', num_ref_dist

           if(num_ref_dist > num_enkf_vars) then
              write(unit = *, fmt = *) '[E21.1] enkf_input.f90, too many reference distributions.'
              stop
           end if
           
           write(unit = *, fmt = *) '[R] : ref_dist_switch'
           ref_dist_switch = .true.
           write(unit = 37, fmt = *) 'ref_dist_switch true '
           write(unit = *, fmt = *) '[R] : ref_dist_var, ref_dist_file_name'

           do i_arg = 1, num_ref_dist

              read(unit = 10, fmt = *) i_dummy, line
              ref_dist_var(i_arg) = i_dummy
              ref_dist_file_name(i_arg) = line
              write(unit = 37, fmt = *) ' ref_dist_var, ref_dist_file_name ',&
                   ref_dist_var(i_arg), ref_dist_file_name(i_arg)
           end do

        else
           write(unit = *, fmt = *) '[NR] : ref_dist_switch, ref_dist_var, ref_dist_file_name'
           write(unit = 37, fmt = *) 'ref_dist_switch false '
        end if

        ! # reference point covariance times
        if(found(10,key_char//' reference point covariance times',line,.FALSE.)) then
           write(unit = *, fmt = *) '[R] : iassout_cov_ref_start, iassout_cov_ref'
           read(unit = 10, fmt = *) iassout_cov_ref_start, iassout_cov_ref
           write(unit = 37, fmt = *) 'iassout_cov_ref_start, iassout_cov_ref ',&
                iassout_cov_ref_start, iassout_cov_ref
        else
           write(unit = *, fmt = *) '[E20] enkf_input.f90 (# reference point covariance times) '
           stop
        end if

        ! # reference point covariance
        switch_cov_ref = .false.
        if(found(10,key_char//' reference point covariance',line,.FALSE.)) then

           ! Get the number of records
           call get_arg('records',line,i_beg_arg,i_end_arg)
           if( i_beg_arg < 1 .or. i_end_arg < i_beg_arg) then
              write(unit = *, fmt = *) '[E21] enkf_input.f90 (# reference point covariance, records=...'
              stop
           else
              read(line(i_beg_arg:i_end_arg),*) num_cov_ref
           end if           
           write(unit = 37, fmt = *) 'num_cov_ref ', num_cov_ref

           if(num_cov_ref > 0) then
              write(unit = *, fmt = *) '[R] : switch_cov_ref'
              switch_cov_ref = .true.
              write(unit = 37, fmt = *) 'switch_cov_ref true '
           else
              write(unit = *, fmt = *) '[R] : switch_cov_ref'
              write(unit = 37, fmt = *) 'switch_cov_ref false '
           end if

           !Put all the values in the matrix mat_cov_ref
           write(unit = *, fmt = *) '[R] : i_cov_ref, j_cov_ref, k_cov_ref, ivar_cov_ref'
           allocate(mat_cov_ref(num_cov_ref,4))
           do i_arg = 1, num_cov_ref
              read(unit = 10, fmt = *) i_cov_ref, j_cov_ref, k_cov_ref, ivar_cov_ref
              write(unit = 37, fmt = *) ' i_cov_ref, j_cov_ref, k_cov_ref, ivar_cov_ref ',&
                   i_cov_ref, j_cov_ref, k_cov_ref, ivar_cov_ref
              mat_cov_ref(i_arg,1) = i_cov_ref
              mat_cov_ref(i_arg,2) = j_cov_ref
              mat_cov_ref(i_arg,3) = k_cov_ref
              mat_cov_ref(i_arg,4) = ivar_cov_ref
           end do

        else
           write(unit = *, fmt = *) '[NR] : switch_cov_ref,',&
                ' i_cov_ref, j_cov_ref, k_cov_ref, ivar_cov_ref'
           write(unit = 37, fmt = *) 'switch_cov_ref false '
        end if

        ! # covariance localisation
        cov_loc_switch = .false.
        if(found(10,key_char//' covariance localisation',line,.FALSE.)) then
           write(unit = *, fmt = *) '[R] : cov_loc_switch'
           write(unit = *, fmt = *) '[R] : cov_loc_len'
           read(unit = 10, fmt = *) cov_loc_switch, cov_loc_lenx, cov_loc_leny, cov_loc_lenz
           write(unit = 37, fmt = *) 'cov_loc_switch ', cov_loc_switch
           write(unit = 37, fmt = *) 'cov_loc_lens  ', cov_loc_lenx,&
                cov_loc_leny, cov_loc_lenz

        else
           write(unit = *, fmt = *) '[NR] : cov_loc_switch, cov_loc_len,'
           write(unit = 37, fmt = *) 'cov_loc_switch false '
        end if
        
        ! # cell centered
        cell_centered = .false.
        if(found(10,key_char//' cell centered',line,.FALSE.)) then
           write(unit = *, fmt = *) '[R] : cell centered'
           read(unit = 10, fmt = '(l1)') cell_centered
           write(unit = 37, fmt = *) 'cell_centered ', cell_centered
        else
           write(unit = *, fmt = *) '[NR] cell centered'
           write(unit = 37, fmt = *) 'cell centered false' 
        end if

        ! # dualenkf
        dual_enkf_switch = .false.
        if(found(10,key_char//' dualenkf',line,.FALSE.)) then
           write(unit = *, fmt = *) '[R] : dual_enkf_switch'
           read(unit = 10, fmt = '(l1)') dual_enkf_switch
           write(unit = 37, fmt = *) 'dual_enkf_switch ', dual_enkf_switch
        end if

        ! # stochbc
        stoch_bc_switch = .false.
        if(found(10,key_char//' stochbc',line,.FALSE.)) then
           
           ! Get the number of records
           call get_arg('records',line,i_beg_arg,i_end_arg)
           if( i_beg_arg < 1 .or. i_end_arg < i_beg_arg) then
              write(unit = *, fmt = *) '[E22] enkf_input.f90 (# stochbc, records=...)'
              stop
           else
              read(line(i_beg_arg:i_end_arg),*) num_stoch_bc
           end if
           write(unit = 37, fmt = *) 'num_stoch_bc ', num_stoch_bc
           

           if(num_stoch_bc > 0) then
              write(unit = *, fmt = *) '[R] : stochbc'
              stoch_bc_switch = .true.
              write(unit = 37, fmt = *) 'stoch_bc_switch true'
              read(unit = 10, fmt = *) stoch_bc_seed_seed(1), stoch_bc_seed_seed(2)
              write(unit = 37, fmt = *) 'stoch_bc_seed_seed ', stoch_bc_seed_seed
              read(unit = 10, fmt = *) (stoch_bc_stddevs(ienkfvar),ienkfvar=1,num_enkf_vars)
              write(unit = 37, fmt = *) (stoch_bc_stddevs(ienkfvar),ienkfvar=1,num_enkf_vars)
           else
              write(unit = *, fmt = *) '[NR]: stoch_bc_switch'
              write(unit = 37, fmt = *) 'stoch_bc_switch false' 
           end if
           
        end if

        ! # generalseed
        general_seed_switch = .false.
        if(found(10,key_char//' generalseed',line,.FALSE.)) then

           write(unit = *, fmt = *) '[R] : general_seed_switch, general_seed_seed'
           read(unit = 10, fmt = *) general_seed_switch, general_seed_seed(1), general_seed_seed(2)
           write(unit = 37, fmt = *) 'general_seed_switch ', general_seed_switch
           write(unit = 37, fmt = *) 'general_seed_seed ', general_seed_seed
        else
           write(unit = *, fmt = *) '[NR]: generalseed'
           write(unit = 37, fmt = *) 'general_seed_switch ', general_seed_switch
        end if

        ! # hybridenkf
        hybrid_switch = .false.
        hybrid_alpha = 0.0d0
        if(found(10,key_char//' hybridenkf',line,.FALSE.)) then
           write(unit = *, fmt = *) '[R] : hybridenkf'
           read(unit = 10, fmt = *) hybrid_switch, hybrid_alpha
           write(unit = 37, fmt = *) 'hybrid_switch true'
           write(unit = 37, fmt = *) 'hybrid_alpha ', hybrid_alpha
        else
           write(unit = *, fmt = *) '[NR]: hybridenkf'
           write(unit = 37, fmt = *) 'hybrid_switch false'
        end if

        ! # assimstp output
        assimstp_switch = .true.
        if(found(10,key_char//' assimstp output',line,.FALSE.)) then
           write(unit = *, fmt = *) '[R] : assimstp_switch'
           read(unit = 10, fmt = *) assimstp_switch
           write(unit = 37, fmt = *) 'assimstp_switch ', assimstp_switch
        else
           write(unit = *, fmt = *) '[NR]: assimstp_switch'
           write(unit = 37, fmt = *) 'assimstp_switch true'
        end if

        ! # vtk output enkf
        vtk_out_enkf = .true.
        if(found(10,key_char//' vtk output enkf',line,.FALSE.)) then
           write(unit = *, fmt = *) '[R] : vtk_out_enkf'
           read(unit = 10, fmt = *) vtk_out_enkf
           write(unit = 37, fmt = *) 'vtk_out_enkf ', vtk_out_enkf
        else
           write(unit = *, fmt = *) '[NR]: vtk_out_enkf'
           write(unit = 37, fmt = *) 'vtk_out_enkf true'
        end if

        ! # vtk output stddev
        vtk_out_stddev = .true.
        if(found(10,key_char//' vtk output stddev',line,.FALSE.)) then
           write(unit = *, fmt = *) '[R] : vtk_out_stddev'
           read(unit = 10, fmt = *) vtk_out_stddev
           write(unit = 37, fmt = *) 'vtk_out_stddev ', vtk_out_stddev
        else
           write(unit = *, fmt = *) '[NR]: vtk_out_stddev'
           write(unit = 37, fmt = *) 'vtk_out_stddev true'
        end if

        ! # vtk output resid
        vtk_out_resid = .true.
        if(found(10,key_char//' vtk output resid',line,.FALSE.)) then
           write(unit = *, fmt = *) '[R] : vtk_out_resid'
           read(unit = 10, fmt = *) vtk_out_resid
           write(unit = 37, fmt = *) 'vtk_out_resid ', vtk_out_resid
        else
           write(unit = *, fmt = *) '[NR]: vtk_out_resid'
           write(unit = 37, fmt = *) 'vtk_out_resid true'
        end if

        ! # vtk output covs
        vtk_out_covs = .true.
        if(found(10,key_char//' vtk output covs',line,.FALSE.)) then
           write(unit = *, fmt = *) '[R] : vtk_out_covs'
           read(unit = 10, fmt = *) vtk_out_covs
           write(unit = 37, fmt = *) 'vtk_out_covs ', vtk_out_covs
        else
           write(unit = *, fmt = *) '[NR]: vtk_out_covs'
           write(unit = 37, fmt = *) 'vtk_out_covs true'
        end if

        ! # vtk output realisations
        vtk_out_realisations = .false.
        if(found(10,key_char//' vtk output realisations',line,.FALSE.)) then
           write(unit = *, fmt = *) '[R] : vtk_out_realisations'
           write(unit = *, fmt = *) '[R] : num_out_realisations'
           read(unit = 10, fmt = *) vtk_out_realisations, num_out_realisations
           write(unit = 37, fmt = *) 'vtk_out_realisations ', vtk_out_realisations
           write(unit = 37, fmt = *) 'num_out_realisations ', num_out_realisations
        else
           write(unit = *, fmt = *) '[NR]: vtk_out_realisations'
           write(unit = 37, fmt = *) 'vtk_out_realisations false'
        end if

        ! # obs standard out
        obs_standard_out = .true.
        if(found(10,key_char//' obs standard out',line,.FALSE.)) then
           write(unit = *, fmt = *) '[R] : obs_standard_out'
           read(unit = 10, fmt = *) obs_standard_out
           write(unit = 37, fmt = *) 'obs_standard_out ', obs_standard_out
        else
           write(unit = *, fmt = *) '[NR]: obs_standard_out'
           write(unit = 37, fmt = *) 'obs_standard_out false'
        end if

        ! # enkf log output
        enkf_log_out = .true.
        if(found(10,key_char//' enkf log output',line,.FALSE.)) then
           write(unit = *, fmt = *) '[R] : enkf_log_out'
           read(unit = 10, fmt = *) enkf_log_out
           write(unit = 37, fmt = *) 'enkf_log_out ', enkf_log_out
        else
           write(unit = *, fmt = *) '[NR]: enkf_log_out'
           write(unit = 37, fmt = *) 'enkf_log_out true'
        end if

        ! # analysis matrix output
        ana_mat_out = .true.
        if(found(10,key_char//' analysis matrix output',line,.FALSE.)) then
           write(unit = *, fmt = *) '[R] : ana_mat_out'
           read(unit = 10, fmt = *) ana_mat_out
           write(unit = 37, fmt = *) 'ana_mat_out ', ana_mat_out
        else
           write(unit = *, fmt = *) '[NR]: ana_mat_out'
           write(unit = 37, fmt = *) 'ana_mat_out true'
        end if

        ! # compute realisation output
        comp_real_out = .true.
        if(found(10,key_char//' compute realisation output',line,.FALSE.)) then
           write(unit = *, fmt = *) '[R] : comp_real_out'
           read(unit = 10, fmt = *) comp_real_out
           write(unit = 37, fmt = *) 'comp_real_out ', comp_real_out
        else
           write(unit = *, fmt = *) '[NR]: comp_real_out'
           write(unit = 37, fmt = *) 'comp_real_out true'
        end if

        ! # iterativeenkf
        iterative_switch = .false.
        iterative_nrobs_int = nrobs_int+1
        iterative_doubleupdate = .false.
        if(found(10,key_char//' iterativeenkf',line,.FALSE.)) then
           write(unit = *, fmt = *) '[R] : iterativeenkf'
           read(unit = 10, fmt = *) iterative_switch, iterative_nrobs_int, iterative_doubleupdate
           write(unit = 37, fmt = *) 'iterative_switch ', iterative_switch
           write(unit = 37, fmt = *) 'iterative_nrobs_int ', iterative_nrobs_int
           write(unit = 37, fmt = *) 'iterative_doubleupdate ', iterative_doubleupdate
        else
           write(unit = *, fmt = *) '[NR]: iterative_switch'
           write(unit = 37, fmt = *) 'iterative_switch false'
        end if

        ! # prescribed velocity
        pres_vel_switch = .false.
        num_pres_vel = 0
        if(found(10,key_char//' prescribed velocity',line,.FALSE.)) then

           ! Get the number of records (num_pres_vel)
           call get_arg('records',line,i_beg_arg,i_end_arg)
           if( i_beg_arg < 1 .or. i_end_arg < i_beg_arg) then
              write(unit = *, fmt = *) '[E24] enkf_input.f90 (# prescribed velocity, records=...)'
              stop
           else
              read(line(i_beg_arg:i_end_arg),*) num_pres_vel
           end if
           write(unit = 37, fmt = *) 'num_pres_vel ', num_pres_vel

           if (num_pres_vel > 0) then
              write(unit = *, fmt = *) '[R] : prescribed velocity'
              read(unit = 10, fmt = *) pres_vel_switch, pres_vel_seed(1), pres_vel_seed(2)
              write(unit = 37, fmt = *) 'pres_vel_switch ', pres_vel_switch
              write(unit = 37, fmt = *) 'pres_vel_seed ', pres_vel_seed
              
              read(unit = 10, fmt= *) vdefault_stddevs(1), vdefault_stddevs(2), vdefault_stddevs(3)
              write(unit = 37, fmt = *) 'vdefault_stddevs ', vdefault_stddevs
              read(unit = 10, fmt= *) vdefault_sysvarmem(1), vdefault_sysvarmem(2), vdefault_sysvarmem(3)
              write(unit = 37, fmt = *) 'vdefault_sysvarmem ', vdefault_sysvarmem
           end if
           
        else
           
           write(unit = *, fmt = *) '[NR]: prescribed velocity'
           write(unit = 37, fmt = *) 'pres_vel_switch false'
           
        end if

        ! # temperature difference
        tempdiff_switch = .false.
        if(found(10,key_char//' temperature difference',line,.FALSE.)) then

           ! Get the number of records (num_tempdiff)
           call get_arg('records',line,i_beg_arg,i_end_arg)
           if( i_beg_arg < 1 .or. i_end_arg < i_beg_arg) then
              write(unit = *, fmt = *) '[E25] enkf_input.f90 (# temperature difference, records=...)'
              stop
           else
              read(line(i_beg_arg:i_end_arg),*) num_tempdiff
           end if
           write(unit = 37, fmt = *) 'num_tempdiff ', num_tempdiff

           if (num_tempdiff > 0) then
              write(unit = *, fmt = *) '[R] : temperature difference'
              read(unit = 10, fmt = *) tempdiff_switch
              write(unit = 37, fmt = *) 'tempdiff_switch ', tempdiff_switch

              
              allocate(tempdiff_inds(num_tempdiff,2))
              do itp = 1, num_tempdiff
                 read(unit = 10, fmt= *) j_dummy, k_dummy
                 tempdiff_inds(itp,1) = j_dummy
                 tempdiff_inds(itp,2) = k_dummy
              end do
              write(unit = 37, fmt = *) 'tempdiff_inds ', tempdiff_inds
           end if
           
        else
           
           write(unit = *, fmt = *) '[NR]: temperature difference'
           write(unit = 37, fmt = *) 'tempdiff_switch false'
           
        end if

        ! # tcon
        tcon_switch = .false.
        if(found(10,key_char//' tcon',line,.FALSE.)) then
           

           ! Get the number of records (num_tcon)
           call get_arg('records',line,i_beg_arg,i_end_arg)
           if( i_beg_arg < 1 .or. i_end_arg < i_beg_arg) then
              write(unit = *, fmt = *) '[E26] enkf_input.f90 (# tcon, records=...)'
              stop
           else
              read(line(i_beg_arg:i_end_arg),*) num_tcon
           end if
           write(unit = 37, fmt = *) 'num_tcon ', num_tcon

           if (num_tcon > 0) then
              write(unit = *, fmt = *) '[R] : thermal conductivity'
              read(unit = 10, fmt = *) tcon_switch, tcon_stddev, tcon_sysvarmem
              write(unit = 37, fmt = *) 'tcon_switch ', tcon_switch
              write(unit = 37, fmt = *) 'tcon_stddev ', tcon_stddev
              write(unit = 37, fmt = *) 'tcon_sysvarmem ', tcon_sysvarmem

              
              allocate(tcon_inds(num_tcon))
              do itc = 1, num_tcon
                 read(unit = 10, fmt= *) tcon_inds(itc)
              end do
              write(unit = 37, fmt = *) 'tcon_inds ', tcon_inds
           end if
           
        else
           
           write(unit = *, fmt = *) '[NR]: tcon'
           write(unit = 37, fmt = *) 'tcon_switch false'
           
        end if

        ! # pilot point
        pp_switch = .false.
        if(found(10,key_char//' pilot point',line,.FALSE.)) then

           ! Get the number of records (num_pp)
           call get_arg('records',line,i_beg_arg,i_end_arg)
           if( i_beg_arg < 1 .or. i_end_arg < i_beg_arg) then
              write(unit = *, fmt = *) '[E27] enkf_input.f90 (# pilot point, records=...)'
              stop
           else
              read(line(i_beg_arg:i_end_arg),*) num_pp
           end if
           write(unit = 37, fmt = *) 'num_pp ', num_pp

           if (num_pp > 0) then
              write(unit = *, fmt = *) '[R] : pilot points'
              read(unit = 10, fmt = *) pp_switch, pp_get_out, pp_ivar
              write(unit = 37, fmt = *) 'pp_switch ', pp_switch
              write(unit = 37, fmt = *) 'pp_get_out ', pp_get_out
              write(unit = 37, fmt = *) 'pp_ivar ', pp_ivar

              if(pp_switch) then
                 allocate(pp_inds(num_pp,3))
                 do ipp = 1, num_pp
                    read(unit = 10, fmt= *) i_pp, j_pp, k_pp
                    pp_inds(ipp,1) = i_pp
                    pp_inds(ipp,2) = j_pp
                    pp_inds(ipp,3) = k_pp
                 end do
                 write(unit = 37, fmt = *) 'pp_inds ', pp_inds
              end if

           end if

        else

           write(unit = *, fmt = *) '[NR]: pilot point'
           write(unit = 37, fmt = *) 'pp_switch false'
        end if

        ! # compute mean
        if(found(10,key_char//' compute mean',line,.FALSE.)) then
           write(unit = *, fmt = *) '[R] : compute mean'
           read(unit = 10, fmt = '(l1)') compute_mean
           write(unit = 37, fmt = *) 'compute_mean ', compute_mean
        else
           write(unit = *, fmt = *) '[E28] enkf_input.f90 (# compute mean) '
           stop 1
        end if

        close(10)

      END
