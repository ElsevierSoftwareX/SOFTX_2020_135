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

!>    @brief declaration of all main variables, arrays and constants
!>    @details
!> definition of global (dynamic) arrays, constants and main descriptions\n
      MODULE arrays
        IMPLICIT NONE

        !> @brief Numerical precision and min/max value.
        !> @details
        !>  Numerical precision and min/max value. \n\n
        !> (1) last digit precision \n
        !> (2) min allowed dble value \n
        !> (3) max allowed dble value \n
        double precision, dimension (3) :: const_dble

        !> @brief Project suffix, filename extension.
        !> @details
        !>  Project suffix, filename extension. \n\n
        !> project name extension/suffix\n
        !> Size: nsmpl.
        character (len=80), allocatable, dimension (:) :: project_sfx

        ! Size of grid
        ! ------------

        !> @brief Cell dimension array, x-direction.
        !> @details
        !> Cell dimension array, x-direction. \n\n
        !> delta size of cell, size i0.
        double precision, allocatable, dimension (:) :: delx

        !> @brief Cell dimension array, y-direction.
        !> @details
        !> Cell dimension array, y-direction. \n\n
        !> delta size of cell, size j0.
        double precision, allocatable, dimension (:) :: dely

        !> @brief Cell dimension array, z-direction.
        !> @details
        !> Cell dimension array, z-direction. \n\n
        !> delta size of cell, size k0.
        double precision, allocatable, dimension (:) :: delz

        !> @brief Absolute cell center positions, x-direction.
        !> @details
        !>  Absolute cell center positions, x-direction. \n\n
        !> absolute position, size  i0.\n\n
        !> With reference to the left, front, bottom corner of the
        !> model.
        double precision, allocatable, dimension (:) :: delxa

        !> @brief Absolute cell center positions, y-direction.
        !> @details
        !>  Absolute cell center positions, y-direction. \n\n
        !> absolute position, size  j0.\n\n
        !> With reference to the left, front, bottom corner of the
        !> model.
        double precision, allocatable, dimension (:) :: delya

        !> @brief Absolute cell center positions, z-direction.
        !> @details
        !>  Absolute cell center positions, z-direction. \n\n
        !> absolute position, size  k0.\n\n
        !> With reference to the left, front, bottom corner of the
        !> model.
        double precision, allocatable, dimension (:) :: delza

        ! Variable arrays
        ! ---------------

        !> @brief Variable array hydraulic potential.
        !> @details
        !>  Variable array hydraulic potential. \n\n
        !> Indices: \n
        !> 1. i-index \n
        !> 2. j-index \n
        !> 3. k-index \n
        !> 4. sample index \n
        double precision, allocatable, dimension (:,:,:,:) :: head

        !> @brief Variable array pressure.
        !> @details
        !>  Variable array pressure. \n\n
        !> Indices: \n
        !> 1. i-index \n
        !> 2. j-index \n
        !> 3. k-index \n
        !> 4. sample index \n
        double precision, allocatable, dimension (:,:,:,:) :: pres

        !> @brief Variable array temperature.
        !> @details
        !>  Variable array temperature. \n\n
        !> Indices: \n
        !> 1. i-index \n
        !> 2. j-index \n
        !> 3. k-index \n
        !> 4. sample index \n
        double precision, allocatable, dimension (:,:,:,:) :: temp

        !> @brief Variable array concentrations.
        !> @details
        !>  Variable array concentrations. \n\n
        !> Indices: \n
        !> 1. i-index \n
        !> 2. j-index \n
        !> 3. k-index \n
        !> 4. species index \n
        !> 5. sample index \n
        double precision, allocatable, dimension (:,:,:,:,:) :: conc

        !> @brief Array of concentration iteration differences.
        !> @details
        !>  Array of concentration iteration differences. \n\n
        !> Used for checking concentration convergence.
        double precision, allocatable, dimension (:,:) :: conc_conv

        !> @brief Array of total salinities.
        !> @details
        !>  Array of total salinities. \n\n
        !> Indices: \n
        !> 1. i-index \n
        !> 2. j-index \n
        !> 3. k-index \n
        !> 4. sample index \n
        double precision, allocatable, dimension (:,:,:,:) :: tsal

        ! Variables for variable time stepping
        ! ------------------------------------

        !> @brief Array with flag for nonlinear iteration maxout.
        !> @details
        !>  Array with flag for nonlinear iteration maxout. \n\n
        !>
        !> - "0": Initialized value \n
        !> - "1": Nonlinear iteration fine, possibly double time
        !>        step. \n
        !> - "-2": Iterative solver or Picard/Newton iteration reached
        !>         maxiter. Half time step. \n
        !>
        !> Variable for variable time stepping.
        integer, allocatable, dimension (:) :: flag_delt

        !> @brief Array with counters for doubling time step.
        !> @details
        !>  Array with counters for doubling time step. \n
        !>
        !> Variable for variable time stepping.
        integer, allocatable, dimension (:) :: delt_count

        !> @brief Array containing information if it is the first timestep.
        !> @details
        !>  Array containing information if it is the first timestep. \n
        !> - 0: First time step \n
        !> - 1: Not the first time step. \n\n
        !>
        !> The information is carried for each sample, the length will
        !> be nsmpl. \n\n
        !>
        !> Variable for variable time stepping.
        integer, allocatable :: flag_1st_timestep(:)

        !> @brief Array of previous time step lengths for variable time step.
        !> @details
        !>  Array of previous time step lengths for variable time step. \n\n
        !>
        !> Variable for variable time stepping.
        double precision, allocatable, dimension (:) :: delt_old

!     maximum of the used units/BC-units ("maxunits"/"bc_maxunits")
        INTEGER nunits
!    number of rock properties + boundary condition types, number of regular rock properties (to load)
        INTEGER nprop, nprop_load
!     first/last index of normal properties in 'nprop'
        INTEGER maxunits, firstidx, lastidx
!     first/last index of bc-units in 'nprop'
        INTEGER bc_maxunits, bc_firstidx, bc_lastidx, nbc

!       string constants for the number of property-units and bc-units
        character (len=2) :: c_npropunit, c_nbcunit, c_npv

        PARAMETER (firstidx=1) ! first rock properties
        PARAMETER (lastidx=17) ! last rock properties
        PARAMETER (bc_firstidx=18) ! first BC unit
        PARAMETER (bc_lastidx=22) ! last BC unit
        PARAMETER (nbc=bc_lastidx-bc_firstidx+1)
        PARAMETER (nprop=lastidx-firstidx+1 +bc_lastidx-bc_firstidx+1)
!       load fewer, except 3*bc-units
        PARAMETER (nprop_load=nprop-nbc)

!       unit index number, unit-cell assignment (rock property for each cell)
        INTEGER, ALLOCATABLE :: uindex(:,:,:)
!       cell index number, no assignment - only for output (grouping)
        INTEGER, ALLOCATABLE :: cindex(:,:,:)

!        def_props = '<name>'
        DOUBLE PRECISION, ALLOCATABLE :: propunit(:,:,:)

!     disable additional hdf5-output
        INTEGER nout_ijk
        PARAMETER (nout_ijk=9)
        LOGICAL out_ijk(nout_ijk)
        INTEGER cout_i, cout_j, cout_k
        PARAMETER (cout_i=1)
        PARAMETER (cout_j=2)
        PARAMETER (cout_k=3)
        INTEGER cout_vx, cout_vy, cout_vz
        PARAMETER (cout_vx=4)
        PARAMETER (cout_vy=5)
        PARAMETER (cout_vz=6)
        INTEGER cout_rhof, cout_visf, cout_uindex
        PARAMETER (cout_rhof=7)
        PARAMETER (cout_visf=8)
        PARAMETER (cout_uindex=9)

        INTEGER ncompress
        PARAMETER (ncompress=4)
        character (len=5) :: compress_suffix(ncompress)
        DATA compress_suffix/'plain', 'bz2', 'gz', 'zip'/
        LOGICAL out_prop(nprop)
        character (len=4) :: properties(nprop+2)
        DATA properties/' por', 'a_kx', 'a_ky', '  kz', 'comp', &
          'a_lx', 'a_ly', '  lz', '   q', '  rc', '  df', '  ec', &
          '  lc', ' bcl', 'bcpd', 's_nr', 's_wr', &
          ' hbc', ' tbc', ' cbc', ' ebc', 'snbc', '  tp', ' sbc'/
        character (len=80) :: doc_properties(nprop+2)
        DATA doc_properties/'porosity', &
          'permeability, anisotropic direction X (ratio of Z), default 1.0', &
          'permeability, anisotropic direction Y (ratio of Z), default 1.0', &
          'permeability, direction Z', &
          'compressibility of rock, default 1.0e-10', &
          'conductivity, anisotropic direction X (ratio of Z), default 1.0', &
          'conductivity, anisotropic direction Y (ratio of Z), default 1.0', &
          'conductivity, direction Z', &
          'heat production, default 0.0', &
          'heat capacity of rock, default 2.06e6', &
          'diffusivity, default 10.0', &
          'electrical conductivity, default 0.0', &
          'coupling coefficient, default 0.0', &
          'Brooks Corey "lambda" - pore size distribution index, default 2.0', &
          'Brooks Corey displacement pressure (PD), default 1.0d6', &
          'residual saturation of the non-wetting phase, default 0.05', &
          'residual saturation of the wetting phase, default 0.2', &
          'flow boundary condition', &
          'temperature boundary condition', &
          'concentration boundary condition', &
          'electric potential boundary condition', &
          'saturation (non-wetting) boundary condition', &
          'time periods', &
          'unspecified single cell boundary condition'/
        DOUBLE PRECISION prop_max(nprop_load)
        DATA prop_max/1.0D0, 1.0D+30, 1.0D+30, 1.0D+30, 1.0D+50, &
          1.0D+30, 1.0D+30, 1.0D+30, 1.0D+30, 1.0D+30, 1.0D+30, 1.0D+30, &
          1.0D+50, 10.0d0, 1.0D+30, 1.0D0, 1.0D0/
        DOUBLE PRECISION prop_min(nprop_load)
        DATA prop_min/1.0D-30, 1.0D-3, 1.0D-3, 1.0D-30, 1.0D-50, &
          1.0D-3, 1.0D-3, 1.0D-3, 1.0D-20, 1.0D-10, 0.0D0, 1.0D-10, &
          1.0D-50, 1.0D-2, 0.0D0, 0.0D0, 0.0D0/
        DOUBLE PRECISION prop_default(nprop_load)
        DATA prop_default/1.0D-7, 1.0D0, 1.0D0, 1.0D-20, 1.0D-10, &
          1.0D0, 1.0D0, 2.0D0, 0.0D0, 2.0D6, 10.0D0, 0.0D0, &
          0.0D0, 2.0D0, 1.0D3, 0.05D0, 0.2D0/

! ------------------------------------------
        INTEGER idx_por
        INTEGER idx_an_kx
        INTEGER idx_an_ky
        INTEGER idx_kz
        INTEGER idx_comp
        INTEGER idx_an_lx
        INTEGER idx_an_ly
        INTEGER idx_lz
        INTEGER idx_q
        INTEGER idx_rc
        INTEGER idx_df
        INTEGER idx_ec
        INTEGER idx_lc


        INTEGER idx_s_nr
        INTEGER idx_s_wr
        INTEGER idx_hbc
        INTEGER idx_tbc
        INTEGER idx_cbc
        INTEGER idx_ebc
        INTEGER idx_snbc
!     integer idx_pbc

        PARAMETER (idx_por=1)
        PARAMETER (idx_an_kx=2)
        PARAMETER (idx_an_ky=3)
        PARAMETER (idx_kz=4)
        PARAMETER (idx_comp=5)
        PARAMETER (idx_an_lx=6)
        PARAMETER (idx_an_ly=7)
        PARAMETER (idx_lz=8)
        PARAMETER (idx_q=9)
        PARAMETER (idx_rc=10)
        PARAMETER (idx_df=11)
        PARAMETER (idx_ec=12)
        PARAMETER (idx_lc=13)


        PARAMETER (idx_s_nr=16)
        PARAMETER (idx_s_wr=17)
        PARAMETER (idx_hbc=18)
        PARAMETER (idx_tbc=19)
        PARAMETER (idx_cbc=20)
        PARAMETER (idx_ebc=21)
        PARAMETER (idx_snbc=22)
!       parameter (idx_pbc = 18)?? like idx_hbc
! ------------------------------------------
        INTEGER idx_tp
        INTEGER idx_sbc
        PARAMETER (idx_tp=nprop+1)
        PARAMETER (idx_sbc=nprop+2)
! ------------------------------------------

! new boundary-condition structures
!     i,j,k - position
        INTEGER cbc_i, cbc_j, cbc_k
        PARAMETER (cbc_i=1)
        PARAMETER (cbc_j=2)
        PARAMETER (cbc_k=3)
!     bc-unit, bc time dependend
        INTEGER cbc_bcu, cbc_bctp
        PARAMETER (cbc_bcu=4)
        PARAMETER (cbc_bctp=5)
!     physical value, boundary type (neuman, dirichlet), sub index (species)
        INTEGER cbc_pv, cbc_bt, cbc_si, cbc_dir
        PARAMETER (cbc_pv=6)
        PARAMETER (cbc_bt=7)
        PARAMETER (cbc_si=8)
        PARAMETER (cbc_dir=9)
!     bc-type max index
        INTEGER nibc, ndbc
        PARAMETER (nibc=9)
!     1:value, 2:bcmy, 3:additional value (e.g. for well function pressure)
        PARAMETER (ndbc=3)

        ! Physical value (pv) indices and names
        ! -------------------------------------

        !> @brief Physical value index: head.
        !> @details
        !>  Physical value index: head. \n\n
        integer, parameter :: pv_head = 1

        !> @brief Physical value index: temp.
        !> @details
        !>  Physical value index: temp. \n\n
        integer, parameter :: pv_temp = 2

        !> @brief Physical value index: conc.
        !> @details
        !>  Physical value index: conc. \n\n
        integer, parameter :: pv_conc = 3

        !> @brief Physical value index: pres.
        !> @details
        !>  Physical value index: pres. \n\n
        integer, parameter :: pv_pres = 5

        !> @brief Physical value index: bhpr.
        !> @details
        !>  Physical value index: bhpr. \n\n
        !> Special case for borehole output.
        integer, parameter :: pv_bhpr = 6

        !> @brief Number of physical value indices.
        !> @details
        !>  Number of physical value indices. \n\n
        !> Or: pv max index.
        integer, parameter :: npv = 6

        !> @brief Array of physical value names.
        !> @details
        !> Array of physical value names. \n\n
        character (len=4), parameter, dimension (npv) :: pv_name = &
            (/'head', 'temp', 'conc', 'nova', 'pres', 'bhpr'/)

        !> @brief Array of physical value output switches.
        !> @details
        !> Array of physical output switches. \n\n
        !> Switches are set according to activity of variables. Output
        !> (vtk/hdf5) for specific variables can be suppressed by
        !> specifying the physical value name in `# disable output`.
        logical out_pv(npv)

!     bc-type
        INTEGER bt_diri, bt_neum, bt_neuw
        PARAMETER (bt_diri=1)
        PARAMETER (bt_neum=2)
        PARAMETER (bt_neuw=3)
        character (len=3) :: bc_name(3)
        DATA bc_name/'bcd', 'bcn', 'bcw'/
!     boundary-condition structures
        INTEGER, ALLOCATABLE :: ibc_data(:,:)
        DOUBLE PRECISION, ALLOCATABLE :: dbc_data(:,:,:)
        DOUBLE PRECISION, ALLOCATABLE :: dbc_dataold(:)

!       borehole logs
        INTEGER nbh_logs
        INTEGER, ALLOCATABLE :: ibh_pos(:,:)
        character (len=256), dimension (:), allocatable :: cbh_name

!     Begin and End of the pv-blocks after sorting
!      reuse *_flow for head and pres !!!
        INTEGER first_flow, last_flow
        INTEGER first_temp, last_temp
        INTEGER first_conc, last_conc
        INTEGER nbc_data

!     SM/simulate
        DOUBLE PRECISION, ALLOCATABLE :: propunitold(:,:)
        DOUBLE PRECISION, ALLOCATABLE :: bcperiodold(:,:,:)

!     -------
!     inverse (AD) & simulate (SIMUL/ENKF)
        INTEGER mpara

!     -------
!     inverse
        INTEGER, ALLOCATABLE :: opti_props(:,:)
        INTEGER, ALLOCATABLE :: opti_bc(:,:)
        INTEGER, ALLOCATABLE :: opti_tp(:,:)
!       1: alpha; 2: beta modificators
        INTEGER mopti_tp
        PARAMETER (mopti_tp=2)
!
!     apriori (inverse)
        DOUBLE PRECISION, ALLOCATABLE :: a_propunit(:,:)

!     error (apriori)
        DOUBLE PRECISION, ALLOCATABLE :: d_propunit(:,:)

!     error (aposteriori)
        DOUBLE PRECISION, ALLOCATABLE :: e_propunit(:,:)

!     main parameter input vector for forward simulation (see "forward_compute")
        DOUBLE PRECISION, ALLOCATABLE :: main_input(:,:)
!     main data output vector for forward simulation (see "forward_compute")
        DOUBLE PRECISION, ALLOCATABLE :: main_output(:,:)
!
!     for jacobi computation (INVERSE) / ENKF /SIMUL
!       parameter master copy
        DOUBLE PRECISION, ALLOCATABLE :: main_input_master(:)

        ! Observed data
        ! -------------

        !> @brief Array of integer data specifications.
        !> @details
        !> Array of integer data specifications. \n
        !> Index 2: [i,j,k, type, sub-index, obs]
        integer, allocatable, dimension (:,:) :: idata

        !> @brief Array of double precision data specifications.
        !> @details
        !> Array of double precision data specifications. \n
        !> Index2: [value, weighting, time, px,py,pz]
        double precision, allocatable, dimension (:,:) :: ddata

        !> @brief Array of simulated data values.
        !> @details
        !> Array of simulated data values. \n
        !> Save the computed values to compare it with
        !> 'ddata(:,cid_pv)' \n
        !> Index 2: Sample index.
        double precision, allocatable, dimension (:,:) :: sdata

        !> @brief Index position in idata for i.
        !> @details
        !> Index position in idata for i. \n
        integer, parameter :: cid_i = 1

        !> @brief Index position in idata for j.
        !> @details
        !> Index position in idata for j. \n
        integer, parameter :: cid_j = 2

        !> @brief Index position in idata for k.
        !> @details
        !> Index position in idata for k. \n
        integer, parameter :: cid_k = 3

        !> @brief Index position in idata for type.
        !> @details
        !> Index position in idata for type. \n
        integer, parameter :: cid_pv = 4

        !> @brief Index position in idata for sub-index.
        !> @details
        !> Index position in idata for sub-index. \n
        integer, parameter :: cid_si = 5

        !> @brief Index position in idata for obs.
        !> @details
        !> Index position in idata for obs. \n
        integer, parameter :: cid_obs = 6

        !> @brief Index-2 dimension of idata.
        !> @details
        !> Index-2 dimension of idata. \n
        !> Number of different integer-parameters arrays in idata.
        integer, parameter :: n_idata = 6

        !> @brief Index position in ddata for value.
        !> @details
        !> Index position in ddata for value. \n
        integer, parameter :: cdd_pv = 1

        !> @brief Index position in ddata for weighting.
        !> @details
        !> Index position in ddata for weighting. \n
        integer, parameter :: cdd_w = 2

        !> @brief Index position in ddata for time.
        !> @details
        !> Index position in ddata for time. \n
        integer, parameter :: cdd_time = 3

        !> @brief Index position in ddata for px.
        !> @details
        !> Index position in ddata for px. \n
        integer, parameter :: cdd_i = 4

        !> @brief Index position in ddata for py.
        !> @details
        !> Index position in ddata for py. \n
        integer, parameter :: cdd_j = 5

        !> @brief Index position in ddata for pz.
        !> @details
        !> Index position in ddata for pz. \n
        integer, parameter :: cdd_k = 6

        !> @brief Index-2 dimension of ddata.
        !> @details
        !> Index-2 dimension of ddata. \n
        !> Number of different double-precision-parameters arrays in
        !> ddata.
        integer, parameter :: n_ddata = 6

!     jump table between parameter index and seeding
        INTEGER, ALLOCATABLE :: seed_para(:,:)
        INTEGER, ALLOCATABLE :: gpara(:)

        ! Coefficients for linear system solver
        ! -------------------------------------

        !> @brief Linear system solver coefficent array a.
        !> @details
        !> Linear system solver coefficent array a. \n
        DOUBLE PRECISION, ALLOCATABLE :: a(:,:,:,:)

        !> @brief Linear system solver coefficent array b.
        !> @details
        !> Linear system solver coefficent array b. \n
        DOUBLE PRECISION, ALLOCATABLE :: b(:,:,:,:)

        !> @brief Linear system solver coefficent array c.
        !> @details
        !> Linear system solver coefficent array c. \n
        DOUBLE PRECISION, ALLOCATABLE :: c(:,:,:,:)

        !> @brief Linear system solver coefficent array d.
        !> @details
        !> Linear system solver coefficent array d. \n
        DOUBLE PRECISION, ALLOCATABLE :: d(:,:,:,:)

        !> @brief Linear system solver coefficent array e.
        !> @details
        !> Linear system solver coefficent array e. \n
        DOUBLE PRECISION, ALLOCATABLE :: e(:,:,:,:)

        !> @brief Linear system solver coefficent array f.
        !> @details
        !> Linear system solver coefficent array f. \n
        DOUBLE PRECISION, ALLOCATABLE :: f(:,:,:,:)

        !> @brief Linear system solver coefficent array g.
        !> @details
        !> Linear system solver coefficent array g. \n
        DOUBLE PRECISION, ALLOCATABLE :: g(:,:,:,:)

        !> @brief Linear system solver coefficent array w.
        !> @details
        !> Linear system solver coefficent array w. \n
        DOUBLE PRECISION, ALLOCATABLE :: w(:,:,:,:)

        !> @brief Linear system solver coefficent array x.
        !> @details
        !> Linear system solver coefficent array x. \n
        DOUBLE PRECISION, ALLOCATABLE :: x(:,:,:,:)

!     [r0] random number vector, used for BiCGStab algorithm
        DOUBLE PRECISION, ALLOCATABLE :: r(:,:,:)
!     openmp-private/locale vectors for linear system solver (lss_*)
!       team-global buffer for boundary exchange (+ismpl)
        DOUBLE PRECISION, ALLOCATABLE :: lss_bound_block(:,:,:,:,:,:)
        DOUBLE PRECISION, ALLOCATABLE :: lss_dnrm(:,:)
        DOUBLE PRECISION, ALLOCATABLE :: lss_tmp(:,:)
!       private copy for preconditioning (+Tlevel_1 +ismpl)
        DOUBLE PRECISION, ALLOCATABLE :: lss_lma(:,:,:), lss_lmb(:,:,:)
        DOUBLE PRECISION, ALLOCATABLE :: lss_lmc(:,:,:), lss_lmd(:,:,:)
        DOUBLE PRECISION, ALLOCATABLE :: lss_lme(:,:,:), lss_lmf(:,:,:)
        DOUBLE PRECISION, ALLOCATABLE :: lss_lmg(:,:,:), lss_lud(:,:,:)
        DOUBLE PRECISION, ALLOCATABLE :: lss_lx(:,:,:), lss_lb(:,:,:)
        DOUBLE PRECISION, ALLOCATABLE :: lss_lloctmp(:,:,:,:)
        DOUBLE PRECISION, ALLOCATABLE :: lss_ldnrm(:,:,:)
        DOUBLE PRECISION, ALLOCATABLE :: lss_ud_block(:,:,:,:)
        DOUBLE PRECISION, ALLOCATABLE :: lss_lr0_hat(:,:,:)

!     for ilu-precond. (shadow vectors)
        DOUBLE PRECISION, ALLOCATABLE :: ud(:,:,:,:)
!     marker of current boundary conditions elements
        CHARACTER, ALLOCATABLE :: bc_mask(:,:)

!     boundary condition info, values 'n','d',' ' on positions [head,temp,conc,pres]
        character (len=npv), dimension (:,:,:), allocatable :: node_info

!     spez. ILU arrays
!       proc. index of block
        INTEGER, ALLOCATABLE :: proza(:,:,:)

        !> @brief Array for storing old head for iteration.
        !> @details
        !>  Array for storing old head for iteration. \n\n
        !> Indices: \n
        !> 1. linear cell-index \n
        !> 2. cgen-level index \n
        !> 3. sample index\n
        !> Size: [I0*J0*K0,ncgen,nsmpl], used for forward newton iteration
        double precision, allocatable, dimension (:,:,:) :: headold

        !> @brief Array for storing old temp for iteration.
        !> @details
        !>  Array for storing old temp for iteration. \n\n
        !> Indices: \n
        !> 1. linear cell-index \n
        !> 2. cgen-level index \n
        !> 3. sample index\n
        !> Size: [I0*J0*K0,ncgen,nsmpl], used for forward newton iteration
        double precision, allocatable, dimension (:,:,:) :: tempold

        !> @brief Array for storing old pressure for iteration.
        !> @details
        !>  Array for storing old pressure for iteration. \n\n
        !> Indices: \n
        !> 1. linear cell-index \n
        !> 2. cgen-level index \n
        !> 3. sample index\n
        !> Size: [I0*J0*K0,ncgen,nsmpl], used for forward newton iteration
        double precision, allocatable, dimension (:,:,:) :: presold

        !> @brief Array for storing old conc for iteration.
        !> @details
        !>  Array for storing old conc for iteration. \n\n
        !> Indices: \n
        !> 1. linear cell-index \n
        !> 2. species index \n
        !> 3. cgen-level index \n
        !> 4. sample index\n
        !> Size: [I0*J0*K0,max(ntrans,1),ncgen,nsmpl], used for forward newton iteration
        double precision, allocatable, dimension (:,:,:,:) :: concold

!     BC time period: (period-index,value-type,TP-ID,sample)
!     - value-type: time, BC-value
        DOUBLE PRECISION, ALLOCATABLE :: bcperiod(:,:,:,:)
!     BC time period - number of periods: (TP-ID)
        INTEGER, ALLOCATABLE :: ibcperiod(:)
!     BC time period - on/off-switch: (period-index,TP-ID)
        LOGICAL, ALLOCATABLE :: lbcperiod(:,:)

        !> @brief Array of output times.
        !> @details
        !> Array of output times. \n
        !> Define output times (time dependend) \n
        !> Read under `# output times`. \n
        double precision, allocatable, dimension (:) :: outt

!     define monitor index, used for "simulate"
        INTEGER, ALLOCATABLE :: smon_idx(:)

        !> @brief Array of time step durations/lengths/values.
        !> @details
        !> Array of time step durations/lengths/values. \n
        !> The time step lengths are computed in `calc_deltatime`
        !> according to the input under `# time periods`: \n
        !> - start and end times of the time periods \n
        !> - number of time steps per period \n
        !> - step type (f.e. distributed linearly, logarithmically) \n
        double precision, allocatable, dimension (:) :: delta_time

!     for inverse computation
        DOUBLE PRECISION, ALLOCATABLE :: a_bcperiod(:,:,:)
        DOUBLE PRECISION, ALLOCATABLE :: d_bcperiod(:,:,:)
        DOUBLE PRECISION, ALLOCATABLE :: e_bcperiod(:,:,:)

        !> @brief Simulation time
        !> @details
        !> Simulation time. \n
        !> Simulation time for each sample are set in
        !> `forward/forward_iter.f90`.
        double precision, allocatable, dimension (:) :: simtime

        !> @brief Samples array for transient switches
        !> @details
        !> Samples array for transient switches. \n\n
        !>
        !> transient execption, to toggle it off
        logical, allocatable, dimension (:) :: tr_switch

        INTEGER c_fhandler, c_foffset
!     [max # of files per thread]
        PARAMETER (c_fhandler=100)
!     file-handler offset
        PARAMETER (c_foffset=30)
!     file handler table
        INTEGER, ALLOCATABLE :: fh_table(:,:)

!     - convergency history buffer -
!     max history length
        INTEGER conv_hlen
        PARAMETER (conv_hlen=50)
!     number of histories
        INTEGER conv_hmax
!     history buffer
        DOUBLE PRECISION, ALLOCATABLE :: conv_history(:,:,:)
!     current history length
        INTEGER, ALLOCATABLE :: conv_chlen(:,:)
!     current position index
        INTEGER, ALLOCATABLE :: conv_ipos(:,:)
!     convergency, when this vector is true
        LOGICAL, ALLOCATABLE :: lcon(:,:)

!     default flow values from read_model.f90
        DOUBLE PRECISION, ALLOCATABLE :: vdefault(:,:)
        LOGICAL :: vdefaultswitch

!     transport
        DOUBLE PRECISION, ALLOCATABLE :: mmas_c(:)
        DOUBLE PRECISION, ALLOCATABLE :: diff_c(:)
        DOUBLE PRECISION, ALLOCATABLE :: beta_c(:)


!     OpenMP REDUCTION arrays
        DOUBLE PRECISION, ALLOCATABLE :: omp_dglobal(:,:,:)
        INTEGER, ALLOCATABLE :: omp_iglobal(:,:,:)

!---- only for DEBUG !!! ----
        INTEGER n_debugout
!     [2,n_debugout]: 1: time step, 2: nl iteration
        INTEGER, ALLOCATABLE :: debugout(:,:)

      END MODULE arrays
