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

!> @brief Module for EnKF-related variables
!> @details
!> Module containing all EnKF-related variables.
module mod_enkf

  !> @brief Number of EnKF variables
  !> @details
  !> Number of EnKF variables \n
  !> The number of variables which could be part
  !> of the EnKF state vector. \n
  !> head, temp, conc, kz, lz, por
  integer, parameter :: num_enkf_vars = 6


  !---------------------------------------------------
  !---------------------------------------------------
  !               ENKF_ITER ARRAYS
  !---------------------------------------------------
  !---------------------------------------------------

  !> @brief Variables of the synthetic true model.
  !> @details
  !> Variables of the synthetic true model. \n
  !> Arrays that will contain the variables of the
  !> synthetic true model. \n
  !> True variable values.
  double precision, allocatable, dimension(:,:) :: var_true

  !> @brief Observation times
  !> @details
  !> Observation times \n
  !> Observation times array.\n
  !> obst(i): Simulation time i \n
  !> Vector containing observation times.
  double precision, allocatable, dimension(:) :: obst

  !> @brief Number of observation locations
  !> @details
  !> Number of observation locations \n
  !> Number of observation locations for each
  !> observation time. \n
  !> nrobs_loc(i): Number of observation grid-cells at obs-time i \n
  !> Contains the number of measurement-points at each
  !> measurement-time.
  integer, allocatable, dimension(:) :: nrobs_loc

  !> @brief Observation activity
  !> @details
  !> Observation activity \n
  !> Array containing the information for each state
  !> vector variable if it is active in a observation.
  !> This is read from the observation file. \n
  !> act_o(iob,i): Signifies if
  !> variable/rock-property i is active at a
  !> certain measurement time \n
  !> Activity of h,t,c,kz,lz,por in measurement
  integer, allocatable, dimension(:,:) :: act_o

  !> @brief i-indices of observations
  !> @details
  !> i-indices of observations \n
  !> Array of i-indices in x-direction of the
  !> observation locations. \n
  !> i_obs(irobs,iloc): i-index of observation iloc
  !> at obs-time irobs
  integer, allocatable, dimension(:,:) :: i_obs
  !> @brief j-indices of observations
  !> @details
  !> j-indices of observations \n
  !> Array of j-indices in y-direction of the
  !> observation locations. \n
  !> j_obs(i,j): j-index of observation j at obs-time i
  integer, allocatable, dimension(:,:) :: j_obs
  !> @brief k-indices of observations
  !> @details
  !> k-indices of observations \n
  !> Array of k-indices in z-direction of the
  !> observation locations. \n
  !> k_obs(i,j): k-index of observation j at obs-time i
  integer, allocatable, dimension(:,:) :: k_obs

  !> @brief Measurement values at observations.
  !> @details
  !> Measurement values at observations. \n
  !> Array of measurement values at certain observation
  !> time, location and for a certain observation variable. \n
  !> var_obs(irobs,iloc,ienkfvar): value of var ienkfvar
  !> for observation point iloc at observation number irobs.
  double precision, allocatable, dimension(:,:,:) :: var_obs

  !> @brief Mean standard deviations at all update steps
  !> @details
  !> Mean standard deviations at all update steps \n
  !> Arrays containing the mean ensemble standard
  !> deviations at all update steps and for all variables. \n
  !> Array holding integrated stddevs for all
  !> obstimes/variables. \n
  !>  The mean variance across the model locations is calculated and
  !>  then the square root is taken. \n
  double precision, allocatable, dimension(:,:) :: stdvar

  !> @brief Root Mean Squared Error
  !> @details
  !> Root Mean Squared Error \n
  !> Array containing the RMSE at different times and for
  !> different variables. \n
  !> Array holding residuals for all measurement times and variables.
  double precision, allocatable, dimension(:,:) :: resvar

  !> @brief State vector index
  !> @details
  !> State vector index \n
  !> Linear index for the cells contained in the
  !> state vector. \n
  !> 1D-array, length number of grid cells \n
  !> contains information if LOCATION is inside state vector \n
  !> array of linear node-indices counted without excluded nodes. \n
  !> Linear index for nodes without non-active nodes. \n
  !> sysindx(l) = l if l in state vector, otherwise sysindx(l) = 0. \n
  integer, allocatable, dimension(:) :: sysindx

  !---------------------------------------------------------
  !---------------------------------------------------------
  !               EVENSEN ENKF ARRAYS
  !---------------------------------------------------------
  !--------------------------------------------------------  

  !> @brief EnKF state vector data matrix
  !> @details
  !> EnKF state vector data matrix \n
  !> The data matrix of realizations with a state vector index and an
  !> ensemble index. \n
  !> The matrix containing the ensemble members (realizations) as columns.
  double precision, allocatable :: mem(:,:)

  !> @brief EnKF Mean vector
  !> @details
  !> EnKF Mean vector \n
  !> Sample mean of the data matrix mem. \n
  !> Vector containing the ensemble averages for the state vector.
  double precision, allocatable :: ave(:)

  !> @brief EnKF Variance vector
  !> @details
  !> EnKF Variance vector \n
  !> Sample variance of the data matrix
  !> mem. \n
  !> Vector containing the ensemble variances for the state vector. \n
  double precision, allocatable :: var(:)

  !> @brief Correlation matrix for cov_ref output
  !> @details
  !> Correlation matrix for cov_ref output \n
  !> Correlation matrices with respect to one point.
  double precision, allocatable :: corr_matrix(:,:)

  !> @brief EnKF system variances
  !> @details
  !> EnKF system variances \n
  !> The system variances from sysvar, where they are associated to
  !> variables, are associated to the state vector variables. Thus
  !> there will be many equal system variances.
  double precision, allocatable :: sysvarmem(:)

  !---------------------------------------------------------
  !---------------------------------------------------------
  !               ENKF INPUT FILE
  !---------------------------------------------------------
  !--------------------------------------------------------  

  !> @brief EnKF Error covariance model
  !> @details
  !> EnKF Error covariance model \n
  !> Possible values: 'gaussian' or 'diagonal' \n
  !> diagonal: no covariances (only variances given below) \n
  !> gaussian: covariances get smaller according distances and
  !> correlation length \n
  character (len=100) :: covmodel

  !> @brief Error covariance matrix approximation switch
  !> @details
  !> Error covariance matrix approximation switch \n
  !> If .true., calculate the exact error covariance
  !> matrix R, if .false., calculate lowrank matrix R.
  !> See: m_enkf_enkf.f90 and analysis.f90 \n
  !> Exact (true) or lowrank (false) R matrix
  logical :: rexact

  !> @brief Observation random noise correction
  !> @details
  !> Observation random noise correction \n
  !> If .true. any random noise is corrected
  !> for mean=0. \n
  !> Random noise corrected for zero-mean and unit variance
  logical :: fixsamp

  !> @brief Mode of the EnKF-analysis
  !> @details
  !> Mode of the EnKF-analysis \n
  !> Possible cases: 11, 21, 12, 22, 13, 23.
  !> The first number determines how the matrix X5
  !> is generated during the analysis: EnKF (1) or SQRT (2). \n
  !> The Second number determines how the inversion
  !> of the matrix R in the EnKF is carried
  !> out (analysis.f90). \n
  !> Choice of EnKF/SquareRoot Filters in Evensen Code. \n
  !> Tested input: 11. \n
  integer :: mode_analysis

  !> @brief Truncation parameter for analysis functions
  !> @details
  !> Truncation parameter for analysis functions \n
  !> Used for eigenvalue decomposition and singular value
  !> decomposition.  See mod_anafunc.f90 \n
  !> Numerical eigenvalue cut-off for pseudo-inversion in
  !> Evensen-Code.
  double precision :: truncation

  !> @brief True Model Check Switch
  !> @details
  !> True Model Check Switch \n
  !> If .true., SHEMAT-Suite reads data
  !> from a true model in the Tecplot file
  !> specified in true_name. \n
  !> Want to compare simulations to a true model
  !> (TrueModelPROJECT.plt)
  logical :: checktrue

  !> @brief Error covariance horizontal correlation length
  !> @details
  !> Error covariance horizontal correlation length \n
  !> Only needed, when covmodel='gaussian' \n
  !> Physical units.
  double precision :: clh

  !> @brief Error covariance vertical correlation length
  !> @details
  !> Error covariance vertical correlation length \n
  !> Only needed, when covmodel='gaussian' \n
  !> Physical units.
  double precision :: clv

  !> @brief Array of observation variances
  !> @details
  !> Array of observation variances \n
  !> One observation variance for each variable
  !> in enkf_variable_names. \n
  !> Variances for the observations (which are perturbed by gaussian
  !> noise)
  double precision, dimension(num_enkf_vars) :: obsvar

  !> @brief Array of system variances
  !> @details
  !> Array of system variances \n
  !> One system variance for each variable
  !> in enkf_variable_names. \n
  !> System variances (variables active in the state vector are
  !> perturbed by gaussian noise before each assimilation)
  double precision, dimension(num_enkf_vars) :: sysvar

  !> @brief Number of observations times
  !> @details
  !> Number of observations times \n
  !> Alternatively: Number of observation intervals. \n
  !> Number of observation times \n
  !> Must be compatible with observation input in the observation file
  !> and with the time discretization in the general input file.
  integer :: nrobs_int

  !> @brief Array of activity in state vector
  !> @details
  !> Array of activity in state vector \n
  !> For all variables in enkf_variable_names.
  !> 1: active, 0: not active. \n
  !> Activity of h,t,c,kz,lz,por in simulation. \n
  !> Active variables/parameters in state vector
  !> 1:h, 2:t, 3:c, 4:kz, 5:lz, 6:por \n
  integer, dimension(num_enkf_vars) :: act_s

  !> @brief Assimilation type
  !> @details
  !> Assimilation type \n
  !> Sequential of combined, read in from EnKF
  !> Input File under '# assimod'.
  !> If 'sequ', assitype=True, if 'comb', assitype=False. \n
  !> Assimilation type for multivariate state vector. \n
  !> combined - all variables are written into one state vector to be
  !> assimilated \n
  !> sequential - all variables are assimilated on their one in one
  !> state vector
  logical :: assitype

  !> @brief Variance normalization switch
  !> @details
  !> Variance normalization switch \n
  !> If combined assimilation is chosen,
  !> variances can be scaled. \n
  !> Scaling of variances to similar values (multivariate)
  logical :: normobs

  !> @brief Damping factor
  !> @details
  !> Damping factor \n
  !> Damping factor for all variables except
  !> the tracer. Values between 0.0d0 and 1.0d0.
  !> Read in the EnKF input file under
  !> '# damping factors' \n
  !> assidamp: damping factor for assimilation (=< 1.d0) (few data
  !> points) \n
  double precision :: assidamp

  !> @brief Damping factor for tracer
  !> @details
  !> Damping factor for tracer \n
  !> Values between 0.0d0 and 1.0d0.
  !> Read in the EnKF input file under
  !> '# damping factors'
  double precision :: tracdamp

  !> @brief Frequency of output
  !> @details
  !> Frequency of output \n
  !> Every iassout assimilation step is printed
  !> for output. EnKF Input File '# iassout'. \n
  !> every 'iassout' assimilation step plot output is printed
  integer :: iassout

  !> @brief Error covariance output switch
  !> @details
  !> Error covariance output switch \n
  !> If .true., the error covariance is
  !> calculated and printed at the end of
  !> the run. \n
  !> true: calculates and prints assimilated error covariance (large,
  !> untested...)
  logical :: err_cov

  !> @brief Unit exclusion threshold
  !> @details
  !> Unit exclusion threshold \n
  !> Exclude all units .le. ex_unit from
  !> the assimilation step. \n
  !> Exclude units with number less than or equal to ex_unit
  integer :: ex_unit

  !> @brief Name of true model Tecplot file
  !> @details
  !> Name of true model Tecplot file \n
  !> Read in from EnKF Input File under '# true_name'
  character (len=100) :: true_name

  !> @brief Name of chemical true model Tecplot file
  !> @details
  !> Name of chemical true model Tecplot file \n
  !> Read in from EnKF Input File under '# true_chem_name'
  character (len=100) :: true_chem_name

  !> @brief Name of observation file
  !> @details
  !> Name of observation file \n
  !> File containing the observations and
  !> observation locations.
  !> EnKF Input File '# observation file names'
  character (len=100) :: obs_name

  !> @brief Assimilation debugging output switch
  !> @details
  !> Assimilation debugging output switch \n
  !> Switch that controls output during assimilation step. \n
  !> Switch that controls assimstp output files. \n
  !> Default is true, because originally, assimstp output was always
  !> produced. \n
  logical :: assimstp_switch

  !> @brief EnKF variables output switch
  !> @details
  !> EnKF variables output switch \n
  !> Switch controlling the output of assim_variable files containing
  !> the variable fields. EnKF Input File under '# vtk output enkf' \n
  !> Switch controlling assim_variables output. \n
  logical :: vtk_out_enkf

  !> @brief EnKF residuals output switch
  !> @details
  !> EnKF residuals output switch \n
  !> Switch controlling residual output.
  logical :: vtk_out_resid

  !> @brief EnKF standard deviations output switch
  !> @details
  !> EnKF standard deviations output switch \n
  !> Switch controlling standard deviation output. \n
  !> Switch controlling assim_variables output
  logical :: vtk_out_stddev

  !> @brief EnKF correlation output switch
  !> @details
  !> EnKF correlation output switch \n
  !> Switch controlling correlation output of
  !> all state vector variables with a single
  !> state vector variable.
  logical :: vtk_out_covs

  !> @brief EnKF single realisations output switch
  !> @details
  !> EnKF single realisations output switch \n
  !> Switch controlling assim_variables output (Single realisations)
  !> \n
  logical :: vtk_out_realisations

  !> @brief Switch controlling standard output every tenth update
  !> @details
  !> Switch controlling standard output every tenth update \n
  logical :: obs_standard_out

  !> @brief Number of realisations for output
  !> @details
  !> Number of realisations for output \n
  !> The first num_out_realisations realisations
  !> are included in assim_variables output.
  integer :: num_out_realisations

  !> @brief Debug output switch
  !> @details
  !> Debug output switch \n
  !> For example, analysis matrices should be dumped. \n
  !> Switch controlling debug output (f.e. dumping analysis matrices)
  !> \n
  logical :: ana_mat_out

  !> @brief EnKF logging output switch
  !> @details
  !> EnKF logging output switch \n
  !> Controls whether to output into enkf.log. \n
  !> Switch controlling debug output (enkf.log)
  logical :: enkf_log_out

  !> @brief Debug output switch
  !> @details
  !> Debug output switch \n
  !> Information about realisations. \n
  !> Switch controlling debug output (info about which realisation is
  !> used ) \n
  logical :: comp_real_out
    
  !---------------------------------------------------------
  !---------------------------------------------------------
  !         LENGTH STATE VECTOR / STDDEV/ RESID
  !---------------------------------------------------------
  !--------------------------------------------------------

  !> @brief Number of realizations
  !> @details
  !> Number of realizations \n
  !> The number of realizations is set to the value of sm_max, which
  !> is set under '# simulate' in the SHEMAT-Suite Input File. Also:
  !> The number of ensemble members.
  integer :: nrens

  !> @brief Number of active variables/parameters in state vector
  !> @details
  !> Number of active variables/parameters in state vector \n Read
  !> from the activity vector act_s. \n
  !> Number of active species in state vector
  integer :: n_act_s

  !> @brief Number of cells in the model
  !> @details
  !> Number of cells in the model \n
  !> Calculated as grid size in nodes: i0*j0*k0 \n
  !> grid size in nodes: i0*j0*k0
  integer :: lstate0

  !> @brief Number of entries of one variable in state vector
  !> @details
  !> Number of entries of one variable in state vector \n
  !> If there is no unit
  !> excluded: lstate == lstate0 \n
  !> Equivalently: Number of active grid cells. \n
  !> Equivalently: Length of one-variable part of the state vector. \n
  !> Equivalently: Number of nodes active in state vector per variable
  integer :: lstate

  !> @brief Length of the state vector.
  !> @details
  !> Length of the state vector. \n
  !> This is: lstate*n_act_s \n
  !> Equivalently: Number of active state variables
  !> (active grid cells * active variables) \n
  !> Length of the state vector (all variables)
  integer :: nstate

  !> @brief Length of first dimension of stddev array
  !> @details
  !> Length of first dimension of stddev array \n
  !> 2*nrobs_int, set in enkf_alloc_stddev()
  integer :: lstd

  !> @brief Length of residual array
  !> @details
  !> Length of residual array \n
  !> 2*nrobs_int \n
  !> Length of the residuals vector
  integer :: lres

  !> @brief Rank of variables in state vector
  !> @details
  !> Rank of variables in state vector \n
  !> Array containing the rank/place of each
  !> variable in the state vector; 0 for the first
  !> variable, 1 for the second and so on. Variables
  !> that are not in the state vector get a -1. \n
  !> Example for head, conc and kz in mem:
  !> rank_mem = [0,-1,1,2,-1,-1]
  integer, dimension(num_enkf_vars) :: rank_mem

  !---------------------------------------------------------
  !---------------------------------------------------------
  !                      OBSERVATIONS
  !---------------------------------------------------------
  !--------------------------------------------------------  

  !> @brief Observation interval index
  !> @details
  !> Observation interval index \n
  !> Index of current observation time. \n
  !> The number of the current measurement.
  integer :: irobs

  !> @brief Number of observation locations
  !> @details
  !> Number of observation locations \n
  !> This is calculated specifically for each update.
  integer :: nrobs

  !> @brief Maximum number of observation locations
  !> @details
  !> Maximum number of observation locations \n
  !> Mostly needed to allocate observation arrays
  !> such that all information can fit for each
  !> update. Set in enkf_read_obs_len(). \n
  !> Maximum number of observation points per observation time
  integer :: max_obs_loc

  !---------------------------------------------------------
  !---------------------------------------------------------
  !               OUTPUT RELATED
  !---------------------------------------------------------
  !--------------------------------------------------------  

  !> @brief Temporary saving point for project name
  !> @details
  !> Temporary saving point for project name \n
  !> The project name is saved in project, which is
  !> changed and used at several times around the
  !> EnKF algorithm. During this time the original
  !> project content is saved in project_org
  character (len=256) :: project_org

  !> @brief Names of the variables/parameters in the EnKF
  !> @details
  !> Names of the variables/parameters in the EnKF \n
  !> Mostly used for output.
  character (len=4), dimension(6) :: &
       enkf_variable_names = (/ "head", "temp", "conc", "kz  ", "lz  ", "por " /)

  !> @brief Matrix of single cell input
  !> @details
  !> Matrix of single cell input \n
  !> mat_single_out(:,1) -> Single cell output i-indices, i_so \n
  !> mat_single_out(:,2) -> Single cell output j-indices, j_so \n
  !> mat_single_out(:,3) -> Single cell output k-indices, k_so \n
  !> mat_single_out(:,4) -> Single cell
  !> output variable indices, ivar_so \n
  !> EnKF Input File under '# single cell output' \n
  !> Location and Variable for single-cell-output 
  !> and whether to take the logarithm (log) or not (lin) \n
  integer, allocatable, dimension(:,:) :: mat_single_out

  !> @brief Single cell output number
  !> @details
  !> Single cell output number \n
  integer :: num_single_out

  !> @brief Starting interval for single cell output
  !> @details
  !> Starting interval for single cell output \n
  !> EnKF Input File under '# single cell output times'
  integer :: iassout_single_start

  !> @brief Interval difference for single cell output
  !> @details
  !> Interval difference for single cell output \n
  !> EnKF Input File under '# single cell output times'
  integer :: iassout_single

  !> @brief Switch for initial single cell output
  !> @details
  !> Switch for initial single cell output \n
  logical :: switch_so_ini

  !> @brief Tecplot or Paraview output switch for cell-centered values
  !> @details
  !> Tecplot or Paraview output switch for cell-centered values \n
  !> Switch to produce cell_centered output.
  logical :: cell_centered

  !> @brief Starting interval for reference point covariance output
  !> @details
  !> Starting interval for reference point covariance output \n
  integer :: iassout_cov_ref_start

  !> @brief Interval difference for reference point covariance output
  !> @details
  !> Interval difference for reference point covariance output \n
  integer :: iassout_cov_ref

  !> @brief Reference point covariance output switch
  !> @details
  !> Reference point covariance output switch \n
  logical :: switch_cov_ref

  !> @brief Reference point covariance output information
  !> @details
  !> Matrix with reference point covariance output information \n
  !> Information from i_cov_ref, j_cov_ref, k_cov_ref
  !> and ivar_cov_ref. \n
  !> mat_cov_ref(:,1) -> Reference point covariance i-indices, i_cov_ref \n
  !> mat_cov_ref(:,2) -> Reference point covariance j-indices, j_cov_ref \n
  !> mat_cov_ref(:,3) -> Reference point covariance k-indices, k_cov_ref \n
  !> mat_cov_ref(:,4) -> Reference point covariance variable indices \n
  !> Reference point covariance switch, location and variable \n
  integer, allocatable, dimension(:,:) :: mat_cov_ref

  !> @brief Number of reference point covariance outputs
  !> @details
  !> Number of reference point covariance outputs \n
  integer :: num_cov_ref

  !> @brief Switch to control computation of mean
  !> @details
  !> Switch to control computation of mean \n
  !> If yes, then the mean of parameters and variables is computed and
  !> then saved and output as ensemble member no. nsmpl. \n
  !> Switch to control whether a final mean is computed and saved to
  !> ensemble member nsmpl \n
  logical :: compute_mean

  !---------------------------------------------------------
  !---------------------------------------------------------
  !               SIMTIME RELATED
  !---------------------------------------------------------
  !--------------------------------------------------------  

  !> @brief Temporary Sample Timestep Index
  !> @details
  !> Temporary Sample Timestep Index \n
  integer :: itimestep_ismpl

  !> @brief Iteration counter
  !> @details
  !> Iteration counter \n
  !> Used in EnKF, simul
  integer :: iter_out

  !> @brief Starting time of forward computation
  !> @details
  !> Starting time of forward computation \n
  double precision :: ttstart

  !> @brief Ending time of forward computation
  !> @details
  !> Ending time of forward computation \n
  double precision :: ttend

  !---------------------------------------------------------
  !---------------------------------------------------------
  !               PARAMETERS
  !---------------------------------------------------------
  !--------------------------------------------------------   

  !> @brief Parameter for double-precision zero
  !> @details
  !> Parameter for double-precision zero \n
  integer, parameter :: zero = 0.d0

  !> @brief Flags for forward computation
  !> @details
  !> Flags for forward computation \n
  character (len=64) :: sflags

  !---------------------------------------------------------
  !---------------------------------------------------------
  !               NORMAL SCORE/VARIABLE DISTRIBUTIONS
  !---------------------------------------------------------
  !-------------------------------------------------------- 

  !> @brief Normal-Score EnKF switch
  !> @details
  !> Normal-Score EnKF switch \n
  !> EnKF Input File under '# normal score' \n
  !> Switch for using normal_score transformation.
  logical :: ns_switch

  !> @brief Normal-Score EnKF backtransformation factor
  !> @details
  !> Normal-Score EnKF backtransformation factor \n
  !> Factor controlling the width of the backtransformation
  !> distribution.
  double precision :: ns_backfactor

  !> @brief Ordered original state vector
  !> @details
  !> Ordered original state vector \n
  !> Original mem-values, but sorted. \n
  !> Saved during Normal Score EnKF and
  !> used for backtransform.
  double precision, allocatable, dimension(:,:) :: mem_ns_original

  !> @brief Ordered Gaussian Ensemble
  !> @details
  !> Ordered Gaussian Ensemble \n
  !> Sorted gaussian, Zero mean, variance one
  double precision, allocatable, dimension(:,:) :: mem_ns_gauss

  !> @brief Normalized observation variances for NS-EnKF
  !> @details
  !> Normalized observation variances for NS-EnKF \n
  !> The original observation variances are taken as normalization.
  double precision, allocatable, dimension(:,:) :: obsvar_ns

  !> @brief Gaussian ensembles at observation locations
  !> @details
  !> Gaussian ensembles at observation locations \n
  !> This array is used in the updated, it is sorted in the same way
  !> as the variable ensembles at the observation locations were
  !> originally. Mean zero, variance one.
  double precision, allocatable, dimension(:,:,:) :: obs_ns

  !> @brief Sorted original variable ensembles at observations locations
  !> @details
  !> Sorted original variable ensembles at observations locations \n
  !> Variable ensembles are set at the observatoin locations and
  !> subsequently sorted.
  double precision, allocatable, dimension(:,:,:) :: obs_ns_original

  !> @brief Sorted gaussian ensembles at observation locations
  !> @details
  !> Sorted gaussian ensembles at observation locations \n
  !> This array is saved for the backtransform. Mean zero, variance
  !> one.
  double precision, allocatable, dimension(:,:,:) :: obs_ns_gauss

  !> @brief Variances of original state variables at observation locations
  !> @details
  !> Variances of original state variables at observation locations \n
  !> Depends on observation location and variable index.
  double precision, allocatable, dimension(:,:) :: obs_ns_orig_vars

  !> @brief Reference distribution
  !> @details
  !> Reference distribution for  \n
  double precision, allocatable, dimension(:) :: h_dist_in

  !> @brief Length of reference distribution
  !> @details
  !> Length of reference distribution \n
  integer :: h_dist_len

  !> @brief Number of reference distributions
  !> @details
  !> Number of reference distributions \n
  integer :: num_ref_dist

  !> @brief Reference distribution switch
  !> @details
  !> Reference distribution switch \n
  logical :: ref_dist_switch

  !> @brief Reference distribution variable
  !> @details
  !> Reference distribution variable \n
  integer, dimension(num_enkf_vars) :: ref_dist_var

  !> @brief File name for reference distributions
  !> @details
  !> File name for reference distributions \n
  character (len=80), dimension(num_enkf_vars) :: ref_dist_file_name

  !---------------------------------------------------------
  !---------------------------------------------------------
  !               COVARIANCE LOCALISATION
  !---------------------------------------------------------
  !-------------------------------------------------------- 

  !> @brief Correlation functions
  !> @details
  !> Correlation functions \n
  !> Correlations between state vector entries and observation
  !> locations.
  double precision, allocatable :: rhos(:,:,:)

  !> @brief Localization switch
  !> @details
  !> Localization switch \n
  !> Switch for Local EnKF. \n
  !> Switch for covariance localisation
  logical :: cov_loc_switch

  !> @brief Correlation length in x direction
  !> @details
  !> Correlation length in x direction \n
  double precision :: cov_loc_lenx

  !> @brief Correlation length in y direction
  !> @details
  !> Correlation length in y direction \n
  double precision :: cov_loc_leny

  !> @brief Correlation length in z direction
  !> @details
  !> Correlation length in z direction \n
  double precision :: cov_loc_lenz

  !> @brief Index for multiple covariance localizations
  !> @details
  !> Index for multiple covariance localizations \n
  !> Not yet implemented: At the moment
  !> this is always set to 1.q
  integer :: cov_loc_assim_id

  !---------------------------------------------------------
  !---------------------------------------------------------
  !               DUAL ENKF
  !---------------------------------------------------------
  !-------------------------------------------------------- 

  !> @brief Dual EnKF switch
  !> @details
  !> Dual EnKF switch \n
  !> Switch, which controls,whether the dual Enkf (Moradkhani 2004)
  !> should be used
  logical :: dual_enkf_switch

  !> @brief Number of active parameters
  !> @details
  !> Number of active parameters \n
  !> Dual EnKF.
  integer :: n_act_s_param

  !> @brief Number of active state variables
  !> @details
  !> Number of active state variables \n
  !> Dual EnKF.
  integer :: n_act_s_state

  !> @brief Length of parameter part of state vector
  !> @details
  !> Length of parameter part of state vector \n
  !> Dual EnKF.
  integer :: nstate_param

  !> @brief Length of state variable part of state vector
  !> @details
  !> Length of state variable part of state vector \n
  !> Dual EnKF.
  integer :: nstate_state

  !> @brief Temporary number of active state vector parameters/variables
  !> @details
  !> Temporary number of active state vector parameters/variables \n
  !> Dual EnKF.
  integer :: n_act_s_tmp

  !> @brief Temporary length of state vector
  !> @details
  !> Temporary length of state vector \n
  !> Dual EnKF.
  integer :: nstate_tmp

  !> @brief Temporary head
  !> @details
  !> Temporary head \n
  !> Dual EnKF.
  double precision, allocatable :: head_tmp(:,:,:,:)

  !> @brief Temporary temp
  !> @details
  !> Temporary temp \n
  !> Dual EnKF.
  double precision, allocatable :: temp_tmp(:,:,:,:)

  !> @brief Temporary conc
  !> @details
  !> Temporary conc \n
  !> Dual EnKF.
  double precision, allocatable :: conc_tmp(:,:,:,:,:)

  !> @brief Temporary simulation time
  !> @details
  !> Temporary simulation time \n
  !> Dual EnKF.
  double precision, allocatable :: simtime_tmp(:)

  !> @brief Temporary number of active state variables/parameters in state vector
  !> @details
  !> Temporary number of active state variables/parameters in state vector \n
  !> Dual EnKF.
  integer, dimension(num_enkf_vars) :: act_s_tmp
  
  !---------------------------------------------------------
  !---------------------------------------------------------
  !               STOCHASTIC BOUNDARY CONDITIONS
  !---------------------------------------------------------
  !-------------------------------------------------------- 

  !> @brief Stochastic boundary condition switch
  !> @details
  !> Stochastic boundary condition switch \n
  !> Switch, which controls,whether the boundary conditions should be
  !> perturbed
  logical :: stoch_bc_switch

  !> @brief Number of stochastic boundary conditions
  !> @details
  !> Number of stochastic boundary conditions \n
  integer :: num_stoch_bc

  !> @brief Seed for stochastic boundary conditions
  !> @details
  !> Seed for stochastic boundary conditions \n
  integer, dimension(2) :: stoch_bc_seed_seed

  !> @brief Stochastic boundary conditions standard deviations
  !> @details
  !> Stochastic boundary conditions standard deviations \n
  double precision, dimension(num_enkf_vars) :: stoch_bc_stddevs

  !> @brief Stochastic boundary condition array
  !> @details
  !> Stochastic boundary condition array \n
  !> Function not clear
  integer, allocatable :: stoch_bc_nbc(:)

  !> @brief Stochastic boundary condition array
  !> @details
  !> Stochastic boundary condition array \n
  !> Function not clear
  integer, allocatable :: stoch_bc_ibc(:)

  !> @brief General random number seed switch
  !> @details
  !> General random number seed switch \n
  !> If .true., the general seed will be
  !> set in enkf_iter.f90
  logical :: general_seed_switch

  !> @brief General random number seed
  !> @details
  !> General random number seed \n
  !> Two number seed.
  integer, dimension(2) :: general_seed_seed 

  !---------------------------------------------------------
  !---------------------------------------------------------
  !               PRESCRIBED VELOCITY ESTIMATION
  !---------------------------------------------------------
  !-------------------------------------------------------- 

  !> @brief Prescribed velocity switch
  !> @details
  !> Prescribed velocity switch \n
  !> Switch, which controls whether Prescribed velocity is
  !> estimated. \n
  logical :: pres_vel_switch

    !> @brief Number of prescribed velocity components
  !> @details
  !> Number of prescribed velocity components \n
  !> Should be 2 or 3.
  integer :: num_pres_vel

  !> @brief Prescribed velocity random seed
  !> @details
  !> Prescribed velocity random seed \n
  !> Seed for Prescribed velocity. \n
  integer, dimension(2) :: pres_vel_seed

  !> @brief Prescribed velocity standard deviations
  !> @details
  !> Prescribed velocity standard deviations \n
  double precision, dimension(3) :: vdefault_stddevs

  !> @brief Prescribed velocity system variances
  !> @details
  !> Prescribed velocity system variances \n
  double precision, dimension(3) :: vdefault_sysvarmem

  !> @brief Thermal conductivity switch
  !> @details
  !> Thermal conductivity switch \n
  !> Switch, whether to estimate single thermal conductivity value \n
  logical :: tcon_switch

  !> @brief Number of thermal conductivities
  !> @details
  !> Number of thermal conductivities \n
  integer :: num_tcon

  !> @brief Indices of thermal conductivities
  !> @details
  !> Indices of thermal conductivities \n
  integer, allocatable :: tcon_inds(:)

  !> @brief Thermal conductivity standard deviation
  !> @details
  !> Thermal conductivity standard deviation \n
  double precision :: tcon_stddev

  !> @brief Thermal conductivity system variance
  !> @details
  !> Thermal conductivity system variance \n
  double precision :: tcon_sysvarmem
  
  !---------------------------------------------------------
  !---------------------------------------------------------
  !            TEMPERATURE DIFFERENCE OBSERVATION
  !---------------------------------------------------------
  !-------------------------------------------------------- 

  !> @brief Temperature difference switch
  !> @details
  !> Temperature difference switch \n
  !> Switch, whether to use temperature difference observations \n
  logical :: tempdiff_switch

  !> @brief Number of temperature difference observations
  !> @details
  !> Number of temperature difference observations \n
  integer :: num_tempdiff

  !> @brief Indices of temperature difference observations
  !> @details
  !> Indices of temperature difference observations \n
  integer, allocatable :: tempdiff_inds(:,:)
  
  !---------------------------------------------------------
  !---------------------------------------------------------
  !               HYBRID ENKF
  !---------------------------------------------------------
  !--------------------------------------------------------

  !> @brief Hybrid EnKF switch
  !> @details
  !> Hybrid EnKF switch \n
  !> Switch, which controls whether a hybrid EnKF is used \n
  logical :: hybrid_switch

  !> @brief Hybrid EnKF covariance matrix
  !> @details
  !> Hybrid EnKF covariance matrix \n
  double precision, allocatable :: pb(:,:)

  !> @brief Hybrid EnKF covariance submatrix
  !> @details
  !> Hybrid EnKF covariance submatrix \n
  double precision, allocatable :: pb_r(:,:)

  !> @brief Hybrid EnKF covariance submatrix
  !> @details
  !> Hybrid EnKF covariance submatrix \n
  double precision, allocatable :: pb_reps(:,:)

  !> @brief Hybrid EnKF mixing parameter
  !> @details
  !> Hybrid EnKF mixing parameter \n
  !> Real determing the weight of static/dynamic parts of the cov \n
  double precision :: hybrid_alpha

  !> @brief Hybrid EnKF covariance kind
  !> @details
  !> Hybrid EnKF covariance kind \n
  !> Not yet implemented: Possibility to
  !> specify different kinds of background
  !> covariances.
  integer :: hybrid_cov_kind

  !---------------------------------------------------------
  !---------------------------------------------------------
  !               ITERATIVE ENKF
  !---------------------------------------------------------
  !--------------------------------------------------------

  !> @brief Original number of observation intervals
  !> @details
  !> Original number of observation intervals \n
  !> This saves the original number of observation
  !> intervals for Dual EnkF and Iterative EnKF,
  !> in which this number is changed.
  integer :: nrobs_int_pure

  !> @brief Iterative EnKF switch
  !> @details
  !> Iterative EnKF switch \n
  !> Switch, which controls whether Iterative EnKF is used \n
  logical :: iterative_switch

  !> @brief Number of observation intervals skipped in Iterative EnKF
  !> @details
  !> Number of observation intervals skipped in Iterative EnKF \n
  !> After this number of intervals Iterative EnKF restarts. \n
  !> Interval after which to return to start variables \n
  integer :: iterative_nrobs_int

  !> @brief Number of starts in Iterative EnKF
  !> @details
  !> Number of starts in Iterative EnKF \n
  !> For the first update, this is set to 1, then it is used to find
  !> the update, when the first restart is done, the number is set to
  !> two and so on.
  integer :: iterative_irobs

  !> @brief Iterative EnKF multiple updates
  !> @details
  !> Iterative EnKF multiple updates \n
  !> Multiple update of the same measurement information.
  logical :: iterative_doubleupdate

  !---------------------------------------------------------
  !---------------------------------------------------------
  !               PILOT POINT ENKF
  !---------------------------------------------------------
  !--------------------------------------------------------

  !> @brief Number of Pilot Points
  !> @details
  !> Number of Pilot Points \n
  integer :: num_pp

  !> @brief State vector before update with Pilot Point restriction
  !> @details
  !> State vector before update with Pilot Point restriction \n State
  !> vector containing all variables and parameter. The parameter, for
  !> which pilot points are chosen is restricted to the pilot
  !> points. This state vector saves the values before the EnKF
  !> update.
  double precision, allocatable :: mem_pp(:,:)

  !> @brief Non-pilot point locations of pilot point parameter
  !> @details
  !> Non-pilot point locations of pilot point parameter \n The
  !> parameter, for which pilot points are chosen is restricted to the
  !> pilot points. This state vector saves the values outside the
  !> pilot point locations. They are needed in the update.
  double precision, allocatable :: mem_rp(:,:)

  !> @brief Temporary nstate in Pilot Points
  !> @details
  !> Temporary nstate in Pilot Points \n
  !> Holds the original nstate for the
  !> time of the update.
  integer :: nstate_pp_temp

  !> @brief Pilot-Point EnKF Switch
  !> @details
  !> Pilot-Point EnKF Switch \n
  !> Switch, whether to use pilot points \n
  logical :: pp_switch

  !> @brief Pilot-Point Background covariance
  !> @details
  !> Pilot-Point Background covariance \n
  !> Set in subroutine enkf_cov_pp.
  double precision, allocatable :: gss(:,:)

  !> @brief Pilot-Point cell indices
  !> @details
  !> Pilot-Point cell indices \n
  !> Indices of the Pilot Point locations read
  !> under '# pilot point'
  integer, allocatable :: pp_inds(:,:)

  !> @brief Pilot-Point state vector indices
  !> @details
  !> Pilot-Point state vector indices \n
  integer, allocatable :: pp_state_inds(:)

  !> @brief Pilot-Point lstate vector indices
  !> @details
  !> Pilot-Point lstate vector indices \n
  integer, allocatable :: pp_lstate_inds(:)

  !> @brief Pilot-Point variable index
  !> @details
  !> Pilot-Point variable index \n
  integer :: pp_ivar

  !> @brief Pilot-Point EnKF Covariance Switch
  !> @details
  !> Pilot-Point EnKF Covariance Switch \n
  !> Possible values:
  !> - 'def': Default, calculate covariance from mem
  !> - 'get': Read gss from file gss.out
  !> - 'out': Write gss to file gss.out
  character (len=3) :: pp_get_out

  !> @brief First original state vector index of pilot point variable
  !> @details
  !> First original state vector index of pilot point variable \n
  integer :: ipp_start

  !> @brief Last original state vector index of pilot point variable
  !> @details
  !> Last original state vector index of pilot point variable \n
  integer :: ipp_end

end module mod_enkf
