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

!>    @brief commonly used global variables and constants
module mod_genrl
      double precision tslocal
!

      !> @brief Memory of data in 8Bytes (std. double precision size).
      !> @details
      !> Memory of data in 8Bytes (std. double precision size). \n
      double precision memory

!     file output compression tool
      integer compress_out

      !> @brief Number of cells in x-direction.
      !> @details
      !> Number of cells in x-direction. \n
      !> Read in under `# grid`, first entry.
      integer :: i0

      !> @brief Number of cells in y-direction.
      !> @details
      !> Number of cells in y-direction. \n
      !> Read in under `# grid`, second entry.
      integer :: j0

      !> @brief Number of cells in z-direction.
      !> @details
      !> Number of cells in z-direction. \n
      !> Read in under `# grid`, third entry.
      integer :: k0

      !> @brief Output orientation of monitoring output.
      !> @details
      !> Output orientation of monitoring output. \n \n
      !> `out_orientation` values: \n
      !>   0: normal (pv column wise), each time step has its own block
      !>      (mp row wise) \n
      !>   1: transposed (pv row wise), each time step has its own
      !>      block (mp column wise) \n
      !>   2: normal (pv column wise), new file for each time step (mp
      !>      row wise) \n
      !>   3: transposed (pv row wise), new file for each time step
      !>      (mp column wise) \n
      !>   4: normal (pv column wise), new file for each monitor point
      !>      (time row wise) \n
      integer :: out_orientation

      !> @brief Runmode integer.
      !> @details
      !> Runmode integer. \n
      !> Read in under `# runmode`. \n\n
      !>
      !> 0: forward modeling \n
      !> 1: forward modeling and data fit \n
      !> 2: inverse modeling or extra steady state \n
      !> 3: inverse modeling \n
      integer :: runmode

      !> @brief Maximum number of nonlinear iterations.
      !> @details
      !> Maximum number of nonlinear iterations. \n
      !> Important for the Picard iterations.\n\n
      !>
      !> Read under `# nlsolve`, first entry.
      integer :: maxiter_nl

      integer nladapt

      !> @brief Switch for checking non-linear iteration convergence.
      !> @details
      !> Switch for checking non-linear iteration convergence. \n
      !> - 0: check convergence (default).\n
      !> - other: do not check convergence.\n
      !>
      !> Read under `# nlsolve`, second line, only entry.
      integer nlconverge
!
      logical head_active, temp_active
      logical trac_active,chem_active,trans_active, pres_active
      logical vtk_out,tec_out,hdf_out,txt_out

      !> @brief Switch informing if an initial flow transformation is needed.
      !> @details
      !>  Switch informing if an initial flow transformation is needed. \n
      !> This variable is .false. if the not-computed flow-variable
      !> (f.e. pres in head-based computation) is given in the input
      !> file. Otherwise `head2pres` or `pres2head` will be called
      !> inside /forward/forward_init.f90.
      logical is_init_flow_trafo_needed

      !> @brief Switch for disabling multiple outputs.
      !> @details
      !> Switch for disabling multiple outputs. \n
      !> Disables output in write_logs, write_data, write_outt,
      !> write_tecdiff, write_monitor.
      logical write_disable
      logical write_smonitor, write_param, write_eoutt

      !> @brief Switch for disabling iteration output.
      !> @details
      !> Switch for disabling iteration output. \n
      !> Disables output from forward iteration.
      logical write_iter_disable

      !> @brief Switch for reading external input files.
      !> @details
      !> Switch for reading external input files. \n
      !> If and only if read_external_input is .true., the input file
      !> is searched for external input files names in the subroutine
      !> read_control.
      logical read_external_input

      !> @brief Hashtag character.
      !> @details
      !> Hashtag character. \n
      !> Used as input keyword signifier.
      character (len=1), parameter :: key_char = '#'

      !> @brief Number of samples.
      !> @details
      !> Number of samples. \n
      !> For example, the number of samples in a Monte Carlo
      !> simulation. Also used in OpenMP parallelization.
      integer nsmpl

      ! Useful constants
      ! ----------------

      !> @brief Very big number.
      !> @details
      !> Very big number. \n\n
      double precision, parameter :: big = 1.0d30

      !> @brief Very small number.
      !> @details
      !> Very small number. \n\n
      double precision, parameter :: small = 1.0d-300

      !> @brief Circle number pi.
      !> @details
      !> Circle number pi. \n\n
      double precision, parameter :: pi = 3.141592653589793d0

      !> @brief Flag for log permeability output.
      !> @details
      !> Flag for log permeability output. \n\n
      logical, parameter :: klogflag = .true.

      !> @brief Number of possible old saved arrays.
      !> @details
      !> Number of possible old saved arrays. \n\n
      !>
      !> max number of "cgen_*"
      integer, parameter :: ncgen = 5

      !> @brief Index constant for old saved arrays.
      !> @details
      !> Index constant for old saved arrays. \n\n
      !>
      !> index constants for HEAD, TEMP and PRES
      !> forward iteration [x_i-1,j];
      !> do not change it without a change of any
      !> "headold" into "headold(1,cgen_fw)" ccc
      integer, parameter :: cgen_fw = 1

      !> @brief Index constant for old saved arrays.
      !> @details
      !> Index constant for old saved arrays. \n\n
      !>
      !> Used around nonlinear forward iteration in `forward_iter`.
      !>
      !> index constants for HEAD, TEMP and PRES
      !> forward iteration [x_i-1,j];
      !> do not change it without a change of any
      !> "headold" into "headold(1,cgen_fw)" ccc\n\n
      !>
      !>       time iteration
      integer, parameter :: cgen_time = 2

      !> @brief Index constant for old saved arrays.
      !> @details
      !> Index constant for old saved arrays. \n\n
      !>
      !> index constants for HEAD, TEMP and PRES
      !> forward iteration [x_i-1,j];
      !> do not change it without a change of any
      !> "headold" into "headold(1,cgen_fw)" ccc\n\n
      !>
      !>       nonlinear-forward iteration (use of ad)
      integer, parameter :: cgen_fwad = 3
!
!       optimisation/realisation iteration [x_*,j]
!         index constants for g_HEAD, g_TEMP and g_PRES
        integer, parameter :: cgen_opti = 4
!       special DD use [f(x_*,j)]
        integer, parameter :: cgen_dd = 5

      !       index of the <master> sample
      integer, parameter :: idx_master=1

        !> @brief Doubling threshold for variable time step size
        !> @details
        !> Doubling threshold for variable time step size. \n\n
        !>
        !> Once delt_count reaches delt_double, the time step size is
        !> doubled. \n\n
        !>
        !> Convergence flag, doubling value for variable time step
        !> size
        integer delt_double

        !> @brief Variable step size nonlinear iteration counter.
        !> @details
        !> Variable step size nonlinear iteration counter. \n\n
        !>
        !> This is a helper-counter that together with the current
        !> number of iterations (in `forward_picard`) decides whether
        !> the doubling counter of the variable time step size utitily
        !> should be moved towards or away from doubling the time step
        !> size.
        integer :: iter_nlold

        !> @brief Flag for Variable step size
        !> @details
        !> Flag for Variable step size. \n\n
        !>
        !> True if variable step size is used.
        logical :: delt_vary

        !> @brief Initial time step size for variable step size
        !> @details
        !> Initial time step size for variable step size. \n\n
        !>
        !> First input under `# variable step size`.
        double precision :: delt_start

        !> @brief Minimum time step size for variable step size
        !> @details
        !> Minimum time step size for variable step size. \n\n
        !>
        !> Second input under `# variable step size`.
        double precision :: delt_min

        !> @brief Maximum time step size for variable step size
        !> @details
        !> Maximum time step size for variable step size. \n\n
        !>
        !> Third input under `# variable step size`.
        double precision :: delt_max

end module mod_genrl
