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

!> @brief global variables for SIMUL (SGSIM, VISIM, ENKF)
module mod_simul

  !> @brief SIMUL(samples) output directory name
  !> @details
  !> SIMUL(samples) output directory name \n
  character (len=15) :: ssample_outdir
  parameter(ssample_outdir = 'samples_output/')

  !> @brief ENKF output directory name
  !> @details
  !> ENKF output directory name \n
  character (len=12) :: senkf_outdir
  parameter(senkf_outdir = 'enkf_output/')

  !> @brief ENKF-05 output directory name
  !> @details
  !> ENKF-05 output directory name \n
  character (len=10) :: s05_outdir
  parameter(s05_outdir = '05_output/')

  !> @brief ENKF- single cell output directory name
  !> @details
  !> ENKF- single cell output directory name \n
  character (len=19) :: ssingle_cell_outdir
  parameter(ssingle_cell_outdir = 'single_cell_output/')

  !> @brief global values constant for each thread
  !> @details
  !> global values constant for each thread \n
  integer, parameter :: MAXVOLS=5

  !> @brief Number of time period units
  !> @details
  !> Number of time period units \n
  !> Connected to /inverse/, not set inside /simul/
  integer mpara_tp

  !> @brief Number of SGSim parameter input files
  !> @details
  !> Number of SGSim parameter input files \n
  !> Read in under "# parameter group" or "# bc group" in the
  !> SHEMAT-Suite input file. max_gpara saves the sum of the
  !> corresponding "records" entries.
  integer max_gpara

  !> @brief Maximum number of nodes in X
  !> @details
  !> Maximum number of nodes in X \n
  !> Connected to SGSim of VISim. Point of setting not found in
  !> /sgsim/.
  integer MAXX

  !> @brief Maximum number of nodes in Y
  !> @details
  !> Maximum number of nodes in Y \n
  !> Connected to SGSim of VISim. Point of setting not found in
  !> /sgsim/.
  integer MAXY

  !> @brief Maximum number of nodes in Z
  !> @details
  !> Maximum number of nodes in Z \n
  !> Connected to SGSim of VISim. Point of setting not found in
  !> /sgsim/.
  integer MAXZ

  !> @brief (Pseudo-)Number of realizations
  !> @details
  !> (Pseudo-)Number of realizations \n

  !> Read in the SHEMAT-Suite input file under "# simulate". \n
  !> If runmode>=2, then sm_max-nsmpl=1. If runmode<2, then
  !> sm_max=nsmpl. \n
  !> For estimation runs (runmode>=2), nsmpl includes a last
  !> realization, in which the mean is saved.
  integer sm_max

  !> @brief Number for sample average approximation
  !> @details
  !> Number for sample average approximation \n
  !> Read in under "# subsample", but seems not to be used.
  integer saa_sample

  !> @brief EnKF postcomputation switch
  !> @details
  !> EnKF postcomputation switch \n
  !> Read in under "# enkf postcompute" in SHEMAT-Suite input file. \n
  !> enkf_post=0: No postcomputation ('none') \n
  !> enkf_post=1: Postcomputation of ensemble ('samples') \n
  !> enkf_post=2: Postcomputation of mean ('mean') \n
  integer enkf_post

  !> @brief Number of EnKF iterations
  !> @details
  !> Number of EnKF iterations \n
  !> Read in under "# enkf iter" in SHEMAT-Suite input file. \n
  !> Multiple EnKF iterations are carried out in shem_sm.f90.
  integer maxiter_enkf

end module mod_simul
