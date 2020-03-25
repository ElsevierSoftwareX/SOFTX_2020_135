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

!> @brief global (dynamic) arrays for data exchange with (SGSIM, VISIM)
      MODULE simul_arrays
        DOUBLE PRECISION, ALLOCATABLE :: simout(:)
        DOUBLE PRECISION, ALLOCATABLE :: sim(:), lvm(:), avepor(:,:)
        DOUBLE PRECISION, ALLOCATABLE :: tmp(:), order(:), dvr(:)
        DOUBLE PRECISION, ALLOCATABLE :: porvar(:,:), krgvar(:)
        DOUBLE PRECISION, ALLOCATABLE :: cd2v(:,:)
        integer, ALLOCATABLE :: novar(:)
!     info to finalise grid adapting
        integer :: sm_i0, sm_j0, sm_k0
        DOUBLE PRECISION sm_delx, sm_dely, sm_delz
!
!     for OpenMP parallelisation, thread-save
!$OMP threadprivate(simout,sim,lvm,avepor,tmp,order)
!$OMP threadprivate(dvr,porvar,krgvar,cd2v,novar)
!$OMP threadprivate(sm_I0, sm_J0, sm_K0, sm_delx, sm_dely, sm_delz)
!
!     parameter file names
        CHARACTER (len=80), ALLOCATABLE :: fnpara(:)
!
      END MODULE simul_arrays
