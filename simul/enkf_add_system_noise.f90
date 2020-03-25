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

!> @brief Adding system noise to state vector
!> @details
!> The samples for the state-vector variables (which may all be
!> identical) are subjected to gaussian system noise. The system noise
!> variances are given in the EnKF-input file.
subroutine enkf_add_system_noise()

  use m_random

  use mod_enkf, only: &
       nstate,&
       nrens,&
       sysvarmem,&
       mem

  implicit none

  double precision, allocatable :: syswork(:)

  integer :: irens, l

  ALLOCATE(syswork(nstate))
  syswork(:) = 0.D0
  
  DO irens = 1, nrens

     ! Normal distribution with mean one and variance zero
     CALL random(syswork,nstate)

     ! Add Gaussian noise scaled by sysvarmem
     DO l = 1, nstate
        mem(l,irens) = mem(l,irens) + &
             dsqrt(sysvarmem(l))*syswork(l)
     END DO

  END DO
  
  DEALLOCATE(syswork)

end subroutine enkf_add_system_noise
