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

!> @brief Backtransform of the Normal-Score EnKF
!> @details
!> First, the Tails for the backtransform are defined using the
!> minimum and maximum of the original state vector, as well as
!> ns_backfactor from the EnKF Input File.
!>
!> The backtransform is carried out by the function backtr() from
!> simul.
subroutine enkf_normal_score_back()

  use mod_enkf, only: &
       nstate,&
       nrens,&
       mem,&
       mem_ns_original,&
       mem_ns_gauss,&
       ns_backfactor

  implicit none

  integer :: istate, i

  double precision :: backtr
  double precision :: max_orig, min_orig
  
  do istate = 1, nstate
     
     max_orig = 0.0d0
     min_orig = 0.0d0
     call enkf_vector_min(mem_ns_original(:,istate),nrens,min_orig)
     call enkf_vector_max(mem_ns_original(:,istate),nrens,max_orig)
     
     ! Length of the tails in backtr (ns_backfactor, default = 3.0d0)
     max_orig = max_orig + ns_backfactor*(max_orig-min_orig)/(real(nrens,8)-1)
     min_orig = min_orig - ns_backfactor*(max_orig-min_orig)/(real(nrens,8)-1)
     
     do i = 1, nrens
        !Function call: backtr from simul/gs/backtr.f
        ! This is the backtransform
        mem(istate,i) = backtr(&
             mem(istate,i),&
             nrens,&
             mem_ns_original(:,istate),&
             mem_ns_gauss(:,istate),&
             min_orig,&
             max_orig,&
             1,&
             1.0d0,&
             1,&
             1.0d0)
     end do

  end do

end subroutine enkf_normal_score_back
