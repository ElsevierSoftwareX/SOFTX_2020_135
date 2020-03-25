! MIT License
!
! Copyright (c) 2019 Geir Evensen
!
! Permission is hereby granted, free of charge, to any person
! obtaining a copy of this software and associated documentation
! files (the "Software"), to deal in the Software without
! restriction, including without limitation the rights to use, copy,
! modify, merge, publish, distribute, sublicense, and/or sell copies
! of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be
! included in all copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
! EXPRESS OR IMPLIED, INCLUDING BUT NOT NO EVENT SHALL THE AUTHORS OR
! COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
! ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE
! OR OTHER DEALINGS IN THE SOFTWARE.

      MODULE m_obs_pert_enkf

      CONTAINS
        SUBROUTINE obs_pert(e,nrens,nrobs,fixsamp)
!  construct observation pertubation
! using a 'diagonal' covariance model

          USE m_random
          USE m_fixsample

          IMPLICIT NONE
          integer, INTENT (IN) :: nrobs
          integer, INTENT (IN) :: nrens
          LOGICAL, INTENT (IN) :: fixsamp
          double precision, INTENT (OUT), dimension (nrobs,nrens) :: e

          CALL random(e,nrobs*nrens)
          IF (fixsamp) CALL fixsample(e,nrobs,nrens)

        END SUBROUTINE obs_pert
      END MODULE m_obs_pert_enkf
