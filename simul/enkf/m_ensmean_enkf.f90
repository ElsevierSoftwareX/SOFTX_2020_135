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

      MODULE m_ensmean_enkf

      CONTAINS
        SUBROUTINE ensmean(a,ave,nx,nrens)
!-----------
! part of EnKF code from Geir Evensen   http://enkf.nersc.no/
!-----------
          IMPLICIT NONE
          integer, INTENT (IN) :: nx
          integer, INTENT (IN) :: nrens
          double precision, INTENT (IN), dimension (nx,nrens) :: a
          double precision, INTENT (OUT), dimension (nx) :: ave
          integer :: j

          ave(:) = a(:,1)
          DO j = 2, nrens
            ave(:) = ave(:) + a(:,j)
          END DO
          ave = (1.0D0/dble(nrens))*ave

        END SUBROUTINE ensmean
      END MODULE m_ensmean_enkf
