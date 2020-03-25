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

      MODULE m_ensvar_enkf

      CONTAINS
        SUBROUTINE ensvar(a,ave,var,nx,nrens)
          IMPLICIT NONE
          integer, INTENT (IN) :: nx
          integer, INTENT (IN) :: nrens
          double precision, INTENT (IN), dimension (nx,nrens) :: a
          double precision, INTENT (IN), dimension (nx) :: ave
          double precision, INTENT (OUT), dimension (nx) :: var
          integer :: j

          var = 0.0D0
          DO j = 1, nrens
            var(:) = var(:) + (a(:,j)-ave(:))*(a(:,j)-ave(:))
          END DO
          var = (1.0D0/dble(real(nrens-1)))*var

        END SUBROUTINE ensvar
      END MODULE m_ensvar_enkf


