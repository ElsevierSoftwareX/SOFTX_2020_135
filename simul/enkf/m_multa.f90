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

      MODULE m_multa

      CONTAINS
        SUBROUTINE multa(a,x,ndim,nrens,iblkmax)
!-----------
! part of EnKF code from Geir Evensen   http://enkf.nersc.no/
!-----------
          IMPLICIT NONE
          integer, INTENT (IN) :: ndim
          integer, INTENT (IN) :: nrens
          integer, INTENT (IN) :: iblkmax
          double precision, INTENT (IN) :: x(nrens,nrens)
          double precision, INTENT (INOUT) :: a(ndim,nrens)
          double precision v(iblkmax,nrens) 
! Automatic work array              
          integer :: ifirst, ilast

          DO ifirst = 1, ndim, iblkmax
            ilast = min(ifirst+iblkmax-1,ndim)
            v(1:ilast-ifirst+1,1:nrens) = a(ifirst:ilast,1:nrens)
            CALL dgemm('n','n',ilast-ifirst+1,nrens,nrens,1.0D0, &
              v(1,1),iblkmax,x(1,1),nrens,0.0D0,a(ifirst,1),ndim)
          END DO

        END SUBROUTINE multa
      END MODULE m_multa
