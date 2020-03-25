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

      MODULE m_randrot

      CONTAINS
        SUBROUTINE randrot(q,nrens)
          IMPLICIT NONE
          integer, INTENT (IN) :: nrens
          double precision, INTENT (OUT), dimension (nrens,nrens) :: q

          double precision, DIMENSION (nrens,nrens) :: a, b
          double precision, dimension (nrens) :: sigma
          double precision, dimension (10*nrens) :: work
          double precision, PARAMETER :: pi = 3.14159253589D0
          integer :: ierr

          CALL random_number(b)
          CALL random_number(a)
          q = dsqrt(-2.D0*dlog(a+dble(tiny(a))))*dcos(2.D0*pi*b)

!$OMP CRITICAL (randrot_critical)
! QR factorization
          CALL dgeqrf(nrens,nrens,q,nrens,sigma,work,10*nrens,ierr)
          IF (ierr/=0) PRINT *, 'randrot: dgeqrf ierr=', ierr

! Construction of Q
          CALL dorgqr(nrens,nrens,nrens,q,nrens,sigma,work,10*nrens, &
            ierr)
          IF (ierr/=0) PRINT *, 'randrot: dorgqr ierr=', ierr
!$OMP END CRITICAL (randrot_critical)


        END SUBROUTINE randrot
      END MODULE m_randrot
