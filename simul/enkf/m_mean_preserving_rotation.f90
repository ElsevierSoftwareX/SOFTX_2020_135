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

      MODULE m_mean_preserving_rotation

      CONTAINS
        SUBROUTINE mean_preserving_rotation(up,nrens)
! Generates the mean preserving random rotation for the EnKF SQRT algorithm
! using the algorithm from Sakov 2006-07.  I.e, generate rotation Up suceh that
! Up*Up^T=I and Up*1=1 (all rows have sum = 1)  see eq 17.
! From eq 18,    Up=B * Upb * B^T
! B is a random orthonormal basis with the elements in the first column equals 1/sqrt(nrens)
! Upb = | 1  0 |
!       | 0  U |
! where U is an arbitrary orthonormal matrix of dim nrens-1 x nrens-1  (eq. 19)

          USE m_randrot

          IMPLICIT NONE

          integer, INTENT (IN) :: nrens
          double precision, INTENT (OUT) :: up(nrens,nrens)

          double precision b(nrens,nrens), q(nrens,nrens), r(nrens,nrens)
          double precision u(nrens-1,nrens-1), upb(nrens,nrens)

          integer :: j, k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Generating the B matrix
! Starting with a random matrix with the correct 1st column
          CALL random_number(b)
          b(:,1) = 1.0D0/dsqrt(dble(real(nrens)))

! modified_gram_schmidt is used to create the orthonormal basis
!   do k=1,nrens
!      R(k,k)=sqrt(dot_product(B(:,k),B(:,k)))
!      Q(:,k)=B(:,k)/R(k,k)
!      do j=k+1,nrens
!         R(k,j)=dot_product(Q(:,k),B(:,j))
!         B(:,j)=B(:,j)- Q(:,k)*R(k,j)
!      enddo
!   enddo
!   B=Q

! with overwriting of B
          DO k = 1, nrens
            r(k,k) = dsqrt(dot_product(b(:,k),b(:,k)))
            b(:,k) = b(:,k)/r(k,k)
            DO j = k + 1, nrens
              r(k,j) = dot_product(b(:,k),b(:,j))
              b(:,j) = b(:,j) - b(:,k)*r(k,j)
            END DO
          END DO


! Check on orthonormality of B
!   do k=1,nrens
!      do j=k,min(k+14,nrens)
!         write(*,'(15f10.4)',advance='no')dot_product(B(:,j),B(:,k))
!      enddo
!      write(*,*)' '
!   enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Creating the orthonormal nrens-1 x nrens-1 U matrix
          CALL randrot(u,nrens-1)
!! Check on orthonormality of U
!   do k=1,nrens-1
!      do j=k,min(k+14,nrens-1)
!         write(*,'(15f10.4)',advance='no')dot_product(U(:,j),U(:,k))
!      enddo
!      write(*,*)' '
!   enddo
!   stop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Creating the orthonormal nrens x nrens Upb matrix
          upb(2:nrens,2:nrens) = u(1:nrens-1,1:nrens-1)
          upb(1,1) = 1.0
          upb(2:nrens,1) = 0.0
          upb(1,2:nrens) = 0.0

! Check on orthonormality of Upb
!   do k=1,nrens
!      do j=k,min(k+14,nrens)
!         write(*,'(15f10.4)',advance='no')dot_product(Upb(:,j),Upb(:,k))
!      enddo
!      write(*,*)' '
!   enddo
!   stop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Creating the random orthonormal mean preserving nrens x nrens Upb matrix: Up=B^T Upb B
          CALL dgemm('n','n',nrens,nrens,nrens,1.0d0,b,nrens,upb,nrens, &
            0.0d0,q,nrens)
          CALL dgemm('n','t',nrens,nrens,nrens,1.0d0,q,nrens,b,nrens, &
            0.0d0,up,nrens)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Checks
!   do k=1,nrens
!      print *,'Up row sum: ',k,sum(Up(k,1:nrens))
!   enddo

! Check on orthonormality of Up
!   do k=1,nrens
!      do j=k,min(k+14,nrens)
!         write(*,'(15f10.4)',advance='no')dot_product(Up(:,j),Up(:,k))
!      enddo
!      write(*,*)' '
!   enddo
!   stop

        END SUBROUTINE
      END MODULE m_mean_preserving_rotation
