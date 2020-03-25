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

!>    @brief direct linear solver call: PLU from LAPACK, solve of : [M] x [x] = [b]
!>    @param[in] I0 lengths of I-dimension of local matrix [M]
!>    @param[in] J0 lengths of J-dimension of local matrix [M]
!>    @param[in] K0 lengths of K-dimension of local matrix [M]
!>    @param[out] x solution vector [x]
!>    @param[in] w right side, vector [b]
!>    @param[in] A 1. diagonal of the system matrix [M]
!>    @param[in] B 2. diagonal of the system matrix [M]
!>    @param[in] C 3. diagonal of the system matrix [M]
!>    @param[in] D 4. diagonal of the system matrix [M]
!>    @param[in] E 5. diagonal of the system matrix [M]
!>    @param[in] F 6. diagonal of the system matrix [M]
!>    @param[in] G 7. diagonal of the system matrix [M]
      SUBROUTINE direct_solve(i0,j0,k0,x,w,a,b,c,d,e,f,g)
        use mod_linfos
        IMPLICIT NONE

        INTEGER i0, j0, k0, kl, ku, ldm

        DOUBLE PRECISION, ALLOCATABLE :: matrix(:,:)
        INTEGER, ALLOCATABLE :: pivots(:)
        INTEGER error, i, j, k, ind, nijk

!     need to detect 64Bit enviroment
        INTEGER test32
        INTEGER (kind=8) test64

        DOUBLE PRECISION x(i0,j0,k0), w(i0,j0,k0), enough, enough2
        DOUBLE PRECISION a(i0,j0,k0), b(i0,j0,k0), c(i0,j0,k0), &
          d(i0,j0,k0), e(i0,j0,k0), f(i0,j0,k0), g(i0,j0,k0)


        nijk = i0*j0*k0
        kl = i0*j0
        IF (i0==1) kl = i0*j0
        IF (j0==1) kl = i0*j0
        IF (k0==1) kl = i0
        ku = kl
        ldm = 2*kl + ku + 1

!     detecting 64Bit-Integer enviroment (t32=t64 for 64Bit)
        test32 = 3000000
        test64 = 3000000
        test32 = test32*test32
        test64 = test64*test64

!     protecting against 32 pointers
        IF ((8.0D0*dble(ldm)*dble(nijk))>2.0D9) THEN
          IF (test32==test64) THEN
            WRITE(*,'(A,A)') 'Warning: try to use 64Bit integer', &
              ' enviroment (needs 64Bit system)'
          ELSE
            WRITE(*,'(A,A,F9.4,A)') &
              'Error: array to large for direct band', &
              ' matrix solver (', 8.0D0*dble(ldm)*dble(nijk)/1024.0D0 &
              **3, ' GBytes).'
            STOP
          END IF
        END IF

        ALLOCATE(matrix(ldm,nijk))
        ALLOCATE(pivots(nijk))
        error = 0

!     right side to [X]
        CALL dcopy(nijk,w,1,x,1)
        CALL set_dval(ldm*nijk,0.D0,matrix)

! --------------------------------------------------------------

        ind = 0
        DO k = 1, k0
          DO j = 1, j0
            DO i = 1, i0
!aw     ind=i +I0*(j-1) +I0*J0*(k-1)
              ind = ind + 1
!              if (i<20 .and. j==1 .and. k==1) write(*,*) "d(",i,",",j,",",k,")=",d(i,j,k)
              IF (d(i,j,k)==0.0D0) THEN
                WRITE(*,'(A,3I5,A,2e15.8)') &
                  'error in direct_solve.f90: main diagonal element equal to zero at ', i, &
                  j, k, '=', d(i,j,k)
                STOP
              ELSE
!         copy diagonals into band strukture
                matrix(ldm-kl,ind) = d(i,j,k)
                IF (ind>1 .AND. i0>1) matrix(ldm-kl+1,ind-1) = c(i,j, &
                  k)
                IF (ind<=nijk-1 .AND. i0>1) matrix(ldm-kl-1,ind+1) &
                  = e(i,j,k)

                IF (ind>i0 .AND. j0>1) matrix(ldm-kl+i0,ind-i0) = b(i, &
                  j,k)
                IF (ind<=nijk-i0 .AND. j0>1) matrix(ldm-kl-i0,ind+i0) &
                  = f(i,j,k)

                IF (ind>i0*j0 .AND. k0>1) matrix(ldm-kl+i0*j0, &
                  ind-i0*j0) = a(i,j,k)
                IF (ind<=nijk-i0*j0 .AND. k0>1) matrix(ldm-kl-i0*j0, &
                  ind+i0*j0) = g(i,j,k)
              END IF

            END DO
          END DO
        END DO

!aw   call dgesv(NIJK,1,matrix,NIJK,pivots,X,NIJK,error)
        CALL dgbsv(nijk,kl,ku,1,matrix,ldm,pivots,x,nijk,error)

! --------------------------------------------------------------

        IF ((linfos(4)>=2) .OR. (error>0)) THEN
          WRITE(*,*) ' PLU: ', error, 'errors was found.'
        END IF

!     print latest, "matrix" is used as tmp. vector
        IF (linfos(4)>=1) THEN
!       MVP
!       [r=matrix]:=[A]x[x]
          CALL s_mvp(i0,j0,k0,x,matrix,a,b,c,d,e,f,g)
!       [r]:=[w]-[A]x[x]
          CALL dscal(nijk,-1.0D0,matrix,1)
          CALL daxpy(nijk,1.0D0,w,1,matrix,1)

!       Ueberpruefung, ob Abbruch
          enough = 0.0D0
          enough2 = 0.0D0

          CALL s_damax(nijk,matrix,enough)
          CALL s_ddot(nijk,matrix,matrix,enough2)
          IF (linfos(4)>=2) THEN
            WRITE(*,'(2(1a,1d20.13))') '  damax(R) =', enough, &
              ', nrm2(R)', dsqrt(enough2)
          ELSE
            WRITE(*,'(1a,1d20.13,1a)') '  nrm2(R) =', dsqrt(enough2), &
              ' (direct)'
          END IF
        END IF

        DEALLOCATE(matrix)
        DEALLOCATE(pivots)

        RETURN
      END
