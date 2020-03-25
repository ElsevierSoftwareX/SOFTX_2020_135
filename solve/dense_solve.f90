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

!>    @brief linear system solver [x]:=[A^-1]x[b], dense matrix (quadratic)
!>    @param[in] NIJK system size
!>    @param[in] NRHS number of right sides/solutions
!>    @param[in] A system matrix
!>    @param[in] B rigth side
!>    @param[out] X solution
!#>   @details
      SUBROUTINE dense_solve(nijk,nrhs,a,b,x)
        use mod_linfos
        IMPLICIT NONE
        INTEGER nijk, nrhs, isnull
        INTEGER, ALLOCATABLE :: pivots(:)
        DOUBLE PRECISION, ALLOCATABLE :: tmp(:), tmpm(:,:)
        INTEGER error
        DOUBLE PRECISION a(nijk,nijk), b(nijk,nrhs)
        DOUBLE PRECISION x(nijk,nrhs), enough, enough2


        IF ((8.0D0*dble(2*nijk+nijk*nijk))>2.0D9) THEN
          WRITE(*,'(A,A,F9.4,A)') &
            'Error: array to large for direct dense', &
            ' matrix solver (', 8.0D0*dble(2*nijk+nijk*nijk)/1024.0D0 &
            **3, ' GBytes).'
          STOP
        END IF

        ALLOCATE(pivots(nijk))
        ALLOCATE(tmp(nijk))
        ALLOCATE(tmpm(nijk,nijk))
        error = 0
        CALL set_ival(nijk,0,pivots)

!     right side to [X]
        CALL dcopy(nijk*nrhs,b,1,x,1)
!     copy matrix
        CALL dcopy(nijk*nijk,a,1,tmpm,1)

! --------------------------------------------------------------

!     solve dense matrix
        CALL dgesv(nijk,nrhs,a,nijk,pivots,x,nijk,error)

        IF (error>0) THEN
          WRITE(*,*) ' dense PLU: ', error, 'errors were found.'
        ELSE IF (linfos(4)>=1) THEN
!aw        write(*,*) ' dense PLU finished.'
        END IF

        IF ((linfos(4)>=1) .AND. (nrhs==1)) THEN
!       MVP
!       [tmp]==[r]:=[A]x[x]
          CALL dgemv('N',nijk,nijk,1.0D0,tmpm,nijk,x,1,0.0D0,tmp,1)

          enough = 0.0D0
          isnull = 0
          enough2 = 0.0D0
          CALL s_ddot(nijk,b,b,enough)
          CALL test_zero(enough,1,isnull)
          IF (isnull/=1) THEN
            CALL s_ddot(nijk,tmp,tmp,enough2)
            enough2 = abs(enough-enough2)/(enough)
          END IF

!       [r]:=[w]-[A]x[x]
          CALL dscal(nijk,-1.0D0,tmp,1)
          CALL daxpy(nijk,1.0D0,b,1,tmp,1)

          enough = 0.0D0
          CALL s_ddot(nijk,tmp,tmp,enough)
          WRITE(*,'(2(a,e12.5))') '   PLU: nrm2(R)', dsqrt(enough), &
            ' rel=', enough2
        END IF

!     restore matrix (only needed for debug)
        CALL dcopy(nijk*nijk,tmpm,1,a,1)

        DEALLOCATE(pivots)
        DEALLOCATE(tmp)
        DEALLOCATE(tmpm)

! --------------------------------------------------------------

        RETURN
      END
