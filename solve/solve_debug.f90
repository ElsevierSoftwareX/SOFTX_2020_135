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

!>    @brief output routine to debug "solve.f90" input, use this instead of "call solve(...)"
!>    @param[in] csf suffix for the output file name
!>    @param[in] x_cur current solution vector [x]
!>    @param[in] ismpl local sample index
      SUBROUTINE solve_debug(csf,x_cur,ismpl)
        use arrays
        use mod_genrl
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k, l
!     residuum
        DOUBLE PRECISION, ALLOCATABLE :: dresid(:)
!     max. normalisation value
        DOUBLE PRECISION, ALLOCATABLE :: dmaxnrm(:)
        DOUBLE PRECISION x_cur(i0,j0,k0,1)
        ! DOUBLE PRECISION dvor
        ! CHARACTER cnach
        character (len=*) :: csf
        LOGICAL test_null
        EXTERNAL test_null
        INTRINSIC max, abs

!     compute "residuum"
        ALLOCATE(dresid(i0*j0*k0))
        ALLOCATE(dmaxnrm(i0*j0*k0))
        CALL s_mvp(i0,j0,k0,x_cur(1,1,1,ismpl),dresid,a(1,1,1,ismpl), &
          b(1,1,1,ismpl),c(1,1,1,ismpl),d(1,1,1,ismpl),e(1,1,1,ismpl), &
          f(1,1,1,ismpl),g(1,1,1,ismpl))
        CALL daxpy(i0*j0*k0,-1.D0,w(1,1,1,ismpl),1,dresid,1)

        l = 0
        DO k = 1, k0
          DO j = 1, j0
            DO i = 1, i0
              l = l + 1
!         dmaxnrm(l) = abs(w(i,j,k,ismpl))
!debug         cnach = 'w'
!debug         dvor = dmaxnrm(l)
!         if (K0.gt.1.and.k.gt.1) dmaxnrm(l) = max(dmaxnrm(l),
!     &     abs(a(i,j,k,ismpl) *x_Cur(i,j,k-1,ismpl)))
!debug         if (dvor.ne.dmaxnrm(l) ) cnach = 'a'
!debug         dvor = dmaxnrm(l)
!         if (J0.gt.1.and.j.gt.1) dmaxnrm(l) = max(dmaxnrm(l),
!     &     abs(b(i,j,k,ismpl) *x_Cur(i,j-1,k,ismpl)))
!debug         if (dvor.ne.dmaxnrm(l) ) cnach = 'b'
!debug         dvor = dmaxnrm(l)
!         if (I0.gt.1.and.i.gt.1) dmaxnrm(l) = max(dmaxnrm(l),
!     &     abs(C(i,j,k,ismpl) *x_Cur(i-1,j,k,ismpl)))
!debug         if (dvor.ne.dmaxnrm(l) ) cnach = 'c'
!debug         dvor = dmaxnrm(l)
!         dmaxnrm(l) = max(dmaxnrm(l),
!     &     abs(d(i,j,k,ismpl) *x_Cur(i,j,k,ismpl)))
!debug         if (dvor.ne.dmaxnrm(l) ) cnach = 'd'
!debug         dvor = dmaxnrm(l)
!         if (I0.gt.1.and.i.lt.I0) dmaxnrm(l) = max(dmaxnrm(l),
!     &     abs(e(i,j,k,ismpl) *x_Cur(i+1,j,k,ismpl)))
!debug         if (dvor.ne.dmaxnrm(l) ) cnach = 'e'
!debug         dvor = dmaxnrm(l)
!         if (J0.gt.1.and.j.lt.J0) dmaxnrm(l) = max(dmaxnrm(l),
!     &     abs(f(i,j,k,ismpl) *x_Cur(i,j+1,k,ismpl)))
!debug         if (dvor.ne.dmaxnrm(l) ) cnach = 'f'
!debug         dvor = dmaxnrm(l)
!         if (K0.gt.1.and.k.lt.K0) dmaxnrm(l) = max(dmaxnrm(l),
!     &     abs(g(i,j,k,ismpl) *x_Cur(i,j,k+1,ismpl)))
!debug         if (dvor.ne.dmaxnrm(l) ) cnach = 'g'
!debug         if (cnach.ne.'d'.and.cnach.ne.'w')
!debug     &     write(*,*) '   ['//cnach//']',l
              IF (bc_mask(l,ismpl)=='+') THEN
!           diagonal dominance
                dmaxnrm(l) = abs(d(i,j,k,ismpl))
                IF ( .NOT. test_null(x_cur(i,j,k, &
                  ismpl))) dmaxnrm(l) = abs(d(i,j,k,ismpl)* &
                  x_cur(i,j,k,ismpl))
              ELSE
                dmaxnrm(l) = 1.D0
              END IF
            END DO
          END DO
        END DO

        OPEN(999,file='system_'//csf//'.dat',status='REPLACE')
        WRITE(999,*) '% w x res err base | a b c d e f g'
        l = 0
        DO k = 1, k0
          DO j = 1, j0
            DO i = 1, i0
              l = l + 1
              WRITE(999,'(5d20.8,2A,7d20.8)') w(i,j,k,ismpl), &
                x_cur(i,j,k,ismpl), dresid(l), dresid(l)/dmaxnrm(l), &
                dmaxnrm(l), ' ' // bc_mask(l,ismpl), ' | ', &
                a(i,j,k,ismpl), b(i,j,k,ismpl), c(i,j,k,ismpl), &
                d(i,j,k,ismpl), e(i,j,k,ismpl), f(i,j,k,ismpl), &
                g(i,j,k,ismpl)
            END DO
          END DO
        END DO
        WRITE(999,*) '% w x res err base | a b c d e f g'
        CLOSE(999)
        DEALLOCATE(dmaxnrm)
        DEALLOCATE(dresid)
        RETURN
      END
