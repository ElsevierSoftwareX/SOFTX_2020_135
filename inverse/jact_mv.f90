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

!>    @brief transposed Jacobi matrix vector product (full Jacobi)
!>    @param[in] xvec vector [x]
!>    @param[out] yvec vector [y], result of [y] = Jacobi^T*[x]
!>    @param[in] ismpl local sample index
      SUBROUTINE jact_mv(xvec,yvec,ismpl)
        use arrays
        use mod_time
        use mod_data
        use mod_inverse
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j
        DOUBLE PRECISION xvec(ndata), yvec(mpara)
        DOUBLE PRECISION, allocatable :: lxvec(:), lyvec(:), ryvec(:)

!
!       zeros needed in "yvec"
        call set_dval(mpara,0.0d0,yvec)
!
#ifndef JACOBI_FREE
! ------------- full matrix
        DO j = 1, mpara
          DO i = 1, ndata
            yvec(j) = yvec(j) +jac(i,j)*xvec(i)
          END DO
        END DO
! -------------
#else
! ------------- matrix free
#ifdef DEBUG
!       compute adjoint by forward mode
        WRITE(*,'(1A)') '<!> :   -> Jacobi^T vector product via AD forward mode'
        allocate(lxvec(mpara), lyvec(ndata), ryvec(mpara))
        do i = 1, mpara
!         seeding
          call set_dval(mpara,0.0d0,lxvec)
          lxvec(i) = 1.0d0
!         compute via forward mode
          call jac_mv(lxvec,lyvec,ismpl)
!         assimilate adjoint
          ryvec(i) = 0.d0
          do j = 1, ndata
            ryvec(i) = ryvec(i) +lyvec(j)*xvec(j)
          enddo
        enddo
#endif
!
!       init all AD variables with zeros
        CALL initzero_ad(ismpl)
!       sample init
        CALL prepare_realisation(ismpl)
!
!       one derived function step
        CALL forward_compute_ad(main_input(1,ismpl), yvec, &
          main_output(1,ismpl), xvec, &
          simtime_0,max_simtime,iter_inv,1,ismpl)
        WRITE(*,*) ' [I] : -> Jacobi^T vector product via AD reverse mode'
!
#ifdef DEBUG
        WRITE(*,'(1A,99G24.16)') '<!> :   forward mode =', ryvec
        deallocate(lxvec, lyvec, ryvec)
        WRITE(*,'(1A,99G24.16)') '<!> :   reverse mode =', yvec
#endif
! -------------
#endif
!
        RETURN
      END
