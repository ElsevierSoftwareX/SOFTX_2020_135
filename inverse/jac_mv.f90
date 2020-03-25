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

!>    @brief Jacobi matrix vector product (full Jacobi)
!>    @param[in] xvec vector [x]
!>    @param[out] yvec vector [y], result of [y] = Jacobi*[x]
!>    @param[in] ismpl local sample index
      SUBROUTINE jac_mv(xvec,yvec,ismpl)
        use arrays
        use mod_time
        use mod_inverse
        use mod_data
        use mod_OMP_TOOLS
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j
        DOUBLE PRECISION xvec(mpara), yvec(ndata)

!
        CALL omp_ordered_create(1)
!
!       zeros needed in "yvec"
        call set_dval(ndata,0.0d0,yvec)
!
#ifndef JACOBI_FREE
! ------------- full matrix
        DO j = 1, ndata
          DO i = 1, mpara
            yvec(j) = yvec(j) +jac(j,i)*xvec(i)
          END DO
        END DO
! -------------
#else
! ------------- matrix free
!       init all AD variables with zeros
        CALL g_initzero(ismpl)
!       sample init
        CALL prepare_realisation(ismpl)
!
!       one derived function step
        CALL g_forward_compute(main_input(1,ismpl), xvec, &
          main_output(1,ismpl), yvec, &
          simtime_0,max_simtime,iter_inv,1,ismpl)
! -------------
#endif
!
        CALL omp_ordered_delete()
!
        RETURN
      END
