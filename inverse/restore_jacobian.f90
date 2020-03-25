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

!>    @brief restore (log) apriori and bayes values
!>    @param[in] p_k parameter vector for update
!>    @param[in] ismpl local sample index
      SUBROUTINE restore_jacobian(p_k,ismpl)
#ifndef JACOBI_FREE
         use arrays
#ifdef AD         
        use g_arrays
#endif
#ifdef AD_RM
        use arrays_ad
#endif
        use mod_genrl
        use mod_inverse
        use mod_linfos
        IMPLICIT NONE
        integer :: ismpl
        integer :: k
        DOUBLE PRECISION p_k(mpara)
        INTEGER pe
        INTEGER param_enabled
        EXTERNAL param_enabled

        IF (linfos(2)>=1) WRITE(*,*) ' ... calculate (exp) values'
!     update parameter in physiCal spaCe
        DO k = 1, mpara
          pe = param_enabled(k,ismpl)
          IF (pe==2) THEN
!           log
            main_input(k,ismpl) = exp(p_k(k))
          ELSE IF (pe==1) THEN
!           lin
            main_input(k,ismpl) = p_k(k)
          END IF
        END DO
!
#endif
        RETURN
      END
