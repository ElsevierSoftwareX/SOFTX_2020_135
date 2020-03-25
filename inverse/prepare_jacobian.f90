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

!>    @brief prepare (log) apriori and bayes values
!>    @param[in] ismpl local sample index
      SUBROUTINE prepare_jacobian(ismpl)
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
        INTEGER pe
        DOUBLE PRECISION get_optia
        INTEGER param_enabled
        EXTERNAL get_optia, param_enabled


        IF (linfos(2)>=1) WRITE(*,*) ' ... calculate (lin/log) values'
!       opti_*    : 0,1,(2) = off, lin, log
        DO k = 1, mpara
          pe = param_enabled(k,ismpl)
          IF (pe==2) THEN
            apri_input(k) = log(get_optia(k,ismpl))
            linlog_input(k) = log(main_input(k,ismpl))
          ELSE IF (pe==1) THEN
            apri_input(k) = get_optia(k,ismpl)
            linlog_input(k) = main_input(k,ismpl)
          END IF
        END DO
        RETURN
      END
