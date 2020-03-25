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

!>    @brief updates the step from p_k to p_k+1
!>    @param[in] p_k parameter vector from current step k
!>    @param[out] p_k parameter vector for step k+1
!>    @param[in] p_delta parameter delta to the next step k+1
!>    @param[in] alpha weigthing-factor for the delta-update
!>    @param[in] ismpl local sample index
      SUBROUTINE update_param(p_k,p_delta,alpha,ismpl)
        use arrays
        use mod_genrl
        use mod_inverse
        use mod_data
        use mod_linfos
        IMPLICIT NONE
        integer :: ismpl
        integer :: i
        INTEGER p_l, param_enabled
        CHARACTER*16 pname
        EXTERNAL param_enabled

        DOUBLE PRECISION alpha
        DOUBLE PRECISION p_k(mpara), p_delta(mpara)
!     max. factor
        DOUBLE PRECISION maxf, e_maxf, apdi, mdelta
!      ### 'log'    [ 1/10pk < p_k < 10pk ]       ###
!      ### 'linear' [ -pk < 0 < p_k < 2pk < 3pk ] ###
!       -> 'maxf' must be greater than 1.0d0 !!!
        PARAMETER (maxf=2.0D0,e_maxf=3.0D0)
        INTRINSIC log, abs
        LOGICAL test_null
        EXTERNAL test_null

!     p_k+1 = p_k + alpha*p_delta
        DO i = 1, mpara
          p_l = param_enabled(i,ismpl)
          apdi = alpha*p_delta(i)
!       is it enabled (level>0) ?
          IF (p_l==2) THEN
!      'log'    [ 1/10pk < p_k < 10pk ]
            mdelta = log(e_maxf)
          ELSE IF (p_l==1) THEN
!      'linear' [ -pk < 0 < p_k < 2pk < 3pk ]
            mdelta = abs(maxf*p_k(i))
            IF (test_null(p_k(i))) mdelta = 1.D0
          ELSE
!     'off', the value should not be changed
!     p_k(i) = p_k(i)
            GO TO 100
          END IF
          IF (apdi>mdelta) THEN
            CALL param_name(i,pname,ismpl)
            WRITE(*,'(3A,2(1e16.8,1A))') 'warning: (+) cut ', pname, &
              '=', p_k(i), ' to ', p_k(i) + mdelta, ' !'
            p_k(i) = p_k(i) + mdelta
          ELSE IF (apdi<-mdelta) THEN
            CALL param_name(i,pname,ismpl)
            WRITE(*,'(3A,2(1e16.8,1A))') 'warning: (-) cut ', pname, &
              '=', p_k(i), ' to ', p_k(i) - mdelta, ' !'
            p_k(i) = p_k(i) - mdelta
          ELSE
            p_k(i) = p_k(i) + apdi
          END IF
100       CONTINUE
        END DO
        RETURN
      END
