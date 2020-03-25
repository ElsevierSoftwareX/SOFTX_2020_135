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

!>    @brief set parameter value with index "seedi"
!>    @param[in] seedi parameter component/index
!>    @param[in] pvalue value to set
!>    @param[in] ismpl local sample index
      SUBROUTINE set_optip(seedi,pvalue,ismpl)
        use arrays
        use mod_genrl
        IMPLICIT NONE
        integer :: ismpl
        INTEGER seedi, s_k, s_u
        DOUBLE PRECISION pvalue

        IF ((seedi<1) .OR. (seedi>mpara)) THEN
          WRITE(*,*) 'error: in "set_optip", no valid index=',seedi, ' out of range'
          STOP
        END IF

        s_k = seed_para(1,seedi)
        s_u = seed_para(2,seedi)

        IF ((s_k<=lastidx) .AND. (s_k>=firstidx)) THEN
!         parameter units
          propunit(s_u,s_k,ismpl) = pvalue
          call stab_param(propunit(s_u,s_k,ismpl),s_k,s_u)
        ELSE IF ((s_k<=bc_lastidx) .AND. (s_k>=bc_firstidx)) THEN
!        bc units
          propunit(s_u,s_k,ismpl) = pvalue
        ELSE IF (s_k==-1) THEN
!        tp units
!        warning, the second parameter need an addition of 1,
!          because of restricted range [1..2] instead of [1..3]
          bcperiod(opti_tp(1,s_u),opti_tp(2,s_u)+1,opti_tp(3,s_u),ismpl) = pvalue
        ELSE IF (s_k==-2) THEN
!        single cell bc
          dbc_data(s_u,1,ismpl) = pvalue
        ELSE
          WRITE(*,*) 'error: "set_optip" has an incorrect index=',seedi, ' value !'
          STOP
        END IF

        RETURN
      END
