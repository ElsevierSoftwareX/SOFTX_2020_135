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

!>    @brief gives the level-mode (off/lin/log) of the parameter with index "seedi"
!>    @param[in] seedi parameter component/index
!>    @param[in] ismpl local sample index
!>    @return 0=disabled, 1,2..=enabled
      INTEGER FUNCTION param_enabled(seedi,ismpl)
        use arrays
        use mod_genrl
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j
        INTEGER seedi, s_k, s_u

!     default disabled
        param_enabled = 0

        s_k = seed_para(1,seedi)
        s_u = seed_para(2,seedi)
        IF ((s_k<=lastidx) .AND. (s_k>=firstidx)) THEN
!        parameter units
          param_enabled = opti_props(s_k-firstidx+1,s_u)
        ELSE IF ((s_k<=bc_lastidx) .AND. (s_k>=bc_firstidx)) THEN
!        bc units
          param_enabled = opti_bc(s_k-bc_firstidx+1,s_u)
        ELSE IF (s_k==-1) THEN
!        tp units
!        param_enabled = opti_tp(4, s_u)
          param_enabled = 1
        ELSE IF (s_k==-2) THEN
!        single cell bc
          i = ibc_data(s_u,cbc_bcu)
          j = ibc_data(s_u,cbc_pv)
          IF (i<=0) THEN
            WRITE(*,'(1A,1I5,1A)') &
              'error: active boundary condition (', seedi, &
              '=index) needs to be assigned to a bc-unit!'
            STOP
          END IF
!        take the value from the bc-unit
          param_enabled = opti_bc(j,i)
        ELSE
          WRITE(*,*) &
            'error: "param_enabled" has an incorrect index=', seedi, &
            ' value !'
          STOP
        END IF

        RETURN
      END
