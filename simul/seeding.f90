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

!> @brief gives the parameter value for the seeding index "seedi"
!> @param[in] seedi seeding component/index
!> @param[in] ismpl local sample index
!> @return parameter value
      DOUBLE PRECISION FUNCTION get_optia(seedi,ismpl)
        use arrays
        use mod_genrl
        use mod_simul
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j
        integer :: seedi, s_k, s_u


        IF ((seedi<1) .OR. (seedi>mpara)) THEN
          WRITE(*,*) 'error : in "get_optia", no valid "seeding=', &
            seedi, '"'
          STOP
        END IF

        s_k = seed_para(1,seedi)
        s_u = seed_para(2,seedi)

        IF ((s_k<=lastidx) .AND. (s_k>=firstidx)) THEN
!        parameter units
          get_optia = a_propunit(s_u,s_k)
        ELSE IF ((s_k<=bc_lastidx) .AND. (s_k>=bc_firstidx)) THEN
!        bc units
          get_optia = a_propunit(s_u,s_k)
        ELSE IF (s_k==-1) THEN
!        tp units
          get_optia = a_bcperiod(opti_tp(1,s_u),opti_tp(2,s_u), &
            opti_tp(3,s_u))
        ELSE IF (s_k==-2) THEN
!        single cell bc
          i = ibc_data(s_u,cbc_bcu)
          j = ibc_data(s_u,cbc_pv)
          IF (i<=0) THEN
            WRITE(*,'(1A,1I5,1A)') &
              'error: active boundary condition (', seedi, &
              '=seeding) needs to be assigned to a bc-unit!'
            STOP
          END IF
!        take the value from the bc-unit
          get_optia = a_propunit(i,j+bc_firstidx-1)
        ELSE
          WRITE(*,*) 'error: "get_optia" has an incorrect seeding=', &
            seedi, ' value !'
          STOP
        END IF

        RETURN
      END

!> @brief gives the standard deviation for the seeding index "seedi"
!> @param[in] seedi seeding component/index
!> @param[in] ismpl local sample index
!> @return standard deviation
      DOUBLE PRECISION FUNCTION get_optid(seedi,ismpl)
        use arrays
        use mod_genrl
        use mod_simul
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j
        integer :: seedi, s_k, s_u


        IF ((seedi<1) .OR. (seedi>mpara)) THEN
          WRITE(*,*) 'error : in "get_optid", no valid "seeding=', &
            seedi, '"'
          STOP
        END IF

        s_k = seed_para(1,seedi)
        s_u = seed_para(2,seedi)

        IF ((s_k<=lastidx) .AND. (s_k>=firstidx)) THEN
!        parameter units
          get_optid = d_propunit(s_u,s_k)
        ELSE IF ((s_k<=bc_lastidx) .AND. (s_k>=bc_firstidx)) THEN
!        bc units
          get_optid = d_propunit(s_u,s_k)
        ELSE IF (s_k==-1) THEN
!        tp units
          get_optid = d_bcperiod(opti_tp(1,s_u),opti_tp(2,s_u), &
            opti_tp(3,s_u))
        ELSE IF (s_k==-2) THEN
!        single cell bc
          i = ibc_data(s_u,cbc_bcu)
          j = ibc_data(s_u,cbc_pv)
          IF (i<=0) THEN
            WRITE(*,'(1A,1I5,1A)') &
              'error: active boundary condition (', seedi, &
              '=seeding) needs to be assigned to a bc-unit!'
            STOP
          END IF
!        take the value from the bc-unit
          get_optid = d_propunit(i,j+bc_firstidx-1)
        ELSE
          WRITE(*,*) 'error: "get_optid" has an incorrect seeding=', &
            seedi, ' value !'
          STOP
        END IF

        RETURN
      END

