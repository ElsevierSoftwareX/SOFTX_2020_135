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

!>    @brief initialisation routine for reinjection
!>    @param[in] ismpl local sample index
      SUBROUTINE user_init(ismpl)
        use arrays
        use mod_genrl
        use mod_linfos
        IMPLICIT NONE
        integer :: ismpl
        integer :: i

        IF (linfos(3)>=2) WRITE(*,*) ' ... user_init (REINJECT-3D)'
!     avoid side effects
        DO i = 1, nbc_data
!        init from the first set
          dbc_dataold(i) = dbc_data(i,1,1)
        END DO
        RETURN
      END

!>    @brief reinjection routine
!>    @param[in] ismpl local sample index
!>    @details
!>injection 1   13       8         5\n
!>injection 2   13       8        13\n
!>injection 3   13       8        21\n
!>production 1  13       18        5\n
!>production 2  13       18       13\n
!>production 3  13       18       21\n
      SUBROUTINE calc_user(ismpl)
        use arrays
        use mod_genrl
        use mod_time
        use mod_linfos
        use mod_wells3d
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j


        DOUBLE PRECISION tempc, pumpc, tempt, deltt
        LOGICAL lbc_found, lbc_found1, lbc_found2
        DOUBLE PRECISION deltat
        EXTERNAL deltat
        INTRINSIC dabs


        IF (linfos(3)>=2) WRITE(*,*) ' ... calc_user (REINJECT-3D)'

!VR --- special case, danger !!! ---

        deltt = deltat(simtime(ismpl),ismpl)
        IF ((simtime(ismpl)+deltt)/tunit>=stoptime) THEN
          tempc = 0.D0
          pumpc = 0.D0
          tempt = 0.D0
!       for all productions
          DO j = 1, num_pro
            lbc_found = .FALSE.
!          search in all boundary condition
            DO i = 1, nbc_data
              IF (ibc_data(i,cbc_i)==ipro(j) .AND. &
                  ibc_data(i,cbc_j)==jpro(j) .AND. &
                  ibc_data(i,cbc_k)==kpro(j) .AND. &
                  ibc_data(i,cbc_pv)==pv_head .AND. &
                  ibc_data(i,cbc_bt)==bt_neum) THEN
                IF (dbc_data(i,1,ismpl)>=0.D0) THEN
                  WRITE(*,'(2A,3I8,1A)') 'error: HEAD Neumann-boundary&
                    & point for production', &
                    ' has a value >= 0.0 at [', ipro(j), jpro(j), &
                    kpro(j), '] !'
                  STOP
                END IF
                lbc_found = .TRUE.
                tempc = tempc + conc(ipro(j),jpro(j),kpro(j),1,ismpl)* &
                  dabs(dbc_data(i,1,ismpl))
                tempt = tempt + temp(ipro(j),jpro(j),kpro(j),ismpl)* &
                  dabs(dbc_data(i,1,ismpl))
                pumpc = pumpc + dabs(dbc_data(i,1,ismpl))
                WRITE(*,'(1A,3(e12.4))') '!!!!! GPK2 center ', &
                  conc(ipro(j),jpro(j),kpro(j),1,ismpl), &
                  temp(ipro(j),jpro(j),kpro(j),ismpl), &
                  dbc_data(i,1,ismpl)
              END IF
            END DO
            IF ( .NOT. lbc_found) THEN
              WRITE(*,'(1A,3I8,1A)') 'error: no HEAD Neumann-boundary &
                &point found for production [', ipro(j), jpro(j), &
                kpro(j), '] !'
              STOP
            END IF
          END DO

          IF (pumpc/=0.D0) THEN
            tempc = tempc/pumpc
            tempt = tempt/pumpc
          END IF

          WRITE(*,'(4(1A,1e12.4))') '!!!!! reinjection rate = ', &
            pumpc, ' l/s,   concentration = ', tempc, &
            ' mmol/l, temperature = ', tempt, ' C, time =', &
            simtime(ismpl)

!       for all sources
          DO j = 1, num_in
            lbc_found = .FALSE.
            lbc_found1 = .FALSE.
            lbc_found2 = .FALSE.
!          search in all boundary condition
            DO i = 1, nbc_data
              IF (ibc_data(i,cbc_i)==iin(j) .AND. &
                  ibc_data(i,cbc_j)==jin(j) .AND. &
                  ibc_data(i,cbc_k)==kin(j) .AND. &
                  ibc_data(i,cbc_pv)==pv_conc .AND. &
                  ibc_data(i,cbc_bt)==bt_diri .AND. &
                  ibc_data(i,cbc_si)==1) THEN
                dbc_data(i,1,ismpl) = tempc
                lbc_found = .TRUE.
              END IF
              IF (ibc_data(i,cbc_i)==iin(j) .AND. &
                  ibc_data(i,cbc_j)==jin(j) .AND. &
                  ibc_data(i,cbc_k)==kin(j) .AND. &
                  ibc_data(i,cbc_pv)==pv_temp .AND. &
                  ibc_data(i,cbc_bt)==bt_diri) THEN
                dbc_data(i,1,ismpl) = tempt
                lbc_found1 = .TRUE.
              END IF
              IF (ibc_data(i,cbc_i)==iin(j) .AND. &
                  ibc_data(i,cbc_j)==jin(j) .AND. &
                  ibc_data(i,cbc_k)==kin(j) .AND. &
                  ibc_data(i,cbc_pv)==pv_head .AND. &
                  ibc_data(i,cbc_bt)==bt_neum) THEN
                dbc_data(i,1,ismpl) = pumpc
                lbc_found2 = .TRUE.
              END IF
            END DO
            IF ( .NOT. lbc_found) THEN
              WRITE(*,'(1A,3I8,1A)') 'error: no CONC Dirichlet-boundar&
                &y point found for injection [', iin(j), jin(j), &
                kin(j), '] !'
              STOP
            END IF
            IF ( .NOT. lbc_found1) THEN
              WRITE(*,'(1A,3I8,1A)') 'error: no TEMP Dirichlet-boundar&
                &y point found for injection [', iin(j), jin(j), &
                kin(j), '] !'
              STOP
            END IF
            IF ( .NOT. lbc_found2) THEN
              WRITE(*,'(1A,3I8,1A)') 'error: no HEAD Dirichlet-boundar&
                &y point found for injection [', iin(j), jin(j), &
                kin(j), '] !'
              STOP
            END IF
          END DO
        ELSE
!       reset the initial values
          DO i = 1, nbc_data
            dbc_data(i,1,ismpl) = dbc_dataold(i)
          END DO
        END IF

!VR --- special case, danger !!! ---

        RETURN
      END
