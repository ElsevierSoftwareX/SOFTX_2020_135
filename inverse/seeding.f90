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

!>    @brief generate a seeding vector (column of jacobi matrix) for the seeding index "seedi"
!>    @param[in] seedi seeding component/index
!>    @param[out] dseed seeding vector
!>    @param[in] ismpl local sample index
!>    @details
!>      kx-,ky-,kz-,lx-,ly-,lz-,q-(,rc-)units and others are seeded\n
      SUBROUTINE seeding(seedi, dseed, ismpl)
        use arrays
#ifdef AD 
        use g_arrays
#endif
#ifdef AD_RM 
        use arrays_ad
#endif
        use mod_data
        use mod_genrl
        use mod_time
        use mod_inverse
        use mod_linfos
        IMPLICIT NONE
        integer :: ismpl
        INTEGER seedi, s_k, s_u, s_level, param_enabled
        CHARACTER*16 pname
        DOUBLE PRECISION dseed(mpara)
        EXTERNAL param_enabled

!
#ifndef AD_RM
        IF ((seedi<1) .OR. (seedi>mpara)) THEN
          WRITE(*,*) 'error : in "seeding", no valid "seeding=', &
            seedi, '"'
          STOP
        END IF
!     clear/zero-init of all variables needed before, supported by "g_initzero"
        s_k = seed_para(1,seedi)
        s_u = seed_para(2,seedi)
        s_level = param_enabled(seedi,ismpl)
        CALL param_name(seedi,pname,ismpl)
        IF ((s_k<=lastidx) .AND. (s_k>=firstidx)) THEN
!        parameter units
          IF (linfos(2)>=1) THEN
            WRITE(*,'(3A,I2,A,I2,A,I5,A,I5)') '  Seeding: ', pname, &
              ', component= ', s_k - firstidx + 1, '/', &
              lastidx - firstidx + 1, ', unit= ', s_u, '/', maxunits
          END IF
        ELSE IF ((s_k<=bc_lastidx) .AND. (s_k>=bc_firstidx)) THEN
!        bc units
          IF (linfos(2)>=1) THEN
            WRITE(*,'(3A,I5,A,I5)') '  Seeding: ', pname, ', unit= ', &
              s_u, '/', bc_maxunits
          END IF
        ELSE IF (s_k==-1) THEN
!        tp units
          IF (linfos(2)>=1) THEN
            WRITE(*,'(3A,I5,A,I5)') '  Seeding: ', pname, ', unit= ', &
              s_u, '/', mpara_tp
          END IF
        ELSE IF (s_k==-2) THEN
!        single cell bc
          IF (linfos(2)>=1) THEN
            WRITE(*,'(3A,I5,A,I5)') '  Seeding: ', pname, ', cell= ', &
              s_u, '/', nbc_data
          END IF
        ELSE
          WRITE(*,*) 'error: incorrect seeding=', seedi, &
            ' value in "seeding.f"!'
          STOP
        END IF
!
        CALL set_dval(mpara,0.0d0,dseed)
!
        IF (s_level==2) THEN
!        logarithm
          dseed(seedi) = main_input(seedi,ismpl)
        ELSE IF (s_level==1) THEN
!        linear
          dseed(seedi) = 1.0d0
        END IF
#else
! AD-Reverse seeding
        if ((seedi<1) .or. (seedi>ndata)) then
           write(*,*) 'error: in "seeding_ad", no valid "seeding=',seedi,'"'
           stop
        end if
        call set_dval(ndata,0.0d0,dseed)
        dseed(seedi)=1.0d0
#endif
!
!        write(*,*) "Seed=",dseed
        RETURN
      END

!>    @brief gives the parameter apriori value for the seeding index "seedi"
!>    @param[in] seedi seeding component/index
!>    @param[in] ismpl local sample index
!>    @return parameter value
      DOUBLE PRECISION FUNCTION get_optia(seedi,ismpl)
        use arrays
        use mod_genrl
        use mod_inverse
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j
        INTEGER seedi, s_k, s_u

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

!>    @brief gives the parameter measurement error for the seeding index "seedi"
!>    @param[in] seedi seeding component/index
!>    @param[in] ismpl local sample index
!>    @return parameter value
      DOUBLE PRECISION FUNCTION get_optid(seedi,ismpl)
        use arrays
        use mod_genrl
        use mod_inverse
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j
        INTEGER seedi, s_k, s_u

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

!>    @brief set parameter numeric error for the seeding index "seedi"
!>    @param[in] seedi seeding component/index
!>    @param[in] value value to set
!>    @param[in] ismpl local sample index
      SUBROUTINE set_optie(seedi,value,ismpl)
        use arrays
        use mod_genrl
        use mod_inverse
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k
        INTEGER seedi, s_k, s_u
        DOUBLE PRECISION value

        IF ((seedi<1) .OR. (seedi>mpara)) THEN
          WRITE(*,*) 'error : in "set_optie", no valid "seeding=', &
            seedi, '"'
          STOP
        END IF

        s_k = seed_para(1,seedi)
        s_u = seed_para(2,seedi)

        IF ((s_k<=lastidx) .AND. (s_k>=firstidx)) THEN
!        parameter units
          e_propunit(s_u,s_k) = value
        ELSE IF ((s_k<=bc_lastidx) .AND. (s_k>=bc_firstidx)) THEN
!        bc units
          e_propunit(s_u,s_k) = value
        ELSE IF (s_k==-1) THEN
!        tp units
          e_bcperiod(opti_tp(1,s_u),opti_tp(2,s_u),opti_tp(3,s_u)) &
            = value
        ELSE IF (s_k==-2) THEN
!        single cell bc, not supported !!!
!        e_dbC_data(s_u,1) = value
        ELSE
          WRITE(*,*) 'error: "set_optie" has an incorrect seeding=', &
            seedi, ' value !'
          STOP
        END IF

        RETURN
      END
