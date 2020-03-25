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

!> @brief wrapps the *SIM-library call and all init and finalising
!> @param[in,out] dv_param parameter vector; initial values at begin, but returns stochastic modification
!> @param[in] realz realisation/ensemble index number
!> @param[in] ismpl local sample index
!> @details
!> allocate temporary variables, calls the *SIM-library (SGSIM, VISIM)\n
!> and setup the parameters with the stochastic errors\n
      SUBROUTINE simul_wrapper(dv_param,ismpl,realz)
        use arrays
        USE simul_arrays
        use mod_genrl
        use mod_genrlc
        use mod_linfos
        use mod_simul
        IMPLICIT NONE
        integer  i,j,k,ismpl
        DOUBLE PRECISION val
!     s_u & s_k: seeding index
        integer :: s_k, s_u, ipara
!     realisation
        integer :: realz, ijk
        integer, dimension (14) :: fileh
!     converts mpara-idx into ijk-idx
        integer, ALLOCATABLE :: ijk_idx(:,:)
        integer, ALLOCATABLE :: ijk_bcidx(:,:)
!       parameter vector
        DOUBLE PRECISION dv_param(mpara)
!
        DOUBLE PRECISION sig, get_optid
        integer :: ll, param_enabled
        LOGICAL test_option
        character (len=255) :: sfile_name
        character (len=8) :: snumber
        integer :: i1,i2,i3,i4
        EXTERNAL get_optid, param_enabled, test_option
        INTRINSIC trim


        IF (max_gpara==0) RETURN
        IF (test_option('-param=postcomp')) THEN
!         generate file name
          CALL chln(project,i1,i2)
          IF (realz>0) THEN
            WRITE(snumber,'(1I7)') realz
          ELSE
            WRITE(*,'(1A)') 'error: realisation number out of range in "simul_wrapper"!'
            STOP
          END IF
          CALL chln(snumber,i3,i4)
          sfile_name = project(i1:i2) // '_postcomp' // &
            '_' // snumber(i3:i4) // '.h5'
          WRITE(*,'(3A)') ' [R] : reading initial properties from "',trim(sfile_name),'"'
!$OMP     CRITICAL
            CALL open_hdf5(' ')
            CALL read2_hdf5(nunits,nprop,propunit(1,1,ismpl),'props_full',sfile_name)
            CALL close_hdf5()
!$OMP     END CRITICAL
          RETURN
        END IF
!     ## needs only to be done by one thread,
!        but there is no need for a barrier since
!        is it constant for all threads !!!
        maxx = i0
        maxy = j0
        maxz = k0
        ALLOCATE(ijk_idx(3,nunits))
        ALLOCATE(ijk_bcidx(3,nunits))
        CALL set_ival(3*nunits,0,ijk_idx)
        CALL set_ival(3*nunits,0,ijk_bcidx)

        DO i3 = 1, k0
          DO i2 = 1, j0
            DO i1 = 1, i0
!              [nunits]-index
              i = uindex(i1,i2,i3)
!              if (ijk_idx(1,i).eq.0.or.ijk_idx(2,i).eq.0.or.ijk_idx(3,i
              IF (ijk_idx(1,i)==0) THEN
!                 setup first ijk-index for each [nunits]
                ijk_idx(1,i) = i1
                ijk_idx(2,i) = i2
                ijk_idx(3,i) = i3
              END IF
            END DO
          END DO
        END DO
        DO i = 1, bc_maxunits
          DO j = 1, nbc_data
            IF (ibc_data(j,cbc_bcu)==i) THEN
              IF (ijk_bcidx(1,i)==0) THEN
!                 setup first ijk-index for each [*bc_data]
                ijk_bcidx(1,i) = ibc_data(j,cbc_i)
                ijk_bcidx(2,i) = ibc_data(j,cbc_j)
                ijk_bcidx(3,i) = ibc_data(j,cbc_k)
                GO TO 100
              END IF
            END IF
          END DO
100       CONTINUE
        END DO

!     init dummy file handler
        DO i = 1, 14
          CALL omp_new_file_handler(fileh(i),i)
        END DO

!     for each parameter/unit group of interesst
        DO ipara = 1, max_gpara
          sm_i0 = i0
          sm_j0 = j0
          sm_k0 = k0
!        start parameter-initialisation/random-seed
          CALL simul(ismpl,realz,ipara)
!        do sanity check only for the first simulation, others are similar
          IF (realz==1) CALL grid_adapt_check(fnpara(ipara))
!        update the property
          DO i = 1, mpara
!          in the group ?
            IF (gpara(i)==ipara) THEN
              s_k = seed_para(1,i)
              s_u = seed_para(2,i)
              IF ((s_k<=lastidx) .AND. (s_k>=firstidx)) THEN
                i1 = ijk_idx(1,s_u)
                i2 = ijk_idx(2,s_u)
                i3 = ijk_idx(3,s_u)
              ELSE IF ((s_k<=bc_lastidx) .AND. (s_k>=bc_firstidx)) &
                  THEN
                i1 = ijk_bcidx(1,s_u)
                i2 = ijk_bcidx(2,s_u)
                i3 = ijk_bcidx(3,s_u)
              ELSE IF (s_k==-1) THEN
                WRITE(*,'(1A)') &
                  'error: active time periods not supported yet!'
                STOP
              ELSE IF (s_k==-2) THEN
                i1 = ibc_data(s_u,cbc_i)
                i2 = ibc_data(s_u,cbc_j)
                i3 = ibc_data(s_u,cbc_k)
              ELSE
                WRITE(*,'(1A)') &
                  'error: "seed_para" out of range ...!'
                STOP
              END IF
!            sanity check
              IF (i1==0 .OR. i2==0 .OR. i3==0) THEN
                WRITE(*,'(2A,2I8,1A)') &
                  'error: no "uindex" connected with the active', &
                  ' parameter: [', s_u, s_k, '], see section "# parameter group" !'
                STOP
              END IF
!AW          ijk = i1 +(i2-1)*I0 +(i3-1)*I0*J0
              CALL grid_adapt_ijk(i1,i2,i3,ijk)
!            modify the value
!CV          val = dv_param(i)
              sig = get_optid(i,ismpl)
              ll = param_enabled(i,ismpl)
              IF (ll==1) THEN
!               lin
!CV             val = val +sig*simout(ijk)
                val = sig*simout(ijk)
              ELSE IF (ll==2) THEN
!               log
!CV             val = exp(log(val) +sig*simout(ijk))
                val = exp(sig*simout(ijk))
              ELSE
                WRITE(*,'(1A)') 'error: wrong lin/log specification for parameters!'
                STOP
              END IF
!             update new value -> main vector variable
              dv_param(i) = val
!             update new value -> into current parameter space for later output
              CALL set_optip(i,val,ismpl)
              IF ((mpara<=1000 .OR. mod(i,1000)==0) .AND. linfos(2)>=1) WRITE(*,'(1A,4I8,1e16.7)') &
                '    # seed-idx,max idx,seed-comp.,seed-unit,value =', i, mpara, s_k, s_u, val
            END IF
          END DO
        END DO

!     clear dummy file handler
        DO i = 1, 14
          CALL omp_del_file_handler(fileh(i))
        END DO
!     free local array
        DEALLOCATE(ijk_bcidx)
        DEALLOCATE(ijk_idx)

        RETURN
      END
