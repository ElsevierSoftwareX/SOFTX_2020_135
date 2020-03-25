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

!>    @brief parameter output for each inverse iteration
!>    @param[in] index_i inverse iteration counter
!>    @param[in] ismpl local sample index
      SUBROUTINE write_parameter(index_i,ismpl)
        use arrays
        use mod_genrl
        use mod_genrlc
        use mod_data
        use mod_inverse
        use mod_linfos
        use mod_OMP_TOOLS

        INCLUDE 'OMP_TOOLS.inc'
	
        INTEGER i, j, k, ismpl
        INTEGER i1, i2, i3, i4, type, index_i, lblank, s_u, lout

        CHARACTER filename*80, snumber*8
        DOUBLE PRECISION tmp1, tmp2, reg_func
        EXTERNAL lblank, reg_func
        CHARACTER*4 tp_albe(2)
        DATA tp_albe/'alfa', 'beta'/


        IF ((mpara<=0) .OR. write_disable .OR. ( .NOT. write_param)) &
          RETURN
#ifdef NOOUT
        RETURN
#endif

        CALL omp_new_file_handler(lout,1)

        CALL chln(project,i1,i2)
        IF (index_i>=0) THEN
          WRITE(snumber,'(1I7)') index_i
        ELSE IF (index_i==-1) THEN
          WRITE(snumber,'(A8)') 'final'
        ELSE IF (index_i==-2) THEN
          WRITE(snumber,'(A8)') 'debug'
        ELSE IF (index_i==-3) THEN
          WRITE(snumber,'(A8)') 'ens_mean'
        ELSE IF (index_i==-4) THEN
          WRITE(snumber,'(A8)') 'mean'
        ELSE IF (index_i==-5) THEN
          WRITE(snumber,'(A8)') 'ens_mean'
        END IF
        CALL chln(snumber,i3,i4)

        filename = project(i1:i2) // '_P_' // snumber(i3:i4) // '.dat'

        OPEN(lout,file=filename,status='unknown',blank='null')

        IF (linfos(3)>=1) THEN
          WRITE(*,'(3A)') '  [W] : Parameter to "', &
            filename(1:lblank(filename)), '"'
        END IF


        WRITE(lout,'(A,i4)')    '% regularization type = ', &
          reg_type
        WRITE(lout,'(A,g14.5)') '% regularization parameter (input) = ', &
          reg_func()
        IF(reg_type > 1) THEN
           WRITE(lout,'(A,i4)') '% regularization parameter scheduled correponding to ', &
          reg_type
        END IF
!     parameter units
        WRITE(lout,'(2A)') key_char//' units (aposteriori), iter =', snumber
        DO i = 1, maxunits
          WRITE(lout,'(1X,'//c_npropunit//'(e14.5,1X),1A,1I8)') &
            (propunit(i,j,ismpl),j=firstidx,lastidx), ' unit = ', i
        END DO

        WRITE(lout,'(a)') key_char//' apriori (units apriori)'
        DO i = 1, maxunits
          WRITE(lout,'(1X,'//c_npropunit//'(e14.5,1X),1A,1I8)') (a_propunit(i,j), &
            j=firstidx,lastidx), ' unit = ', i
        END DO

        WRITE(lout,'(a)') '% errors (variances aposteriori)'
        DO i = 1, maxunits
          WRITE(lout,'(1X,'//c_npropunit//'(e14.5,1X),1A,1I8)') (e_propunit(i,j), &
            j=firstidx,lastidx), ' unit = ', i
        END DO

        WRITE(lout,'(2A)') &
          key_char//' errors (variances apriori) (no reg_par)'
        DO i = 1, maxunits
          WRITE(lout,'(1X,'//c_npropunit//'(e14.5,1X),1A,1I8)') (d_propunit(i,j), &
            j=firstidx,lastidx), ' unit = ', i
        END DO

!     boundary conditions
        IF (bc_maxunits>0) THEN
          WRITE(lout,'(2A)') key_char//' bcunits (aposteriori), iter =', &
            snumber
          DO i = 1, bc_maxunits
            WRITE(lout,'('//c_nbcunit//'(1e14.5,4x),1A,1I5)') (propunit(i,j,ismpl), &
              j=bc_firstidx,bc_lastidx), ' unit = ', i
          END DO

          WRITE(lout,'(a)') key_char//' bcapriori (bcunits apriori)'
          DO i = 1, bc_maxunits
            WRITE(lout,'('//c_nbcunit//'(1e14.5,4x),1A,1I5)') (a_propunit(i,j), &
              j=bc_firstidx,bc_lastidx), ' unit = ', i
          END DO

          WRITE(lout,'(a)') '% bcerrors (variances aposteriori)'
          DO i = 1, bc_maxunits
            WRITE(lout,'('//c_nbcunit//'(1e14.5,4x),1A,1I5)') (e_propunit(i,j), &
              j=bc_firstidx,bc_lastidx), ' unit = ', i
          END DO

          WRITE(lout,'(a)') key_char//' bcerrors (variances apriori)'
          DO i = 1, bc_maxunits
            WRITE(lout,'('//c_nbcunit//'(1e14.5,4x),1A,1I5)') (d_propunit(i,j), &
              j=bc_firstidx,bc_lastidx), ' unit = ', i
          END DO
        END IF

!     time periods
        IF (mpara_tp>0) THEN
          WRITE(lout,'(1A,1I4,2A)') key_char//' tpunits, records=', mpara_tp, &
            ' (aposteriori), iter =', snumber
          DO i = mpara - mpara_tp + 1, mpara
            s_u = seed_para(2,i)
            WRITE(lout,'(2I8,1e14.5,2x,1A)') opti_tp(1,s_u), &
              opti_tp(3,s_u), bcperiod(opti_tp(1,s_u),opti_tp(2,s_u)+1 &
              ,opti_tp(3,s_u),ismpl), tp_albe(opti_tp(2,s_u))
          END DO

          WRITE(lout,'(1A,1I2,1A)') key_char//' tpapriori, records=', &
            mpara_tp, ' (tp apriori)'
          DO i = mpara - mpara_tp + 1, mpara
            s_u = seed_para(2,i)
            WRITE(lout,'(2I8,1e14.5,2x,1A)') opti_tp(1,s_u), &
              opti_tp(3,s_u), a_bcperiod(opti_tp(1,s_u),opti_tp(2,s_u) &
              ,opti_tp(3,s_u)), tp_albe(opti_tp(2,s_u))
          END DO

          WRITE(lout,'(1A,1I2,1A)') '% tperrors, records=', mpara_tp, &
            ' (variances aposteriori)'
          DO i = mpara - mpara_tp + 1, mpara
            s_u = seed_para(2,i)
            WRITE(lout,'(2I8,1e14.5,2x,1A)') opti_tp(1,s_u), &
              opti_tp(3,s_u), e_bcperiod(opti_tp(1,s_u),opti_tp(2,s_u) &
              ,opti_tp(3,s_u)), tp_albe(opti_tp(2,s_u))
          END DO

          WRITE(lout,'(1A,1I4,A)') key_char//' tperrors, records=', mpara_tp, &
            ' (variances apriori)'
          DO i = mpara - mpara_tp + 1, mpara
            s_u = seed_para(2,i)
            WRITE(lout,'(2I8,1e14.5,2x,1A)') opti_tp(1,s_u), &
              opti_tp(3,s_u), d_bcperiod(opti_tp(1,s_u),opti_tp(2,s_u) &
              ,opti_tp(3,s_u)), tp_albe(opti_tp(2,s_u))
          END DO
        END IF


        CLOSE(lout)
        CALL omp_del_file_handler(lout)

        RETURN
      END
