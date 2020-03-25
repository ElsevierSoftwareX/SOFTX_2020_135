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

!> @brief parameter output for each SIMUL/ENKF iteration
!> @param[in] ident index number for file name
!> @param[in] ismpl local sample index
      SUBROUTINE write_parameter(ident,ismpl)
        use arrays
        use mod_genrl
        use mod_genrlc
        use mod_data
        use mod_simul
        use mod_linfos
        use mod_OMP_TOOLS
        IMPLICIT NONE
        
	INCLUDE 'OMP_TOOLS.inc'
        integer :: i, j, k, ident, ismpl
        integer :: i1, i2, i3, i4, type, lblank, s_u, lout
        EXTERNAL lblank
        character (len=80) :: filename
        character (len=20) :: snumber
        DOUBLE PRECISION tmp1, tmp2
        character (len=4), dimension (2) :: tp_albe
        DATA tp_albe/'alfa', 'beta'/


        IF ((mpara<=0) .OR. ( .NOT. write_param)) RETURN
#ifdef NOOUT
        RETURN
#endif
!
!     get his own file discriptor index
        CALL omp_new_file_handler(lout,15)
!
        CALL chln(project,i1,i2)
!
        IF (ident>=0) THEN
          WRITE(snumber,'(1I7)') ident
        ELSE IF (ident==-1) THEN
          WRITE(snumber,'(A20)') 'final'
        ELSE IF (ident==-2) THEN
          WRITE(snumber,'(A20)') 'debug'
        ELSE IF (ident==-3) THEN
          WRITE(snumber,'(A20)') 'ens_mean'
        ELSE IF (ident==-4) THEN
          WRITE(snumber,'(A20)') 'mean'
        ELSE IF (ident==-5) THEN
          WRITE(snumber,'(A20)') 'ens_mean'
        END IF
        CALL chln(snumber,i3,i4)
!
        filename = project(i1:i2) // '_P_' // snumber(i3:i4) // '.dat'
        OPEN(lout,file=filename,status='unknown',blank='null')
        IF (linfos(3)>=1) THEN
          WRITE(*,'(3A)') '  [W] : Parameter to "', &
            filename(1:lblank(filename)), '"'
        END IF

!     parameter units
        WRITE(lout,'(2A)') key_char//' units (aposteriori), iter =', snumber
        DO i = 1, maxunits
          WRITE(lout,'(1X,'//c_npropunit//'(e14.5,1X),1A,1I8)') &
            (propunit(i,j,ismpl),j=firstidx,lastidx), ' unit = ', i
        END DO	

!     boundary conditions
        IF (bc_maxunits>0) THEN
          WRITE(lout,'(2A)') key_char//' bcunits (aposteriori), iter =', &
            snumber
          DO i = 1, bc_maxunits
            WRITE(lout,'(1X,'//c_nbcunit//'(1e14.5,1X),1A,1I7)') &
              (propunit(i,j,ismpl),j=bc_firstidx,bc_lastidx), &
              ' unit = ', i
          END DO
          WRITE(lout,'(2A)') key_char//' splitted bcunits, iter =', snumber
          WRITE(lout,'(1A)') ' <i-idx> <j-idx> <k-idx> <bc-unit> &
            &<pv-idx>       <value>'
          DO k = 1, mpara
            IF (seed_para(1,k)==-2) WRITE(lout, &
              '(3(1I7,1X),1I9,1X,1I8,1X,1e14.5)') &
              ibc_data(seed_para(2,k),cbc_i), ibc_data(seed_para(2,k), &
              cbc_j), ibc_data(seed_para(2,k),cbc_k), &
              ibc_data(seed_para(2,k),cbc_bcu), &
              ibc_data(seed_para(2,k),cbc_pv), &
              dbc_data(seed_para(2,k),1,ismpl)
          END DO
        END IF

!     time periods
        IF (mpara_tp>0 .AND. mpara>0) THEN
          WRITE(lout,'(1A,1I4,2A)') key_char//' tpunits, records=', mpara_tp, &
            ' (aposteriori), iter =', snumber
          DO i = mpara - mpara_tp + 1, mpara
            s_u = seed_para(2,i)
            WRITE(lout,'(2I8,1e14.5,2x,1A)') opti_tp(1,s_u), &
              opti_tp(3,s_u), bcperiod(opti_tp(1,s_u),opti_tp(2,s_u)+1 &
              ,opti_tp(3,s_u),ismpl), tp_albe(opti_tp(2,s_u))
          END DO
        END IF
 
        CLOSE(lout)
        CALL omp_del_file_handler(lout)

        RETURN
      END
