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

!>    @brief final divided difference computation of the Jacobian matrix (h-step)
!>    @param[in] seedi seeding component/index
!>    @param[in] hh current step size h
!>    @param[in] hh_bak old step size h
!>    @param[in] ismpl local sample index
!>    @details
!>      temp,head,conc,pres : filled out (jacoby values DD)\n
      SUBROUTINE dd_jacobian(seedi,hh,hh_bak,ismpl)
        use arrays
#ifdef AD 
        use g_arrays
#endif
#ifdef AD_RM 
        use arrays_ad
#endif
        use mod_genrl
        use mod_genrlc
        use mod_conc
        use mod_inverse
        use mod_linfos
        IMPLICIT NONE
        integer :: ismpl
        INTEGER seedi, s_k, s_u, p_offs, p_end
        DOUBLE PRECISION hh, hh_bak, hh_log, s_diff, s_diff2, ddot
        LOGICAL test_null
        INTEGER lblank, param_enabled
        EXTERNAL ddot, lblank, test_null, param_enabled
        INTRINSIC dsqrt


        IF ((seedi<1) .OR. (seedi>mpara)) THEN
          WRITE(*,*) 'error: in "dd_jacobian", no valid "seeding=', &
            seedi, '"'
          STOP
        END IF

!     setup for log optimization
        hh_log = hh_bak
!     restore parameter
        CALL set_optip(seedi,hh_bak,ismpl)

!     in case of linear optimization
        IF (param_enabled(seedi,ismpl)==1) hh_log = 1.0D0

!      [f(x+h)-f(x)]/h
!   -> hh_bak*[f(x+h)-f(x)]/h !!! hh_bak<>1.0d0 in case of log optimization
        CALL daxpy(i0*j0*k0,-1.0D0,headold(1,cgen_dd,idx_master),1, &
          head(1,1,1,ismpl),1)
        CALL daxpy(i0*j0*k0,-1.0D0,tempold(1,cgen_dd,idx_master),1, &
          temp(1,1,1,ismpl),1)
        CALL daxpy(i0*j0*k0*ntrans,-1.0D0,concold(1,1,cgen_dd,idx_master), &
          1,conc(1,1,1,1,ismpl),1)
        CALL daxpy(i0*j0*k0,-1.0D0,presold(1,cgen_dd,idx_master),1, &
          pres(1,1,1,ismpl),1)

        IF ( .NOT. test_null(hh_bak)) THEN
          CALL dscal(i0*j0*k0,hh_log/(hh_bak*hh),head(1,1,1,ismpl),1)
          CALL dscal(i0*j0*k0,hh_log/(hh_bak*hh),temp(1,1,1,ismpl),1)
          CALL dscal(i0*j0*k0*ntrans,hh_log/(hh_bak*hh), &
            conc(1,1,1,1,ismpl),1)
          CALL dscal(i0*j0*k0,hh_log/(hh_bak*hh),pres(1,1,1,ismpl),1)
        ELSE
          CALL dscal(i0*j0*k0,hh_log/hh,head(1,1,1,ismpl),1)
          CALL dscal(i0*j0*k0,hh_log/hh,temp(1,1,1,ismpl),1)
          CALL dscal(i0*j0*k0*ntrans,hh_log/hh,conc(1,1,1,1,ismpl),1)
          CALL dscal(i0*j0*k0,hh_log/hh,pres(1,1,1,ismpl),1)
        END IF

        RETURN
      END
