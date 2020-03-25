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

!>    @brief wrapper for the parameter output
!>    @param[in] r "runmode"
!>    @param[in] i inverse iteration counter
!>    @param[in] s (start) inverse iteration counter
!>    @param[in] ismpl local sample index
      SUBROUTINE para_write(r,i,s,ismpl)
        IMPLICIT NONE
        INTEGER r, i, s, ismpl, ll

        IF (r>1) CALL write_parameter(i,ismpl)
        IF (r>1) CALL write_hdfparameter(s,ismpl)
        RETURN
      END

!>    @brief wrapper for the restart output (currently disabled)
!>    @param[in] restart_name restart file name
!>    @param[in] iter_inv inverse iteration counter
!>    @param[in] ismpl local sample index
      SUBROUTINE resinverse_write(restart_name,iter_inv,ismpl)
        IMPLICIT NONE
        CHARACTER restart_name*(*)
        INTEGER iter_inv, ismpl

        CALL write_restartinv(restart_name,iter_inv,ismpl)
!     against problems during writing
        CALL write_restartinv('_'//restart_name,iter_inv,ismpl)
        RETURN
      END

!>    @brief writes the header for combinded parameter & data output
!>    @param[in] ndata number of observed data
!>    @param[in] mpara number of active parameters
!>    @param[in] iter_inv inverse iteration counter
!>    @param[in] rms_data data fit (RMS)
!>    @param[in] rms_para parameter fit (RMS)
!>    @param[in] ismpl local sample index
      SUBROUTINE datapara_output(ndata,mpara,iter_inv,rms_data, &
          rms_para,ismpl)
        IMPLICIT NONE
        INTEGER ndata, mpara, iter_inv, ismpl
        DOUBLE PRECISION rms_data, rms_para
        INTRINSIC dsqrt, dble

        IF ((ndata/=0) .AND. (mpara/=0)) THEN
          WRITE(*,'(1A,I16)') '  objective function, iteration=', &
            iter_inv
          WRITE(*,'(1A,2e24.16)') '   data     = ', rms_data, &
            dsqrt(rms_data)/dble(ndata)
          WRITE(*,'(1A,1e24.16)') '   parameter= ', rms_para
        END IF
        RETURN
      END

!>    @brief writes the header for data output
!>    @param[in] ndata number of observed data
!>    @param[in] iter_inv inverse iteration counter
!>    @param[in] maxiter_inv max inverse iteration number
!>    @param[in] tol_inv tolerance criteria
!>    @param[in] rms_data data fit (RMS)
!>    @param[in] ismpl local sample index
      SUBROUTINE data_output(ndata,iter_inv,maxiter_inv,tol_inv, &
          rms_data,ismpl)
        IMPLICIT NONE
        INTEGER ndata, iter_inv, maxiter_inv, ismpl
        DOUBLE PRECISION tol_inv, rms_data
        INTRINSIC dsqrt, dble

        IF (ndata>0) THEN
          IF (dsqrt(rms_data)<=tol_inv*dble(ndata)) WRITE(*, &
            '(A,1e20.12,A,1e20.6)') '  data rms ', &
            dsqrt(rms_data)/dble(ndata), ' reached tolerance ', &
            tol_inv
          IF (iter_inv>=maxiter_inv) WRITE(*,'(A,1e20.12,A,I6)') &
            '  data rms ', dsqrt(rms_data)/dble(ndata), &
            ' reached at max number of iterations, ', maxiter_inv
        END IF
        RETURN
      END
