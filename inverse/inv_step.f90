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

!>    @brief encapsulated inversion step
!>    @param[in] qval_old old quality function value
!>    @param[in,out] qval_new new quality function value
!>    @param[in] p_k current parameter vector
!>    @param[in] p_apr apriori parameter vector
!>    @param[in] pgrad current derivative values
!>    @param[in] pgrad_sec old derivative values
!>    @param[in,out] ptmp_vec temporary vectors (unspecified values)
!>    @param[in,out] ptmp_mat temporary matrix (unspecified values)
!>    @param[in] ismpl local sample index
      SUBROUTINE inv_step(qval_new,qval_old,p_k,p_apr,pgrad,pgrad_sec, &
        ptmp_vec,ptmp_mat,ismpl)
      use arrays
      use mod_genrl
      use mod_data
      use mod_inverse
      IMPLICIT NONE
      integer :: ismpl
      DOUBLE PRECISION qval_old,qval_new,val,sum
      DOUBLE PRECISION p_k(*),p_apr(*),ptmp_vec(*),pgrad(*),pgrad_sec(*)
      DOUBLE PRECISION ptmp_mat(*)
      DOUBLE PRECISION, allocatable :: tmp_dmat(:,:), p_delta(:)

      ALLOCATE(p_delta(mpara))
!
#ifdef STBAY
!     std. Step-Bayes
      CALL step_bayes(p_k,p_apr,ptmp_vec,ptmp_mat,p_delta,ismpl)
#elif PSPAC
!     std. Parameter-Space
      CALL step_p_space(p_k,p_apr,ptmp_vec,ptmp_mat,p_delta,ismpl)
#elif DSPAC
!     std. Data-Space
      WRITE(*,'(2a,1e12.5,a)') ' warning: allocate extra memory ', &
        '[ndata x ndata]: ',dble(ndata*ndata/ (1024**2)),'MB !'
      ALLOCATE(tmp_dmat(ndata,ndata))
      CALL step_d_space(p_k,p_apr,ptmp_vec,tmp_dmat,p_delta,ismpl)
      DEALLOCATE(tmp_dmat)
#elif MF_STBAY
!     matrix free Step-Bayes (for reverse mode)
      CALL step_bayes_mf(p_k,p_apr,ptmp_vec,p_delta,ismpl)
#else
      WRITE(*,'(1A)') 'error: no inversion step method specified !!!'
      STOP
#endif
!
      CALL update_model(qval_new,qval_old,p_k,p_apr,p_delta,ismpl)
!
      DEALLOCATE(p_delta)
!
      RETURN
      END

!>    @brief initialisation for different solver (optimisation methods)
!>    @param[in] irelax max iteration number for relaxiation (ignored)
      SUBROUTINE init_step()
        USE mod_inverse
      IMPLICIT NONE
      INTEGER iopti_method

!
      iopti_method = -1
      abort_inv = .FALSE.
      hdiag_provide = .FALSE.
!
#ifdef STBAY
      iopti_method = 0
      relax_max = 3
      relax_out = 1
#elif  PSPAC
      iopti_method = 1
      relax_max = 3
      relax_out = 1
#elif  DSPAC
      iopti_method = 2
      relax_max = 3
      relax_out = 1
#elif  MF_STBAY
      iopti_method = 5
      relax_max = 3
      relax_out = 1
#else
      WRITE(*,'(1A)') 'error: no inversion step method specified !!!'
      STOP
#endif
!
#ifdef JACOBI_FREE
      IF (iopti_method==0 .OR. iopti_method==1 .OR. iopti_method==2) THEN
        WRITE(*,'(1A)') 'error: (Jacobi) matrix free computation not supported for this optimisation method!!!'
        STOP
      END IF
#endif
!
      RETURN
      END
