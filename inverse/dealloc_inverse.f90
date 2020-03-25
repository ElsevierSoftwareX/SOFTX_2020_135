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

!>    @brief free memory of objects for inverse use
!>    @param[in] ismpl local sample index
      SUBROUTINE dealloc_inverse(ismpl)
        use arrays
#ifdef AD
        use g_arrays
#endif
#ifdef AD_RM 
        use arrays_ad
#endif
        use mod_genrl
        use mod_genrlc
        use mod_time
        use mod_inverse
        use mod_data
        use mod_linfos
        IMPLICIT NONE
        integer :: ismpl
        integer :: i

        INTEGER i_max


        IF (linfos(1)>=2) WRITE(*,*) ' [I] : ... dealloc_inverse'

        DEALLOCATE(seed_para)
        memory = memory - 2*mpara

        DEALLOCATE(apri_input)
        memory = memory - mpara

        DEALLOCATE(linlog_input)
        memory = memory - mpara

        DEALLOCATE(main_input_master)
        memory = memory - mpara
#ifdef AD
        DEALLOCATE(g_main_input)
        memory = memory - mpara*nsmpl
#endif
#ifdef AD_RM
        DEALLOCATE(main_input_ad)
        memory = memory - ndata*nsmpl
#endif

#ifndef JACOBI_FREE
        DEALLOCATE(jac)
        memory = memory - max(ndata,1)*max(mpara,1)
#ifdef AD_RM
        DEALLOCATE(jacT)
        memory = memory - max(ndata,1)*max(mpara,1)
#endif 
!
        DEALLOCATE(tmp_vec)
        memory = memory - mpara*max(ndata,mpara)*5
        DEALLOCATE(tmp_mat)
        memory = memory - mpara*mpara
#else
        DEALLOCATE(tmp_vec)
        memory = memory - max(ndata,mpara)*4
        DEALLOCATE(tmp_mat)
        memory = memory - 1*1
#endif

        DEALLOCATE(grad)
        DEALLOCATE(grad_sec)

#ifndef JACOBI_FREE
        IF ((covar/=0) .OR. (resmat/=0)) THEN
          DEALLOCATE(covar_p)
          memory = memory - mpara*mpara
        END IF
#endif

        DEALLOCATE(wdt)
        memory = memory - 2

#ifndef JACOBI_FREE
        IF (resmat/=0) THEN
          DEALLOCATE(resmat_p)
          memory = memory - mpara*mpara
          IF (ndata>0) THEN
            memory = memory - mpara*ndata
          END IF
        END IF
#endif

        IF (para_weight>=1) THEN
          DEALLOCATE(covar_prior_p)
          memory = memory - mpara*mpara
        END IF
        IF (data_weight>=1) THEN
          DEALLOCATE(covar_prior_d)
          memory = memory - ndata*ndata
        END IF

!     sample initialisation
        i_max = max(maxunits,bc_maxunits)
        DEALLOCATE(propunitold)
        memory = memory - i_max*nprop
        DEALLOCATE(bcperiodold)
        memory = memory - ngsmax*3*max(nbctp,1)

        IF (lread_joutt) DEALLOCATE(seed_index)

        RETURN
      END SUBROUTINE dealloc_inverse
