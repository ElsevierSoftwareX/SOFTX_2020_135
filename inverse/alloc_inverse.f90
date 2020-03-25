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

!>    @brief allocate objects for inverse use
!>    @param[in] ismpl local sample index
      SUBROUTINE alloc_inverse(ismpl)
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
        integer :: i, j, k

        INTEGER i_max,si,sj


        IF (linfos(1)>=2) WRITE(*,*) ' [I] : ... alloc_inverse'

        ALLOCATE(seed_para(2,mpara))
        memory = memory + 2*mpara

        ALLOCATE(apri_input(mpara))
        memory = memory + mpara
        ALLOCATE(linlog_input(mpara))
        memory = memory + mpara

!     init seed index
        k = 0
        DO j = 1, maxunits
          DO i = firstidx, lastidx
            IF (opti_props(i-firstidx+1,j)>0) THEN
              k = k + 1
              seed_para(1,k) = i
              seed_para(2,k) = j
            END IF
          END DO
        END DO
        DO j = 1, bc_maxunits
          DO i = bc_firstidx, bc_lastidx
            IF (opti_bc(i-bc_firstidx+1,j)>0) THEN
              k = k + 1
              seed_para(1,k) = i
              seed_para(2,k) = j
            END IF
          END DO
        END DO
        DO i = 1, mpara_tp
          k = k + 1
          seed_para(1,k) = -1
          seed_para(2,k) = i
        END DO
!     "mpara" should be equal to "k"
        IF (mpara/=k) THEN
          WRITE(*,'(2A)') 'error: lost some optimization', &
            ' paramaters in "read_invers.f" !!!'
          STOP
        END IF

        ALLOCATE(main_input_master(mpara))
        memory = memory + mpara
#ifdef AD
        ALLOCATE(g_main_input(mpara,nsmpl))
        memory = memory + mpara*nsmpl
#endif
#ifdef AD_RM
        ALLOCATE(main_input_ad(max(ndata,1),nsmpl))
        memory = memory + max(ndata,1)+nsmpl
#endif

#ifndef JACOBI_FREE
        ALLOCATE(jac(max(ndata,1),max(mpara,1)))
#ifdef AD_RM
        ALLOCATE(jacT(max(mpara,1),max(ndata,1)))
#endif
        memory = memory + max(ndata,1)*max(mpara,1)
        CALL set_dval(ndata*mpara,0.D0,jac)

        ALLOCATE(tmp_vec(mpara,max(ndata,mpara),5))
        memory = memory + mpara*max(ndata,mpara)*5
        ALLOCATE(tmp_mat(mpara,mpara))
        memory = memory + mpara*mpara
#else
        ALLOCATE(tmp_vec(max(ndata,mpara),1,4))
        memory = memory + max(ndata,mpara)*4
        ALLOCATE(tmp_mat(1,1))
        memory = memory + 1*1
#endif

!     dummy allocation against zero pointers
        ALLOCATE(grad_sec(1))
        ALLOCATE(grad(1))
!     Resmat need Covar (see calc_covariances)
#ifndef JACOBI_FREE
        IF ((covar/=0) .OR. (resmat/=0)) THEN
          ALLOCATE(covar_p(mpara,mpara))
          memory = memory + mpara*mpara
        END IF
#endif

        ALLOCATE(wdt(2))
        memory = memory + 2
        wdt(1) = 1.D0
        wdt(2) = 1.D0

#ifndef JACOBI_FREE
        IF (resmat/=0) THEN
          ALLOCATE(resmat_p(mpara,mpara))
          memory = memory + mpara*mpara
          IF (ndata>0) THEN
            ALLOCATE(ginv(mpara,ndata))
            memory = memory + mpara*ndata
          END IF
        END IF
#endif

!     sample initialisation
        i_max = max(maxunits,bc_maxunits)
        ALLOCATE(propunitold(i_max,nprop))
        memory = memory + i_max*nprop
        ALLOCATE(bcperiodold(ngsmax,3,max(nbctp,1)))
        memory = memory + ngsmax*3*max(nbctp,1)

        RETURN
      END SUBROUTINE alloc_inverse
