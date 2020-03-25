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

!>    @brief evaluate bayesian quality function for iterate p_k
!>    @param[in] p_k parameter (iteration k)
!>    @param[in] p_apr apriori parameter value
!>    @param[in] ismpl local sample index
!>    @return scalar quality function value
!>    @details
!> evaluate bayesian quality function for iterate p_k\n
!> corresponding to the definition of Tarantola (2004) p.68\n
!>    fkt = [(d -forw)^T *C_dd^-1 *(d -forw) \n
!>          +reg_func() *(p_k -p_apr)^T *C_pp^-1 *(p_k -p_apr)]\n
!>    [1] = {([ndata]-[ndata])^T *[ndata^2] *([ndata] -[ndata])\n
!>          +[1] *([mpara] -[mpara])^T *[mpara^2] *([mpara] -[mpara]}\n
!>    tmp1 = d -forw\n
!>    [ndata] = [ndata] -[ndata]\n
!>    tmp2 = C_dd^-1 *tmp1\n
!>    [ndata] = [ndata^2] *[ndata]\n
!>    tmp3 = p_k -p_apr\n
!>    [mpara] = [mpara] -[mpara]\n
!>    tmp4 = C_pp^-1 *tmp3\n
!>    [mpara] = [mpara^2] *[mpara]\n
!>    fkt = (tmp1^T *tmp2 +reg_func() *tmp3^T *tmp4)\n
!>    [1] = ([ndata]^T *[ndata] +[1] *[mpara]^T *[mpara])\n
      DOUBLE PRECISION FUNCTION qfunc(p_k,p_apr,ismpl)
        use arrays
        use mod_genrl
        use mod_genrlc
        use mod_inverse
        use mod_data
        use mod_linfos
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j

        DOUBLE PRECISION p_k(mpara), p_apr(mpara), t_d, t_p
        DOUBLE PRECISION, ALLOCATABLE :: tmp(:,:)

        DOUBLE PRECISION wp2_ps, wd2_ps, reg_func
!              [mpara^2], [ndata^2]
        EXTERNAL wp2_ps, wd2_ps, reg_func


        ALLOCATE(tmp(max(mpara,ndata),4))

!     tmp1 = d -forw
!     [ndata] = [ndata] -[ndata]
        DO i = 1, ndata
          tmp(i,1) = ddata(i,cdd_pv) - sdata(i,ismpl)
        END DO

!     tmp2 = C_dd^-1 *tmp1
!     [ndata] = [ndata^2] *[ndata]
        DO j = 1, ndata
          tmp(j,2) = 0.0D0
          DO i = 1, ndata
            tmp(j,2) = tmp(j,2) + wd2_ps(j,i,ismpl)*tmp(i,1)
          END DO
        END DO

!     tmp3 = p_k -p_apr
!     [mpara] = [mpara] -[mpara]
        DO i = 1, mpara
          tmp(i,3) = p_k(i) - p_apr(i)
        END DO

!     tmp4 = C_pp^-1 *tmp3
!     [mpara] = [mpara^2] *[mpara]
        DO j = 1, mpara
          tmp(j,4) = 0.0D0
          DO i = 1, mpara
            tmp(j,4) = tmp(j,4) + wp2_ps(j,i,ismpl)*tmp(i,3)
          END DO
        END DO

!     fkt = (tmp1^T *tmp2 + tmp3^T *tmp4)
!     [1] = ([ndata]^T *[ndata] +[mpara]^T *[mpara])
        CALL s_ddot(ndata,tmp(1,1),tmp(1,2),t_d)
        CALL s_ddot(mpara,tmp(1,3),tmp(1,4),t_p)
        qfunc = 0.5D0*(t_d+reg_func()*t_p)

        IF (linfos(2)>=2) WRITE(*,'(A,3(e15.9,A))') &
          '  objective function=', qfunc, ', ( d:', t_d, ', p:', t_p, ' )'

        rms_data = t_d
        rms_para = t_p
        DEALLOCATE(tmp)

        RETURN
      END
