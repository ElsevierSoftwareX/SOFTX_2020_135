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

!>    @brief (weighted) jacoby-matrix output for each physical value
!>    @param[in] ismpl local sample index
      SUBROUTINE write_jacw(ismpl)
        use arrays
#ifndef JACOBI_FREE
#ifndef AD_RM
        use g_arrays
#else
        use arrays_ad
#endif
        use mod_genrl
        use mod_genrlc
        use mod_inverse
        use mod_data
        use mod_linfos
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k, l

        DOUBLE PRECISION dx, dy, dz, errd
        DOUBLE PRECISION, ALLOCATABLE :: jac_row(:)
        INTEGER i1, i2, i3, ii, kk, type, ll_pv, lout
        INTEGER i1_, i2_, i3_, i4_

        CHARACTER filename*80, snumber*8
        INTEGER lblank
        EXTERNAL lblank

        IF (write_disable) RETURN
#ifdef NOWJAC
        RETURN
#endif

        IF (mpara>=9999) THEN
!       needs modifying on the write-format
          WRITE(*,*) 'error : in "write_jacw" need "mpara"<9999 !'
          STOP
        END IF

        CALL omp_new_file_handler(lout,16)

        ALLOCATE(jac_row(mpara))
        CALL chln(project,i1_,i2_)
        WRITE(snumber,'(1I7)') iter_inv
        CALL chln(snumber,i3_,i4_)

        DO ll_pv = 1, npv
          kk = 0
          IF (ll_pv==pv_head) kk = ndata_h
          IF (ll_pv==pv_temp) kk = ndata_t
          IF (ll_pv==pv_conc) kk = ndata_c
          IF (ll_pv==pv_pres) kk = ndata_p

          IF (kk>0) THEN
            filename = project(i1_:i2_) // '_JW' // pv_name(ll_pv) // &
              '_' // snumber(i3_:i4_) // '.plt'
            OPEN(lout,file=filename,status='unknown',blank='null')

            IF (linfos(2)>=0) THEN
              WRITE(*,'(5A)') '  [W] : ', pv_name(ll_pv), &
                '-JacW to "', filename(1:lblank(filename)), '"'
            END IF

! 'x y z uid pv-value por1 akx1 aky1 kz1 alx1 aly1 lz1 q1 por2 akx2 aky2 kz2 alx2 aly2 lz2 q2'
            CALL head_tecpl(lout,kk,1,1,pv_name(ll_pv),ismpl)

            DO l = 1, ndata
              type = idata(l,cid_pv)
              IF (type==ll_pv) THEN
                i = idata(l,cid_i)
                j = idata(l,cid_j)
                k = idata(l,cid_k)
                errd = 1.D0/ddata(l,cdd_w)
                dx = delxa(i)
                dy = delya(j)
                dz = delza(k)
                DO ii = 1, mpara
                  jac_row(ii) = jac(l,ii)
                END DO
                WRITE(lout,'(3(e16.7),1I8,9999(e16.7,1X))') dx, dy, &
                  dz, uindex(i,j,k), sdata(l,ismpl), &
                  (errd*jac_row(ii),ii=1,mpara)
              END IF
            END DO
            CLOSE(lout)
          END IF
        END DO
        DEALLOCATE(jac_row)

        CALL omp_del_file_handler(lout)
#endif
        RETURN
      END
