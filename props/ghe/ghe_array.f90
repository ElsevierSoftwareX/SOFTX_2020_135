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

      MODULE work_array_ghe



!        IMPLICIT NONE



        DOUBLE PRECISION, ALLOCATABLE :: Tu(:,:,:), Td(:,:,:), &
          Tout(:,:,:), Tin_new(:,:), Tin(:,:), dTsghe(:), Hghe(:,:),&
          dTs(:), Tsghe(:,:), qshe(:,:,:), dTsgheprime(:), &
          Tsghe_old(:,:), pghe(:,:), Temp_flip(:,:,:), H_flip(:,:,:), &
          pres_flip(:,:,:), Ts(:,:), dTsdepth(:), ypp(:), press(:,:), &
          Tsspline(:), pressspline(:), lagrange(:,:,:), yppq(:), &
          dpres(:,:)

        INTEGER nt, nghe, tWRITE2, time, endtime, ng, it, k_p, iper, &
          tWRITE, np, nr, modswitch
        INTEGER, ALLOCATABLE :: ighe(:),jghe(:),fghe(:),k_end(:), &
          k_start(:)
        DOUBLE PRECISION  QF, kw, Cpw, Dw, ru, dz, ruo, Bu, ks, kfill, &
          r1, fluid_type, deltatime, ntdouble, Cpwu, Dwu, Cpwd, Dwd, &
          depthup, kwd, kwu, dy, Cpwi, Dwi, kwi, kinvisci, &
          time_count, dghe

        DOUBLE PRECISION, ALLOCATABLE :: rdu(:), rd(:), rgr(:), ku(:), &
          kd(:), kgr(:), Rdb(:), Rudu(:), Rdown(:), Rup(:), Ru_du_cd(:),&
          Rb_d_cd(:), depth(:),rb(:)

        INTEGER mp,nl,p, fup, kdepth, ghetype, kl  ! l

      END MODULE work_array_ghe
