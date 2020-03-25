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

!>    @brief apply 7point-star matrix multiply [as]:=[M]x[s], (OpenMP version)
!>    @param[in] N_I lengths of I dimension of local matrix [M]
!>    @param[in] N_J lengths of J dimension of local matrix [M]
!>    @param[in] N_K lengths of K dimension of local matrix [M]
!>    @param[in] s vector [s]
!>    @param[out] as vector [as]
!>    @param[in] A 1. diagonal of the system matrix [M]
!>    @param[in] B 2. diagonal of the system matrix [M]
!>    @param[in] C 3. diagonal of the system matrix [M]
!>    @param[in] D 4. (main) diagonal of the system matrix [M]
!>    @param[in] E 5. diagonal of the system matrix [M]
!>    @param[in] F 6. diagonal of the system matrix [M]
!>    @param[in] G 7. diagonal of the system matrix [M]
!>    @details
!>    OpenMP parallelised, general version - no special blocking\n
!>    apply 7point-star matrix multiply\n
!>    compute [as]:=[M]x[s], [s],[as],[M] given in 3-D-structure\n
!>    Data-Cube :\n
!>    @image html cube.png
!       k     * * * *
!     /     *     * *
!    0 -j * * * *   *
!    |    *     *   *
!    i    *     * *
!         * * * *
      SUBROUTINE omp_mvp(n_i,n_j,n_k,s,as,a,b,c,d,e,f,g)
        use mod_OMP_TOOLS
        use mod_blocking_size
        IMPLICIT NONE
        INCLUDE 'OMP_TOOLS.inc'

        INTEGER n_i, n_j, n_k
!      use mod_blocking_size

        DOUBLE PRECISION s(*), as(*)
        DOUBLE PRECISION a(*), b(*), c(*), d(*), e(*), f(*), g(*)

        INTEGER pim, pip, pjm, pjp, pkm, pkp
        INTEGER aim, aip, ajm, ajp, akm, akp
        INTEGER bm, pm, pmb, am, pmmax

!     thread stuff
        INTEGER tpos, tanz


        tpos = 1
        tanz = n_i*n_j*n_k
!$      call OMP_PART(N_I*N_J*N_K,tpos,tanz)

        bm = tpos + tanz - 1
        pm = tpos
        pmmax = n_i*n_j*n_k + 1


!     not ready yet, should modify for better performance !
!     +2.5% for '/27'
        pmb = min(pm+int(bl_size/bldiv_mvp),bm+1)
        am = pmb - pm

!     I-dim
        pim = max(pm-1,1)
        aim = pmb - 1 - pim
        pip = min(pmb+1,pmmax)
        aip = pip - 1 - pm

!     J-dim
        pjm = max(pm-n_i,1)
        ajm = pmb - n_i - pjm
        pjp = min(pmb+n_i,pmmax)
        ajp = pjp - n_i - pm

!     K-dim
        pkm = max(pm-n_i*n_j,1)
        akm = pmb - n_i*n_j - pkm
        pkp = min(pmb+n_i*n_j,pmmax)
        akp = pkp - n_i*n_j - pm

100     CONTINUE

        CALL dxyz(am,s(pm),d(pm),as(pm))
!AW C$OMP barrier

        IF (n_i>1) THEN
          CALL dxypz(aip,s(pm+1),e(pm),as(pm))
!AW C$OMP barrier
        END IF
        IF (n_j>1) THEN
          CALL dxypz(ajp,s(pm+n_i),f(pm),as(pm))
!AW C$OMP barrier
        END IF
        IF (n_k>1) THEN
          CALL dxypz(akp,s(pm+n_i*n_j),g(pm),as(pm))
!AW C$OMP barrier
        END IF

        IF (n_i>1) THEN
          CALL dxypz(aim,s(pim),c(pim+1),as(pim+1))
!AW C$OMP barrier
        END IF
        IF (n_j>1) THEN
          CALL dxypz(ajm,s(pjm),b(pjm+n_i),as(pjm+n_i))
!AW C$OMP barrier
        END IF
        IF (n_k>1) THEN
          CALL dxypz(akm,s(pkm),a(pkm+n_i*n_j),as(pkm+n_i*n_j))
!AW C$OMP barrier
        END IF

!aw      write(*,*)' ',OMP_GET_HIS_THREAD_NUM(),pm,pm+am-1,pmb-1
!aw      write(*,*)'iC',OMP_GET_HIS_THREAD_NUM(),pim+1,-1,aim
!aw      write(*,*)'iE',OMP_GET_HIS_THREAD_NUM(),pm,1,aip
!aw      write(*,*)'jB',OMP_GET_HIS_THREAD_NUM(),pjm+N_I,-N_I,ajm
!aw      write(*,*)'jF',OMP_GET_HIS_THREAD_NUM(),pm,N_I,ajp
!aw      write(*,*)'kA',OMP_GET_HIS_THREAD_NUM(),pkm+N_I*N_J,-N_I*N_J,akm
!aw      write(*,*)'kG',OMP_GET_HIS_THREAD_NUM(),pm,N_I*N_J,akp

        pm = pmb
        pmb = min(pm+int(bl_size/bldiv_mvp),bm+1)
        am = pmb - pm

!     I-dim
        pim = max(pm-1,1)
        aim = pmb - 1 - pim
        pip = min(pmb+1,pmmax)
        aip = pip - 1 - pm

!     J-dim
        pjm = max(pm-n_i,1)
        ajm = pmb - n_i - pjm
        pjp = min(pmb+n_i,pmmax)
        ajp = pjp - n_i - pm

!     K-dim
        pkm = max(pm-n_i*n_j,1)
        akm = pmb - n_i*n_j - pkm
        pkp = min(pmb+n_i*n_j,pmmax)
        akp = pkp - n_i*n_j - pm

        IF (am>0) GO TO 100

!$OMP  barrier
!      need barrier here

        RETURN
      END

!>    @brief apply 7point-star matrix multiply [as]:=[M]x[s], serial (no OpenMP) implementation
!>    @param[in] N_I lengths of I dimension of local matrix [M]
!>    @param[in] N_J lengths of J dimension of local matrix [M]
!>    @param[in] N_K lengths of K dimension of local matrix [M]
!>    @param[in] s vector [s]
!>    @param[out] as vector [as]
!>    @param[in] A 1. diagonal of the system matrix [M]
!>    @param[in] B 2. diagonal of the system matrix [M]
!>    @param[in] C 3. diagonal of the system matrix [M]
!>    @param[in] D 4. (main) diagonal of the system matrix [M]
!>    @param[in] E 5. diagonal of the system matrix [M]
!>    @param[in] F 6. diagonal of the system matrix [M]
!>    @param[in] G 7. diagonal of the system matrix [M]
!>    @details
!>    serial (not OpenMP parallelised), general version - no special blocking\n
!>    apply 7point-star matrix multiply\n
!>    compute [as]:=[M]x[s], [s],[as],[M] given in 3-D-structure\n
!>    Data-Cube :\n
!>    @image html cube.png
!       k     * * * *
!     /     *     * *
!    0 -j * * * *   *
!    |    *     *   *
!    i    *     * *
!         * * * *
      SUBROUTINE s_mvp(n_i,n_j,n_k,s,as,a,b,c,d,e,f,g)
        use mod_blocking_size
        IMPLICIT NONE

        INTEGER n_i, n_j, n_k
!      use mod_blocking_size

        DOUBLE PRECISION s(*), as(*)
        DOUBLE PRECISION a(*), b(*), c(*), d(*), e(*), f(*), g(*)

        INTEGER pim, pip, pjm, pjp, pkm, pkp
        INTEGER aim, aip, ajm, ajp, akm, akp
        INTEGER bm, pm, pmb, am, pmmax

!     thread stuff (disabled here)
        INTEGER tpos, tanz


        tpos = 1
        tanz = n_i*n_j*n_k

        bm = tpos + tanz - 1
        pm = tpos
        pmmax = n_i*n_j*n_k + 1

!     not ready yet, should modify for better performance !
!     +2.5% for '/27'
        pmb = min(pm+int(bl_size/bldiv_mvp),bm+1)
        am = pmb - pm

!     I-dim
        pim = max(pm-1,1)
        aim = pmb - 1 - pim
        pip = min(pmb+1,pmmax)
        aip = pip - 1 - pm

!     J-dim
        pjm = max(pm-n_i,1)
        ajm = pmb - n_i - pjm
        pjp = min(pmb+n_i,pmmax)
        ajp = pjp - n_i - pm

!     K-dim
        pkm = max(pm-n_i*n_j,1)
        akm = pmb - n_i*n_j - pkm
        pkp = min(pmb+n_i*n_j,pmmax)
        akp = pkp - n_i*n_j - pm

100     CONTINUE

        CALL dxyz(am,s(pm),d(pm),as(pm))

        IF (n_i>1) THEN
          CALL dxypz(aip,s(pm+1),e(pm),as(pm))
        END IF
        IF (n_j>1) THEN
          CALL dxypz(ajp,s(pm+n_i),f(pm),as(pm))
        END IF
        IF (n_k>1) THEN
          CALL dxypz(akp,s(pm+n_i*n_j),g(pm),as(pm))
        END IF

        IF (n_i>1) THEN
          CALL dxypz(aim,s(pim),c(pim+1),as(pim+1))
        END IF
        IF (n_j>1) THEN
          CALL dxypz(ajm,s(pjm),b(pjm+n_i),as(pjm+n_i))
        END IF
        IF (n_k>1) THEN
          CALL dxypz(akm,s(pkm),a(pkm+n_i*n_j),as(pkm+n_i*n_j))
        END IF

        pm = pmb
        pmb = min(pm+int(bl_size/bldiv_mvp),bm+1)
        am = pmb - pm

!     I-dim
        pim = max(pm-1,1)
        aim = pmb - 1 - pim
        pip = min(pmb+1,pmmax)
        aip = pip - 1 - pm

!     J-dim
        pjm = max(pm-n_i,1)
        ajm = pmb - n_i - pjm
        pjp = min(pmb+n_i,pmmax)
        ajp = pjp - n_i - pm

!     K-dim
        pkm = max(pm-n_i*n_j,1)
        akm = pmb - n_i*n_j - pkm
        pkp = min(pmb+n_i*n_j,pmmax)
        akp = pkp - n_i*n_j - pm

        IF (am>0) GO TO 100

        RETURN
      END

!>    @brief BLAS adapting : z=z+x*y
!>    @param[in] N length of vectors [z],[x],[y]
!>    @param[in] x vector [x]
!>    @param[in] y vector [y]
!>    @param[in] z vector [z]
!>    @param[out] z vector [z]
      SUBROUTINE dxypz(n,x,y,z)
        IMPLICIT NONE
        INTEGER n, i, m
        DOUBLE PRECISION x(n), y(n), z(n)

        m = mod(n,4)
        DO i = 1, m
          z(i) = z(i) + x(i)*y(i)
        END DO

        m = m + 1
        DO i = m, n, 4
          z(i+0) = z(i+0) + x(i+0)*y(i+0)
          z(i+1) = z(i+1) + x(i+1)*y(i+1)
          z(i+2) = z(i+2) + x(i+2)*y(i+2)
          z(i+3) = z(i+3) + x(i+3)*y(i+3)
        END DO

        RETURN
      END


!>    @brief BLAS adapting : z=x*y
!>    @param[in] N length of vectors [z],[x],[y]
!>    @param[in] x vector [x]
!>    @param[in] y vector [y]
!>    @param[out] z vector [z]
      SUBROUTINE dxyz(n,x,y,z)
        IMPLICIT NONE
        INTEGER n, i, m
        DOUBLE PRECISION x(n), y(n), z(n)

        m = mod(n,4)
        DO i = 1, m
          z(i) = x(i)*y(i)
        END DO

        m = m + 1
        DO i = m, n, 4
          z(i+0) = x(i+0)*y(i+0)
          z(i+1) = x(i+1)*y(i+1)
          z(i+2) = x(i+2)*y(i+2)
          z(i+3) = x(i+3)*y(i+3)
        END DO

        RETURN
      END
