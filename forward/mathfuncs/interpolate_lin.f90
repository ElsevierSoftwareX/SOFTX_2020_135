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

!>    @brief inter cell interpolation
!>    @param[in] I0 i-dimension
!>    @param[in] J0 j-dimension
!>    @param[in] K0 k-dimension
!>    @param[in] i cell index, direction I0
!>    @param[in] j cell index, direction J0
!>    @param[in] k cell index, direction K0
!>    @param[in] vals state variables (values)
!>    @param[in] px I0-direction interpolation position
!>    @param[in] py J0-direction interpolation position
!>    @param[in] pz K0-direction interpolation position
!>    @param[in] delx I0-direction cell dimensions (delta size)
!>    @param[in] dely J0-direction cell dimensions (delta size)
!>    @param[in] delz K0-direction cell dimensions (delta size)
!>    @param[in] delxa I0-direction absolute cell positions
!>    @param[in] delya J0-direction absolute cell positions
!>    @param[in] delza K0-direction absolute cell positions
!>    @return interpolated value
      DOUBLE PRECISION FUNCTION interpolatelin(i0,j0,k0, i,j,k, vals, px,py,pz, delx,dely,delz, delxa,delya,delza)
        IMPLICIT NONE
!       i,j,k cell index (first corner - 3D)
        INTEGER i, j, k, i0, j0, k0
!       neighbours i,j,k cell index (second corner - 3D)
        INTEGER ni, nj, nk
!       state variables (values)
        DOUBLE PRECISION vals(i0,j0,k0)
!       x,y,z cell dimensions (delta size)
        DOUBLE PRECISION delx(i0), dely(j0), delz(k0)
!       x,y,z absolute cell position
        DOUBLE PRECISION delxa(i0), delya(j0), delza(k0)
!       x,y,z interpolation position
        DOUBLE PRECISION px, py, pz
!
        DOUBLE PRECISION d_a0, d_b0, d_a1, d_b1, d_a2, d_b2
        DOUBLE PRECISION lin_interpol
        EXTERNAL lin_interpol
        INTRINSIC min, max, abs

!       get neighbour i-index
        ni = i+1
        IF (px<delxa(i)) ni = i-1
        ni = min(max(ni,1),i0)
!       get neighbour j-index
        nj = j+1
        IF (py<delya(j)) nj = j-1
        nj = min(max(nj,1),j0)
!       get neighbour k-index
        nk = k+1
        IF (pz<delza(k)) nk = k-1
        nk = min(max(nk,1),k0)
!
!       x1
        d_a0 = vals(i,j,k)
!       x2
        d_b0 = vals(ni,j,k)
!       y1
        d_a1 = lin_interpol(d_a0,d_b0,abs(delxa(i)-px),(delx(i)+delx(ni))*0.5d0)
!
!       x1
        d_a0 = vals(i,nj,k)
!       x2
        d_b0 = vals(ni,nj,k)
!       y2
        d_b1 = lin_interpol(d_a0,d_b0,abs(delxa(i)-px),(delx(i)+delx(ni))*0.5d0)
!
!       z1
        d_a2 = lin_interpol(d_a1,d_b1,abs(delya(j)-py),(dely(j)+dely(nj))*0.5d0)
!
!       x1
        d_a0 = vals(i,j,nk)
!       x2
        d_b0 = vals(ni,j,nk)
!       y1
        d_a1 = lin_interpol(d_a0,d_b0,abs(delxa(i)-px),(delx(i)+delx(ni))*0.5d0)
!
!       x1
        d_a0 = vals(i,nj,nk)
!       x2
        d_b0 = vals(ni,nj,nk)
!       y2
        d_b1 = lin_interpol(d_a0,d_b0,abs(delxa(i)-px),(delx(i)+delx(ni))*0.5d0)
!
!       z2
        d_b2 = lin_interpol(d_a1,d_b1,abs(delya(j)-py),(dely(j)+dely(nj))*0.5d0)
!
        interpolatelin = lin_interpol(d_a2,d_b2,abs(delza(k)-pz),(delz(k)+delz(nk))*0.5d0)
!
        RETURN
      END

!>    @brief linear interpolation
!>    @param[in] da value A
!>    @param[in] db value B
!>    @param[in] dp relative position
!>    @param[in] dd (delta) cell size
!>    @return linear interpolated value
      DOUBLE PRECISION FUNCTION lin_interpol(da, db, dp, dd)
        IMPLICIT NONE
!       da: value A; db: value B; dp: relative position; dd: (delta) cell size
        DOUBLE PRECISION da, db, dp, dd
        lin_interpol = (da*(dd-dp) + db*dp)/dd
        RETURN
      END
