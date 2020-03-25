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
!>    @param[in] delxa I0-direction absolute cell position
!>    @param[in] delya J0-direction absolute cell position
!>    @param[in] delza K0-direction absolute cell position
!>    @return interpolated value
      DOUBLE PRECISION FUNCTION interpolatepol(i0,j0,k0, i,j,k, vals, px,py,pz, delx,dely,delz, delxa,delya,delza)
        IMPLICIT NONE
!       i,j,k cell index (first corner - 3D)
        INTEGER i, j, k, i0, j0, k0, mode3d
!       corrected i,j,k, can be eliminated later - px,py,pz is inside the [i,j,k] cell ([0,0,0] corner)
        INTEGER c_i, c_j, c_k
!       state variables (values)
        DOUBLE PRECISION vals(i0,j0,k0)
!       x,y,z cell dimensions (delta size)
        DOUBLE PRECISION delx(i0), dely(j0), delz(k0)
!       x,y,z absolute cell position
        DOUBLE PRECISION delxa(i0), delya(j0), delza(k0)
!       x,y,z interpolation position
        DOUBLE PRECISION px, py, pz
!
        DOUBLE PRECISION interpolate_3d, interpolate_am
        EXTERNAL interpolate_3d, interpolate_am
        INTRINSIC max, min

!       correct i-index
        c_i = i
        IF (px<delxa(i)) c_i = i-1
        IF (px>delxa(min(i+1,i0))) c_i = i+1
        c_i = min(max(c_i,1),i0)
!       correct j-index
        c_j = j
        IF (py<delya(j)) c_j = j-1
        IF (py>delya(min(j+1,j0))) c_j = j+1
        c_j = min(max(c_j,1),j0)
!       correct k-index
        c_k = k
        IF (pz<delza(k)) c_k = k-1
        IF (pz>delza(min(k+1,k0))) c_k = k+1
        c_k = min(max(c_k,1),k0)
!
        mode3d = 3
        interpolatepol = interpolate_3d(i0,j0,k0, c_i,c_j,c_k, mode3d, vals, px,py,pz, delx,dely,delz, delxa,delya,delza)
!
        RETURN
      END

!>    @brief inter cell interpolation (piece wise ???-methode), mode version
!>    @param[in] I0 i-dimension
!>    @param[in] J0 j-dimension
!>    @param[in] K0 k-dimension
!>    @param[in] i cell index, direction I0
!>    @param[in] j cell index, direction J0
!>    @param[in] k cell index, direction K0
!>    @param[in] mode3d direction mode, 1==I0, 2==J0, 3==K0
!>    @param[in] vals state variables (values)
!>    @param[in] px I0-direction interpolation position
!>    @param[in] py J0-direction interpolation position
!>    @param[in] pz K0-direction interpolation position
!>    @param[in] delx I0-direction cell dimensions (delta size)
!>    @param[in] dely J0-direction cell dimensions (delta size)
!>    @param[in] delz K0-direction cell dimensions (delta size)
!>    @param[in] delxa I0-direction absolute cell position
!>    @param[in] delya J0-direction absolute cell position
!>    @param[in] delza K0-direction absolute cell position
!>    @return piece wise ???-methode interpolated value
      DOUBLE PRECISION RECURSIVE FUNCTION interpolate_3d(i0,j0,k0, i,j,k, mode3d, vals, px,py,pz, &
        delx,dely,delz, delxa,delya,delza)
        IMPLICIT NONE
!       i,j,k cell index (first corner - 3D)
        INTEGER i, j, k, i0, j0, k0, mode3d
!       state variables (values)
        DOUBLE PRECISION vals(i0,j0,k0)
!       x,y,z cell dimensions (delta size)
        DOUBLE PRECISION delx(i0), dely(j0), delz(k0)
!       x,y,z absolute cell position
        DOUBLE PRECISION delxa(i0), delya(j0), delza(k0)
!       x,y,z interpolation position
        DOUBLE PRECISION px, py, pz
!
        DOUBLE PRECISION d_a2, d_a3, d_a4, d_a5, d_d23, d_d34, d_d45, dp
        DOUBLE PRECISION pw_uk_interpol, pw_cr_interpol, interpolate_3d_wrapper, lin_interpol
        EXTERNAL pw_uk_interpol, pw_cr_interpol, interpolate_3d_wrapper, lin_interpol
        INTRINSIC min, max

        IF (mode3d==0) THEN
          interpolate_3d = vals(i,j,k)
          RETURN
        END IF
        IF (mode3d==1) THEN
!         prepare X-direction interpolation, compute directly
!         cell size deltas
          d_d23 = (delx(max(i-1,1))+delx(i))*0.5d0
          d_d34 = (delx(i)+delx(min(i+1,I0)))*0.5d0
          d_d45 = (delx(min(i+1,I0))+delx(min(i+2,I0)))*0.5d0
!         (x) delta
          dp = px -delxa(i)
!         function values
          d_a2 = interpolate_3d_wrapper(i0,j0,k0, max(i-1,1),j,k, mode3d-1, vals, px,py,pz, delx,dely,delz, delxa,delya,delza)
          d_a3 = interpolate_3d_wrapper(i0,j0,k0, i,j,k, mode3d-1, vals, px,py,pz, delx,dely,delz, delxa,delya,delza)
          d_a4 = interpolate_3d_wrapper(i0,j0,k0, min(i+1,I0),j,k, mode3d-1, vals, px,py,pz, delx,dely,delz, delxa,delya,delza)
          d_a5 = interpolate_3d_wrapper(i0,j0,k0, min(i+2,I0),j,k, mode3d-1, vals, px,py,pz, delx,dely,delz, delxa,delya,delza)
!         sanity check
          IF (px<delxa(i).OR.(px>delxa(min(i+1,I0)).AND.i<I0)) THEN
            WRITE(*,'(1a,1i7,1a,1g16.8)') 'error: i-cell index (',i, &
              ') incorrect for position :',px
            STOP
          END IF
        END IF
        IF (mode3d==2) THEN
!         prepare Y-direction interpolation, consist of X-direction interpolations
!         cell size deltas
          d_d23 = (dely(max(j-1,1))+dely(j))*0.5d0
          d_d34 = (dely(j)+dely(min(j+1,J0)))*0.5d0
          d_d45 = (dely(min(j+1,J0))+dely(min(j+2,J0)))*0.5d0
!         (y) delta
          dp = py -delya(j)
!         function values
          d_a2 = interpolate_3d_wrapper(i0,j0,k0, i,max(j-1,1),k, mode3d-1, vals, px,py,pz, delx,dely,delz, delxa,delya,delza)
          d_a3 = interpolate_3d_wrapper(i0,j0,k0, i,j,k, mode3d-1, vals, px,py,pz, delx,dely,delz, delxa,delya,delza)
          d_a4 = interpolate_3d_wrapper(i0,j0,k0, i,min(j+1,J0),k, mode3d-1, vals, px,py,pz, delx,dely,delz, delxa,delya,delza)
          d_a5 = interpolate_3d_wrapper(i0,j0,k0, i,min(j+2,J0),k, mode3d-1, vals, px,py,pz, delx,dely,delz, delxa,delya,delza)
!         sanity check
          IF (py<delya(j).OR.(py>delya(min(j+1,J0)).AND.j<J0)) THEN
            WRITE(*,'(1a,1i7,1a,1g16.8)') 'error: j-cell index (',j, &
              ') incorrect for position :',py
            STOP
          END IF
        END IF
        IF (mode3d==3) THEN
!         prepare Z-direction interpolation, consist of Y-direction interpolations
!         cell size deltas
          d_d23 = (delz(max(k-1,1))+delz(k))*0.5d0
          d_d34 = (delz(k)+delz(min(k+1,K0)))*0.5d0
          d_d45 = (delz(min(k+1,K0))+delz(min(k+2,K0)))*0.5d0
!         (z) delta
          dp = pz -delza(k)
!         function values
          d_a2 = interpolate_3d_wrapper(i0,j0,k0, i,j,max(k-1,1), mode3d-1, vals, px,py,pz, delx,dely,delz, delxa,delya,delza)
          d_a3 = interpolate_3d_wrapper(i0,j0,k0, i,j,k, mode3d-1, vals, px,py,pz, delx,dely,delz, delxa,delya,delza)
          d_a4 = interpolate_3d_wrapper(i0,j0,k0, i,j,min(k+1,K0), mode3d-1, vals, px,py,pz, delx,dely,delz, delxa,delya,delza)
          d_a5 = interpolate_3d_wrapper(i0,j0,k0, i,j,min(k+2,K0), mode3d-1, vals, px,py,pz, delx,dely,delz, delxa,delya,delza)
!         sanity check
          IF (pz<delza(k).OR.(pz>delza(min(k+1,K0)).AND.k<K0)) THEN
            WRITE(*,'(1a,1i7,1a,1g16.8)') 'error: k-cell index (',k, &
              ') incorrect for position :',pz
            STOP
          END IF
        END IF
!
!       compute unknown-methode interpolation
        interpolate_3d = pw_uk_interpol(d_a2, d_a3, d_a4, d_a5, dp, d_d23, d_d34, d_d45)
!.UNUSED/ compute Catmull-Rom interpolation
!        interpolate_3d = pw_cr_interpol(d_a2, d_a3, d_a4, d_a5, dp, d_d23, d_d34, d_d45)
!.UNUSED/ compute linear interpolation
!        interpolate_3d = lin_interpol(d_a3, d_a4, dp, d_d34)
!
        RETURN
      END

!>    @brief wrapper routine for the inter cell interpolation (piece wise ???-methode), mode version
!>    @param[in] I0 i-dimension
!>    @param[in] J0 j-dimension
!>    @param[in] K0 k-dimension
!>    @param[in] i cell index, direction I0
!>    @param[in] j cell index, direction J0
!>    @param[in] k cell index, direction K0
!>    @param[in] mode3d direction mode, 1==I0, 2==J0, 3==K0
!>    @param[in] vals state variables (values)
!>    @param[in] px I0-direction interpolation position
!>    @param[in] py J0-direction interpolation position
!>    @param[in] pz K0-direction interpolation position
!>    @param[in] delx I0-direction cell dimensions (delta size)
!>    @param[in] dely J0-direction cell dimensions (delta size)
!>    @param[in] delz K0-direction cell dimensions (delta size)
!>    @param[in] delxa I0-direction absolute cell position
!>    @param[in] delya J0-direction absolute cell position
!>    @param[in] delza K0-direction absolute cell position
!>    @return piece wise ???-methode interpolated value
      DOUBLE PRECISION RECURSIVE FUNCTION interpolate_3d_wrapper(i0,j0,k0, i,j,k, mode3d, vals, px,py,pz, &
        delx,dely,delz, delxa,delya,delza)
        IMPLICIT NONE
!       i,j,k cell index (first corner - 3D)
        INTEGER i, j, k, i0, j0, k0, mode3d
!       state variables (values)
        DOUBLE PRECISION vals(i0,j0,k0)
!       x,y,z cell dimensions (delta size)
        DOUBLE PRECISION delx(i0), dely(j0), delz(k0)
!       x,y,z absolute cell position
        DOUBLE PRECISION delxa(i0), delya(j0), delza(k0)
!       x,y,z interpolation position
        DOUBLE PRECISION px, py, pz
!
        DOUBLE PRECISION interpolate_3d
        EXTERNAL interpolate_3d
        interpolate_3d_wrapper = interpolate_3d(i0,j0,k0, i,j,k, mode3d, vals, px,py,pz, delx,dely,delz, delxa,delya,delza)
        RETURN
      END

!>    @brief piece wise unknown-methode interpolation
!>    @param[in] f_am value x_a-1
!>    @param[in] f_a value x_a
!>    @param[in] f_b value x_b
!>    @param[in] f_bp value x_b+1
!>    @param[in] dp relative position >= 0, for [x_a x_b] intervall, difference from x_a
!>    @param[in] d_aam (delta) cell size, for [x_a-1 x_a] interval
!>    @param[in] d_ba (delta) cell size, for [x_a x_b] interval
!>    @param[in] d_bpb (delta) cell size, for [x_b x_b+1] interval
!>    @return interpolated value
      DOUBLE PRECISION FUNCTION pw_uk_interpol(f_am, f_a, f_b, f_bp, dp, d_aam, d_ba, d_bpb)
        IMPLICIT NONE
!       f_* function values for all nodes
!       dp: relative position >= 0
!       d_* (delta) cell size, for all x_* intervals
        DOUBLE PRECISION f_am, f_a, f_b, f_bp, dp, d_aam, d_ba, d_bpb
        DOUBLE PRECISION ds, gf_sm, gf_xm, gf_am, gf_bm, i_s, gf_s
        DOUBLE PRECISION gf_b_a
        INTRINSIC max, min, dabs
        LOGICAL test_null
        EXTERNAL test_null

! g_s = (2.0d0*f_b -2.0d0*f_a -gf_am*(ds-a) -gf_bm*(b-ds))/(b-a)
! ds = (g_s*(b-a) -2.0d0*f_b +2.0d0*f_a -gf_am*a +gf_bm*b) /(gf_bm -gf_am)
! mit a:=0, b:=d_ba ->
! g_s = (f_b+f_b -f_a-f_a -gf_am*ds -gf_bm*(d_ba-ds))/d_ba
! ds = (g_s*d_ba -f_b-f_b +f_a+f_a +gf_bm*d_ba) /(gf_bm -gf_am)

!       mean gradient between [a b]
        gf_s = (f_b -f_a) /d_ba
!       gradient mean for (a)
        gf_am = 0.5d0 *((f_a -f_am) /d_aam +gf_s)
!       gradient mean for (b)
        gf_bm = 0.5d0 *(gf_s +(f_bp -f_b) /d_bpb)
!
        gf_b_a = gf_bm -gf_am
!       delta (x) sattle point
        ds = 0.5d0*d_ba
!       gradient for the sattle point (s)
        gf_sm = 2.0d0*gf_s -0.5d0 *(gf_am +gf_bm)
!
        IF (dp<ds) THEN
!         gradient for (x)
          gf_xm = gf_am +dp *(gf_sm -gf_am) /ds
!         f_x, result for (x)
          pw_uk_interpol = f_a +(gf_am +gf_xm) *dp *0.5d0
        ELSE
!         integral area [a s]
          i_s = (gf_sm +gf_am) *ds *0.5d0
!         gradient for (x)
          gf_xm = gf_sm +(dp -ds) *(gf_bm -gf_sm) /(d_ba -ds)
!         f_x, result for (x)
          pw_uk_interpol = f_a +(gf_sm +gf_xm) *(dp -ds) *0.5d0 +i_s
        ENDIF
!
        RETURN
      END
