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

!#>   @brief <description>
!#>   @param[in] p <description>
!#>   @param[in] t <description>
!#>   @param[in] ismpl local sample index
!#>   @return rho <description>
!#>   @details
      DOUBLE PRECISION FUNCTION rhow(p,t)
        IMPLICIT NONE
        DOUBLE PRECISION t, p, p2, t2, t4, t6
        DOUBLE PRECISION a, b, c, d, e, f, g
        DOUBLE PRECISION a0, a1, a2, b0, b1, b2, c0, c1, c2, d0, d1, &
          d2, e0, e1, e2, f0, f1, f2, g0, g1, g2
        DATA a0, a1, a2, b0, b1, b2, c0, c1, c2, d0, d1, d2, e0, e1, &
          e2, f0, f1, f2, g0, g1, g2/ + 9.99792877961606D+02, &
          + 5.07605113140940D-04, -5.2842547816413D-10, &
          + 5.13864847162196D-02, -3.61991396354483D-06, &
          + 7.97204102509724D-12, -7.53557031774437D-03, &
          + 6.37212093275576D-05, -1.66203631393248D-13, &
          + 4.60380647957350D-05, -5.61299059722121D-10, &
          + 1.80924436489400D-15, -2.26651454175013D-07, &
          + 3.36874416675978D-12, -1.30352149261326D-17, &
          + 6.14889851856743D-10, -1.06165223196756D-14, &
          + 4.75014903737416D-20, -7.39221950969522D-13, &
          + 1.42790422913922D-17, -7.13130230531541D-23/


        p2 = p*p
        t2 = t*t
        t4 = t2*t2
        t6 = t4*t2

! liquid density
        a = a0 + a1*p + a2*p2
        b = b0 + b1*p + b2*p2
        c = c0 + c1*p + c2*p2
        d = d0 + d1*p + d2*p2
        e = e0 + e1*p + e2*p2
        f = f0 + f1*p + f2*p2
        g = g0 + g1*p + g2*p2

        rhow = a + b*t + c*t2 + d*t2*t + e*t4 + f*t4*t + g*t6

        RETURN
      END
      

!#>   @brief <description>
!#>   @param[in] p <description>
!#>   @param[in] t <description>
!#>   @param[in] ismpl local sample index
!#>   @return rho <description>
!#>   @details
      DOUBLE PRECISION FUNCTION drhowdt(p,t)
        IMPLICIT NONE
        DOUBLE PRECISION t, p
        DOUBLE PRECISION a0, a1, a2, b0, b1, b2, c0, c1, c2, d0, d1, &
          d2, e0, e1, e2, f0, f1, f2, g0, g1, g2
        DATA a0, a1, a2, b0, b1, b2, c0, c1, c2, d0, d1, d2, e0, e1, &
          e2, f0, f1, f2, g0, g1, g2/ + 9.99792877961606D+02, &
          + 5.07605113140940D-04, -5.2842547816413D-10, &
          + 5.13864847162196D-02, -3.61991396354483D-06, &
          + 7.97204102509724D-12, -7.53557031774437D-03, &
          + 6.37212093275576D-05, -1.66203631393248D-13, &
          + 4.60380647957350D-05, -5.61299059722121D-10, &
          + 1.80924436489400D-15, -2.26651454175013D-07, &
          + 3.36874416675978D-12, -1.30352149261326D-17, &
          + 6.14889851856743D-10, -1.06165223196756D-14, &
          + 4.75014903737416D-20, -7.39221950969522D-13, &
          + 1.42790422913922D-17, -7.13130230531541D-23/

! pure water density temperature derivative

        drhowdt = b0+b1*p+t*(c0+c1*p+c2*p**2)*2d0+b2*p**2+t**2.d0*(d0+d1*p+d2*p**2)*3.d0+ &
                  t**3*(e0+e1*p+e2*p**2)*4.d0+t**4*(f0+f1*p+f2*p**2)*5.d0+ &
                  t**5*(g0+g1*p+g2*p**2)*6.d0

        RETURN
      END
      
!#>   @brief <description>
!#>   @param[in] p <description>
!#>   @param[in] t <description>
!#>   @param[in] ismpl local sample index
!#>   @return rho <description>
!#>   @details
      DOUBLE PRECISION FUNCTION drhowdp(p,t)
        IMPLICIT NONE
        DOUBLE PRECISION t, p
        DOUBLE PRECISION a0, a1, a2, b0, b1, b2, c0, c1, c2, d0, d1, &
          d2, e0, e1, e2, f0, f1, f2, g0, g1, g2
        DATA a0, a1, a2, b0, b1, b2, c0, c1, c2, d0, d1, d2, e0, e1, &
          e2, f0, f1, f2, g0, g1, g2/ + 9.99792877961606D+02, &
          + 5.07605113140940D-04, -5.2842547816413D-10, &
          + 5.13864847162196D-02, -3.61991396354483D-06, &
          + 7.97204102509724D-12, -7.53557031774437D-03, &
          + 6.37212093275576D-05, -1.66203631393248D-13, &
          + 4.60380647957350D-05, -5.61299059722121D-10, &
          + 1.80924436489400D-15, -2.26651454175013D-07, &
          + 3.36874416675978D-12, -1.30352149261326D-17, &
          + 6.14889851856743D-10, -1.06165223196756D-14, &
          + 4.75014903737416D-20, -7.39221950969522D-13, &
          + 1.42790422913922D-17, -7.13130230531541D-23/

! pure water density pressure derivative

        drhowdp = b0+b1*p+t*(c0+c1*p+c2*p**2)*2d0+b2*p**2+ &
                  t**2*(d0+d1*p+d2*p**2)*3d0+t**3*(e0+e1*p+e2*p**2)*4d0+ &
                  t**4*(f0+f1*p+f2*p**2)*5d0+t**5*(g0+g1*p+g2*p**2)*6d0

        RETURN
      END

