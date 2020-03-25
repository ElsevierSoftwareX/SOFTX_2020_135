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

!>    @brief inititialise cache/block size
!>    @details
!> initialise cache size and block parameters\n
!> and numerical basis constants (zero values, precision ...)\n
      SUBROUTINE model_init(ismpl)
        USE arrays
        use mod_genrl
        use mod_blocking_size
        IMPLICIT NONE
        integer :: ismpl
        DOUBLE PRECISION dlamch
        EXTERNAL dlamch
        INTRINSIC dsqrt

!       !! defines lin. solver constants [part 1] !!
        bldiv_cg = 7.0D0
        bldiv_bicg(1) = 10.5D0
        bldiv_bicg(2) = 5.5D0
        bldiv_mvp = 11.5D0
        bldiv_dot(1) = 2.0D0
        bldiv_dot(2) = 6.5D0
        bldiv_dot(3) = 8.0D0
!
!       !! defines lin. solver constants [part 2] !!
#ifdef CLsun
!       cache size for SUN Ultra-Sparc-III (8 MByte) cpu
        cache_size = 384000
#elif CLopt
!       cache size for AMD Opteron (1 MByte) cpu
        cache_size = 1024000
#elif CLpen
!       cache size for Intel P4 (512 KByte) cpu
        cache_size = 512000
#else
!       cache size for Intel P4/Celeron (128 KByte) cpu
        cache_size = 128000
#endif
!       block size (number of doubles minus 20%)
        bl_size = cache_size / 10

!       initalise some values
        memory = 0.0D0

        const_dble(1) = dlamch('E')
!       3/4 of the min value
        const_dble(2) = dsqrt(dlamch('U'))
        const_dble(2) = const_dble(2)*dsqrt(const_dble(2))
!       3/4 of the max value
        const_dble(3) = dsqrt(dlamch('O'))
        const_dble(3) = const_dble(3)*dsqrt(const_dble(3))

        RETURN
      END
