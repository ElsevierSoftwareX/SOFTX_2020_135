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


!#>   @brief allocates space for variables for permafrost
!     icefrac
!     liq liquidus temperature in degree Celsius (above which ice is molten)
!     sol solidus temperature in degree Celsius (below which water is frozen)
!     i0 number of cells in x direction
!     j0 number of cells in y direction
!     k0 number of cells in z direction
      SUBROUTINE ice_allocate
                use ice
        use mod_genrl
        IMPLICIT NONE

        ALLOCATE(icefrac(i0,j0,k0))
        ALLOCATE(liq(i0,j0,k0))
        ALLOCATE(sol(i0,j0,k0))

        RETURN
      END

!#>   @brief deallocates space for
!     icefrac 
!     liq liquidus temperature in degree Celsius (above which ice is molten)
!     sol solidus temperature in degree Celsius (below which water is frozen)

      SUBROUTINE ice_deallocate
                use ice
        IMPLICIT NONE

        DEALLOCATE(icefrac)
        DEALLOCATE(liq)
        DEALLOCATE(sol)

        RETURN
      END
