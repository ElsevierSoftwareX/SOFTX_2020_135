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

!> @brief allocate memory of objectes used for VISIM
      SUBROUTINE alloc_visim()
        USE simul_arrays
        use mod_genrl
        use mod_simul
        use mod_linfos
        IMPLICIT NONE

        IF (max_gpara==0) RETURN
        ALLOCATE(sim(i0*j0*k0))
        ALLOCATE(lvm(i0*j0*k0))
        ALLOCATE(avepor(i0,j0))
        memory = memory + i0*j0*k0 + i0*j0*k0 + i0*j0
        ALLOCATE(tmp(i0*j0*k0))
        ALLOCATE(order(i0*j0*k0))
        ALLOCATE(dvr(i0*j0*k0))
        ALLOCATE(simout(i0*j0*k0))
        memory = memory + 4*i0*j0*k0
        ALLOCATE(porvar(i0,j0))
        ALLOCATE(krgvar(i0*j0*k0))
        memory = memory + i0*j0 + i0*j0*k0
        ALLOCATE(novar(i0*j0*k0))
        memory = memory + i0*j0*k0
        ALLOCATE(cd2v(i0*j0*k0,maxvols))
        memory = memory + i0*j0*k0*maxvols

        RETURN
      END

!> @brief free memory of objects used for VISIM
      SUBROUTINE dealloc_visim()
        USE simul_arrays
        use mod_genrl
        use mod_simul
        use mod_linfos
        IMPLICIT NONE

        IF (max_gpara==0) RETURN
        DEALLOCATE(sim)
        DEALLOCATE(lvm)
        DEALLOCATE(avepor)
        memory = memory - i0*j0*k0 - i0*j0*k0 - i0*j0
        DEALLOCATE(tmp)
        DEALLOCATE(order)
        DEALLOCATE(dvr)
        DEALLOCATE(simout)
        memory = memory - 4*i0*j0*k0
        DEALLOCATE(porvar)
        DEALLOCATE(krgvar)
        memory = memory - i0*j0 - i0*j0*k0
        DEALLOCATE(novar)
        memory = memory - i0*j0*k0
        DEALLOCATE(cd2v)
        memory = memory - i0*j0*k0*maxvols

        RETURN
      END
