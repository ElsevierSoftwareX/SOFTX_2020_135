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

!> @brief allocate memory of objects used for simulation (SGSIM and VISIM)
!> @param[in] ismpl local sample index
      SUBROUTINE alloc_simul(ismpl)
        use arrays
        USE simul_arrays
        use mod_genrl
        use mod_time
        use mod_linfos
        IMPLICIT NONE
        integer :: ismpl

        integer :: i_max
        INTRINSIC max

        IF (linfos(1)>=2) WRITE(*,*) ' [I] : ... alloc_simul'

        i_max = max(maxunits,bc_maxunits)
        ALLOCATE(propunitold(i_max,nprop))
        memory = memory - i_max*nprop
        ALLOCATE(bcperiodold(ngsmax,3,max(nbctp,1)))
        memory = memory + ngsmax*3*max(nbctp,1)
!
        ALLOCATE(main_input_master(mpara))
        memory = memory + mpara

! --- will be allocated in "read_simul" before ---
!      allocate(fnpara(1))
!      allocate(gpara(1))
!      allocate(seed_para(2,1))

        RETURN
      END

!> @brief free memory of objects used for simulation (SGSIM and VISIM)
!> @param[in] ismpl local sample index
      SUBROUTINE dealloc_simul(ismpl)
        use arrays
        USE simul_arrays
        use mod_genrl
        use mod_time
        use mod_linfos
        IMPLICIT NONE
        integer :: ismpl

        integer :: i_max
        INTRINSIC max

        IF (linfos(1)>=2) WRITE(*,*) ' [I] : ... dealloc_simul'

        i_max = max(maxunits,bc_maxunits)
        DEALLOCATE(propunitold)
        memory = memory - i_max*nprop
        DEALLOCATE(bcperiodold)
        memory = memory - ngsmax*3*max(nbctp,1)
!
        DEALLOCATE(main_input_master)
        memory = memory - mpara

        DEALLOCATE(seed_para)
        DEALLOCATE(gpara)
        DEALLOCATE(fnpara)

        RETURN
      END
