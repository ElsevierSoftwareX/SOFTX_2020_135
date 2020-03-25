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

!> @brief parallelisation wrapper for "omp_old_save"
!> @param[in] level level number
!> @param[in] ismpl local sample index
      subroutine old_save(level,ismpl)

        use mod_OMP_TOOLS

        implicit none

        integer :: ismpl

        include 'OMP_TOOLS.inc'

        integer, intent (in) :: level


#ifdef fOMP
!$OMP parallel num_threads(Tlevel_1)
!$      call omp_binding(abs(ismpl))
#endif

        call omp_old_save(level,ismpl)

#ifdef fOMP
!$OMP end parallel
#endif
!
        return

      end subroutine old_save

!> @brief save current state as an old version
!> @param[in] level level number (which old version)
!> @param[in] ismpl local sample index
      subroutine omp_old_save(level,ismpl)

        use arrays
        use mod_genrl
        use mod_conc

        implicit none

        ! local sample index
        integer :: ismpl

        ! cgen level index
        integer, intent (in) :: level

        ! directional cell-indices
        integer :: i, j, k

        ! counter for concentration
        integer :: l

        ! Start position in array for process
        integer :: tpos

        ! Number of array elements for process
        integer :: tanz

        ! OpenMP partition, get tpos, tanz
        call omp_part(i0*j0*k0,tpos,tanz)

        ! Get i, j, k indices for tpos
        call ijk_m(tpos,i,j,k)

        ! Copy all variable arrays to old: var -> varold
        call dcopy(tanz,head(i,j,k,abs(ismpl)),1, headold(tpos,level,max(1,ismpl)), 1)
        call dcopy(tanz,temp(i,j,k,abs(ismpl)),1, tempold(tpos,level,max(1,ismpl)), 1)
        call dcopy(tanz,pres(i,j,k,abs(ismpl)),1, presold(tpos,level,max(1,ismpl)), 1)
        do l = 1, ntrans
          call dcopy(tanz,conc(i,j,k,l,abs(ismpl)),1, concold(tpos,l,level,max(1,ismpl)),1)
        end do

        return

      end subroutine omp_old_save
