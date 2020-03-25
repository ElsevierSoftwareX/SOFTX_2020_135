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

!>    @brief allocate data arrays (dynamic size)
!>    @param[in] ismpl local sample index
!>    @details
      SUBROUTINE alloc_data(ismpl)
        use arrays
        use mod_genrl
        use mod_genrlc
        use mod_data
        use mod_linfos
        IMPLICIT NONE
        integer :: ismpl
        INTRINSIC max

        IF (linfos(1)>=2) THEN
          WRITE(*,*) ' '
          WRITE(*,*) ' [I] : ... alloc_data'
          WRITE(*,*) ' '
        END IF
!
        IF (ndata>=1) THEN
          ALLOCATE(ddata(ndata,n_ddata))
          memory = memory + ndata*n_ddata
          ALLOCATE(idata(ndata,n_idata))
          memory = memory + ndata*n_idata
          ALLOCATE(resid(ndata))
          memory = memory + ndata
        END IF
!
        DEALLOCATE(sdata)
        memory = memory - nsmpl
        ALLOCATE(sdata(max(ndata,1),nsmpl))
        memory = memory + max(ndata,1)*nsmpl
!
        RETURN
      END SUBROUTINE alloc_data


!>    @brief free data arrays (with dynamic size)
!>    @param[in] ismpl local sample index
!>    @details
      SUBROUTINE dealloc_data(ismpl)
        use arrays
        use mod_genrl
        use mod_genrlc
        use mod_data
        use mod_linfos
        IMPLICIT NONE
        integer :: ismpl
        INTRINSIC max

        IF (linfos(1)>=2) THEN
          WRITE(*,*) ' '
          WRITE(*,*) ' [I] : ... dealloc_data'
          WRITE(*,*) ' '
        END IF
!
        IF (ndata>=1) THEN
          DEALLOCATE(ddata)
          memory = memory - ndata*n_ddata
          DEALLOCATE(idata)
          memory = memory - ndata*n_idata
          DEALLOCATE(resid)
          memory = memory - ndata
        END IF
        DEALLOCATE(sdata)
        memory = memory - max(ndata,1)*nsmpl
!
        RETURN
      END SUBROUTINE dealloc_data
