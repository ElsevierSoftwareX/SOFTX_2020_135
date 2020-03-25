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

!>    @brief gives the parameter name (ascii) with index "seedi"
!>    @param[in] seedi parameter component/index
!>    @param[in] ismpl local sample index
!>    @param[out] pname parameter name
      SUBROUTINE param_name(seedi,pname,ismpl)
        use arrays
        use mod_genrl
        IMPLICIT NONE
        integer :: ismpl
        character (len=16) :: pname
        INTEGER seedi, s_k, s_u

        IF ((seedi<1) .OR. (seedi>mpara)) THEN
          WRITE(*,*) 'error: in "param_name", no valid index=', &
            seedi, ' out of range'
          STOP
        END IF

        s_k = seed_para(1,seedi)
        s_u = seed_para(2,seedi)
!     time period handling
        IF (s_k==-1) s_k = idx_tp
!     single cell bc handling
        IF (s_k==-2) s_k = idx_sbc

        IF (s_k>nprop+2) THEN
          WRITE(*,*) 'error : in "param_name", no valid "s_k=', s_k, &
            '"'
          STOP
        END IF

        pname = ' '

200     FORMAT (1A4,'_unit',1I7.7)
        WRITE(pname,200) properties(s_k), s_u

        RETURN
      END
