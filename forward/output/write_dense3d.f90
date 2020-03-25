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

!>    @brief write values in fortran77 compressed text form
!>    @param[in] A value array to write
!>    @param[in] I0 i-dimension
!>    @param[in] J0 j-dimension
!>    @param[in] K0 k-dimension
!>    @param[in] CHAN file handler
!>    @param[in] ismpl local sample index
      SUBROUTINE write_dense3d(a,i0,j0,k0,chan,ismpl)
        IMPLICIT NONE
        INTEGER i0, j0, k0, i, chan, g, nijk, ismpl
        DOUBLE PRECISION a(i0*j0*k0), va, vn
!     Intel compiler buffer workaround
        INTEGER maxb, l10, lenb
        PARAMETER (lenb=65536)
!     character (len=lenB) :: rowbuffer
        character (len=65536) :: rowbuffer
!     format buffer
        character (len=25) :: element
        INTRINSIC dabs, int, log, dble

!     init
        nijk = i0*j0*k0
        IF (nijk<=0) RETURN

        g = 1
        maxb = 1
        va = a(1)
        IF (dabs(va)<1.0D-99) va = 0.0D0
        rowbuffer = ' '

!     compute minimal place holder for the counter
        l10 = int(log(dble(nijk+1))/log(10.0D0)) + 1
        element = '(I' // achar(l10+48) // ',A1,SP,1e24.17,A1)'

        DO i = 2, nijk
          vn = a(i)
          IF (dabs(vn)<1.0D-99) vn = 0.0D0
          IF (vn/=va) THEN
            IF (g==1) THEN
              WRITE(rowbuffer(maxb:maxb+25),'(A1,SP,e24.17,A1)') ' ', &
                va, ','
              maxb = maxb + 26
              g = 0
            ELSE
              WRITE(rowbuffer(maxb:maxb+25+l10),element) g, '*', va, &
                ','
              maxb = maxb + 26 + l10
              g = 0
            END IF
            IF (maxb>=(lenb-26-l10)) THEN
              WRITE(chan,'(A)') rowbuffer(1:maxb-1)
              maxb = 1
            END IF
            va = vn
          END IF
          g = g + 1
        END DO
        IF (g==1) THEN
          WRITE(rowbuffer(maxb:maxb+24),'(A1,SP,1e24.17)') ' ', va
          maxb = maxb + 25
        ELSE
          WRITE(rowbuffer(maxb:maxb+25+l10),element) g, '*', va, ' '
          maxb = maxb + 26 + l10
        END IF
        WRITE(chan,'(A)') rowbuffer(1:maxb-1)

        RETURN
      END
