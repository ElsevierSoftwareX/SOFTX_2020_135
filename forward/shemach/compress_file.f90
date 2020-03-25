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

!>    @brief compress a file (to *.bz2, *.zip, *.gz)
!>    @param[in] fname file to compress
!>    @param[in] ctool index number of the compression tool
      SUBROUTINE compress_file(ctool,fname)
        use arrays
        IMPLICIT NONE
        character (len=*) :: fname
        INTEGER ctool
        INTRINSIC trim

!     no compression
        IF (compress_suffix(ctool)=='plain') RETURN
!     bzip2
        IF (compress_suffix(ctool)=='bz2') CALL system('bzip2 -f -9 "' &
          //trim(fname)//'"')
!     gnu zip
        IF (compress_suffix(ctool)=='gz') CALL system('gzip -f -9 "'// &
          trim(fname)//'"')
!     std. zip
        IF (compress_suffix(ctool)=='zip') CALL system('zip -m -9 "'// &
          trim(fname)//'.zip" "'//trim(fname)//'"')
        RETURN
      END
