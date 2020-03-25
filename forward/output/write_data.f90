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

!>    @brief shows observation data differences
!>    @param[in] ident index number for file name
!>    @param[in] ismpl local sample index
      SUBROUTINE write_data(ident,ismpl)
        use arrays
        use mod_genrl
        use mod_genrlc
        use mod_flow
        use mod_temp
        use mod_data
        use mod_linfos
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k, l
        INCLUDE 'version.inc'
!
        DOUBLE PRECISION res, val, resn
        INTEGER i1s, i2s, i1, i2, i3, i4, ui, i_type, ozone, ident, &
          lblank, i_si, lout
        EXTERNAL lblank
        character (len=80) :: filename
        character (len=20) :: snumber


        IF (write_disable) RETURN
#ifdef NOOUT
        RETURN
#endif
        IF (ndata<=0) RETURN

        CALL omp_new_file_handler(lout,16)

        CALL chln(project,i1,i2)
        CALL chln(project_sfx(ismpl),i1s,i2s)
        IF (ident>=0) THEN
          WRITE(snumber,'(1I7)') ident
        ELSE IF (ident==-1) THEN
          WRITE(snumber,'(A20)') 'final'
        ELSE IF (ident==-2) THEN
          WRITE(snumber,'(A20)') 'debug'
        ELSE IF (ident==-3) THEN
          WRITE(snumber,'(A20)') 'ens_mean'
        ELSE IF (ident==-4) THEN
          WRITE(snumber,'(A20)') 'mean'
        ELSE IF (ident==-5) THEN
          WRITE(snumber,'(A20)') 'ens_mean'
        END IF
        CALL chln(snumber,i3,i4)

        IF (i1s==0) THEN
          filename = project(i1:i2) // '_' // snumber(i3:i4) // '.dat'
        ELSE
          filename = project(i1:i2) // project_sfx(ismpl) (i1s:i2s) // &
            '_' // snumber(i3:i4) // '.dat'
        END IF

        OPEN(lout,file=filename,status='unknown',blank='null')

        IF (linfos(3)>=1) THEN
          WRITE(*,'(3A)') '  [W] : Data to "', &
            filename(1:lblank(filename)), '"'
        END IF

        WRITE(lout,'(2A)') '% Shemat-Suite version: ', version
        WRITE(lout,'(2A)') '% data fit at iteration ', snumber
        WRITE(lout,'(3A)') '%    i     j     k     unit     type   ', &
          '            calc             obs             err', &
          '            diff             res            time  ozone'

        DO l = 1, ndata
          i = idata(l,cid_i)
          j = idata(l,cid_j)
          k = idata(l,cid_k)
          ui = uindex(i,j,k)
          i_type = idata(l,cid_pv)
          i_si = idata(l,cid_si)
          ozone = idata(l,cid_obs)
          val = sdata(l,ismpl)

          res = val - ddata(l,cdd_pv)
          resn = res/ddata(l,cdd_w)
          WRITE(lout,'(3i6,3x,i6,3x,i6,3x,6e16.8,i6)') i, j, k, ui, &
            i_type, val, ddata(l,cdd_pv), ddata(l,cdd_w), res, resn, &
            ddata(l,cdd_time), ozone
        END DO

        CLOSE(lout)
        CALL omp_del_file_handler(lout)

        RETURN
      END
