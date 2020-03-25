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

!>    @brief compares two different runs
!>    @param[in] ismpl local sample index
!>    @details
!> compare the main state-variables (head,temp,pres),\n
!> only available when configured\n
      SUBROUTINE compare_run(ismpl)
        use arrays
        use mod_genrl
        use mod_genrlc
        use mod_linfos
        IMPLICIT NONE
        integer :: ismpl
        integer :: i

        DOUBLE PRECISION shead, stemp, spres, mhead, mtemp, &
          mpres
        INTEGER lblank, i1, i2
        LOGICAL test_option
        character (len=80) :: filename
        character (len=7) :: sijk
        character (len=8) :: sfield
        EXTERNAL lblank, test_option
        INTRINSIC dsqrt


        CALL chln(project,i1,i2)
        filename = project(i1:i2) // '.compare'
        WRITE(sijk,'(1I7)') i0*j0*k0

        IF (test_option('init_comp')) THEN
          OPEN(76,file=filename,status='unknown',blank='null')
          IF (linfos(3)>=1) THEN
            WRITE(*,'(3A)') '  [W] : Comparison to "', &
              filename(1:lblank(filename)), '"'
          END IF
          WRITE(76,*) 'head:'
          WRITE(76,'('//sijk//'G25.17)') head
          WRITE(76,*) 'temp:'
          WRITE(76,'('//sijk//'G25.17)') temp
          WRITE(76,*) 'pres:'
          WRITE(76,'('//sijk//'G25.17)') pres
          CLOSE(76)
        END IF

        IF (test_option('run_comp')) THEN
          OPEN(76,file=filename)
          IF (linfos(3)>=1) THEN
            WRITE(*,'(3A)') '  [R] : Comparison to "', &
              filename(1:lblank(filename)), '"'
          END IF
          READ(76,'(1A)') sfield
          READ(76,'('//sijk//'G25.17)') (headold(i,cgen_fw,ismpl), &
            i=1,i0*j0*k0)
          READ(76,'(1A)') sfield
          READ(76,'('//sijk//'G25.17)') (tempold(i,cgen_fw,ismpl), &
            i=1,i0*j0*k0)
          READ(76,'(1A)') sfield
          READ(76,'('//sijk//'G25.17)') (presold(i,cgen_fw,ismpl), &
            i=1,i0*j0*k0)
          CALL s_damax(i0*j0*k0,head(1,1,1,ismpl),mhead)
          CALL daxpy(i0*j0*k0,-1.0D0,headold(1,cgen_fw,ismpl),1, &
            head(1,1,1,ismpl),1)
          shead = 0.0D0
          CALL s_ddot(i0*j0*k0,head(1,1,1,ismpl),head(1,1,1,ismpl), &
            shead)
          CALL s_damax(i0*j0*k0,temp(1,1,1,ismpl),mtemp)
          CALL daxpy(i0*j0*k0,-1.0D0,tempold(1,cgen_fw,ismpl),1, &
            temp(1,1,1,ismpl),1)
          stemp = 0.0D0
          CALL s_ddot(i0*j0*k0,temp(1,1,1,ismpl),temp(1,1,1,ismpl), &
            stemp)
          CALL s_damax(i0*j0*k0,pres(1,1,1,ismpl),mpres)
          CALL daxpy(i0*j0*k0,-1.0D0,presold(1,cgen_fw,ismpl),1, &
            pres(1,1,1,ismpl),1)
          spres = 0.0D0
          CALL s_ddot(i0*j0*k0,pres(1,1,1,ismpl),pres(1,1,1,ismpl), &
            spres)
          WRITE(*,'(6X,1A,1e24.16)') '||rel.Head||2=', &
            dsqrt(shead)/mhead
          WRITE(*,'(6X,1A,1e24.16)') '||rel.Temp||2=', &
            dsqrt(stemp)/mtemp
          WRITE(*,'(6X,1A,1e24.16)') '||rel.Pres||2=', &
            dsqrt(spres)/mpres
          CLOSE(76)
        END IF

        RETURN
      END
