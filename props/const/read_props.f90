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

!> @brief read user defined additionally parameter
!> @param[in] ismpl local sample index
      subroutine read_props(ismpl)
        use arrays
        use mod_flow, only: rref
        use mod_genrl
        use mod_genrlc, only: project
        use mod_const, only: npropsf, fprops, pconst_rhof, &
            pconst_compf, pconst_cpf, pconst_lamf, pconst_visf
        use mod_linfos, only: linfos

        implicit none

        ! Sample index
        integer :: ismpl

        ! Filename of input file
        character (len=80) :: filename

        ! String of curretly read line from input file
        character (len=80) :: line

        ! Index
        integer :: i

        ! String utilities
        integer, external :: lblank
        logical, external ::  found

        ! External input utility
        logical, external :: no_ext_link


        ! Set filename
        filename = project(:lblank(project))
        write(*,*) ' '
        write(*,*) '  reading constant fluid properties:'
        write(*,*) '    from file "', filename, '"'
        write(*,*) ' '

        ! open project config file
        open(79,file=filename,status='old')

        ! input from keyword '# fluid props'
        if (found(79,key_char//' fluid props',line,.false.)) then
          read(79,*) (fprops(i),i=1,npropsf)

          if (linfos(3)>=1) then
            write(*,*) ' [R] : constant fluid properties'
            write(*,*) ' '
            write(*,'(a,/a,a)') '     fluid properties: ', &
              '        rhof       compf         cpf        lamf', &
              '        visf'
            write(*,'(6e12.4)') (fprops(i),i=1,npropsf)
            write(*,*) ' '
          end if
          
        else                    

          ! default fluid density
          fprops(pconst_rhof) = 998.0d0
          ! default fluid compressibility
          fprops(pconst_compf) = 5.0d-8
          ! default fluid volumetric heat capacity (rho*c)
          fprops(pconst_cpf) = 4218.0d0
          ! default fluid thermal conductivity
          fprops(pconst_lamf) = 0.65d0
          ! default fluid viscosity
          fprops(pconst_visf) = 1.0d-3

          write(*,*) ' [D] : constant fluid properties assumed'
          write(*,*) ' '

          if (linfos(3)>=1) then
            write(*,'(a,/a,a)') '     fluid properties: ', &
              '        rhof       compf         cpf        lamf', &
              '        visf'
            write(*,'(6e12.4)') (fprops(i),i=1,npropsf)
            write(*,*) ' '
          end if

        end if

        ! Overwriting rref with the constant density read in this file
        write(*,*) ' [?] : WARNING: overwrite "rref" with constant "rhof"'
        rref = fprops(pconst_rhof)


        ! close project config file
        close(79)

        return

      end subroutine read_props
