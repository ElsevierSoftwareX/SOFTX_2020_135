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


!>    @brief read user defined additionally parameter
!>    @param[in] ismpl local sample index
      SUBROUTINE read_props(ismpl)


        use arrays
                use ice
        use mod_genrl
        use mod_genrlc
        IMPLICIT NONE
        integer :: ismpl
        integer :: i, j, k

        character (len=80) :: line

        INTEGER lblank
        LOGICAL found, no_ext_link
        EXTERNAL found, no_ext_link, lblank



! PLEASE comment out the following line if you want to use the additional reading
        RETURN

!     std. screen information
        WRITE(*,*)
        WRITE(*,*) '  reading additionally user parameter:'
        WRITE(*,*) '    from file "', project(:lblank(project)), '"'
!     open model file
        OPEN(79,file=project,status='old')
!     init HDF5 support, when available
        CALL open_hdf5(' ')

! --------- begin reading part ---------

!     !!! [liq] reading example for an existing [I0,J0,K0] array

!     searching for the right data part, special keyword '# liq init'
        IF (found(79,'# liq init',line,.FALSE.)) THEN
!        check for an HDF5 entry 'liq', if available read the array [liq] from the HDF5 file
          IF (no_ext_link(i0,j0,k0,liq,'liq',line)) THEN
!           read the values from the model file instead of the HDF5 file
            READ(79,*) (((liq(i,j,k),i=1,i0),j=1,j0),k=1,k0)
          END IF
!        std. screen information
          WRITE(*,*) ' [R] : liq'
        ELSE
!        std. screen information, when '# liq init' not found
          WRITE(*,*) ' <D> : liq = -0.04'
!        make a default initialisation of the array [liq] with '-0.04d0'
          CALL set_dval(i0*j0*k0,-0.04D0,liq)
        END IF

!     !!! add here additionally array readings ...

! --------- begin reading part ---------

!     finish HDF5 support
        CALL close_hdf5()
!     close model file
        CLOSE(79)
        RETURN
      END
