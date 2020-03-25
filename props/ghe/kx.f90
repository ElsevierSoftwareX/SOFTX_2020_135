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

!> @brief assign permeability in x direction to cell
!> @param[in] i grid indices
!> @param[in] j grid indices
!> @param[in] k grid indices
!> @param[in] ismpl local sample index
!> @return  permeability                        (m^2)
!> @details
!> kx returns the permeability in x-direction [m2] at node(i,j,k) from
!> the input file.\n\n
!>
!> The permeability in x-direction is the product of the permeability
!> in z-direction and the anisotropy factor for the x-direction.
      double precision function kx(i,j,k,ismpl)
        use arrays, only: propunit, uindex, idx_kz, idx_an_kx

        implicit none

        ! Location indices
        integer, intent (in) :: i
        integer, intent (in) :: j
        integer, intent (in) :: k

        ! Sample index
        integer :: ismpl


        kx = propunit(uindex(i,j,k),idx_kz,ismpl)* &
            propunit(uindex(i,j,k),idx_an_kx,ismpl)

!       Wird das Setzen der Werte benoetigt?

!!       ANFANG SCHLEIFE ÜBER SONDEN #######################
!        DO n=1,nghe
!!         hier wird variable Tiefe der Sonden gesetzt
!          k_start(n)=K0-(k_end(n)-depth_hpr/delz(1))
!          k_end(n)=K0-k_end(n)
!
!          WRITE(*,*), 'kdepth', k_start(n), k_end(n)
!          DO l=k_end(n),k_start(n)
!!           setze permeabilitaet an den sonden kl
!            IF (ighe(n).eq.i && jghe(n).eq.j && (K0-l+1).eq.k)
!              kx=1e-25
!            END IF
!          END DO
!        END DO
!!       ENDE SCHLEIFE ÜBER SONDEN. #########################

        return

      end function kx
