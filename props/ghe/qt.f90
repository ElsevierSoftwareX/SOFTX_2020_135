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

!>    @brief assign heat production to cell
!>    @param[in] i grid indices
!>    @param[in] j grid indices
!>    @param[in] k grid indices
!>    @param[in] ismpl local sample index
!>    @return heat production
      DOUBLE PRECISION FUNCTION qt(i,j,k,ismpl)
        USE work_array_ghe
	USE arrays
 !       USE mod_genrl
 !       USE mod_genrlc
!        IMPLICIT NONE
        INTEGER i, j, k, ismpl
 !       DOUBLE PRECISION por

! vr : aCtor (1-por) ?
!      por=propunit(uindex(i,j,k),idx_por,ismpl)
!      qt = (1.d0- por)*(propunit(uindex(i,j,k),idx_q,ismpl)
        qt = propunit(uindex(i,j,k),idx_q,ismpl)

        IF (simtime(ismpl).le.Tin(mp,1)) THEN
          mp=mp
	  
        ELSE 
          mp=mp+1 	
        END IF			   
!        WRITE(*,*),Tin(m,1)

!       Anfang Schleife über Sonden
        DO nl=1,nghe
          IF (Tin(mp,2).ne.0) THEN !  An/aus
!            DO kl=k_start(nl),k_end(nl)
              IF (i.eq.(ighe(nl)).and.(j.eq.jghe(nl)).and.((k.ge.k_start(nl)) &
	         .and.(k.le.k_end(nl))))  THEN
                qt = -Tin(mp,3)*1000/&
                  (delx(ighe(nl))*dely(jghe(nl)))
!	      ELSE
!                qt =propunit(uindex(i,j,k),idx_q,ismpl)	
              END IF
 !           END DO
	  ELSE
	    qt=propunit(uindex(i,j,k),idx_q,ismpl)	
          END IF
        END DO
!       Ende Schleife über Sonden

!        WRITE(*,*), 'Periode: (iper)' , mp, 'Zeit', simtime(ismpl)
		
        RETURN
      END
