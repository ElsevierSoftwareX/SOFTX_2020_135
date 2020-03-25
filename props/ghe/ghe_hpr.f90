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

!************************************************************************
!*****                                                              *****
!*****                                                              *****
!*****                                                              *****
!************************************************************************
!*****                                                              *****
!*****  THIS PROGRAM APPRXIMATES THE EFFECT OF GHE ON THE SUBSURFACE*****
!*****                                                              *****
!*****  BY USING AN EFFECTIVE HEAT PRODUCTION                       *****
!*****                                                              *****
!*****                                                              *****
!*****                dm 270808 modified by jh and rs and dm 250112 *****
!************************************************************************
!*
!**********************************************************************
      SUBROUTINE GHE_HPROD
        USE arrays
        USE work_array_ghe
	USE mod_genrl
        USE mod_genrlc
!        IMPLICIT NONE

        DOUBLE PRECISION depth_hpr

!       Alle Sonden gleich...zunächst
        ALLOCATE(dTs(K0))

!       READ all geometrical borehole parameters************************
        WRITE(*,*) "geometry"

        OPEN(1,FILE='ghe_new.par')

!       Anzahl der GHE
        READ(1,*) nghe

        READ(1,*) depth_hpr                !

!       i Stelle GHEs
        ALLOCATE(ighe(nghe))

!       j Stelle GHEs
        ALLOCATE(jghe(nghe))

!       Sondenkopf an der Erdoberfläche im Modell (k von Shemat)
        ALLOCATE(k_end(nghe))

!       Sondenfuß (k von Shemat)
        ALLOCATE(k_start(nghe))

!       Ort(e) der Sonden
        READ(1,*) (ighe(i), jghe(i),k_end(i),i=1,nghe)
        CLOSE(1)

        WRITE(*,*)  nghe
        DO nl=1,nghe
          WRITE(*,*) ighe(nl),jghe(nl)
        END DO

        OPEN(2,FILE='ghe_ini_new.ini')

!       Anzahl der Perioden
        READ(2,*) iper
        WRITE(*,*)
        WRITE(*,*) '  READING NEW Ground Heat Exchanger model parameters'
        WRITE(*,*)
        WRITE(*,*)  'Anzahl der Perioden', iper

        ALLOCATE(Tin(iper,4))

!       Periode (h), Flow an/aus, Leistung (NEU: Leistung (kW) im Gesamtfeld!)
        READ(2,*) (Tin(i,1),Tin(i,2),Tin(i,3),i=1,iper)
        CLOSE(2)

!       Zeiten in Sekunden
        Tin(1,1)=Tin(1,1)*3600.0d0

        DO i=2,iper
           Tin(i,1)=Tin(i,1)*3600.0d0+Tin(i-1,1)
        END DO
!       auf Sonde und Tiefe verteilt...spaeter Hprod
        DO i=1,iper
          Tin(i,3)=Tin(i,3)/nghe/depth_hpr
        END DO
!----------------% ende einlesen ini------

! Hier muss die Anzahl der Sonden berücksichtigt werden!

!       ANFANG SCHLEIFE ÜBER SONDEN #######################
        DO nl=1,nghe
!         hier wird variable Tiefe der Sonden gesetzt
         k_start(nl)=k_end(nl)-depth_hpr/delz(1)+1
  !        k_end(nl)=K0-k_end(nl)

          WRITE(*,*), 'kdepth', k_start(nl), k_end(nl)
!          DO k=k_end(n),k_start(n)
!!           setze permeabilitaet an den sonden kl
!            kx(ighe(n),jghe(n),K0-k+1,ismpl)=1e-25
!            ky(ighe(n),jghe(n),K0-k+1,ismpl)=1e-25
!            kz(ighe(n),jghe(n),K0-k+1,ismpl)=1e-25
!          END DO
!
        END DO
!!       ENDE SCHLEIFE ÜBER SONDEN. #########################

!       Für nächste Routine.
        mp=1
      END SUBROUTINE


!*******ONE TIME STEP***********************************************************

!      SUBROUTINE GHETIMESTEP_HPR()
!        USE work_array_ghe
!        implicit none

!        IF (sdelt.le.Tin(m,1)) THEN
!          m=m
!        ELSE
!          m=m+1
!        END IF
!        WRITE(*,*),Tin(m,1)

!!       Anfang Schleife über Sonden
!        DO n=1,nghe
!          IF (Tin(m,2).ne.0) THEN !  An/aus
!            DO k=k_end(n),k_start(n)-1
!              qt(ighe(n),jghe(n),K0-k+1)=-Tin(m,3)*1000/&
!                (delx(ighe(n))*dely(jghe(n)))
!            END DO
!	      ELSE
!            DO k=k_end(n),k_start(n)-1
!              qt(ighe(n),jghe(n),K0-k+1)=0
!            END DO
!          END IF
!        END DO
!!       Ende Schleife über Sonden

!        WRITE(*,*), 'Periode: (iper)' , m, 'Zeit', sdelt
!      END SUBROUTINE
