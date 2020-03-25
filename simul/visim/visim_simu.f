!>Copyright (C) 1996, The Board of Trustees of the Leland Stanford     %\n
!>Junior University.  All rights reserved.                             %\n
!>
!>The programs in GSLIB are distributed in the hope that they will be  %\n
!>useful, but WITHOUT ANY WARRANTY.  No author or distributor accepts  %\n
!>responsibility to anyone for the consequences of using them or for   %\n
!>whether they serve any particular purpose or work at all, unless he  %\n
!>says so in writing.  Everyone is granted permission to copy, modify  %\n
!>and redistribute the programs in GSLIB, but only under the condition %\n
!>that this notice and the above copyright notice remain intact.       %\n
      double precision function  simu(cmean1, cstdev1)
      
c     
c     This function draws from the local conditional distribution and return
c     the value simu. The drawing depends on the type of local distribution 
c     specified in idrawopt
c     
      use simul_arrays
        use mod_simul
      implicit none
      include 'visim.inc'
      double precision acorni, drawfrom_condtab
      external acorni, drawfrom_condtab
      double precision p  
      double precision cmean1, cstdev1
      double precision aunif, bunif   
      double precision cvarn, cmn, zt, cvar
      integer :: ierr, idum
c
c
      p = acorni(idum)
      simu = 0.0d0
c
      if(idrawopt.eq.0) then   
         call gauinv(dble(p),zt,ierr)
         simu=zt*cstdev1+cmean1            
      else if(idrawopt.eq.1) then   
c     USE DSSIM HR CODE
c         simu = drawfrom_condtab(cmean1,cstdev1,p)
         simu = drawfrom_condtab(cmean1,cstdev1,p)
      else	
         write(*,*) 'Error: drawing option larger than 1'     
         write(*,*) 'No implementation for this option'
      endif
c
      return
      end
