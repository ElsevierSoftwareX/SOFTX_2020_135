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
      subroutine create_condtab()
c-----------------------------------------------------------------------
c     
c     Builds a lookup table forthe local shape of the 
c     conditional probaility.
c     See Oz et. al 2003 or Deutch 2000, for details.
c     *********************************************
c     
c     INPUT VARIABLES:
c     
c     OUTPUT VARIABLES:
c     
c     condtab : conditoinal prob llokup table 
c     
c     ORIGINAL : Thomas Mejer Hansen                       DATE: August 2005
c
c
c-----------------------------------------------------------------------
      use simul_arrays
        use mod_simul
        use mod_OMP_TOOLS
      implicit none
      include 'visim.inc'
      include 'OMP_TOOLS.inc'
      integer :: i,j,k,im,iv,ierr,ierror,i_monte
      integer :: llout,llout1,llout2,llout3,llout4
      double precision Gmean,gGvar
      double precision sum_sim,sum_sim2
      double precision mean_sim, var_sim
      double precision p,zt,simu
      double precision arr(20)
      double precision target(nbt), target_nscore(nbt)
      double precision target_p(nbt), temp(nbt)
      double precision target_weight(nbt)
      integer :: itarget_quan
c      double precision zmin,zmax
      double precision q_norm, q_back(599) 
      double precision x_cpdf(500) 
      double precision backtrans
      integer :: index_cdf
      character (len=80) :: tmpfl
      double precision te
      double precision dummy1,dummy2

      integer :: GmeanType,GvarType
c
      double precision backtr
      external backtr
c

      GmeanType=0
      GvarType=0

      do i=1,nbt
         target(i)=bootvar(i);
c     NEXT LINE TO MAKE SURE ALLE DATA HAVE WEIGHT 1... 
c     THERE IS PROBABLY A BUG IN THE visim_readpar_uncertainty.f file here
c     AS bootvar is NOT 1 when it has to be..
         target_weight(i)=1;
c         write(*,*) i,target(i),target_weight(i)
      enddo


c
c     Compute Normal Score of TARGET HISTOGRAM
c
      if (idbg.gt.0) then 
         write(*,*) 'zmin=',zmin
         write(*,*) 'zmax=',zmax
         write(*,*) 'ltail=',ltail,ltpar
         write(*,*) 'utail=',utail,utpar
      endif 

      call nscore(nbt,target,zmin,zmax,0,target_weight,temp,-1,
     &     target_nscore,ierror)
CAW-      call nscore(nbt,target,zmin,zmax,0,target_weight,temp,1,
CAW-     &     target_nscore,ierror)

c     WRITING NSCORE TABLE TO FILE
      write(tmpfl,771) 'nscore',outfl
      call omp_new_file_handler(llout,3)
      open(llout, file=tmpfl, status ='unknown', position='APPEND')
      do i=1,nbt
         write(llout,*) target(i),target_nscore(i),target_weight(i),
     &        zmin,zmax
      enddo
      close(llout)

      do i=1,n_q
         x_quan(i)=(1.d0/n_q)/2+(i-1)*(1.d0/n_q)
      enddo


      if (idbg.gt.0) then 
         write(*,*) ' Nscore MEAN range = ',min_Gmean, max_Gmean,n_Gmean
         write(*,*) ' Nscore VAR range = ',min_Gvar, max_Gvar, n_Gvar
         write(*,*) ' Number of quantiles = ',n_q
         write(*,*)' Number of samples drawn in nscore space = ',n_monte
      endif

      write(*,*) 'Calc CondPDF Lookup n_Gmean,n_Gvar=',n_Gmean,n_Gvar 
      do im=1,n_Gmean

         if (GmeanType.eq.1) then 
c     Focus on middle range mean
            if (im.lt.(n_Gmean/2)) then
               Gmean = min_Gmean + 0.5d0*(max_Gmean-min_Gmean)*
     1              (1- exp(-1.0d0*im/(n_Gmean/20)) )
            else
               Gmean = min_Gmean+ 0.5d0*(max_Gmean-min_Gmean) +
     1              0.5d0*(max_Gmean-min_Gmean)*
     1              (1-(1-exp(-1.0d0*(n_Gmean-im)/(n_Gmean/20)) ))
            endif
         else
c     linear mean range
            Gmean=min_Gmean+(im-1)*(max_Gmean-min_Gmean)/(n_Gmean-1)
         endif

         if (idbg.gt.2) write(*,*) 'precalc lookup im,n_Gmean=',
     1        im,n_Gmean

         do iv=1,n_Gvar

            if (GvarType.eq.1) then 
c Cosine, focus on low variances
               gGvar = min_Gvar + 
     1              (max_Gvar-min_Gvar)*(1-cos(.5d0*iv*3.14d0/n_Gvar))
            elseif (GvarType.eq.2) then
c Cosine, focus on high variances
               gGvar = min_Gvar + 
     1              (max_Gvar-min_Gvar)*(cos(.5d0*iv*3.14d0/n_Gvar))
            else
c     Linear 
               gGvar=min_Gvar+(iv-1)*(max_Gvar-min_Gvar)/(n_Gvar-1)
            endif
            if (iv.eq.1) gGvar=min_Gvar
c     BACK TRANSFORM QUANTILES
            dummy1=0
            dummy2=0
            do i=1,n_q
c               x_quan(i)=(1./n_q)/2+(i-1)*(1./n_q)
               call gauinv(dble(x_quan(i)) ,zt,ierr)
               q_norm=zt*sqrt(gGvar)+Gmean            

               x_cpdf(i) = backtr(q_norm,nbt,target,target_nscore,
     +              zmin,zmax,ltail,ltpar,utail,utpar)

               dummy1 = dummy1 + x_cpdf(i)
               dummy2 = dummy2 + x_cpdf(i)*x_cpdf(i)

            enddo

            dummy1 = dummy1 / n_q
            dummy2 = dummy2 / n_q
            dummy2 = dummy2 - dummy1*dummy1


C
C     MONTE CARLO SAMPLING OF GAUSSIAN PDF
C             sum_sim=0
C             sum_sim2=0
C             do imonte=1,n_monte
C                p = acorni2(idum)         
C                do i=1,10
C                   p = acorni2(idum) 
C                   if ((p.gt.x_quan(1)).AND.(p.lt.x_quan(n_q))) then
C                      exit
C                   else
C                      if (i.gt.20) then
C                         write(*,*) 'QUANTILE OUTSIDE RANGE',
C      1                       i,p,x_quan(1),x_quan(n_q)
C                      endif
C                   endif
C                enddo
C c     Just in case quantile still out of range
C                if (p .le. x_quan(1)) then
C                   index_cdf = 1
C                else if (p .gt. x_quan(n_q)) then
C                   index_cdf = n_q
C c                  write(*,*) 'x_cpdf(index_cdf)',x_cpdf(index_cdf),n_q
C                endif
C                call locate(x_quan,n_q,1,n_q,p,index_cdf)
C                sum_sim =  sum_sim + x_cpdf(index_cdf)
C                sum_sim2 = sum_sim2 + x_cpdf(index_cdf)*x_cpdf(index_cdf)
C             enddo               
C             sum_sim=sum_sim/n_monte
C             sum_sim2=sum_sim2/n_monte
C             mean_sim = sum_sim
C             var_sim = sum_sim2 - sum_sim*sum_sim

            mean_sim=dummy1
            var_sim=dummy2
            
            if (var_sim.lt.0) var_sim=0

            condlookup_mean(im,iv)=mean_sim
            condlookup_var(im,iv)=var_sim
            do i=1,n_q
               condlookup_cpdf( im,iv,i) = x_cpdf(i)
            enddo
            
         enddo

         if (idbg.gt.1) write(*,*) 'Gmean, gGvar, SimMean, SimVar= ',
     1        Gmean,gGvar,mean_sim
         
      enddo

c     wirte lookup tables to disk      
      if (idbg.gt.0) then 
         call omp_new_file_handler(llout1,4)
         call omp_new_file_handler(llout2,5)
         call omp_new_file_handler(llout3,6)
         call omp_new_file_handler(llout4,7)
         write(tmpfl,771) 'cond_imean',outfl
         open(llout1, file=tmpfl, status = 'unknown')
         write(tmpfl,771) 'cond_mean',outfl
         open(llout2, file=tmpfl, status = 'unknown')
         write(tmpfl,771) 'cond_var',outfl
         open(llout3, file=tmpfl, status = 'unknown')
         if (idbg.ge.6) then
            write(tmpfl,771) 'cond_cpdf',outfl
            open(llout4, file=tmpfl,status='unknown',position='APPEND')
         endif
 771     format(A,'_',A)
         do im=1,n_Gmean
            do iv=1,n_Gvar
               write(llout1,*) im
               write(llout2,*) condlookup_mean(im,iv)
               write(llout3,*) condlookup_var(im,iv)
               if (idbg.gt.2) then
                  do i=1,n_q
                     write(llout4,*) condlookup_cpdf(im,iv,i)
                  enddo
               endif
            enddo     
         enddo
         close(llout1)
         close(llout2)
         close(llout3)
         close(llout4)
      endif

      return

      end
      


      double precision function drawfrom_condtab(cmean,cvar,p)
c-----------------------------------------------------------------------
c     
c     Draw from a lookup table for the local shape of the 
c     conditional probaility.
c     See Oz et. al 2003 or Deutch 2000, for details.
c     *********************************************
c     
c     INPUT VARIABLES:
c     
c     OUTPUT VARIABLES:
c     
c     condtab : conditoinal prob llokup table 
c     
c     ORIGINAL : Thomas Mejer Hansen                       DATE: August 2005
c     
c
c-----------------------------------------------------------------------
      use simul_arrays
        use mod_simul
      implicit none
      include 'visim.inc'
      double precision cmean, cvar
      integer :: im,iv,iq,ie,is,xid
      double precision cmean_arr(500), i_arr(500)
      double precision i_mean(500), mean(500)
      double precision dist, dm, dv
      double precision mindist
      integer :: im_sel, iv_sel, i, idum
      double precision m_sel, v_sel
      double precision p, acorni2
      external acorni2
      integer :: index_cdf
      double precision Kmean, Kstd, Fmean, Fstd, draw
      character (len=80) :: tmpfl
c
c
      cvar=cvar*cvar

c     NEXT TMH
CAW      dm=xmax-xmin;
c     NEXT OZ
      dm=skgmean
      dv=gvar

c      write(*,*) 'gmean=',gmean
      
      mindist=1d+9
      do im=1,n_Gmean
         do iv=1,n_Gvar


C     TMH STYLE
c     BUT TOO HIGH SILL VALUE
            dist=( (condlookup_mean(im,iv)-cmean)/dm )**2+
     +        ( (condlookup_var(im,iv)-cvar)/dv )**2            
C     OZ STYLE            
C     WORSE MATCH TO HISTOGRAM THAN ABOVE
c            dist=( (condlookup_mean(im,iv)-cmean)/dm )**2+
c     +        abs (condlookup_var(im,iv)-cvar)/dv             

            
            if (dist.lt.mindist) then
               mindist=dist
               im_sel=im
               iv_sel=iv
            endif
         enddo
      enddo

      m_sel = condlookup_mean(im_sel,iv_sel)
      v_sel = condlookup_var(im_sel,iv_sel)

c
      write(tmpfl,871) 'kriging',outfl
      open(lout_krig, file=tmpfl, position='APPEND', status='unknown')
      write(lout_krig,86) cmean, cvar,m_sel,v_sel
      close(lout_krig)
 86   format(f12.6,f15.9,f12.6,f15.9)  
 871  format(A,'_',A)
c

c     NOW DRAW FROM LOCAL CPDF

c     MAKE SURE THAT QUANTILE IS WITHIN BOUNDS
      do i=1,10
c     NEXT LINE COMMENTED OUT ONCE - WHY 
         p = acorni2(idum) 
         if ((p.gt.x_quan(1)).AND.(p.lt.x_quan(n_q))) then
            exit
         else
            if (i.gt.1) then
           write(*,*) 'QUANTILE OUTSIDE RANGE',i,p,x_quan(1),x_quan(n_q)
            endif
         endif
      enddo

      if (p .le. x_quan(1)) then
         index_cdf = 1
         write(*,*)'bad low quantile',p,x_quan(1)
      else if (p .gt. x_quan(n_q)) then
         index_cdf = n_q
         write(*,*)'bad high quantile',p,x_quan(n_q)
      else
         call locate(x_quan,n_q,1,n_q,p,index_cdf)
      endif

c      index_cdf = 1 + int(p*n_q)
      draw = condlookup_cpdf(im_sel,iv_sel,index_cdf) 
c      write(*,*) 'INDEX ME',index_cdf


c     CORRECTION ACCORDING TO Oz et al, 2003

      Fmean = condlookup_mean(im_sel,iv_sel)
      Fstd = sqrt( condlookup_var(im_sel,iv_sel) )

      Kmean= cmean
      Kstd = sqrt(cvar)
      
      drawfrom_condtab = ( draw - Fmean ) * ( Kstd / Fstd) + Kmean

c     NEXT LINE TO NOT USE CORRECTION
c     SOMETIMES THIS CORRECTION MAKES MORE BAD THAN GOOD
      drawfrom_condtab = draw


      if (idbg.gt.14) then
         write(*,*) 'INDEX CDF = ',index_cdf
         write(*,*) 'cmean,cvar 2->',cmean,cvar
         write(*,*) 'im_sel,iv_sel -->',im_sel,iv_sel
         write(*,*) 'm_sel,v_sel -->',m_sel,v_sel
         write(*,*) 'Kmean,Kstd -->',Kmean,Kstd
         write(*,*) 'Fmean,Fstd -->',Fmean,Fstd
         write(*,*) 'Fmean,Fstd -->',drawfrom_condtab,draw
      endif



      return 


      end
      
