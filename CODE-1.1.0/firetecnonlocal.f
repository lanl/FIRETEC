      subroutine firetecnonlocal(it)
      use metryic
      use pres
      use gridsetup
      use fireteca
      use constants
      use turba
      use turbb
      use xvo
      use xve
      use workavg
      use radiation
      use msga
      use ignite
      use nonlocal

      Implicit None


      !JAS 3/7/06 added explicit declarations to comply with implicit none
      integer,intent(in) :: it
      integer :: i,j,k,its,ia,ja  !KOO  icallrad moved to radiation.mod
      real,external :: zcart
      real :: perchydroremaining,ztopcell !,hydrofactor
      real :: urhoa,vrhoa,speedrho2,zla
      real :: rhosold ! ,rhosnph
      real :: siesold,tempsold,rignitionperiod,travelrate,pretime
      real :: finterp,x1,y1,x2,y2,dist,speedfuelrho2


c      call rmaxmin1(temps,'temps',1-ih,np+ih,1-ih,mp+ih,l)
c      call rmaxmin1(tempg,'tempg',1-ih,np+ih,1-ih,mp+ih,l)

      frho=0.
      fhc=0.
      fhcb=0.
      fox=0.
      foxb=0.
      frhof=0.
      frhosies=0.
      frhowater=0.
      if(irhovapor.eq.1)
     .frhovapor=0.
      if(irhovapor.eq.1)
     .frhovaporb=0.

         do k=1,l
            do j=1,mp
               do i=1,np
                  if (rhofinitial(i,j,k).gt.0.0) then
                  perchydroremaining=
     &       max(0.,(rhof(i,j,k)-rhohydrothresh*rhofinitial(i,j,k))
     .                 /(rhofinitial(i,j,k)*(1.-rhohydrothresh)))
C                  cf(i,j,k)=cfhydro*perchydroremaining
C     .                     +cfchar*(1-perchydroremaining)    
                    cf(i,j,k)=0.001
                  else
                    cf(i,j,k)=0.0
                  endif

c               hydrofactor=exp(-rno/rnfuel*psif(i,j,k)
c     &                        *rhof(i,j,k)/xvb(i,j,k,8))
c               thetasolid(i,j,k)=perchydroremaining*.25*hydrofactor
c     &                          +(1-perchydroremaining)*.9
               thetasolid(i,j,k)=perchydroremaining*.25
     &                          +(1-perchydroremaining)*.9


c                 if (actualfueldepth(i,j,k).lt.0.or.k.gt.1) then
                 if (actualfueldepth(i,j,k).lt.0.or.k.lt.1) then
                    xvfuel(i,j,k,1)=xvb(i,j,k,1)
                    xvfuel(i,j,k,2)=xvb(i,j,k,2)
                    xvfuel(i,j,k,3)=xvb(i,j,k,3)
                    xvfuel(i,j,k,4)=xvb(i,j,k,5)
                    xvfuel(i,j,k,5)=xvb(i,j,k,6)
 
                 elseif (k.lt.1) then
                    ztopcell=zcart(zedge(k+1),i,j)-zs(i,j)
                    urhoa=xvb(i,j,k+1,1)
                    xvfuel(i,j,k,1)=xvb(i,j,k,1)*
     &                        sqrt(actualfueldepth(i,j,k)/ztopcell)
                    vrhoa=xvb(i,j,k+1,2)
                    if (abs(xvb(i,j,k+1,2)).le.abs(xvb(i,j,k,2)).or.
     &                 xvb(i,j,k,2)/xvb(i,j,k+1,2).le.0.)
     &                 vrhoa=xvb(i,j,k,2)
                    xvfuel(i,j,k,2)=xvb(i,j,k,2)*
     &                        sqrt(actualfueldepth(i,j,k)/ztopcell)

                    xvfuel(i,j,k,3)=xvb(i,j,k,3)*
     &                        sqrt(actualfueldepth(i,j,k)/ztopcell)
                    speedrho2=(xvb(i,j,k,1)**2+
     &                             xvb(i,j,k,2)**2+
     &                             xvb(i,j,k,3)**2)
                    speedfuelrho2=(xvfuel(i,j,k,1)**2+
     &                             xvfuel(i,j,k,2)**2+
     &                             xvfuel(i,j,k,3)**2)

                    xvfuel(i,j,k,4)=xvb(i,j,k,5)*speedfuelrho2/speedrho2
                    xvfuel(i,j,k,5)=xvb(i,j,k,6)*speedfuelrho2/speedrho2
                 endif
 
               enddo
            enddo
         enddo
c      call rmaxmin1(cf,'cf in firetec',1-ih,np+ih,1-ih,mp+ih,l)
c      call rmaxmin1(rhofinitial,'rhofinitial
c     . in firetec',1-ih,np+ih,1-ih,mp+ih,l)
c234567
      if (ilapdo.eq.1) then

         call lapdo(xvb,xvb(1-ih,1-ih,1,7),foxb,
     .   d11b,d22b,d33b,d12b,d13b,d23b,1-ih,np+ih,1-ih,mp+ih,l,nv)
         call lapdo(xvb,xvb(1-ih,1-ih,1,8),fhcb,
     .   d11b,d22b,d33b,d12b,d13b,d23b,1-ih,np+ih,1-ih,mp+ih,l,nv)     

         do k=1,l
            do j=1,mp
               do i=1,np
                  foxb(i,j,k)=foxb(i,j,k)*dti
                  fhcb(i,j,k)=fhcb(i,j,k)*dti
               enddo
            enddo
         enddo

      endif
         do k=1,l
            do j=1,mp
               do i=1,np
                  pr(i,j,k)=(xvb(i,j,k,4)*rg/prrcp)**(cp/cv)
                  tempg(i,j,k)=xvb(i,j,k,4)/xvb(i,j,k,nv)
     +                 *(pr(i,j,k)*1.e-5)**(rg/cp)
               enddo
            enddo
         enddo
c      if (temps(i,j,k).gt.2000) then
c         write (*,*) 'endo of firetec i=',i,' j=',j,' k=',k
c         write(*,*) 'temps(i,j,k)=',temps(i,j,k)
ccc         write(*,*) 'rhof(i,j,k)=',rhof(i,j,k)
c         write(*,*) 'rhowater(i,j,k)=',rhowater(i,j,k)
c       endif

c       icallrad=10
c  JLWrad10 change    !MMC 5/1/07 ignition/restart fix by JAS
c      if (irad.eq.1.and.((it.eq.1.and.irst.eq.0).or.mod(it,10).eq.0)) then
       if ((it.eq.1.and.irst.eq.0).or.(mod(it,icallrad).eq.0)) then !KOO

      firad=0.
      frhosiesrad=0.
c         call fire_radiation(xvb(1-ih,1-ih,1,nv),xvb(1-ih,1-ih,1,4),
         !call fire_radiation(xvb(1-ih,1-ih,1,nv),tempg,
         !+           firad,xvb(1-ih,1-ih,1,7),rhof,cv,temps,frhosiesrad,
         !+                c13,c23,gi,gmul,x,y,z,zs,
         !+                E,Ef,Es,rnetgas,rnetsol,tambientarray)
         !JAS 3/2/06 changed signature of fire_radiation as cv is not used.
c  KOO begin     
         if (irad.EQ.1) call fire_radiation()
         !if (irad.EQ.1) call fire_radiation(xvb(1-ih,1-ih,1,nv),tempg,
      !+        firad,xvb(1-ih,1-ih,1,7),rhof,sizescale,temps,frhosiesrad,
      !+                Ef,Es,rnetgas,rnetsol,tambientarray)
      !+        tambientarray)
         if (irad.EQ.2) call firerad_MC()
         !if (irad.EQ.2) call firerad_MC(xvb(1-ih,1-ih,1,nv),tempg,
      !+        firad,xvb(1-ih,1-ih,1,7),rhof,sizescale,temps,frhosiesrad,
      !+                 zs,zb,rnetgas,rnetsol,tambientarray,it) ! ma version
c KOO end

      endif   !end if (irad.eq.1.and.(it.eq.1.or.mod(it,10).eq.0))

c      if (temps(i,j,k).gt.2000) then
c         write (*,*) 'after rad  i=',i,' j=',j,' k=',k
c         write(*,*) 'temps(i,j,k)=',temps(i,j,k)
c         write(*,*) 'rhof(i,j,k)=',rhof(i,j,k)
c         write(*,*) 'rhowater(i,j,k)=',rhowater(i,j,k)
c       endif

      do its=1,ntp
      ifuelcount=its        !flag that can be used for debugging standard I/O
         do k=1,l
            do j=1,mp
               do i=1,np
                  pr(i,j,k)=(xvb(i,j,k,4)*rg/prrcp)**(cp/cv)
                  tempg(i,j,k)=xvb(i,j,k,4)/xvb(i,j,k,nv)
     +                 *(pr(i,j,k)*1.e-5)**(rg/cp)
               enddo
            enddo
         enddo

         call fuelnonlocal(xvb,xvfuel,1-ih,np+ih,1-ih,mp+ih,l,nv)

         do k=1,l
            do j=1,mp
               do i=1,np
                  frho(i,j,k)=ff(i,j,k)*dtp
                  frhof(i,j,k)=-ff(i,j,k)*dtp
                  fhc(i,j,k)=ff(i,j,k)*dtp*rnhc-fg(i,j,k)*dtp*rng
     &                       +0.5*fhcb(i,j,k)*dtp                 
                  frhowater(i,j,k)=-fw(i,j,k)*dtp
                  fox(i,j,k)=(-fg(i,j,k)*rnonl+0.5*foxb(i,j,k))*dtp
                enddo
             enddo
         enddo

         call fueltempnonlocal(xvb,xvfuel,1-ih,np+ih,1-ih,mp+ih,l,nv)
         call energysourcenonlocal(xvb,1-ih,np+ih,1-ih,mp+ih,l,nv)

         do k=1,l
            do j=1,mp
               do i=1,np
                  zla=zcart(z(k),i,j)-zs(i,j)
                  xvb(i,j,k,1)=xvb(i,j,k,1)+0.5*f1avg(i,j,k)*dtp
                  xvb(i,j,k,2)=xvb(i,j,k,2)+0.5*f2avg(i,j,k)*dtp
                  xvb(i,j,k,3)=xvb(i,j,k,3)+0.5*f3avg(i,j,k)*dtp
                  xvb(i,j,k,4)=xvb(i,j,k,4)+0.5*fi(i,j,k)
                  xvb(i,j,k,5)=xvb(i,j,k,5)+0.5*fka(i,j,k)*dtp
                  xvb(i,j,k,6)=xvb(i,j,k,6)+0.5*fkb(i,j,k)*dtp
                  xvb(i,j,k,7)=xvb(i,j,k,7)+fox(i,j,k)
                  xvb(i,j,k,8)=xvb(i,j,k,8)+fhc(i,j,k)
                  xvb(i,j,k,nv)=xvb(i,j,k,nv)+frho(i,j,k)

c234567890123456789012345678901234567890123456789023456789012345678
                  qflux(i,j,k)=qflux(i,j,k)+0.5*frhosiesrad(i,j,k)
                  rhof(i,j,k)=rhof(i,j,k)+frhof(i,j,k)
                  rhowater(i,j,k)=rhowater(i,j,k)+frhowater(i,j,k)
c                  if (rhowater(i,j,k).lt.0.0) write (*,*) 'water 
c     &     undershoot',
c     &        i,j,k,rhowater(i,j,k),rhos(i,j,k),
c     &   rhof(i,j,k),psiwmax(i,j,k)
                  if (rhowater(i,j,k).lt.0.0) rhowater(i,j,k)=0.0
                  rhosold=rhos(i,j,k)
                  rhos(i,j,k)=rhof(i,j,k)+rhowater(i,j,k)+rhodirt(i,j,k)
c                  rhosnph=(rhosold+rhos(i,j,k))/2
                  siesold=sies(i,j,k)
                  if(zla.le.fueldepth) then
c                    sies(i,j,k)=
c     &                 ((frhosies(i,j,k)+.5*frhosiesrad(i,j,k)*dtp)
c     &                  +siesold*(rhosnph-.5*(rhos(i,j,k)-rhosold)))/
c     &                   (rhosnph+.5*(rhos(i,j,k)-rhosold))   !rrl  7/9/05
                     sies(i,j,k)=(sies(i,j,k)*rhosold
     .                           +frhosies(i,j,k)
     .                         +0.5*frhosiesrad(i,j,k)*dtp)/rhos(i,j,k)
                     rmoist(i,j,k)=rhowater(i,j,k)/rhof(i,j,k)
                     cpsolid(i,j,k)=(rhodirt(i,j,k)*cpdirt+rhof(i,j,k)*
     &                              (cpwood+cpwater*rmoist(i,j,k)))
     &                              /rhos(i,j,k)
                     tempsold=temps(i,j,k)
                     temps(i,j,k)=sies(i,j,k)/cpsolid(i,j,k)
c      if (temps(i,j,k).gt.1000) then
c         write (*,*) 'after temps calc  i=',i,' j=',j,' k=',k
c         write(*,*) 'temps(i,j,k)=',temps(i,j,k)
c         write(*,*) 'ff(i,j,k)=',ff(i,j,k)
c         write(*,*) 'frhof(i,j,k)=',frhof(i,j,k)
c         write(*,*) 'rhof(i,j,k)=',rhof(i,j,k)
c         write(*,*) 'rhowater(i,j,k)=',rhowater(i,j,k)
c         write(*,*) 'sies(i,j,k)=',sies(i,j,k)
c         write(*,*) 'siesold=',siesold
c         write(*,*) 'cpsolid(i,j,k)=',cpsolid(i,j,k)
c         write(*,*) 'rhosold=',rhosold
c         write(*,*) 'rhos(i,j,k)=',rhos(i,j,k)
c         write(*,*) 'frhosies(i,j,k)=',frhosies(i,j,k)
c         write(*,*) 'frhosiesrad(i,j,k)=',frhosiesrad(i,j,k)
c       endif

c     if((i.eq.20.and.k.eq.1).and.mpi_rank.eq.0) 
c    .print*,i,k,sies(i,j,k),cpsolid(i,j,k)

c            if (temps(i,j,k).lt.270..and.k.lt.8) then
c               write (*,*) 'slipping temps(i,j,k)'
c               write (*,*) i,j,k,temps(i,j,k),tempsold,sies(i,j,k),
c     &         rhos(i,j,k),frhosies(i,j,k),frhosiesrad(i,j,k)*.5
c            write (*,*) 'Data for iproc=',mpi_rank,
c     &                   ' i=',i,' j=',j,' k=',k
c            write (*,*) 'rhof=',rhof(i,j,k)
c            write (*,*) 'rhofinitial=',rhofinitial(i,j,k)
c            write (*,*) 'rhowater=',rhowater(i,j,k)
c            write (*,*) 'rhos=',rhos(i,j,k)
c            write (*,*) 'cpsolid=',cpsolid(i,j,k)
c            write (*,*) 'fw=',fw(i,j,k)
c            write (*,*) 'ff=',ff(i,j,k)
c            write (*,*) 'psiwmax=',psiwmax(i,j,k)
c            write (*,*) 'psiw=',psiw(i,j,k)
c            write (*,*) 'temps=',temps(i,j,k)
c            write (*,*) 'temps(i,j,k)=',temps(i+1,j+1,k)
c            write (*,*) 'tempg(i,j,k)=',tempg(i+1,j+1,k)
c            write (*,*) 'convht=',convht(i,j,k)
c            write (*,*) 'frhosies=',frhosies(i,j,k)
c            write (*,*) ' '
c            write (*,*) 'This run is being terminated!!! '
c            call mpi_finalize
c              stop
c            endif
                  endif
c                  if (temps(i,j,k).ge.350) 
c     &               write (*,*) 'compare firad and frhosiesrad'
c     &                           ,i,j,k,temps(i,j,k),cpsolid(i,j,k),
c     &            frhosiesrad(i,j,k),ff(i,j,k)
cc                  if (temps(i,j,k).ge.1000) 
c                  if (temps(i,j,k).ge.350.and.mpi_rank.eq.1) 
cc     &      write (*,*) 'mpi_rank=',mpi_rank,i,j,k,temps(i,j,k),
cc     &                           tempg(i,j,k)
c     &                  ,frhosies(i,j,k),
cc     &                  frhosies(i,j,k),frhosies(i,j,k)/cpsolid(i,j,k)
c     &            frhosiesrad(i,j,k)*.5*dtp ,thetasolid(i,j,k)
c            if (mpi_rank.eq.1) then 
c            if (i.eq.2.and.j.eq.17.and.k.eq.3) then 
c            write (*,*) 
c     &        'Data for iproc=',iproc,' i=',i,' j=',j,' k=',k
c            write (*,*) 'rhof=',rhof(i,j,k)	
c            write (*,*) 'rhowater=',rhowater(i,j,k)	
c            write (*,*) 'rhos=',rhos(i,j,k)	
c            write (*,*) 'cpsolid=',cpsolid(i,j,k)	
c            write (*,*) 'fw=',fw(i,j,k)	
c            write (*,*) 'ff=',ff(i,j,k)	
c            write (*,*) 'psiwmax=',psiwmax(i,j,k)
c            write (*,*) 'psiw=',psiw(i,j,k)
c            write (*,*) 'temps=',temps(i,j,k)
c            write (*,*) 'temps(40,40,1)=',temps(i+1,j+1,k)
c            write (*,*) 'tempg(40,40,1)=',tempg(i+1,j+1,k)
c            write (*,*) 'convht=',convht(i,j,k)
c            write (*,*) 'frhosies=',frhosies(i,j,k)
c               endif
c              endif

               enddo
            enddo
         enddo
c      call rmaxmin1(frhosies,'frhosiesafter',1-ih,np+ih,1-ih,mp+ih,l)
c      call rmaxmin1(sies,'siesafter',1-ih,np+ih,1-ih,mp+ih,l)
         
!JAS 4/25/07       !if (igntype==2.and.irst==0)then   !if we use terratorch style ignition
       if (igntype==2.and.((time+restarttime).lt.endigntime))then   !if we use terratorch style ignition
c234567
c          startigntime=0.2
c          endigntime=5.0
c          xfirelinelow=20.
c          xfirelinehigh=20.
c          yfirelinelow=15.
c          yfirelinehigh=25.
c          flamedistance=6.
          rignitionperiod=endigntime-startigntime
          travelrate=sqrt((xfirelinelow-xfirelinehigh)**2
     &             +(yfirelinelow-yfirelinehigh)**2)/
     &          rignitionperiod
          pretime=flamedistance/travelrate
          if (pretime.gt.startigntime) then
            startigntime=pretime
            endigntime=pretime+(rignitionperiod)
          endif
   
          if ((time+restarttime).ge.startigntime-pretime
     &       .and.(time+restarttime).le.endigntime+pretime) then
c          write (*,*) time
      
            do k=1,l
               do j=1,mp
                  do i=1,np
                     ia=(npos-1)*np+i
                     ja=(mpos-1)*mp+j
                     x1=(float(ia)-.5)*dx
                     y1=(float(ja)-.5)*dy
c                          write (*,*) 'x',x1,xfirelinehigh,xfirelinelow
c                          write (*,*) 'y',y1,yfirelinehigh,yfirelinelow
                     finterp=((time+restarttime)-startigntime)
     &                      /(endigntime-startigntime)
                     x2=(xfirelinehigh-xfirelinelow)*finterp
     &                     +xfirelinelow
                     y2=(yfirelinehigh-yfirelinelow)*finterp
     &                     +yfirelinelow
                     if (y1.ge.yfirelinelow
     &                     .and.y1.le.yfirelinehigh
     &                     .and.k.lt.4.) then
             
                        dist=sqrt((x2-x1)**2+(y2-y1)**2)
                        tempsold=temps(i,j,k)
                        
                        temps(i,j,k)=max(temps(i,j,k),
     &                      min(tfire,
     &                          tfire-  tfire/flamedistance*dist)
     &                            )
                           if (temps(i,j,k).ne.tempsold)  then
                             cpsolid(i,j,k)=cpwood
                             sies(i,j,k)=cpsolid(i,j,k)*temps(i,j,k)
                             rhowater(i,j,k)=0.0
                             rmoist(i,j,k)=0.0
                             rhos(i,j,k)=rhof(i,j,k)
                           endif
      
                     endif
                  enddo
               enddo
            enddo
          endif !end if (time.ge.startigntime-pretime.and.time.le.endigntime+pretime)

       elseif ((igntype.eq.1).and.((time+restarttime).ge.startigntime).and.
     &  ((time+restarttime).le.(startigntime+(targettemp-tambient)/ramprate)))then
                            !elseif we use rinitfire style ignition
         do k=1,l
            do j=1,mp
               do i=1,np
                  if (ifirestart(i,j,k).eq.1) then
                     temps(i,j,k)=max(temps(i,j,k),
     &                 min(((time+restarttime)-startigntime)*(ramprate)
     &                 +tambient,targettemp))
                     sies(i,j,k)=cpsolid(i,j,k)*temps(i,j,k)
c                     if(mpi_rank.eq.0)then
c                       write (*,*) 'resetting temps in ',
c     &                 mpi_rank,' i=',i,' j=',j,' k=',k,
c     &                 ' temps=',temps(i,j,k),
c     &                 (time-startigntime)*(ramprate) 
c     &                 +tambientarray(i,j,k),targettemp
c                     endif
                  endif
               enddo
            enddo
         enddo
        endif    !end if terratorch or rinitfire style ignition
       
         do k=1,l
            do j=1,mp
               do i=1,np

                  fi(i,j,k)=0.
                  frho(i,j,k)=0.
                  fhc(i,j,k)=0.
                  fox(i,j,k)=0.
                  frhosies(i,j,k)=0.
                  frhof(i,j,k)=0.
                  frhowater(i,j,k)=0.
               enddo
            enddo
         enddo
      
c       call rmaxmin1(frhosiesrad,'frhosiesrad',1-ih,np+ih,1-ih,mp+ih,l)
c       call rmaxmin1(firad,'firad',10,40,4,16,10)
      
      enddo

      return
      end
