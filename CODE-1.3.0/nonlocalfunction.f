!***************************************************
       module nonlocalfunction
       Implicit None

       save

       contains
       
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
                  rhos(i,j,k)=rhof(i,j,k)+rhowater(i,j,k)
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
                     cpsolid(i,j,k)=(rhof(i,j,k)*
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
      end subroutine firetecnonlocal
*******************************************************************************
      subroutine fuelnonlocal(xv,xvfuel,il,iu,jl,ju,lls,nvp)
      use xve
      use metryic
      use fireteca
      use turba
      use turbb
      use constants
      use gridsetup
      use msga
      use pres
      use nonlocal

      Implicit None

      !JAS 3/7/06 added explicit declarationis to comply with implicit
      !none
      integer,intent(in) :: il,iu,jl,ju,lls,nvp
      
      !JAS 3/7/06 added explicit declarationis to comply with implicit
      !none
      integer :: i,j,k
      real,external :: zcart
      real :: rkctempref,zla,sigmac,th,frampmod,trampk,tfstepk !,speed2
      real :: rlcrit,rlramp,rljoint1,tjoint1,scalefactor,rpsi,rltop
      real :: rlramp2,psijoint2,rljoint2,rpsi2,rl,tshift,temptest
      real :: perchydroremaining,taps,tapsw,xmap !,hydrofactor
      real :: slambdaog,psiwmaxold,fglimit ! ,slambdaof
      real xv(il:iu,jl:ju,lls,nvp)
      real xvfuel(il:iu,jl:ju,lls,5)
      
      real er,xerf,c1,c2,c3
      real p,a1,a2,a3,terf,sum
      real firsttime
      real erf


      data c1,c2,c3/0.5,0.0079,1.0/
      data p,a1,a2,a3/0.47047,0.3480242,-0.0958798,0.7478556/
      
      real,allocatable:: rkctemp(:,:,:)
      
      firsttime = 300.0
 
      if(irod.eq.1) allocate (rkctemp(1-ih:np+ih, 1-ih:mp+ih,l))

      psif=0.
      psiw=0.

c      do k=2,l
      do k=1,l
         do j=1,mp
            do i=1,np
               rkctempref=.2*xv(i,j,k,6)/xv(i,j,k,nv)
c it should be noted that xvfuel(i,j,k,5) corrspondes to xv(i,j,k,6)
c the .2 factor is included as an extrapolation for turbulence at the c
c scale.  The .2 factor should be replaced by a function of sc and sb.

               rkctemp(i,j,k)=.2*xvfuel(i,j,k,5)/xv(i,j,k,nv)
c               if (tempg(i,j,k).gt.500)
c     &               write (*,*) rkctempref,rkctemp(i,j,k)
            enddo
         enddo
      enddo
c         do j=1,mp
c            do i=1,np
c               rkctemp(i,j,1)=rkctemp(i,j,2)
c            enddo
c         enddo
      
      do k=1,l
         do j=1,mp
            do i=1,np
               zla=zcart(z(k),i,j)-zs(i,j)
               if(zla.le.fueldepth) then
c                  tamb=xe(i,j,k,4)/xe(i,j,k,nv)*(pr(i,j,k)*1.e-5)
c     .**(rg/cp)
c                   tamb=300.0
                  sigmac=0.
                  if (rkctemp(i,j,k).gt.1.e-06)
c   readdress the .09 s k**.5 issue             !(rrl)
c the traditional .09 has been left out of sigmac, however the effect of
c value is absorbed directly in cf and sigmac is only used in the
c reaction rate
c the .5 was inserted in the conversion to isotropic turb.  It 
c can be thought of as part of cf.
     &               sigmac=sc*0.5*sqrt(rkctemp(i,j,k))
c     &               sigmac=xv(i,j,k,nv)*sc(i,j,k)*0.5*
c     &               ( d11b(i,j,k)+d11b(i+1,j,k)
c     &                +d22b(i,j,k)+d22b(i,j+1,k)
c     &                +d33b(i,j,k)+d33b(i,j,k+1))
c     &               /sqrt(0.2*rkctemp(i,j,k))

c taps is the bottom of the ramp and th is the top of the ramp
                  th=2*tcrit-tfstep
             
c      speed2=(xv(i,j,k,1)/xv(i,j,k,nv))**2
c     &      +(xv(i,j,k,2)/xv(i,j,k,nv))**2
c     &      +(xv(i,j,k,3)/xv(i,j,k,nv))**2

c----------------------------------RRL's
      if(ifuel==1)then            !Use RRL's method of calculating psif!
c      frampmod=1.0
c      frampmod=min(1.0,exp(-1.5*(xv(i,j,k,6)/xv(i,j,k,nv)-.0002)))
c     &        *min(1.0,exp(-4.5*(speed2*.1-.04)))
c      trampk=-(tcrit-tramp)*frampmod+tcrit
c      tfstepk=-(tcrit-tfstep)*frampmod+tcrit
c
c      tjoint=psijoint*2.*(tcrit-trampk)+trampk-tfstepk
c      tjoint2=(2*tcrit-tfstepk)-
c     &           ((1.-psijoint)*2.*(tcrit-trampk)+trampk)
c      bramp=(psijoint-tjoint/6.0/(tcrit-trampk))*3./tjoint**2
c      aramp=(1./(2*(tcrit-trampk))-2.*bramp*tjoint)/(3.*tjoint**2)
c
cc taps is the bottom of the ramp and th is the top of the ramp
c                  th=2*tcrit-tfstepk
c
c                  if (temps(i,j,k).le.tfstepk ) then
c                     psif(i,j,k)=0.0
c                  elseif (temps(i,j,k).le.tjoint+tfstepk) then
c                     tshift=temps(i,j,k)-tfstepk
c                     psif(i,j,k)=aramp*tshift**3+tshift**2*bramp
c                  elseif (temps(i,j,k).le.th-tjoint2) then
c                     psif(i,j,k)=.5*(temps(i,j,k)-tramp)/(tcrit-trampk)
c                  elseif (temps(i,j,k).le.th) then
c                     tshift=th-temps(i,j,k)
c                     psif(i,j,k)=1.-(aramp*tshift**3+tshift**2*bramp)
c                  elseif (temps(i,j,k).gt.th) then
c                     psif(i,j,k)=1.
c                  endif
      frampmod=1.
c      frampmod=exp-(7.*max(0.,xv(15,20,1,6)))
      trampk=-(tcrit-tramp)*frampmod+tcrit
      trampk=tramp
      tfstepk=tfstep

      rlcrit=(tcrit-tfstepk)
      rlramp=(trampk-tfstepk)
c      write (*,*) 'rlcrit is',rlcrit,' rlramp is',rlramp

      rljoint1=psijoint*2.*(rlcrit-rlramp)+rlramp
      tjoint1=rljoint1+tfstepk

      scalefactor=sqrt(rlramp**2-(rljoint1-rlramp)**2)/psijoint
      rpsi=2.0/scalefactor*(rlcrit-rlramp)*rljoint1
     &     +scalefactor*psijoint


      rltop=2.*rlcrit
      rlramp2=rltop-rlramp
      psijoint2=1.-psijoint
c      rljoint2=rlcrit+(psijoint2-.5)*(rlcrit-rlramp)*2.
      rljoint2=rltop-rljoint1
      tjoint2=(rljoint2+tfstepk)
      rpsi2=scalefactor-rpsi
Cc      write (*,*) 'rpsi is ',rpsi,' rpsi2 is ',rpsi2
C
C
Cc      write (*,*) psijoint, tcrit,trampk,tfstepk,tjoint2
Cc      write (*,*) psijoint, rlcrit,rlramp,rljoint1,rljoint2
Cc      write (*,*) 'the bottom of the curve is at ',tfstepk
Cc      write (*,*) 'the first joint is at ',tjoint1
Cc      write (*,*) 'the first joint after shift is at ',rljoint1
Cc      write (*,*) 'the second joint is at ',tjoint2
Cc      write (*,*) 'the second joint after shift is at ',rljoint2
Cc      write (*,*) 'the value at the first joint is',
Cc     &   1./scalefactor*(rpsi-sqrt(rpsi**2-rljoint1**2))
Cc      write (*,*) 'the value at the first joint should be',psijoint
Cc      write (*,*) 'the slope at the first joint is ',
Cc     &    .5/scalefactor*2.*rljoint1/sqrt(rpsi**2-rljoint1**2)
Cc      write (*,*) 'the slope at the first joint should be ',
Cc     &              .5/(rlcrit-rlramp)
Cc      write (*,*) 'the value at the second joint be',
Cc     & 1./scalefactor*(rpsi2+sqrt(rpsi**2-(rljoint2-rltop)**2))
Cc      write (*,*) 'the value at the second joint should be',
Cc     &    psijoint2
Cc      write (*,*) 'the slope at the second joint is ',
Cc     &    -.5/scalefactor*2.*(rljoint2-rltop)
Cc     &          /sqrt(rpsi**2-(rljoint2-rltop)**2)
Cc      write (*,*) 'the slope should be ',1/(tcrit-trampk)*.5
C
                 rl=temps(i,j,k)-tfstepk
                  if (rl.le.0.0 ) then
        psif(i,j,k)=0.0
                  elseif (rl.le.rljoint1) then
        psif(i,j,k)=1./scalefactor*
     &              (rpsi-sqrt(rpsi**2-rl**2))
                  elseif (rl.le.rljoint2) then
        psif(i,j,k)=.5/(rlcrit-rlramp)*(rl-rlramp)
                  elseif (rl.le.rltop) then
                     tshift=th-temptest
        psif(i,j,k)=1./scalefactor*
     &         (rpsi2+sqrt(rpsi**2-(rl-rltop)**2))
                  elseif (rl.gt.rltop) then
        psif(i,j,k)=1.
                  endif
c------------------------------- end RRL's

      elseif(ifuel==2)then            !Use Michael Clark's method of
calculating psif!
                  
c*******************************************************************
c Michael Clark 7/11/05
c 
c A subroutine that calculates the error function and scales it
c
c  erf(x) is calculated using a rational approximation
c    source: Handbook of Mathematical Functions
c            Edited by: Abramowitx & Stegun
c            Dover Publications, Inc. New York
c
c   EQ 7.1.25 (p. 299)
c     erf x = 1-(a1*t + a2*t**2 + a3*t**3)*exp(-x**2) + eps(x)
c
c         t = 1/(1+p*x)
c
c     |eps(x)| < 2.5*10**-5
c
c       p =  0.47047
c      a1 =  0.3480242
c      a2 = -0.0958798
c      a3 =  0.7478556
c-------------------------------------------------------------------
c   scaled_erf = c1*(erf[c2*(t-tcrit)]+c3)
c
c   e.g. If c1=0.5 & c3=1.0, erf(x) will assymptotically approach
c        1 at +infinity and 0 at -infinity
c 
c  implemented by JJC 7/13/05
c******************************************************************* 
          
              xerf = c2*(temps(i,j,k)-tcrit)
              terf = 1./(1.+p*abs(xerf))
              sum  = a1*terf+a2*terf*terf+a3*terf*terf*terf
              er   = 1-sum*exp(-xerf*xerf)
              psif(i,j,k) = c1*(er+c3)
              if (xerf.lt.0.0) psif(i,j,k)=-psif(i,j,k)+c3
c
c-------------------------------------------------------------------i
      elseif(ifuel==3)then            !Use JAS's method of calculating
c  reimplemented by JAS 2/2/06 to use intrinsic erf (much faster!)

              
              !look closely at the integral erf calculates to see why
              !here!
c      psif(i,j,k) =
c1*(er*((temps(i,j,k)-tcrit)/abs(temps(i,j,k)-tcrit))+c3)  

              if(temps(i,j,k).le.tfstep)then 
                psif(i,j,k)=0.0
              elseif(temps(i,j,k).le.(2*(tcrit-tfstep)+tfstep))then
                er = erf(c2*(temps(i,j,k)-tcrit))                         
                psif(i,j,k) = c1*(c3+er)
              else 
                psif(i,j,k) = 1.0
              endif  
           
c---------------------------------------------------------------------
      endif    !end if(ifuel==1,2,or 3)        
                  
                  if (rhofinitial(i,j,k).gt.0.0) then
                  perchydroremaining=
     &       max(0.,(rhof(i,j,k)-rhohydrothresh*rhofinitial(i,j,k))
     .                 /(rhofinitial(i,j,k)*(1.-rhohydrothresh)))
C                cf(i,j,k)=cfhydro*perchydroremaining
C     .                    +cfchar*(1-perchydroremaining)

                cf(i,j,k) = 0.0008
              else
                cf(i,j,k)=0.0
              endif


c               hydrofactor=exp(-1.*rno/rnfuel*psif(i,j,k)
c     &                        *rhof(i,j,k)/xv(i,j,k,7))
c               thetasolid(i,j,k)=perchydroremaining*.25*hydrofactor
c     &                          +(1-perchydroremaining)*.75
              thetasolid(i,j,k)=perchydroremaining*.25
     &                          +(1-perchydroremaining)*.75
c              thetasolid(i,j,k)=0.0
c              if (thetasolid(i,j,k).le.0.d0.or.
c     &            thetasolid(i,j,k).ge.1.d0) then
c                write(*,*) 'thetasolid calculation'
c                write(*,*) mpi_rank,i,j,k,perchydroremaining,
c     &                     hydrofactor, rhof(i,j,k),psif(i,j,k),
c     &                     xv(i,j,k,7)
c              endif     
c               taps=tambient+tstepk+(tcrit-373)
               taps=tramp
c               th=2*tcrit-taps
               xmap=(temps(i,j,k)-taps)/(tcrit-taps)
               if (xmap.lt.0) xmap=0.0
               if (xmap.gt.2) xmap=2.0

c similar process for the water evap. ramp

c      tapsw=tambient+tstep*(twvap-tambent)/(tcrit-tambient)
c      tapsw=tambient+.25*tstep*(twvap-tambient)/(tcrit-tambient)
      tapsw=tambientarray(i,j,k)+tstep
c      tapsw=tamb+80.
      th=2*twvap-tapsw

      if(temps(i,j,k).gt.tapsw.and.temps(i,j,k).lt.th)
     .psiw(i,j,k)=(temps(i,j,k)-tapsw)/(th-tapsw)        !(rrl 6/1/99)
      if(temps(i,j,k).ge.th) psiw(i,j,k)=1.
      if(temps(i,j,k).le.tapsw) psiw(i,j,k)=0.
      psiwmaxold=psiwmax(i,j,k)
      psiwmax(i,j,k)=max(psiwmax(i,j,k),psiw(i,j,k))

c      slambdaof=rhof(i,j,k)*xv(i,j,k,7)/(rhof(i,j,k)/rnfuel+ !(rrl
c     .                                xv(i,j,k,7)/rno)**2
!***********************************
!***** Solid Fuel Forcing Term *****
!***********************************
c      ff(i,j,k)= cf(i,j,k)*rhof(i,j,k)*xv(i,j,k,7)*sigmac*
c     .          psif(i,j,k)*slambdaof/(rhoref*sx*sx)     !(rrl 6/1/99)
      ff(i,j,k)= cf(i,j,k)*rhof(i,j,k)*sigmac*
     .          psif(i,j,k)/(sx*sx)                       !(jjc 5/11/06)     

      fw(i,j,k)=0.0
      if (psiwmax(i,j,k).lt.1.and.rhowater(i,j,k).gt.0.0)
     & fw(i,j,k)=min(rhowater(i,j,k)*
     &           (psiwmax(i,j,k)-psiwmaxold)/
     &          (1.00-psiwmax(i,j,k))/dtp,rhowater(i,j,k)/dtp)

      endif    ! end of fuel depth block
      
c    now set psif term for gas phase combustion
            if(xv(i,j,k,8).gt.1.d-12.and.
     &            tempg(i,j,k).gt.tambient) then
     
       sigmac=sc*0.5*sqrt(rkctemp(i,j,k))
       th=2*tcrit-tfstep
c----------------------------------RRL's
      if(ifuel==1)then            !Use RRL's method of calculating psif

      frampmod=1.
c      frampmod=exp-(7.*max(0.,xv(15,20,1,6)))
      trampk=-(tcrit-tramp)*frampmod+tcrit
      trampk=tramp
      tfstepk=tfstep

      rlcrit=(tcrit-tfstepk)
      rlramp=(trampk-tfstepk)
c      write (*,*) 'rlcrit is',rlcrit,' rlramp is',rlramp

      rljoint1=psijoint*2.*(rlcrit-rlramp)+rlramp
      tjoint1=rljoint1+tfstepk

      scalefactor=sqrt(rlramp**2-(rljoint1-rlramp)**2)/psijoint
      rpsi=2.0/scalefactor*(rlcrit-rlramp)*rljoint1
     &     +scalefactor*psijoint


      rltop=2.*rlcrit
      rlramp2=rltop-rlramp
      psijoint2=1.-psijoint
c      rljoint2=rlcrit+(psijoint2-.5)*(rlcrit-rlramp)*2.
      rljoint2=rltop-rljoint1
      tjoint2=(rljoint2+tfstepk)
      rpsi2=scalefactor-rpsi

                 rl=tempg(i,j,k)-tfstepk
                  if (rl.le.0.0 ) then
        psig(i,j,k)=0.0
                  elseif (rl.le.rljoint1) then
        psig(i,j,k)=1./scalefactor*
     &              (rpsi-sqrt(rpsi**2-rl**2))
                  elseif (rl.le.rljoint2) then
        psig(i,j,k)=.5/(rlcrit-rlramp)*(rl-rlramp)
                  elseif (rl.le.rltop) then
                     tshift=th-temptest
        psig(i,j,k)=1./scalefactor*
     &         (rpsi2+sqrt(rpsi**2-(rl-rltop)**2))
                  elseif (rl.gt.rltop) then
        psig(i,j,k)=1.
                  endif
c------------------------------- end RRL's

      elseif(ifuel==2)then            !Use Michael Clark's method of
calculating psi!
                  
c*******************************************************************
c Michael Clark 7/11/05
c 
c A subroutine that calculates the error function and scales it
c
c  erf(x) is calculated using a rational approximation
c    source: Handbook of Mathematical Functions
c            Edited by: Abramowitx & Stegun
c            Dover Publications, Inc. New York
c
c   EQ 7.1.25 (p. 299)
c     erf x = 1-(a1*t + a2*t**2 + a3*t**3)*exp(-x**2) + eps(x)
c
c         t = 1/(1+p*x)
c
c     |eps(x)| < 2.5*10**-5
c
c       p =  0.47047
c      a1 =  0.3480242
c      a2 = -0.0958798
c      a3 =  0.7478556
c-------------------------------------------------------------------
c   scaled_erf = c1*(erf[c2*(t-tcrit)]+c3)
c
c   e.g. If c1=0.5 & c3=1.0, erf(x) will assymptotically approach
c        1 at +infinity and 0 at -infinity
c 
c  implemented by JJC 7/13/05
c******************************************************************* 
          
              xerf = c2*(tempg(i,j,k)-tcrit)
              terf = 1./(1.+p*abs(xerf))
              sum  = a1*terf+a2*terf*terf+a3*terf*terf*terf
              er   = 1-sum*exp(-xerf*xerf)
              psig(i,j,k) = c1*(er+c3)
              if (xerf.lt.0.0) psig(i,j,k)=-psig(i,j,k)+c3
c
c-------------------------------------------------------------------i
      elseif(ifuel==3)then            !Use JAS's method of calculating
c  reimplemented by JAS 2/2/06 to use intrinsic erf (much faster!)


              if(tempg(i,j,k).le.tfstep)then 
                psig(i,j,k)=0.0
              elseif(tempg(i,j,k).le.(2*(tcrit-tfstep)+tfstep))then
                er = erf(c2*(tempg(i,j,k)-tcrit))                         
                psig(i,j,k) = c1*(c3+er)
              else 
                psig(i,j,k) = 1.0
              endif  
           
c---------------------------------------------------------------------
      endif    !end if(ifuel==1,2,or 3)      

c     now we set slambdaog, the stoichiometry term

      slambdaog=xv(i,j,k,8)*xv(i,j,k,7)/(xv(i,j,k,8)/rng+     !(jjc
     &                                   xv(i,j,k,7)/rnonl)**2
     
c     creating a 'flux limiter' term here

       fglimit = tanh((xv(i,j,k,8)-1.d-6)*(5.d5/pi))

c     Finally we set the 'fuel forcing' term

      fg(i,j,k)=(xv(i,j,k,8)*xv(i,j,k,7)*cg*sigmac*psig(i,j,k)
     &          *slambdaog*fglimit)
     &          /(2*sx*sx*(rnonl*xv(i,j,k,8)+rng*xv(i,j,k,7)))     

c      ia=(npos-1)*np+i  !  This is just for the if statement below.
c      ja=(mpos-1)*mp+j  !  It's the global position.
c      if(fg(i,j,k).ne.0) then
c      if(tempg(i,j,k).gt.350.and.ifuelcount.eq.1.and.k.eq.1) then
c        write(*,*) ia,ja,k
c        write(*,*) 'tg=', tempg(i,j,k),' fg=', fg(i,j,k), 
c     &             'ts=', temps(i,j,k),' ff=', ff(i,j,k)
c        write(*,*) 'sigmac=', sigmac,' psig=',psig(i,j,k),
c     &             ' psif=',psif(i,j,k)
c        write(*,*) 'slambdaog = ',slambdaog,' rhof= ',rhof(i,j,k),
c     &             'xv8 = ',xv(i,j,k,8),'  xv7= ',xv(i,j,k,7)
c      endif
      endif  !end of the xv8, tempg shunt      
   
      enddo
      enddo
      enddo
      !call rmaxmin1(ff,'ff',1-ih,np+ih,1-ih,mp+ih,l)
      !call rmaxmin1(fw,'fw',1-ih,np+ih,1-ih,mp+ih,l)
c      call rmaxmin1(psiwmax,'psiwmax',1-ih,np+ih,1-ih,mp+ih,l)

      if(irod.eq.1) deallocate (rkctemp)
 
      return
      end subroutine fuelnonlocal

*******************************************************************************
      subroutine energysourcenonlocal(xv,il,iu,jl,ju,lls,nvp)
      use turba
      use pres
      use fireteca
      use constants
      use gridsetup
      use xvo
      use xve
      use msga
      use nonlocal
      use radiation, only:irad

      Implicit None

      !JAS 3/7/06 added explicit declarations to comply with implicit
      !none
      integer,intent(in) :: il,iu,jl,ju,lls,nvp

      real xv(il:iu,jl:ju,lls,nvp)
      real,allocatable:: u(:, :,:),
     .                   v(:, :,:), 
     .                   w(:, :,:)
c     .                  ,tbar(:,: ,:)
      
      !JAS 3/7/06 added explicit declarations to comply with implicit
      !none
      integer :: i,j,k ! ,irhohydro
      real :: rrhomicro,sstemp,ufuel,vfuel,wfuel ! ,rhovapor,rhohydro
      real :: rktemp,usp,re,h,av,convhtb,qxt
      real :: gammaterm,rneteng,rwatergainht
      real :: rnetmass,rnetengwater,capqterm,energy,rhogas
     
      if(irod.eq.1) then
      allocate (u(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (v(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (w(1-ih:np+ih, 1-ih:mp+ih,l))
c      allocate (tbar(1-ih:np+ih, 1-ih:mp+ih,l))
      endif

      do k=1,l
      do j=1,mp
      do i=1,np
      u(i,j,k)=xv(i,j,k,1)/xv(i,j,k,nv)
      v(i,j,k)=xv(i,j,k,2)/xv(i,j,k,nv)
      w(i,j,k)=xv(i,j,k,3)/xv(i,j,k,nv)
c     tbar(i,j,k)=xv(i,j,k,4)/xv(i,j,k,nv)*(pr(i,j,k)*1.e-5)**(rg/cp)
      enddo
      enddo
      enddo

      do k=1,l
      do j=1,mp
      do i=1,np
               rrhomicro=1./rhomicro(i,j,k)
               sstemp=sizescale(i,j,k) !FP

      rhogas=xv(i,j,k,nv)

c      usp=sqrt(u(i,j,k)*u(i,j,k)+v(i,j,k)*v(i,j,k)+w(i,j,k)*w(i,j,k))
c       rktemp=xv(i,j,k,6)/xv(i,j,k,nv)+
c     +            xv(i,j,k,5)/xv(i,j,k,nv)
c      re=ss*(usp+sqrt(rktemp))/2.e-05
!rrl
c      h=2*0.683*re**0.466*thermcondair/ss                          !rrl
                  ufuel=xvfuel(i,j,k,1)/xv(i,j,k,nv)
                  vfuel=xvfuel(i,j,k,2)/xv(i,j,k,nv)
                  wfuel=xvfuel(i,j,k,3)/xv(i,j,k,nv)
                  rktemp=xvfuel(i,j,k,4)/xv(i,j,k,nv)+
     +            1.2*(xvfuel(i,j,k,5)/xv(i,j,k,nv))

                  usp=sqrt(rktemp)+sqrt(ufuel*ufuel
     &                         +vfuel*vfuel+wfuel*wfuel)
                  re=sizescale(i,j,k)*usp/2.e-05 !FP
!rrl
c                  h=0.683*re**0.466*thermcondair/ss                !rrl
c                  h=2*0.683*re**0.466*thermcondair/ss        !rrl*2
c                  h=2*0.683*re**0.466*thermcondair/ss  
                  h=.5*0.683*re**0.466*thermcondair/sizescale(i,j,k) !FP
!4/23/03
c      h=0.683*re**0.466*thermcondair/ss                          !rrl
c      av=2.*(rhof(i,j,k)*0.004+rhowater(i,j,k)*.0001)/ss    !rrl
      av=2.*(rhof(i,j,k)*rrhomicro)/sizescale(i,j,k) !FP
      sizescale(i,j,k)=sstemp  !FP
      tambientarray(i,j,k)=xe(i,j,k,4)/xe(i,j,k,nv)
     &      *(pre(i,j,k)*1.0e-5)**(rg/cp)
      convhtb=h*av*(temps(i,j,k)-tempg(i,j,k))
c      if(mpi_rank.eq.0) write (*,*) 'energysource',j,k,temps(i,j,k),
c     .   tempg(i,j,k),xv(i,j,k,4),pr(i,j,k),rg,cp

      
      

c      qxt=rke*sigma/sb(i,j,k)*(tempg(i,j,k)**4-tambient**4)  !wss
      qxt=rke*sigma/sqrt(dx*dy)*(tempg(i,j,k)**4
     &                           -tambientarray(i,j,k)**4)  !wss
      if(irad.GE.1) qxt=0.                                  !rrl !KOO
c      gammaterm= cvvapor*rhovapor/(gammav*cv*rhogas)
c     .          +cvoxygen*xv(i,j,k,7)/(gammo*cv*rhogas)
c     .          +(rhogas-rhovapor-xv(i,j,k,7))
c     .                   /(gamma*rhogas)
c      gammaterm=gammaterm*(convhtb+firad(i,j,k)-qxt)
      gammaterm=convhtb+firad(i,j,k)-qxt
      if (xvb(i,j,k,8).ne.0.d0) then
        rneteng=(1.-(thetasolid(i,j,k)*ff(i,j,k)*dtp*rnhc/
     &         xvb(i,j,k,8)))*hfgas*fg(i,j,k)
      else
        rneteng=0.d0
      endif      
      rwatergainht=fw(i,j,k)*cpwater*twvap      !same as below
      rnetmass=ff(i,j,k)*tcrit*cpwood+rwatergainht   !rwaterg.. term
                                                           !cpwood added
                                                           !jjc rrl
                                                           !09/21/01
      rnetengwater=fw(i,j,k)*cvvapor*twvap
      energy=gammaterm+rneteng+rnetmass
c      if (gammaterm.lt.-100000.or.
c     &     tempg(i,j,k).gt.900)then 
c      if (-1*firad(i,j,k).lt.rneteng+convhtb.or.
c     &     tempg(i,j,k).gt.1500) then
c      pause
c      write (*,*) 'firad no2',i,
c     &     j,k,tempg(i,j,k),gammaterm,rneteng
c     &       ,convhtb,firad(i,j,k),xvb(i,j,k,7)/xvb(i,j,k,8),rnetmass,
c     &       rnetengwater
c      endif
c      endif
      capqterm=energy*(1.e5/pr(i,j,k))**(rg/cp)/cp             !rrl
c      fi(i,j,k)=.5*fib(i,j,k)*dtp+2.*energy*dtp               !rrl
      fi(i,j,k)=fib(i,j,k)*dtp+2.*capqterm*dtp                 !rrl
c      if(mpi_rank.eq.0.and.temps(i,j,k).gt.300) 
c     & write(*,*) 'energysource',i,j,k,temps(i,j,k),
c     &       tempg(i,j,k),energy,capqterm,fi(i,j,k)
c      temp=tbar(i,j,k)+0.5*fi(i,j,k)/(cv*xv(i,j,k,nv))
c      temp=tempg(i,j,k)+0.5*fi(i,j,k)/xv(i,j,k,nv)
c     .     *(pr(i,j,k)*1.e-5)**(rg/cp)
c     if(temp.lt.tambient.and.abs(energy).gt.0.) 
c    .print*,i,k,convhtb,temps(i,j,k),tempg(i,j,k)
c     if((i.eq.23.and.k.eq.3).and.mpi_rank.eq.0) 
c    .print*,fi(22,1,3),'fi',rneteng,ff(i,j,k),hf,thetag
c     tmp=(xv(i,j,k,4)+0.5*fi(i,j,k))/(cv*xv(i,j,k,nv))
c     if(tmp.lt.tbar(i,j,k)) fi(i,j,k)=
c    .2.*(tbar(i,j,k)*cv*xv(i,j,k,nv)-xv(i,j,k,4))
      enddo
      enddo
      enddo
      deallocate (u)
      deallocate (v)
      deallocate (w)



    
      return
      end subroutine energysourcenonlocal
*******************************************************************************
      subroutine fueltempnonlocal(xv,xvfuel,il,iu,jl,ju,lls,nvp)
      use turba
      use pres
      use fireteca
      use constants
      use gridsetup
      use msga
      use nonlocal

      Implicit None

      !JAS 3/7/06 added explicit declarations to comply with implicit
      !none
      integer,intent(in) :: il,iu,jl,ju,lls,nvp

      real xv(il:iu,jl:ju,lls,nvp) 
      real xvfuel(il:iu,jl:ju,lls,5) 

      !JAS 3/7/06 added explicit declarations to comply with implicit
      !none
      integer :: i,j,k
      real :: rrhomicro,sstemp,u,v,w,uref,vref,wref,rktemp,rktempref
      real :: sp,spref,re,reref,h,href,av,reactht,rmassloss,waterevpht
      real :: rwaterlossht,tmp

      do k=1,l
         do j=1,mp
            do i=1,np
               rrhomicro=1./rhomicro(i,j,k)
               !sstemp=ss
               sstemp=sizescale(i,j,k) !FP

               if (rhof(i,j,k).gt.1.e-04) then
                  uref=xv(i,j,k,1)/xv(i,j,k,nv)
                  vref=xv(i,j,k,2)/xv(i,j,k,nv)
                  wref=xv(i,j,k,3)/xv(i,j,k,nv)
                  rktempref=xv(i,j,k,6)/xv(i,j,k,nv)+
     +            xv(i,j,k,5)/xv(i,j,k,nv)
                  spref=sqrt(rktempref)
     &                 +sqrt(uref*uref+vref*vref+wref*wref) 
                  !reref=ss*spref/2.e-05
                  !!rrl
                  reref=sizescale(i,j,k)*spref/2.e-05                !FP
!rrl
                  !href=2*0.683*reref**0.466*thermcondair/ss
                  !!rrl*2 8/15/01
                  href=2*0.683*reref**0.466*thermcondair/sizescale(i,j,k)
!FP               !rrl*2 8/15/01

                  u=xvfuel(i,j,k,1)/xv(i,j,k,nv)
                  v=xvfuel(i,j,k,2)/xv(i,j,k,nv)
                  w=xvfuel(i,j,k,3)/xv(i,j,k,nv)
                  rktemp=xvfuel(i,j,k,4)/xv(i,j,k,nv)+
     +            1.2*(xvfuel(i,j,k,5)/xv(i,j,k,nv))

                  sp=sqrt(rktemp)+sqrt(u*u+v*v+w*w) 
                  !re=ss*sp/2.e-05                              !rrl
                  re=sizescale(i,j,k)*sp/2.e-05
!rrl
                  !h=.5*0.683*re**0.466*thermcondair/ss         !rrl
                  h=.5*0.683*re**0.466*thermcondair/sizescale(i,j,k) !FP
!rrl
c                  h=2*0.683*re**0.466*thermcondair/ss         !rrl*2
c                  h=2*0.683*re**0.466*thermcondair/ss         !rrl*2
c      av=2.*(rhof(i,j,k)*0.004+rhowater(i,j,k)*.0001)/ss      !rrl
c      av=2.*(rhof(i,j,k)*.004)/ss      !rrl
      av=2.*(rhof(i,j,k)*rrhomicro)/sizescale(i,j,k)      !rrl
      sizescale(i,j,k)=sstemp
c      tg=xv(i,j,k,4)/xv(i,j,k,nv)*(pr(i,j,k)*1.0e5)**(rg/cp)  !rrl
      convht(i,j,k)=h*av*(tempg(i,j,k)-temps(i,j,k))    !wss
c      solidqxt=av*sigma*solidemisivity*(temps(i,j,k))
      if (xv(i,j,k,8).ne.0.d0) then
        reactht=(thetasolid(i,j,k)*ff(i,j,k)*dtp*rnhc/xv(i,j,k,8))
     &            *hfgas*fg(i,j,k)+hfsolid*ff(i,j,k)
!jjc
      else
        reactht=hfsolid*ff(i,j,k)
      endif
      rmassloss=-tcrit*cpwood*ff(i,j,k)                 !rrl
      waterevpht=-fw(i,j,k)*hwevap                             !rrl
      rwaterlossht=-fw(i,j,k)*cpwater*twvap                    !rrl

      tmp= convht(i,j,k)                                       !rrl
     .     +reactht+rmassloss
     .     +waterevpht+rwaterlossht
      
c      frhosies(i,j,k)=frhosies(i,j,k)+tmp*dtp
      frhosies(i,j,k)=tmp*dtp
c      qflux(i,j,k)=convht(i,j,k)+reactht+waterevpht
      qflux(i,j,k)=convht(i,j,k)+reactht        !rrl
c      if (i.gt.39.and.i.lt.43.and.j.gt.36.and.j.lt.40..and.k.eq.2) 
c234567
c       if (temps(i,j,k).lt.-11000.or.temps(i,j,k).gt.14000) 
c     &   write (*,*) i,j,k,rhos(i,j,k),rhof(i,j,k),rhowater(i,j,k)
c     &   ,temps(i,j,k),tempg(i,j,k),
c     &fw(i,j,k),qflux(i,j,k),rwaterlossht,
c     &   waterevpht,convht(i,j,k),reactht
      endif
      enddo
      enddo
      enddo
c     call rmaxmin1(convht,'convht',1-ih,np+ih,1-ih,mp+ih,l)
c     call rmaxmin1(frhosies,'frhosies',1-ih,np+ih,1-ih,mp+ih,l)
     
 
      return
      end subroutine fueltempnonlocal
*******************************************************************************
      subroutine lapdf(rk,fk,s,d11,d22,d33,d12,d13,d23,il,iu,jl,ju,lls)
      use gridsetup
      use turba
      use metryic
      use msga
      
      Implicit None

      !JAS 3/7/06 added explicit declarations to comply with implicit
      !none
      integer,intent(in) :: il,iu,jl,ju,lls

      real      rk(il:iu, jl:ju,lls),
     .           s(il:iu, jl:ju,lls),
     .          fk(il:iu, jl:ju,lls)
      real            d13(il:iu+1, jl:ju,lls+1), 
     .                d12(il:iu+1, jl:ju+1,lls),
     .                d23(il:iu, jl:ju+1,lls+1),
     .                d11(il:iu+1, jl:ju,lls),
     .                d22(il:iu, jl:ju+1,lls),
     .                d33(il:iu, jl:ju,lls+1)
      real,allocatable::
     .           r(:, :,:),
     .          hx(:, :,:),
     .          hy(:, :,:),
     .          hz(:, :,:),
     .          pz(:, :,:),
     .        srff(:, :),
     .           u(:, :,:),
     .           v(:, :,:),
     .           w(:, :,:)

      !JAS 3/7/06 added explicit declarations to comply with implicit
      !none
      integer :: i,j,k,illim,iulim,jllim,julim,ip1,im1,jp1,jm1
      real :: hdxi,hdyi,hdzi,g33,g23,g13,gii,px,py,pxa,pya,pza
      real :: saa,sqrtka 
      real :: hxa,hya,hza,hx1,hy1,hx2,hy2,hxp1,hxp2,hyp1,hyp2,hzp1,hzp2
      real :: hxn,hxnp1,hymp1,hzLp1
      real :: coef,coefa,d12a,d13a,d21a,d23a,d31a,d32a


      if(iturb.ge.1) then
      allocate (r(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (hx(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (hy(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (hz(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (pz(1-ih:np+ih+1, 1-ih:mp+ih+1,l+1))
      allocate (srff(1-ih:np+ih, 1-ih:mp+ih))
      allocate (u(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (v(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (w(1-ih:np+ih, 1-ih:mp+ih,l))
      endif
 
      hdxi=0.5*dxi
      hdyi=0.5*dyi
      hdzi=0.5*dzi
      call updated(rk,rk,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)

compute z-derivatives at (i,j,k+-1/2)
      do k=2,L
        do j=1,mp
          do i=1,np
            pz(i,j,k)=dzi*(rk(i,j,k)-rk(i,j,k-1))
          end do
        end do
      end do
      do j=1,mp
        do i=1,np
          pz(i,j,L+1)=-pz(i,j,L)
        end do
      end do

      do j=1,mp
      do i=1,np
        srff(i,j)=0.
      enddo
      enddo

      do j=1,mp
      do i=1,np
      g33=(c13(i,j)*gmul(1))**2+(c23(i,j)*gmul(1))**2+gi(i,j,1)**2
      srff(i,j) = sqrt(g33)*srff(i,j)
      enddo
      enddo
c This section needs to be addressed with new ibcx=0 bc rrl       
      jllim = 1  + (j3-ibcy)*botdedge                 !added d rrl
      julim = mp + (ibcy-j3)*topdedge                 !added d rrl
      do 21 j=jllim,julim
         if (j3.eq.1) then
            if (topdedge.eq.1 .and. j.eq.mp) then                 !added
d rrl
               jp1 = mp + 2
            else
               jp1 = j + 1
            end if
            if (botdedge.eq.1 .and. j.eq.1) then                 !added
d rrl
               jm1 = -1
            else
               jm1 = j - 1
            end if
         else
            jp1=1
            jm1=1
         end if
c        jp1=j+j3-j/m*(m-1)
c        jm1=j-j3+(m-j)/(m-j3)*(m-j3)
         illim = 1  + (1-ibcx)*leftdedge                 !added d rrl
         iulim = np + (ibcx-1)*rightdedge                 !added d rrl
         do 21 i=illim,iulim
            if (rightdedge.eq.1 .and. i.eq.np) then
!added d rrl
               ip1 = np + 2
            else
               ip1 = i + 1
            end if
            if (leftdedge.eq.1 .and. i.eq.1) then                 !added
d rrl
               im1 = -1
            else
               im1 = i - 1
            end if
      g13=c13(i,j)*gmul(1)
      g23=c23(i,j)*gmul(1)
      g33=g13**2+g23**2+gi(i,j,1)**2
      px=hdxi*(rk(ip1,j,1)-rk(im1,j,1))
      py=hdyi*(rk(i,jp1,1)-rk(i,jm1,1))*j3
   21 pz(i,j,1)=-2.*(srff(i,j)+g13*px+g23*py)/g33-pz(i,j,2)
      if(ibcx.eq.0) then
         illim = 1*leftdedge   + np*(1-leftdedge)                 !added
d rrl
         iulim = np*rightdedge +  1*(1-rightdedge)
!added d rrl
         do 211 i=illim,iulim,np-1
            jllim = 1  + (j3-ibcy)*botdedge                 !added d rrl
            julim = mp + (ibcy-j3)*topdedge                 !added d rrl
            do 2111 j=jllim,julim
               if (j3.eq.1) then
                  if (topdedge.eq.1 .and. j.eq.mp) then
!added d rrl
                     jp1 = mp + 2
                  else
                     jp1 = j + 1
                  end if
                  if (botdedge.eq.1 .and. j.eq.1) then
!added d rrl
                     jm1 = -1
                  else
                     jm1 = j - 1
                  end if
               else
                  jp1=1
                  jm1=1
               end if
c        jp1=j+j3-j/m*(m-1)
c        jm1=j-j3+(m-j)/(m-j3)*(m-j3)
               g13=c13(i,j)*gmul(1)
               g23=c23(i,j)*gmul(1)
               g33=g13**2+g23**2+gi(i,j,1)**2
               py=hdyi*(rk(i,jp1,1)-rk(i,jm1,1))*j3
 2111      pz(i,j,1)=-2.*(srff(i,j)+g23*py)/(g33-g13*g13)-pz(i,j,2)
           if(ibcy.eq.0.and.j3.eq.1) then
              jllim = 1*botdedge  + mp*(1-botdedge)
!added d rrl
              julim = mp*topdedge +  1*(1-topdedge)
!added d rrl
              do 2112 j=jllim,julim,mp-j3
 2112         pz(i,j,1)=-2.*srff(i,j)/(gi(i,j,1)**2)-
     .                  pz(i,j,2)
           endif
c what is the value of px?  
 211    continue
      endif
      if(ibcy.eq.0.and.j3.eq.1) then
         jllim = 1*botdedge  + mp*(1-botdedge)                 !added d
         julim = mp*topdedge +  1*(1-topdedge)                 !added d
         do 212 j=jllim,julim,mp-j3
            illim = 1  + (1-ibcx)*leftdedge                 !added d rrl
            iulim = np + (ibcx-1)*rightdedge                 !added d
            do 2121 i=illim,iulim
               if (rightdedge.eq.1 .and. i.eq.np) then
!added d rrl
                  ip1 = np + 2
               else
                  ip1 = i + 1
               end if
               if (leftdedge.eq.1 .and. i.eq.1) then
!added d rrl
                  im1 = -1
               else
                  im1 = i - 1
               end if
c        ip1=i+1-i/n*(n-1)
c        im1=i-1+(n-i)/(n-1)*(n-1)
               g13=c13(i,j)*gmul(1)
               g23=c23(i,j)*gmul(1)
               g33=g13**2+g23**2+gi(i,j,1)**2
               px=hdxi*(rk(ip1,j,1)-rk(im1,j,1))
 2121      pz(i,j,1)=-2.*(srff(i,j)+g13*px)/(g33-g23*g23)-pz(i,j,2)
c  in this case what is the value of py
 212    continue
      endif
c the following lines use rightedge and topedge and need examination.
c They might be ok.  rrl
      if (rightdedge.eq.0 .and. topdedge.eq.0) then
!added d rrl
         call updated(pz,pz,np,mp,l+1,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      else if (rightdedge.eq.1 .and. topdedge.eq.0) then
!added d rrl
         call updated(pz,pz,np+1,mp,l+1,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      else if (rightdedge.eq.0 .and. topdedge.eq.1) then
!added d rrl
         call updated(pz,pz,np,mp+1,l+1,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      else
         call updated(pz,pz,np+1,mp+1,l+1,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      end if
 
compute x-flux at (i+-1/2,j,k)
      do j=1,mp
        do i=1+leftdedge,np                 !added d rrl
          do k=1,L 
            g13=0.5*gmul(k)*(c13(i-1,j)+c13(i,j))
            pza=0.25*(pz(i-1,j,k)+pz(i,j,k)+pz(i-1,j,k+1)+pz(i,j,k+1))
            Pxa=dxi*(rk(i,j,k)-rk(i-1,j,k))
            hx(i,j,k)=( pxa + g13*pza )
          end do
        end do
      end do
create boundary conditions at i=1
      call updated(hx,hx,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      if(leftdedge.eq.1) then                 !added d rrl
      do j=1,mp
        do k=1,L
          hx(1,j,k) = (ibcx-1)*hx(2,j,k) + ibcx*hx(0,j,k)
        end do
      end do
      endif
      call updated(hx,hx,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
c
compute y-flux at (i,j+-1/2,k)
      if (j3.eq.1) then
        do i=1,np
          do j=1+botdedge,mp                 !added d rrl
            do k=1,L
              g23=0.5*gmul(k)*(c23(i,j-j3)+c23(i,j))
              pza=0.25*(  pz(i,j-j3,k  ) + pz(i,j,k  )
     .                  + pz(i,j-j3,k+1) + pz(i,j,k+1) )
              pya=dyi*(rk(i,j,k)-rk(i,j-j3,k))
              hy(i,j,k)=( pya + g23*pza )
            end do
          end do
        end do
create boundary conditions at j=1
        call updated(hy,hy,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
        if(botdedge.eq.1) then                 !added d rrl
        do k=1,l
          do i=1,np
            hy(i,1,k)= (ibcy-1)*hy(i,2,k) + ibcy*hy(i,0,k)
          end do
        end do
        endif
      endif
        call updated(hy,hy,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
 
compute z-flux at (i,j,k+-1/2)
! i) include the dh/dx and dh/dz terms
      do k=2,L
        do j=1,mp
          do i=1,np-1*rightdedge                 !added d rrl
            g13=0.5*(gmul(k)+gmul(k-1))*c13(i,j)
            gii=0.5*(gi(i,j,k)+gi(i,j,k-1))
            hxa=0.25*(hx(i,j,k-1)+hx(i+1,j,k-1)+hx(i,j,k)+hx(i+1,j,k))
            hza=gii*pz(i,j,k)
            hz(i,j,k)= gii*hza + g13*hxa
          end do
          if (rightdedge.eq.1) then                 !added d rrl
corporate b.c. for hx on i=n+1
             g13=0.5*(gmul(k)+gmul(k-1))*c13(np,j)
             gii=0.5*(gi(np,j,k)+gi(np,j,k-1))
             hx1 = (ibcx-1)*hx(np,j,k-1) + ibcx*hx(np+2,j,k-1)
             hx2 = (ibcx-1)*hx(np,j,k) + ibcx*hx(np+2,j,k)
             hxa=0.25*(hx(np,j,k-1)+hx1+hx(np,j,k)+hx2)
             hza=gii*pz(np,j,k)
             hz(np,j,k)= gii*hza + g13*hxa
          end if
        end do
      end do
! ii) include the dh/dy term if 3D
      if(j3.eq.1) then
      call updated(hy,hy,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
        do k=2,L
          do i=1,np
            do j=1,mp-topdedge                 !added d rrl
              g23=0.5*(gmul(k)+gmul(k-1))*c23(i,j)
              hya=0.25*(hy(i,j,k-1)+hy(i,j+j3,k-1)+
     1                  hy(i,j,k)+hy(i,j+j3,k))
              hz(i,j,k)=hz(i,j,k) + g23*hya
            end do
          end do
        end do
      if (topdedge.eq.1) then                 !added d rrl
        do k=2,L
          do i=1,np
corporate b.c. for hy on j=m+1
            g23=0.5*(gmul(k)+gmul(k-1))*c23(i,mp)
            hy1= (ibcy-1)*hy(i,mp,k-1) + ibcy*hy(i,mp+2,k-1)
            hy2= (ibcy-1)*hy(i,mp,k) + ibcy*hy(i,mp+2,k)
            hya=0.25*(hy(i,mp,k-1)+hy1+hy(i,mp,k)+hy2)
            hz(i,mp,k)=hz(i,mp,k) + g23*hya
          end do
        end do
      end if
      endif

create boundary conditions at k=1; for k=L see divergence below
c surface fluxes:
      do j=1,mp
      do i=1,np
      hz(i,j,1)=-hz(i,j,2)
      end do
      end do
      call updated(hz,hz,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)

compute u
      do k=1,l
      do j=1,mp
      do i=1+leftdedge,np                 !added d rrl
      saa=0.5*(s(i-1,j,k)+s(i,j,k))
      sqrtka=0.5*(sqrtk(i-1,j,k)+sqrtk(i,j,k))
      coef=0.
      if(sqrtka.ne.0.) coef=saa/sqrtka
      coefa=coef*dt*(dxi**2+j3*dyi**2+dzi**2)
      d12a=0.
      if(j3.eq.1) d12a=0.5*(d12(i,j,k)+d12(i,j+1,k))
      d13a=0.5*(d13(i,j,k)+d13(i,j,k+1))
      if(j3.eq.1) hyp1=hy(i,j+1,k)
      if(j3.eq.1) hyp2=hy(i-1,j+1,k)
      if((topdedge.eq.1.and.j.eq.mp).and.j3.eq.1)
!added d rrl
     .               hyp1=hy(i,j,k)*(ibcy-1)
     .                   +hy(i,mp+2,k)*ibcy
      if((topdedge.eq.1.and.j.eq.mp).and.j3.eq.1)
!added d rrl
     .               hyp2=hy(i-1,j,k)*(ibcy-1)
     .                   +hy(i-1,mp+2,k)*ibcy
      hya=0.
      if(j3.eq.1) hya=0.25*(hy(i,j,k)+hyp1+hy(i-1,j,k)+hyp2)
      if(k.lt.l) hzp1=hz(i,j,k+1)
      if(k.lt.l) hzp2=hz(i-1,j,k+1)
      if(k.eq.l) hzp1=-hz(i,j,k)
      if(k.eq.l) hzp2=-hz(i-1,j,k)
      hza=0.25*(hz(i,j,k)+hzp1+hz(i-1,j,k)+hzp2)
c     if(abs(d11(i,j,k)*hx(i,j,k))*coefa.gt.0.5) print*,'error'
c     if(abs(d12a*hy(i,j,k))*coefa.gt.0.5) print*,'error'
c     if(abs(d13a*hz(i,j,k))*coefa.gt.0.5) print*,'error'
      u(i,j,k)=coef*(d11(i,j,k)*hx(i,j,k)+d12a*hya+d13a*hza)
      enddo
      enddo
      enddo
      call updated(u,u,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)

      if(leftdedge.eq.1) then                 !added d rrl
      do j=1,mp
        do k=1,L
          u(1,j,k) = (ibcx-1)*u(2,j,k) + ibcx*u(0,j,k)
        end do
      end do
      endif

compute v
      if(j3.eq.1) then
      do k=1,l
      do j=1+botdedge,mp                 !added d rrl
      do i=1,np
      saa=0.5*(s(i,j-1,k)+s(i,j,k))
      sqrtka=0.5*(sqrtk(i,j-1,k)+sqrtk(i,j,k))
      coef=0.
      if(sqrtka.ne.0.) coef=saa/sqrtka
      d21a=0.5*(d12(i,j,k)+d12(i+1,j,k))
      d23a=0.5*(d23(i,j,k)+d23(i,j,k+1))
      hxp1=hx(i+1,j,k)
      hxp2=hx(i+1,j-1,k)
      if(rightdedge.eq.1.and.i.eq.np) hxp1=hx(i,j,k)*(ibcx-1)
!added d rrl
     .                                   +hx(np+2,j,k)*ibcx
      if(rightdedge.eq.1.and.i.eq.np) hxp2=hx(i,j-1,k)*(ibcx-1)
!added d rrl
     .                                   +hx(np+2,j-1,k)*ibcx
      hxa=0.25*(hx(i,j,k)+hxp1+hx(i,j-1,k)+hxp2)
      if(k.ne.l) hzp1=hz(i,j,k+1)
      if(k.ne.l) hzp2=hz(i,j-1,k+1)
      if(k.eq.l) hzp1=-hz(i,j,k)
      if(k.eq.l) hzp2=-hz(i,j-1,k)
      hza=0.25*(hz(i,j,k)+hzp1+hz(i,j-1,k)+hzp2)
      v(i,j,k)=coef*(d21a*hxa+d22(i,j,k)*hy(i,j,k)+d23a*hza)
      enddo
      enddo
      enddo
      call updated(v,v,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      if(botdedge.eq.1) then                 !added d rrl
      do k=1,l
      do i=1,np
      v(i,1,k)= (ibcy-1)*v(i,2,k) + ibcy*v(i,0,k)
      end do
      end do
      endif
      endif
compute w
      do k=2,l
      do j=1,mp
      do i=1,np
      saa=0.5*(s(i,j,k)+s(i,j,k-1))
      sqrtka=0.5*(sqrtk(i,j,k)+sqrtk(i,j,k-1))
      if(sqrtka.eq.0.) coef=0.
      if(sqrtka.ne.0.) coef=saa/sqrtka
      d31a=0.5*(d13(i,j,k)+d13(i+1,j,k))
      d32a=0. 
      if(j3.eq.1) d32a=0.5*(d23(i,j,k)+d23(i,j+1,k))
      coefa=coef*dt*(dxi**2+j3*dyi**2+dzi**2)
      hxp1=hx(i+1,j,k)
      hxp2=hx(i+1,j,k-1)
      if(rightdedge.eq.1.and.i.eq.np)                 !added d rrl 
     .                    hxp1=hx(i,j,k)*(ibcx-1)
     .                        +hx(np+2,j,k)*ibcx
      if(rightdedge.eq.1.and.i.eq.np)                  !added d rrl
     .                    hxp2=hx(i,j,k-1)*(ibcx-1)
     .                        +hx(np+2,j,k-1)*ibcx
      hxa=0.25*(hx(i,j,k)+hx(i,j,k-1)+hxp1+hxp2)
      if(j3.eq.1) hyp1=hy(i,j+1,k)
      if(j3.eq.1) hyp2=hy(i,j+1,k-1)
      if((topdedge.eq.1.and.j.eq.mp).and.j3.eq.1)
!added d rrl
     .                    hyp1=hy(i,j,k)*(ibcy-1)
     .                        +hy(i,mp+2,k)*ibcy
      if((topdedge.eq.1.and.j.eq.mp).and.j3.eq.1)
!added d rrl
     .                    hyp2=hy(i,j,k-1)*(ibcy-1)
     .                        +hy(i,mp+2,k-1)*ibcy
      hya=0.
      if(j3.eq.1) hya=0.25*(hy(i,j,k)+hyp1+hy(i,j,k-1)+hyp2)
c     if(abs(d31a*hx(i,j,k))*coefa.gt.0.5) print*,'error'
c     if(abs(d32a*hy(i,j,k))*coefa.gt.0.5) print*,'error'
c     if(abs(d33(i,j,k)*hz(i,j,k))*coefa.gt.0.5) print*,'error'
      w(i,j,k)=coef*(d31a*hxa+d32a*hya+d33(i,j,k)*hz(i,j,k))
      enddo
      enddo
      enddo

      do j=1,mp
      do i=1,np
      w(i,j,1)=-w(i,j,2)
      enddo
      enddo
       
      do k=1,l
      do j=1,mp
      do i=1,np
c     w(i,j,k)=0.
      hx(i,j,k)=u(i,j,k) 
      hy(i,j,k)=v(i,j,k) 
      hz(i,j,k)=w(i,j,k) 
      enddo
      enddo
      enddo
      call updated(hx,hx,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      call updated(hy,hy,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      call updated(hz,hz,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)

compute Laplacian term by term
      do k=1,l
      do j=1,mp
      do i=1,np
      r(i,j,k)=0.
      end do
      end do
      end do
compute d/dx(dh/dx)
      do k=1,L
        do j=1,mp
          do i=1,np-rightdedge                 !added d rrl
            r(i,j,k) = dxi*(hx(i+1,j,k)-hx(i,j,k))
          end do
        end do
      end do
create boundary conditions on hx at i=n+1
      if (rightdedge.eq.1) then                 !added d rrl
         do k=1,L
            do j=1,mp
               hxn=hx(np,j,k)
               hxnp1 = (ibcx-1)*hxn + ibcx*hx(np+2,j,k)
               r(np,j,k) = dxi*(hxnp1-hxn)
            end do
         end do
      end if
c for 3D problem compute d/dy(dh/dy) and add to r
      if(j3.eq.1) then
        do k=1,L
          do i=1,np
            do j=1,mp-topdedge                 !added d rrl
              r(i,j,k) = r(i,j,k) + dyi*(hy(i,j+1,k)-hy(i,j,k))
            end do
          end do
        end do
       if (topdedge.eq.1) then                 !added d rrl
          do k=1,L
             do i=1,np
create boundary conditions at j=m+1
                hymp1= (ibcy-1)*hy(i,mp,k) + ibcy*hy(i,mp+2,k)
                r(i,mp,k) = r(i,mp,k) + dyi*(hymp1-hy(i,mp,k))
             end do
          end do
       end if
      endif
compute d/dz(dh/dz) and add to r
      do j=1,mp
        do i=1,np
          do k=1,L-1
            r(i,j,k) = r(i,j,k) + dzi*(hz(i,j,k+1)-hz(i,j,k))
          end do
corporate b.c. hz(i,j,L+1)=-hz(i,j,L) at k=L
          hzLp1 =-hz(i,j,L)
          r(i,j,L) = r(i,j,L) + dzi*(hzLp1-hz(i,j,L))
        end do
      end do
 
      do k=1,l
      do j=1,mp
      do i=1,np
      r(i,j,k)=0.09*2.*r(i,j,k)*gi(i,j,k)/3.
      fk(i,j,k)=fk(i,j,k)+2.*r(i,j,k)*dt
      end do
      end do
      end do
      call updated(fk,fk,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)

      deallocate (r)
      deallocate (hx)
      deallocate (hy)
      deallocate (hz)
      deallocate (pz)
      deallocate (srff)
      deallocate (u)
      deallocate (v)
      deallocate (w)
 
      return
      end subroutine lapdf
!lapdf*****************************!
      !subroutine lapdfs(rk,fk,s,
      !.d11,d22,d33,d12,d13,d23,il,iu,jl,ju,lls,iflg)
      !JAS 3/7/06 changed routine signature here 
      ! as "s" and "d-arrays" are not used!
      subroutine lapdfs(rk,fk,il,iu,jl,ju,lls,iflg)
      use gridsetup
      use turba
      use xvo
      use xve, only:xe
      use metryic
      use msga

      Implicit None

      !JAS 3/7/06 added explicit declarations to comply with implicit
      !none
      integer,intent(in) :: il,iu,jl,ju,lls

      real      rk(il:iu, jl:ju,lls),
      !.           s(il:iu, jl:ju,lls),
     .          fk(il:iu, jl:ju,lls)
      !real            d13(il:iu+1, jl:ju,lls+1), 
      !.                d12(il:iu+1, jl:ju+1,lls),
      !.                d23(il:iu, jl:ju+1,lls+1),
      !.                d11(il:iu+1, jl:ju,lls),
      !.                d22(il:iu, jl:ju+1,lls),
      !.                d33(il:iu, jl:ju,lls+1)
      real,allocatable::
     .           r(:, :,:),
     .          hx(:, :,:),
     .          hy(:, :,:),
     .          hz(:, :,:),
     .          pz(:, :,:),
     .        srff(:, :),
     .           u(:, :,:),
     .           v(:, :,:),
     .           w(:, :,:)!,
!     .         rkb(:,:,:)

      !JAS 3/7/06 added explicit declarations to comply with implicit
      !none
      integer :: i,j,k,illim,iulim,jllim,julim,ip1,im1,jp1,jm1,iflg
      real :: hdxi,hdyi,hdzi,g33,g23,g13,gii,px,py,pxa,pya,pza
      real :: hxa,hya,hza,hx1,hy1,hx2,hy2
      real :: hxn,hxnp1,hymp1
      real :: coefx,coefy,coefz,aswitch ! ,dcr1

      if(iturb.ge.1) then
      allocate (r(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (hx(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (hy(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (hz(1-ih:np+ih, 1-ih:mp+ih,l+1))
      allocate (pz(1-ih:np+ih+1, 1-ih:mp+ih+1,l+1))
      allocate (srff(1-ih:np+ih, 1-ih:mp+ih))
      allocate (u(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (v(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (w(1-ih:np+ih, 1-ih:mp+ih,l))
!      if (iturb.eq.2) then
!        allocate (rkb(1-ih:np+ih, 1-ih:mp+ih,l))
!        rkb(:,:,:)=xvb(:,:,:,6)/xvb(:,:,:,nv)
!      endif
      endif
 
c      aswitch=-1.*real(iflg-6)
      aswitch=0.
c      write (*,*) 'iflg=',iflg,' aswitch= ',aswitch
      hdxi=0.5*dxi
      hdyi=0.5*dyi
      hdzi=0.5*dzi
      call updated(rk,rk,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
c      write (*,*) 'inside lapdfs',xe(5,5,l,iflg),xe(5,5,l,nv)

compute z-derivatives at (i,j,k+-1/2)
      do k=2,L
        do j=1,mp
          do i=1,np
            pz(i,j,k)=dzi*(rk(i,j,k)-rk(i,j,k-1))
          end do
        end do
      end do
      do j=1,mp
        do i=1,np
          pz(i,j,L+1)=dzi*(xe(i,j,l,iflg)/xe(i,j,l,nv)-rk(i,j,l))
c      write (*,*) 'after pz top',i,j,dzi,xe(i,j,l,iflg),
c     &                xe(i,j,l,nv),
c     &                rk(i,j,l)
c23456
        end do
      end do

      do j=1,mp
      do i=1,np
        srff(i,j)=0.
      enddo
      enddo

      do j=1,mp
      do i=1,np
      g33=(c13(i,j)*gmul(1))**2+(c23(i,j)*gmul(1))**2+gi(i,j,1)**2
      srff(i,j) = sqrt(g33)*srff(i,j)
      enddo
      enddo
       
      jllim = 1  + (j3-ibcy)*botdedge                 !added d rrl
      julim = mp + (ibcy-j3)*topdedge                 !added d rrl
      do 21 j=jllim,julim
         if (j3.eq.1) then
            if (topdedge.eq.1 .and. j.eq.mp) then                 !added
               jp1 = mp + 2
            else
               jp1 = j + 1
            end if
            if (botdedge.eq.1 .and. j.eq.1) then                 !added
               jm1 = -1
            else
               jm1 = j - 1
            end if
         else
            jp1=1
            jm1=1
         end if
c        jp1=j+j3-j/m*(m-1)
c        jm1=j-j3+(m-j)/(m-j3)*(m-j3)
         illim = 1  + (1-ibcx)*leftdedge                 !added d rrl
         iulim = np + (ibcx-1)*rightdedge                 !added d rrl
         do 21 i=illim,iulim
            if (rightdedge.eq.1 .and. i.eq.np) then
!added d rrl
               ip1 = np + 2
            else
               ip1 = i + 1
            end if
            if (leftdedge.eq.1 .and. i.eq.1) then                 !added
               im1 = -1
            else
               im1 = i - 1
            end if
      g13=c13(i,j)*gmul(1)
      g23=c23(i,j)*gmul(1)
      g33=g13**2+g23**2+gi(i,j,1)**2
      px=hdxi*(rk(ip1,j,1)-rk(im1,j,1))
      py=hdyi*(rk(i,jp1,1)-rk(i,jm1,1))*j3
c   21 pz(i,j,1)=-2.*(srff(i,j)+g13*px+g23*py)/g33-pz(i,j,2)
   21 pz(i,j,1)=-1.*(srff(i,j)+g13*px+g23*py)/g33
     &          +aswitch*2.*dzi*rk(i,j,1)
c   21 pz(i,j,1)=-1.*(srff(i,j)+g13*px+g23*py)/g33
      if(ibcx.eq.0) then
         illim = 1*leftdedge   + np*(1-leftdedge)                 !added
         iulim = np*rightdedge +  1*(1-rightdedge)
!added d rrl
         do 211 i=illim,iulim,np-1
            jllim = 1  + (j3-ibcy)*botdedge                 !added d rrl
            julim = mp + (ibcy-j3)*topdedge                 !added d rrl
            do 2111 j=jllim,julim
               if (j3.eq.1) then
                  if (topdedge.eq.1 .and. j.eq.mp) then
!added d rrl
                     jp1 = mp + 2
                  else
                     jp1 = j + 1
                  end if
                  if (botdedge.eq.1 .and. j.eq.1) then
!added d rrl
                     jm1 = -1
                  else
                     jm1 = j - 1
                  end if
               else
                  jp1=1
                  jm1=1
               end if
c        jp1=j+j3-j/m*(m-1)
c        jm1=j-j3+(m-j)/(m-j3)*(m-j3)
               g13=c13(i,j)*gmul(1)
               g23=c23(i,j)*gmul(1)
               g33=g13**2+g23**2+gi(i,j,1)**2
               py=hdyi*(rk(i,jp1,1)-rk(i,jm1,1))*j3
c 2111  pz(i,j,1)=-2.*(srff(i,j)+g23*py)/(g33-g13*g13)-pz(i,j,2)
c234567
 2111  pz(i,j,1)=-1.*(srff(i,j)+g23*py)/(g33-g13*g13)
     &           +aswitch*2.*dzi*rk(i,j,1)
c 2111  pz(i,j,1)=-1.*(srff(i,j)+g23*py)/(g33-g13*g13)
           if(ibcy.eq.0.and.j3.eq.1) then
              jllim = 1*botdedge  + mp*(1-botdedge)
!added d rrl
              julim = mp*topdedge +  1*(1-topdedge)
!added d rrl
              do 2112 j=jllim,julim,mp-j3
c 2112         pz(i,j,1)=-2.*srff(i,j)/(gi(i,j,1)**2)-
c     .                  pz(i,j,2)
 2112         pz(i,j,1)=-1.*srff(i,j)/(gi(i,j,1)**2)
     &                  +aswitch*2.*dzi*rk(i,j,1)
c 2112         pz(i,j,1)=-1.*srff(i,j)/(gi(i,j,1)**2)
           endif
 211    continue
      endif
      if(ibcy.eq.0.and.j3.eq.1) then
         jllim = 1*botdedge  + mp*(1-botdedge)                 !added d
         julim = mp*topdedge +  1*(1-topdedge)                 !added d
         do 212 j=jllim,julim,mp-j3
            illim = 1  + (1-ibcx)*leftdedge                 !added d rrl
            iulim = np + (ibcx-1)*rightdedge                 !added d
            do 2121 i=illim,iulim
               if (rightdedge.eq.1 .and. i.eq.np) then
!added d rrl
                  ip1 = np + 2
               else
                  ip1 = i + 1
               end if
               if (leftdedge.eq.1 .and. i.eq.1) then
!added d rrl
                  im1 = -1
               else
                  im1 = i - 1
               end if
c        ip1=i+1-i/n*(n-1)
c        im1=i-1+(n-i)/(n-1)*(n-1)
               g13=c13(i,j)*gmul(1)
               g23=c23(i,j)*gmul(1)
               g33=g13**2+g23**2+gi(i,j,1)**2
               px=hdxi*(rk(ip1,j,1)-rk(im1,j,1))
c 2121      pz(i,j,1)=-2.*(srff(i,j)+g13*px)/(g33-g23*g23)-pz(i,j,2)
 2121      pz(i,j,1)=-1.*(srff(i,j)+g13*px)/(g33-g23*g23)
     &                  +aswitch*2.*dzi*rk(i,j,1)
c 2121      pz(i,j,1)=-1.*(srff(i,j)+g13*px)/(g33-g23*g23)
 212    continue
      endif

      if (rightdedge.eq.0 .and. topdedge.eq.0) then
!added d rrl
         call updated(pz,pz,np,mp,l+1,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      else if (rightdedge.eq.1 .and. topdedge.eq.0) then
!added d rrl
         call updated(pz,pz,np+1,mp,l+1,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      else if (rightdedge.eq.0 .and. topdedge.eq.1) then
!added d rrl
         call updated(pz,pz,np,mp+1,l+1,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      else
         call updated(pz,pz,np+1,mp+1,l+1,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      end if
 
compute x-flux at (i+-1/2,j,k)
      do j=1,mp
        do i=1+leftdedge,np
          do k=1,L 
c            coef=0.66*0.09*0.5*(sb(i,j,k)+sb(i-1,j,k))*   !rrl
c     .           sqrt(1.2*0.5*(rkb(i,j,k)+rkb(i-1,j,k)))  !rrl 
            coefx=0.09*0.5*(saxy(i,j,k)+saxy(i-1,j,k))*
     .           0.5*(sqrtk(i,j,k)+sqrtk(i-1,j,k))
     .      *0.5*(xvb(i-1,j,k,nv)+xvb(i,j,k,nv))
            coefz=0.09*0.5*(saz(i,j,k)+saz(i-1,j,k))*
     .           0.5*(sqrtk(i,j,k)+sqrtk(i-1,j,k))
     .      *0.5*(xvb(i-1,j,k,nv)+xvb(i,j,k,nv))
            !dcr1=3.*coef*dt*(dxi**2+j3*dyi**2+dzi**2)
c            if(abs(dcr1).gt.0.20)
c     .coef=0.20/(3.*dt*(dxi**2+j3*dyi**2+dzi**2))
            g13=0.5*gmul(k)*(c13(i-1,j)+c13(i,j))
            pza=0.25*(pz(i-1,j,k)+pz(i,j,k)+pz(i-1,j,k+1)+pz(i,j,k+1))
            Pxa=dxi*(rk(i,j,k)-rk(i-1,j,k))
            hx(i,j,k)=(coefx* pxa + coefz*g13*pza )
          end do
        end do
      end do
create boundary conditions at i=1
      call updated(hx,hx,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      if(leftdedge.eq.1) then
      do j=1,mp
        do k=1,L
          hx(1,j,k) = (ibcx-1)*hx(2,j,k) + ibcx*hx(0,j,k)
        end do
      end do
      endif
      call updated(hx,hx,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
c
compute y-flux at (i,j+-1/2,k)
      if (j3.eq.1) then
        do i=1,np
          do j=1+botdedge,mp                 !added d rrl
            do k=1,L
c           coef=0.66*0.09*0.5*(sb(i,j,k)+sb(i,j-1,k))*
c     .           sqrt(1.2*0.5*(rkb(i,j,k)+rkb(i,j-1,k)))
            coefy=0.09*0.5*(saxy(i,j,k)+saxy(i,j-1,k))*
     .           0.5*(sqrtk(i,j,k)+sqrtk(i,j-1,k))
     .     *0.5*(xvb(i,j-1,k,nv)+xvb(i,j,k,nv))
            coefz=0.09*0.5*(saz(i,j,k)+saz(i,j-1,k))*
     .           0.5*(sqrtk(i,j,k)+sqrtk(i,j-1,k))
     .     *0.5*(xvb(i,j-1,k,nv)+xvb(i,j,k,nv))
c            dcr1=3.*coef*dt*(dxi**2+j3*dyi**2+dzi**2)
c            if(abs(dcr1).gt.0.20)
c     .coef=0.20/(3.*dt*(dxi**2+j3*dyi**2+dzi**2))
              g23=0.5*gmul(k)*(c23(i,j-j3)+c23(i,j))
              pza=0.25*(  pz(i,j-j3,k  ) + pz(i,j,k  )
     .                  + pz(i,j-j3,k+1) + pz(i,j,k+1) )
              pya=dyi*(rk(i,j,k)-rk(i,j-j3,k))
              hy(i,j,k)=(coefy*pya + g23*pza*coefz )
            end do
          end do
        end do
create boundary conditions at j=1
        call updated(hy,hy,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
        if(botdedge.eq.1) then                 !added d rrl
        do k=1,l
          do i=1,np
            hy(i,1,k)= (ibcy-1)*hy(i,2,k) + ibcy*hy(i,0,k)
          end do
        end do
        endif
      endif
        call updated(hy,hy,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
 
compute z-flux at (i,j,k+-1/2)
! i) include the dh/dx and dh/dz terms
      do k=2,L
        do j=1,mp
          do i=1,np-1*rightdedge                 !added d rrl
c             coef=0.66*0.09*0.5*(sb(i,j,k)+sb(i,j,k-1))*
c     .           sqrt(1.2*0.5*(rkb(i,j,k)+rkb(i,j,k-1)))
            coefz=0.09*0.5*(saz(i,j,k)+saz(i,j,k-1))*
     .           0.5*(sqrtk(i,j,k)+sqrtk(i,j,k-1))
     .      *0.5*(xvb(i,j,k-1,nv)+xvb(i,j,k,nv))
c            dcr1=3.*coef*dt*(dxi**2+j3*dyi**2+dzi**2)
c            if(abs(dcr1).gt.0.20)
c     .coef=0.20/(3.*dt*(dxi**2+j3*dyi**2+dzi**2))
            g13=0.5*(gmul(k)+gmul(k-1))*c13(i,j)
            gii=0.5*(gi(i,j,k)+gi(i,j,k-1))
            hxa=0.25*(hx(i,j,k-1)+hx(i+1,j,k-1)+hx(i,j,k)+hx(i+1,j,k))
            hza=gii*pz(i,j,k)
            hz(i,j,k)= coefz*gii*hza + g13*hxa
          end do
          if (rightdedge.eq.1) then                 !added d rrl
corporate b.c. for hx on i=n+1
c             coef=0.66*0.09*0.5*(sb(np,j,k)+sb(np,j,k-1))*
c     .           sqrt(1.2*0.5*(rkb(np,j,k)+rkb(np,j,k-1)))
            coefz=0.09*0.5*(saz(np,j,k)+saz(np,j,k-1))*
     .           0.5*(sqrtk(np,j,k)+sqrtk(np,j,k-1))
     .      *0.5*(xvb(np,j,k-1,nv)+xvb(np,j,k,nv))
c            dcr1=3.*coef*dt*(dxi**2+j3*dyi**2+dzi**2)
c            if(abs(dcr1).gt.0.20)
c     .coef=0.20/(3.*dt*(dxi**2+j3*dyi**2+dzi**2))
             g13=0.5*(gmul(k)+gmul(k-1))*c13(np,j)
             gii=0.5*(gi(np,j,k)+gi(np,j,k-1))
             hx1 = (ibcx-1)*hx(np,j,k-1) + ibcx*hx(np+2,j,k-1)
             hx2 = (ibcx-1)*hx(np,j,k) + ibcx*hx(np+2,j,k)
             hxa=0.25*(hx(np,j,k-1)+hx1+hx(np,j,k)+hx2)
             hza=gii*pz(np,j,k)
             hz(np,j,k)= coefz*gii*hza + g13*hxa
          end if
        end do
      end do
      k=1
        do j=1,mp
          do i=1,np-1*rightdedge                 !added d rrl
c             coef=0.66*0.09*0.5*(sb(i,j,k)+sb(i,j,k-1))*
c     .           sqrt(1.2*0.5*(rkb(i,j,k)+rkb(i,j,k-1)))
            coefz=0.09*saz(i,j,k)*
     .           sqrtk(i,j,k)
     .          *xvb(i,j,k,nv)
c            dcr1=3.*coef*dt*(dxi**2+j3*dyi**2+dzi**2)
c            if(abs(dcr1).gt.0.20)
c     .coef=0.20/(3.*dt*(dxi**2+j3*dyi**2+dzi**2))
            g13=gmul(k)*c13(i,j)
            gii=gi(i,j,k)
            hxa=0.5*(hx(i,j,k)+hx(i+1,j,k))
            hza=gii*pz(i,j,k)
            hz(i,j,k)= coefz*gii*hza + g13*hxa
          end do
          if (rightdedge.eq.1) then                 !added d rrl
corporate b.c. for hx on i=n+1
c             coef=0.66*0.09*.5*(sb(np,j,k)+sb(np,j,k-1))*
c     .           sqrt(1.2*0.5*(rkb(np,j,k)+rkb(np,j,k-1)))
            coefz=0.09*saz(np,j,k)*
     .          sqrtk(np,j,k)
     .         *xvb(np,j,k,nv)
c            dcr1=3.*coef*dt*(dxi**2+j3*dyi**2+dzi**2)
c            if(abs(dcr1).gt.0.20)
c     .coef=0.20/(3.*dt*(dxi**2+j3*dyi**2+dzi**2))
             g13=gmul(k)*c13(np,j)
             gii=gi(np,j,k)
c             hx1 = (ibcx-1)*hx(np,j,k-1) + ibcx*hx(np+2,j,k-1)
             hx2 = (ibcx-1)*hx(np,j,k) + ibcx*hx(np+2,j,k)
             hxa=0.5*(hx(np,j,k)+hx2)
             hza=gii*pz(np,j,k)
             hz(np,j,k)= coefz*gii*hza + g13*hxa
          end if
        end do
        k=l+1
        do j=1,mp
          do i=1,np-1*rightdedge                 !added d rrl
c             coef=0.66*0.09*sb(i,j,k-1)*
c     .           sqrt(1.2*rkb(i,j,k-1))
            coefz=0.09*saz(i,j,k-1)*
     .           sqrtk(i,j,k-1)
     .      *xvb(i,j,k-1,nv)
c            dcr1=3.*coef*dt*(dxi**2+j3*dyi**2+dzi**2)
c            if(abs(dcr1).gt.0.20)
c     .coef=0.20/(3.*dt*(dxi**2+j3*dyi**2+dzi**2))
            g13=gmul(k-1)*c13(i,j)
            gii=gi(i,j,k-1)
            hxa=0.5*(hx(i,j,k-1)+hx(i+1,j,k-1))
            hza=gii*pz(i,j,k)
            hz(i,j,k)= coefz*gii*hza + g13*hxa
          end do
          if (rightdedge.eq.1) then                 !added d rrl
corporate b.c. for hx on i=n+1
c             coef=0.66*0.09*0.5*(sb(np,j,k)+sb(np,j,k-1))*
c     .           sqrt(1.2*0.5*(rkb(np,j,k)+rkb(np,j,k-1)))
            coefz=0.09*saz(np,j,k-1)*
     .           sqrtk(np,j,k-1)
     .      *xvb(np,j,k-1,nv)
c            dcr1=3.*coef*dt*(dxi**2+j3*dyi**2+dzi**2)
c            if(abs(dcr1).gt.0.20)
c     .coef=0.20/(3.*dt*(dxi**2+j3*dyi**2+dzi**2))
             g13=gmul(k-1)*c13(np,j)
             gii=gi(np,j,k-1)
             hx1 = (ibcx-1)*hx(np,j,k-1) + ibcx*hx(np+2,j,k-1)
c             hx2 = (ibcx-1)*hx(np,j,k) + ibcx*hx(np+2,j,k)
             hxa=0.5*(hx(np,j,k-1)+hx1)
             hza=gii*pz(np,j,k)
             hz(np,j,k)= coefz*gii*hza + g13*hxa
          end if
        end do

! ii) include the dh/dy term if 3D
      if(j3.eq.1) then
      call updated(hy,hy,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
        do k=2,L
          do i=1,np
            do j=1,mp-topdedge                 !added d rrl
              g23=0.5*(gmul(k)+gmul(k-1))*c23(i,j)
              hya=0.25*(hy(i,j,k-1)+hy(i,j+j3,k-1)+
     1                  hy(i,j,k)+hy(i,j+j3,k))
              hz(i,j,k)=hz(i,j,k) + g23*hya
            end do
          end do
        end do
      if (topdedge.eq.1) then                 !added d rrl
        do k=2,L
          do i=1,np
corporate b.c. for hy on j=m+1
            g23=0.5*(gmul(k)+gmul(k-1))*c23(i,mp)
            hy1= (ibcy-1)*hy(i,mp,k-1) + ibcy*hy(i,mp+2,k-1)
            hy2= (ibcy-1)*hy(i,mp,k) + ibcy*hy(i,mp+2,k)
            hya=0.25*(hy(i,mp,k-1)+hy1+hy(i,mp,k)+hy2)
            hz(i,mp,k)=hz(i,mp,k) + g23*hya
          end do
        end do
      end if
        k=1
          do i=1,np
            do j=1,mp-topdedge                 !added d rrl
              g23=gmul(k)*c23(i,j)
              hya=0.5*(
     1                  hy(i,j,k)+hy(i,j+j3,k))
              hz(i,j,k)=hz(i,j,k) + g23*hya
            end do
        end do
      if (topdedge.eq.1) then                 !added d rrl
          do i=1,np
corporate b.c. for hy on j=m+1
            g23=gmul(k)*c23(i,mp)
c            hy1= (ibcy-1)*hy(i,mp,k-1) + ibcy*hy(i,mp+2,k-1)
            hy2= (ibcy-1)*hy(i,mp,k) + ibcy*hy(i,mp+2,k)
            hya=0.5*(hy(i,mp,k)+hy2)
            hz(i,mp,k)=hz(i,mp,k) + g23*hya
          end do
      end if
       k=l+1
          do i=1,np
            do j=1,mp-topdedge                 !added d rrl
              g23=gmul(k-1)*c23(i,j)
              hya=0.5*(hy(i,j,k-1)+hy(i,j+j3,k-1))
              hz(i,j,k)=hz(i,j,k) + g23*hya
            end do
          end do
      if (topdedge.eq.1) then                 !added d rrl
          do i=1,np
corporate b.c. for hy on j=m+1
            g23=gmul(k-1)*c23(i,mp)
            hy1= (ibcy-1)*hy(i,mp,k-1) + ibcy*hy(i,mp+2,k-1)
c            hy2= (ibcy-1)*hy(i,mp,k) + ibcy*hy(i,mp+2,k)
            hya=0.5*(hy(i,mp,k-1)+hy1)
            hz(i,mp,k)=hz(i,mp,k) + g23*hya
          end do
      end if

      endif

create boundary conditions at k=1; for k=L see divergence below
c surface fluxes:
c     do j=1,mp
c     do i=1,np
c     hz(i,j,1)=-hz(i,j,2)
c     end do
c     end do
      call updated(hz,hz,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)

compute Laplacian term by term
      do k=1,l
      do j=1,mp
      do i=1,np
      r(i,j,k)=0.
      end do
      end do
      end do
compute d/dx(dh/dx)
      do k=1,L
        do j=1,mp
          do i=1,np-rightdedge                 !added d rrl
            r(i,j,k) = dxi*(hx(i+1,j,k)-hx(i,j,k))
          end do
        end do
      end do
create boundary conditions on hx at i=n+1
      if (rightdedge.eq.1) then                 !added d rrl
         do k=1,L
            do j=1,mp
               hxn=hx(np,j,k)
               hxnp1 = (ibcx-1)*hxn + ibcx*hx(np+2,j,k)
               r(np,j,k) = dxi*(hxnp1-hxn)
            end do
         end do
      end if
c for 3D problem compute d/dy(dh/dy) and add to r
      if(j3.eq.1) then
        do k=1,L
          do i=1,np
            do j=1,mp-topdedge                 !added d rrl
              r(i,j,k) = r(i,j,k) + dyi*(hy(i,j+1,k)-hy(i,j,k))
            end do
          end do
        end do
       if (topdedge.eq.1) then                 !added d rrl
          do k=1,L
             do i=1,np
create boundary conditions at j=m+1
                hymp1= (ibcy-1)*hy(i,mp,k) + ibcy*hy(i,mp+2,k)
                r(i,mp,k) = r(i,mp,k) + dyi*(hymp1-hy(i,mp,k))
             end do
          end do
       end if
      endif
compute d/dz(dh/dz) and add to r
      do j=1,mp
        do i=1,np
          do k=1,L
            r(i,j,k) = r(i,j,k) + dzi*(hz(i,j,k+1)-hz(i,j,k))
          end do
corporate b.c. hz(i,j,L+1)=-hz(i,j,L) at k=L
c         hzLp1 =-hz(i,j,L)
c         r(i,j,L) = r(i,j,L) + dzi*(hzLp1-hz(i,j,L))
        end do
      end do
 
      do k=1,l
      do j=1,mp
      do i=1,np
      r(i,j,k)=r(i,j,k)*gi(i,j,k)
      fk(i,j,k)=fk(i,j,k)+2.*r(i,j,k)*dt
      end do
      end do
      end do
      call updated(fk,fk,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)

      deallocate (r)
      deallocate (hx)
      deallocate (hy)
      deallocate (hz)
      deallocate (pz)
      deallocate (srff)
      deallocate (u)
      deallocate (v)
      deallocate (w)
      !deallocate (rkb)
 
      return
      end subroutine lapdfs

***********************************************************************
      subroutine lapdi(xv,fi,s,d11,d22,d33,d12,d13,d23,il,iu,jl,ju,lls,nvp)
      use gridsetup
      use metryic
      use turba
      use msga

      Implicit None

      !JAS 3/7/06 added explicit declarations to comply with implicit
      !none
      integer,intent(in) :: il,iu,jl,ju,lls,nvp

      real      xv(il:iu, jl:ju,lls,nvp),
     .          s(il:iu, jl:ju,lls),
     .          fi(il:iu, jl:ju,lls)
      real            d13(il:iu+1, jl:ju,lls+1),
     .                d12(il:iu+1, jl:ju+1,lls),
     .                d23(il:iu, jl:ju+1,lls+1),
     .                d11(il:iu+1, jl:ju,lls),
     .                d22(il:iu, jl:ju+1,lls),
     .                d33(il:iu, jl:ju,lls+1)
            real,allocatable::
     .           r(:, :,:),
     .          ri(:, :,:),
     .         rkb(:, :,:),
     .          hx(:, :,:),
     .          hy(:, :,:),
     .          hz(:, :,:),
     .          pz(:, :,:),
     .        srff(:, :),
     .           u(:, :,:),
     .           v(:, :,:),
     .           w(:, :,:)

      !JAS 3/7/06 added explicit declarations to comply with implicit
      !none
      integer :: i,j,k,illim,iulim,jllim,julim,ip1,im1,jp1,jm1,kmz
      real :: hdxi,hdyi,hdzi,g33,g23,g13,gii,px,py,pxa,pya,pza
      real :: rkba
      real :: hxa,hya,hza,hx1,hy1,hx2,hy2,hxp1,hxp2,hyp1,hyp2,hzp1,hzp2
      real :: hxn,hxnp1,hymp1,hzLp1
      real :: coef,d12a,d13a,d21a,d23a,d31a,d32a

      if(iturb.ge.1) then
      allocate (ri(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (r(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (rkb(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (hx(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (hy(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (hz(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (pz(1-ih:np+ih+1, 1-ih:mp+ih+1,l+1))
      allocate (srff(1-ih:np+ih, 1-ih:mp+ih))
      allocate (u(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (v(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (w(1-ih:np+ih, 1-ih:mp+ih,l))
      endif


 
      hdxi=0.5*dxi
      hdyi=0.5*dyi
      hdzi=0.5*dzi

      do k=1,l 
      do j=1,mp 
      do i=1,np 
      ri(i,j,k)=xv(i,j,k,4)/xv(i,j,k,nv)
      rkb(i,j,k)=xv(i,j,k,5)/xv(i,j,k,nv)
      enddo
      enddo
      enddo
      call updated(ri,ri,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      call updated(rkb,rkb,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)

compute z-derivatives at (i,j,k+-1/2)
      do k=2,L
        do j=1,mp
          do i=1,np
            pz(i,j,k)=dzi*(ri(i,j,k)-ri(i,j,k-1))
          end do
        end do
      end do
      do j=1,mp
        do i=1,np
          pz(i,j,L+1)=-pz(i,j,L)
        end do
      end do

      do j=1,mp
      do i=1,np
        srff(i,j)=0.
      enddo
      enddo

      do j=1,mp
      do i=1,np
      g33=(c13(i,j)*gmul(1))**2+(c23(i,j)*gmul(1))**2+gi(i,j,1)**2
      srff(i,j) = sqrt(g33)*srff(i,j)
      enddo
      enddo
      jllim = 1  + (j3-ibcy)*botdedge                       !add d rrl
      julim = mp + (ibcy-j3)*topdedge                       !add d rrl
      do 21 j=jllim,julim
         if (j3.eq.1) then
            if (topdedge.eq.1 .and. j.eq.mp) then
!add d rrl
               jp1 = mp + 2
            else
               jp1 = j + 1
            end if
            if (botdedge.eq.1 .and. j.eq.1) then
!add d rrl
               jm1 = -1
            else
               jm1 = j - 1
            end if
         else
            jp1=1
            jm1=1
         end if
c        jp1=j+j3-j/m*(m-1)
c        jm1=j-j3+(m-j)/(m-j3)*(m-j3)
         illim = 1  + (1-ibcx)*leftdedge                       !add d
         iulim = np + (ibcx-1)*rightdedge                       !add d
         do 21 i=illim,iulim
            if (rightdedge.eq.1 .and. i.eq.np) then
!add d rrl
               ip1 = np + 2
            else
               ip1 = i + 1
            end if
            if (leftdedge.eq.1 .and. i.eq.1) then
!add d rrl
               im1 = -1
            else
               im1 = i - 1
            end if
      g13=c13(i,j)*gmul(1)
      g23=c23(i,j)*gmul(1)
      g33=g13**2+g23**2+gi(i,j,1)**2
      px=hdxi*(ri(ip1,j,1)-ri(im1,j,1))
      py=hdyi*(ri(i,jp1,1)-ri(i,jm1,1))*j3
   21 pz(i,j,1)=-2.*(srff(i,j)+g13*px+g23*py)/g33-pz(i,j,2)
      if(ibcx.eq.0) then
         illim = 1*leftdedge   + np*(1-leftdedge)
!add d rrl
         iulim = np*rightdedge +  1*(1-rightdedge)
!add d rrl
         do 211 i=illim,iulim,np-1
            jllim = 1  + (j3-ibcy)*botdedge                       !add d
            julim = mp + (ibcy-j3)*topdedge                       !add d
            do 2111 j=jllim,julim
               if (j3.eq.1) then
                  if (topdedge.eq.1 .and. j.eq.mp) then
!add d rrl
                     jp1 = mp + 2
                  else
                     jp1 = j + 1
                  end if
                  if (botdedge.eq.1 .and. j.eq.1) then
!add d rrl
                     jm1 = -1
                  else
                     jm1 = j - 1
                  end if
               else
                  jp1=1
                  jm1=1
               end if
c        jp1=j+j3-j/m*(m-1)
c        jm1=j-j3+(m-j)/(m-j3)*(m-j3)
               g13=c13(i,j)*gmul(1)
               g23=c23(i,j)*gmul(1)
               g33=g13**2+g23**2+gi(i,j,1)**2
               py=hdyi*(ri(i,jp1,1)-ri(i,jm1,1))*j3
 2111      pz(i,j,1)=-2.*(srff(i,j)+g23*py)/(g33-g13*g13)-pz(i,j,2)
           if(ibcy.eq.0.and.j3.eq.1) then
              jllim = 1*botdedge  + mp*(1-botdedge)
!add d rrl
              julim = mp*topdedge +  1*(1-topdedge)
!add d rrl
              do 2112 j=jllim,julim,mp-j3
 2112         pz(i,j,1)=-2.*srff(i,j)/(gi(i,j,1)**2)-
     .                  pz(i,j,2)
           endif
 211    continue
      endif
      if(ibcy.eq.0.and.j3.eq.1) then
         jllim = 1*botdedge  + mp*(1-botdedge)
!add d rrl
         julim = mp*topdedge +  1*(1-topdedge)
!add d rrl
         do 212 j=jllim,julim,mp-j3
            illim = 1  + (1-ibcx)*leftdedge                       !add d

            iulim = np + (ibcx-1)*rightdedge                       !add
            do 2121 i=illim,iulim
               if (rightdedge.eq.1 .and. i.eq.np) then
!add d rrl
                  ip1 = np + 2
               else
                  ip1 = i + 1
               end if
               if (leftdedge.eq.1 .and. i.eq.1) then
!add d rrl
                  im1 = -1
               else
                  im1 = i - 1
               end if
c        ip1=i+1-i/n*(n-1)
c        im1=i-1+(n-i)/(n-1)*(n-1)
               g13=c13(i,j)*gmul(1)
               g23=c23(i,j)*gmul(1)
               g33=g13**2+g23**2+gi(i,j,1)**2
               px=hdxi*(ri(ip1,j,1)-ri(im1,j,1))
 2121      pz(i,j,1)=-2.*(srff(i,j)+g13*px)/(g33-g23*g23)-pz(i,j,2)
 212    continue
      endif

      if (rightdedge.eq.0 .and. topdedge.eq.0) then
!add d rrl
         call updated(pz,pz,np,mp,l+1,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      else if (rightdedge.eq.1 .and. topdedge.eq.0) then
!add d rrl
         call updated(pz,pz,np+1,mp,l+1,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      else if (rightdedge.eq.0 .and. topdedge.eq.1) then
!add d rrl
         call updated(pz,pz,np,mp+1,l+1,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      else
         call updated(pz,pz,np+1,mp+1,l+1,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      end if

       
compute x-flux at (i+-1/2,j,k)
      do j=1,mp
        do i=1+leftdedge,np                       !add d rrl
          do k=1,L
            g13=0.5*gmul(k)*(c13(i-1,j)+c13(i,j))
            pza=0.25*(pz(i-1,j,k)+pz(i,j,k)+pz(i-1,j,k+1)+pz(i,j,k+1))
            Pxa=dxi*(ri(i,j,k)-ri(i-1,j,k))
            hx(i,j,k)=( pxa + g13*pza )
          end do
        end do
      end do
create boundary conditions at i=1
      call updated(hx,hx,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      if(leftdedge.eq.1) then                       !add d rrl
      do j=1,mp
        do k=1,L
          hx(1,j,k) = (ibcx-1)*hx(2,j,k) + ibcx*hx(0,j,k)

        end do
      end do
      endif
      call updated(hx,hx,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
c
compute y-flux at (i,j+-1/2,k)
      if (j3.eq.1) then
        do i=1,np
          do j=1+botdedge,mp                       !add d rrl
            do k=1,L
              g23=0.5*gmul(k)*(c23(i,j-j3)+c23(i,j))
              pza=0.25*(  pz(i,j-j3,k  ) + pz(i,j,k  )
     .                  + pz(i,j-j3,k+1) + pz(i,j,k+1) )
              pya=dyi*(ri(i,j,k)-ri(i,j-j3,k))
              hy(i,j,k)=( pya + g23*pza )
            end do
          end do
        end do
create boundary conditions at j=1
        call updated(hy,hy,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
        if(botdedge.eq.1) then                       !add d rrl
        do k=1,l
          do i=1,np
            hy(i,1,k)= (ibcy-1)*hy(i,2,k) + ibcy*hy(i,0,k)
          end do
        end do
        endif
        call updated(hy,hy,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      endif

compute z-flux at (i,j,k+-1/2)
! i) include the dh/dx and dh/dz terms
      do k=2,L
        do j=1,mp
          do i=1,np-1*rightdedge                       !add d rrl
            g13=0.5*(gmul(k)+gmul(k-1))*c13(i,j)
            gii=0.5*(gi(i,j,k)+gi(i,j,k-1))
            hxa=0.25*(hx(i,j,k-1)+hx(i+1,j,k-1)+hx(i,j,k)+hx(i+1,j,k))
            hza=gii*pz(i,j,k)
            hz(i,j,k)= gii*hza + g13*hxa
          end do
          if (rightdedge.eq.1) then                       !add d rrl
corporate b.c. for hx on i=n+1
             g13=0.5*(gmul(k)+gmul(k-1))*c13(np,j)
             gii=0.5*(gi(np,j,k)+gi(np,j,k-1))
             hx1 = (ibcx-1)*hx(np,j,k-1) + ibcx*hx(np+2,j,k-1)
             hx2 = (ibcx-1)*hx(np,j,k) + ibcx*hx(np+2,j,k)
             hxa=0.25*(hx(np,j,k-1)+hx1+hx(np,j,k)+hx2)
             hza=gii*pz(np,j,k)
             hz(np,j,k)= gii*hza + g13*hxa
          end if
        end do
      end do
! ii) include the dh/dy term if 3D
      if(j3.eq.1) then
      call updated(hy,hy,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
        do k=2,L
          do i=1,np
            do j=1,mp-topdedge                       !add d rrl
              g23=0.5*(gmul(k)+gmul(k-1))*c23(i,j)
              hya=0.25*(hy(i,j,k-1)+hy(i,j+j3,k-1)+

     1                  hy(i,j,k)+hy(i,j+j3,k))
              hz(i,j,k)=hz(i,j,k) + g23*hya
            end do
          end do
        end do
      if (topdedge.eq.1) then                       !add d rrl
        do k=2,L
          do i=1,np
corporate b.c. for hy on j=m+1
            g23=0.5*(gmul(k)+gmul(k-1))*c23(i,mp)
            hy1= (ibcy-1)*hy(i,mp,k-1) + ibcy*hy(i,mp+2,k-1)
            hy2= (ibcy-1)*hy(i,mp,k) + ibcy*hy(i,mp+2,k)
            hya=0.25*(hy(i,mp,k-1)+hy1+hy(i,mp,k)+hy2)
            hz(i,mp,k)=hz(i,mp,k) + g23*hya
          end do
        end do
      end if
      endif
create boundary conditions at k=1; for k=L see divergence below
c surface fluxes:
      do j=1,mp
      do i=1,np
      hz(i,j,1)=-hz(i,j,2)
      end do
      end do
      call updated(hz,hz,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
  
compute u
      do k=1,l
      do j=1,mp
      do i=1+leftdedge,np                       !add d rrl
      rkba=0.5*(rkb(i,j,k)+rkb(i-1,j,k))
c     coef=0.09*0.5*(sc(i,j,k)+sc(i-1,j,k))/sqrt(0.2*rkba)
      coef=0.09*0.5*(s(i,j,k)+s(i-1,j,k))/sqrt(rkba)
      d12a=0.
      if(j3.eq.1) d12a=0.5*(d12(i,j,k)+d12(i,j+1,k))
      d13a=0.5*(d13(i,j,k)+d13(i,j,k+1))
      if(j3.eq.1) hyp1=hy(i,j+1,k)
      if(j3.eq.1) hyp2=hy(i-1,j+1,k)
      if((topdedge.eq.1.and.j.eq.mp).and.j3.eq.1)
!add d rrl
     .                        hyp1=hy(i,j,k)*(ibcy-1)
     .                            +hy(i,mp+2,k)*ibcy
      if((topdedge.eq.1.and.j.eq.mp).and.j3.eq.1)
!add d rrl
     .                        hyp2=hy(i-1,j,k)*(ibcy-1)
     .                            +hy(i-1,mp+2,k)*ibcy
      hya=0.
      if(j3.eq.1) hya=0.25*(hy(i,j,k)+hyp1+hy(i-1,j,k)+hyp2)
      if(k.ne.l) hzp1=hz(i,j,k+1)
      if(k.ne.l) hzp2=hz(i-1,j,k+1)
      if(k.eq.l) hzp1=-hz(i,j,k)
      if(k.eq.l) hzp2=-hz(i-1,j,k)
      hza=0.25*(hz(i,j,k)+hzp1+hz(i-1,j,k)+hzp2)
      u(i,j,k)=0.2*coef*(d11(i,j,k)*hx(i,j,k)+d12a*hya+
     .               d13a*hza)
      enddo
      enddo
      enddo
      call updated(u,u,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)

      if(leftdedge.eq.1) then                       !add d rrl
      do j=1,mp
        do k=1,L
          u(1,j,k) = (ibcx-1)*u(2,j,k) + ibcx*u(0,j,k)
        end do
      end do
      endif

compute v
      if(j3.eq.1) then
      do k=1,l
      do j=1+botdedge,mp                       !add d rrl
      do i=1,np
      rkba=0.5*(rkb(i,j,k)+rkb(i,j-1,k))
c     coef=0.09*0.5*(sc(i,j,k)+sc(i,j-1,k))/sqrt(0.2*rkba)
      coef=0.09*0.5*(s(i,j,k)+s(i,j-1,k))/sqrt(rkba)
      d21a=0.5*(d12(i,j,k)+d12(i+1,j,k))
      d23a=0.5*(d23(i,j,k)+d23(i,j,k+1))
      hxp1=hx(i+1,j,k)
      hxp2=hx(i+1,j-1,k)
      if(rightdedge.eq.1.and.i.eq.np) hxp1=hx(i,j,k)*(ibcx-1)
!add d rrl
     .                                   +hx(np+2,j,k)*ibcx
      if(rightdedge.eq.1.and.i.eq.np) hxp2=hx(i,j-1,k)*(ibcx-1)
!add d rrl
     .                                   +hx(np+2,j-1,k)*ibcx
      hxa=0.25*(hx(i,j,k)+hxp1+hx(i,j-1,k)+hxp2)
      if(k.ne.l) hzp1=hz(i,j,k+1)
      if(k.ne.l) hzp2=hz(i,j-1,k+1)
      if(k.eq.l) hzp1=-hz(i,j,k)
      if(k.eq.l) hzp2=-hz(i,j-1,k)
      hza=0.25*(hz(i,j,k)+hzp1+hz(i,j-1,k)+hzp2)
      v(i,j,k)=0.2*coef*(d21a*hxa+d22(i,j,k)*hy(i,j,k)+
     .                   d23a*hza)
      enddo
      enddo
      enddo
      call updated(v,v,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      if(botdedge.eq.1) then                       !add d rrl
      do k=1,l
      do i=1,np
      v(i,1,k)= (ibcy-1)*v(i,2,k) + ibcy*v(i,0,k)
      end do
      end do
      endif
      call updated(v,v,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      endif
compute w
      do k=2,l
      do j=1,mp
      do i=1,np
      kmz=k-1
      rkba=0.5*(rkb(i,j,k)+rkb(i,j,k-1))
c     coef=0.09*0.5*(sc(i,j,k)+sc(i,j,kmz))/sqrt(0.2*rkba)
      coef=0.09*0.5*(s(i,j,k)+s(i,j,kmz))/sqrt(rkba)
      d31a=0.5*(d13(i,j,k)+d13(i+1,j,k))
      d32a=0.
      if(j3.eq.1) d32a=0.5*(d23(i,j,k)+d23(i,j+1,k))
      hxp1=hx(i+1,j,k)
      hxp2=hx(i+1,j,k-1)
      if(rightdedge.eq.1.and.i.eq.np) hxp1=hx(i,j,k)*(ibcx-1)
!add d rrl
     .                                   +hx(np+2,j,k)*ibcx
      if(rightdedge.eq.1.and.i.eq.np) hxp2=hx(i,j,k-1)*(ibcx-1)
!add d rrl
     .                                   +hx(np+2,j,k-1)*ibcx
      hxa=0.25*(hx(i,j,k)+hx(i,j,k-1)+hxp1+hxp2)
      if(j3.eq.1) hyp1=hy(i,j+1,k)
      if(j3.eq.1) hyp2=hy(i,j+1,k-1)
      if((topdedge.eq.1.and.j.eq.mp).and.j3.eq.1)
!add d rrl
     .                             hyp1=hy(i,j,k)*(ibcy-1)
     .                                 +hy(i,mp+2,k)*ibcy
      if((topdedge.eq.1.and.j.eq.mp).and.j3.eq.1)
!add d rrl
     .                             hyp2=hy(i,j,k-1)*(ibcy-1)
     .                                 +hy(i,mp+2,k-1)*ibcy
      hya=0.
      if(j3.eq.1) hya=0.25*(hy(i,j,k)+hyp1+hy(i,j,k-1)+hyp2)
      w(i,j,k)=0.2*coef*(d31a*hxa+d32a*hya+
     .             d33(i,j,k)*hz(i,j,k))
      enddo
      enddo
      enddo

      do j=1,mp
      do i=1,np
      w(i,j,1)=-w(i,j,2)
      enddo
      enddo
       
      do k=1,l
      do j=1,mp
      do i=1,np
      hx(i,j,k)=u(i,j,k)
      hy(i,j,k)=v(i,j,k)
      hz(i,j,k)=w(i,j,k)
      enddo
      enddo
      enddo
      call updated(hx,hx,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      call updated(hy,hy,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      call updated(hz,hz,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)

compute Laplacian term by term
      do k=1,l
      do j=1,mp
      do i=1,np
      r(i,j,k)=0.
      end do
      end do
      end do
compute d/dx(dh/dx)
      do k=1,L
        do j=1,mp
          do i=1,np-rightdedge                       !add d rrl
            r(i,j,k) = dxi*(hx(i+1,j,k)-hx(i,j,k))
          end do
        end do
      end do
create boundary conditions on hx at i=n+1
      if (rightdedge.eq.1) then                       !add d rrl
         do k=1,L
            do j=1,mp
               hxn=hx(np,j,k)
               hxnp1 = (ibcx-1)*hxn + ibcx*hx(np+2,j,k)
               r(np,j,k) = dxi*(hxnp1-hxn)
            end do
         end do
      end if
c for 3D problem compute d/dy(dh/dy) and add to r
      if(j3.eq.1) then
        do k=1,L
          do i=1,np
            do j=1,mp-topdedge                       !add d rrl
              r(i,j,k) = r(i,j,k) + dyi*(hy(i,j+1,k)-hy(i,j,k))
            end do
          end do
        end do
       if (topdedge.eq.1) then                       !add d rrl
          do k=1,L
             do i=1,np
create boundary conditions at j=m+1
                hymp1= (ibcy-1)*hy(i,mp,k) + ibcy*hy(i,mp+2,k)
                r(i,mp,k) = r(i,mp,k) + dyi*(hymp1-hy(i,mp,k))
             end do
          end do
       end if
      endif
compute d/dz(dh/dz) and add to r
      do j=1,mp
        do i=1,np
          do k=1,L-1
            r(i,j,k) = r(i,j,k) + dzi*(hz(i,j,k+1)-hz(i,j,k))
          end do
corporate b.c. hz(i,j,L+1)=-hz(i,j,L) at k=L
          hzLp1 =-hz(i,j,L)
          r(i,j,L) = r(i,j,L) + dzi*(hzLp1-hz(i,j,L))
        end do
      end do

      call updated(fi,fi,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      do k=1,l
      do j=1,mp
      do i=1,np
       r(i,j,k)=r(i,j,k)*gi(i,j,k)
      fi(i,j,k)=fi(i,j,k)+2.*r(i,j,k)*dt
      enddo
      enddo
      enddo
      call updated(fi,fi,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
150   format(5x,e30.15)
151   format(e30.15,2x,e30.15)

      deallocate (ri)
      deallocate (r)
      deallocate (rkb)
      deallocate (hx)
      deallocate (hy)
      deallocate (hz)
      deallocate (pz)
      deallocate (srff)
      deallocate (u)
      deallocate (v)
      deallocate (w)
 
      return
      end subroutine lapdi
!**************************end subroutine
!lapdi*****************************!
      !subroutine lapdis(xv,fi,s,
      !.d11,d22,d33,d12,d13,d23,il,iu,jl,ju,lls,nvp)
      !JAS 3/7/06 changed routine signature here and in calling routine
      !(turb) as "s" and "d-arrays" are not used!
      subroutine lapdis(xv,fi,il,iu,jl,ju,lls,nvp)
      use gridsetup
      use metryic
      use xve, only:xe
      use turba
      use msga

      Implicit None

      !JAS 3/7/06 added explicit declarations to comply with implicit
      !none
      integer,intent(in) :: il,iu,jl,ju,lls,nvp

      real      xv(il:iu, jl:ju,lls,nvp),
      !.          s(il:iu, jl:ju,lls),
     .          fi(il:iu, jl:ju,lls)
      !real            d13(il:iu+1, jl:ju,lls+1),
      !.                d12(il:iu+1, jl:ju+1,lls),
      !.                d23(il:iu, jl:ju+1,lls+1),
      !.                d11(il:iu+1, jl:ju,lls),
      !.                d22(il:iu, jl:ju+1,lls),
      !.                d33(il:iu, jl:ju,lls+1)
            real,allocatable::
     .           r(:, :,:),
     .          ri(:, :,:),
     .         rkb(:, :,:),
     .          hx(:, :,:),
     .          hy(:, :,:),
     .          hz(:, :,:),
     .          pz(:, :,:),
     .        srff(:, :),
     .           u(:, :,:),
     .           v(:, :,:),
     .           w(:, :,:)

      !JAS 3/7/06 added explicit declarations to comply with implicit
      !none
      integer :: i,j,k,illim,iulim,jllim,julim,ip1,im1,jp1,jm1
      real :: hdxi,hdyi,hdzi,g33,g23,g13,gii,px,py,pxa,pya,pza
      real :: hxa,hya,hza,hx1,hy1,hx2,hy2
      real :: hxn,hxnp1,hymp1,hzLp1
      real :: coef,dcr1

      if(iturb.ge.1) then
      allocate (ri(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (r(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (rkb(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (hx(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (hy(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (hz(1-ih:np+ih, 1-ih:mp+ih,l+1))
      allocate (pz(1-ih:np+ih+1, 1-ih:mp+ih+1,l+1))
      allocate (srff(1-ih:np+ih, 1-ih:mp+ih))
      allocate (u(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (v(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (w(1-ih:np+ih, 1-ih:mp+ih,l))
      endif


 
      hdxi=0.5*dxi
      hdyi=0.5*dyi
      hdzi=0.5*dzi

      do k=1,l 
      do j=1,mp 
      do i=1,np 
      ri(i,j,k)=xv(i,j,k,4)/xv(i,j,k,nv)
      rkb(i,j,k)=xv(i,j,k,5)/xv(i,j,k,nv)  !FP
      enddo
      enddo
      enddo
      call updated(ri,ri,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      call updated(rkb,rkb,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)

compute z-derivatives at (i,j,k+-1/2)
      do k=2,L
        do j=1,mp
          do i=1,np
            pz(i,j,k)=dzi*(ri(i,j,k)-ri(i,j,k-1))
          end do
        end do
      end do
      do j=1,mp
        do i=1,np
c          pz(i,j,L+1)=-pz(i,j,L)
          pz(i,j,L+1)=dzi*(xe(i,j,l,4)/xe(i,j,l,nv)-ri(i,j,l))   !rrl
        end do
      end do

      do j=1,mp
      do i=1,np
        srff(i,j)=0.
      enddo
      enddo

      do j=1,mp
      do i=1,np
      g33=(c13(i,j)*gmul(1))**2+(c23(i,j)*gmul(1))**2+gi(i,j,1)**2
      srff(i,j) = sqrt(g33)*srff(i,j)
      enddo
      enddo
      jllim = 1  + (j3-ibcy)*botdedge                       !add d rrl
      julim = mp + (ibcy-j3)*topdedge                       !add d rrl
      do 21 j=jllim,julim
         if (j3.eq.1) then
            if (topdedge.eq.1 .and. j.eq.mp) then
!add d rrl
               jp1 = mp + 2
            else
               jp1 = j + 1
            end if
            if (botdedge.eq.1 .and. j.eq.1) then
!add d rrl
               jm1 = -1
            else
               jm1 = j - 1
            end if
         else
            jp1=1
            jm1=1
         end if
c        jp1=j+j3-j/m*(m-1)
c        jm1=j-j3+(m-j)/(m-j3)*(m-j3)
         illim = 1  + (1-ibcx)*leftdedge                       !add d
         iulim = np + (ibcx-1)*rightdedge                       !add d
         do 21 i=illim,iulim
            if (rightdedge.eq.1 .and. i.eq.np) then
!add d rrl
               ip1 = np + 2
            else
               ip1 = i + 1
            end if
            if (leftdedge.eq.1 .and. i.eq.1) then
!add d rrl
               im1 = -1
            else
               im1 = i - 1
            end if
      g13=c13(i,j)*gmul(1)
      g23=c23(i,j)*gmul(1)
      g33=g13**2+g23**2+gi(i,j,1)**2
      px=hdxi*(ri(ip1,j,1)-ri(im1,j,1))
      py=hdyi*(ri(i,jp1,1)-ri(i,jm1,1))*j3
c   21 pz(i,j,1)=-2.*(srff(i,j)+g13*px+g23*py)/g33-pz(i,j,2)
   21 pz(i,j,1)=-1.*(srff(i,j)+g13*px+g23*py)/g33   !rrl (insulated
c ground)
      if(ibcx.eq.0) then
         illim = 1*leftdedge   + np*(1-leftdedge)
!add d rrl
         iulim = np*rightdedge +  1*(1-rightdedge)
!add d rrl
         do 211 i=illim,iulim,np-1
            jllim = 1  + (j3-ibcy)*botdedge                       !add d
            julim = mp + (ibcy-j3)*topdedge                       !add d
            do 2111 j=jllim,julim
               if (j3.eq.1) then
                  if (topdedge.eq.1 .and. j.eq.mp) then
!add d rrl
                     jp1 = mp + 2
                  else
                     jp1 = j + 1
                  end if
                  if (botdedge.eq.1 .and. j.eq.1) then
!add d rrl
                     jm1 = -1
                  else
                     jm1 = j - 1
                  end if
               else
                  jp1=1
                  jm1=1
               end if
c        jp1=j+j3-j/m*(m-1)
c        jm1=j-j3+(m-j)/(m-j3)*(m-j3)
               g13=c13(i,j)*gmul(1)
               g23=c23(i,j)*gmul(1)
               g33=g13**2+g23**2+gi(i,j,1)**2
               py=hdyi*(ri(i,jp1,1)-ri(i,jm1,1))*j3
c 2111      pz(i,j,1)=-2.*(srff(i,j)+g23*py)/(g33-g13*g13)-pz(i,j,2)
 2111      pz(i,j,1)=-1.*(srff(i,j)+g23*py)/(g33-g13*g13)   !rrl
c insulated ground)
           if(ibcy.eq.0.and.j3.eq.1) then
              jllim = 1*botdedge  + mp*(1-botdedge)
!add d rrl
              julim = mp*topdedge +  1*(1-topdedge)
!add d rrl
              do 2112 j=jllim,julim,mp-j3
c 2112         pz(i,j,1)=-2.*srff(i,j)/(gi(i,j,1)**2)-
c     .                  pz(i,j,2)
 2112         pz(i,j,1)=-1.*srff(i,j)/(gi(i,j,1)**2)   !rrl (insulated
c ground)
           endif
 211    continue
      endif
      if(ibcy.eq.0.and.j3.eq.1) then
         jllim = 1*botdedge  + mp*(1-botdedge)
!add d rrl
         julim = mp*topdedge +  1*(1-topdedge)
!add d rrl
         do 212 j=jllim,julim,mp-j3
            illim = 1  + (1-ibcx)*leftdedge                       !add d

            iulim = np + (ibcx-1)*rightdedge                       !add
            do 2121 i=illim,iulim
               if (rightdedge.eq.1 .and. i.eq.np) then
!add d rrl
                  ip1 = np + 2
               else
                  ip1 = i + 1
               end if
               if (leftdedge.eq.1 .and. i.eq.1) then
!add d rrl
                  im1 = -1
               else
                  im1 = i - 1
               end if
c        ip1=i+1-i/n*(n-1)
c        im1=i-1+(n-i)/(n-1)*(n-1)
               g13=c13(i,j)*gmul(1)
               g23=c23(i,j)*gmul(1)
               g33=g13**2+g23**2+gi(i,j,1)**2
               px=hdxi*(ri(ip1,j,1)-ri(im1,j,1))
 2121      pz(i,j,1)=-1.*(srff(i,j)+g13*px)/(g33-g23*g23)  !rrl
c (insulated ground)
 212    continue
      endif

      if (rightdedge.eq.0 .and. topdedge.eq.0) then
!add d rrl
         call updated(pz,pz,np,mp,l+1,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      else if (rightdedge.eq.1 .and. topdedge.eq.0) then
!add d rrl
         call updated(pz,pz,np+1,mp,l+1,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      else if (rightdedge.eq.0 .and. topdedge.eq.1) then
!add d rrl
         call updated(pz,pz,np,mp+1,l+1,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      else
         call updated(pz,pz,np+1,mp+1,l+1,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      end if

       
compute x-flux at (i+-1/2,j,k)
      do j=1,mp
        do i=1+leftdedge,np                       !add d rrl
          do k=1,L
c            coef=0.66*0.09*0.5*(sb(i,j,k)+sb(i-1,j,k))*   
c     .           sqrt(1.2*0.5*(rkb(i,j,k)+rkb(i-1,j,k)))
c     .      *0.5*(xv(i-1,j,k,nv)+xv(i,j,k,nv))
              coef=0.09*rturbprandtl*0.5*(sb(i,j,k)+sb(i-1,j,k))*   
     .           sqrt(0.5*(rkb(i,j,k)+rkb(i-1,j,k)))
     .        *0.5*(xv(i-1,j,k,nv)+xv(i,j,k,nv))
            dcr1=3.*coef*dt*(dxi**2+j3*dyi**2+dzi**2)
            if(abs(dcr1).gt.0.20)
     .      write(6,*) 'diffusion too strong in theta',i,j,k
c     .coef=0.20/(3.*dt*(dxi**2+j3*dyi**2+dzi**2))
            g13=0.5*gmul(k)*(c13(i-1,j)+c13(i,j))
            pza=0.25*(pz(i-1,j,k)+pz(i,j,k)+pz(i-1,j,k+1)+pz(i,j,k+1))
            Pxa=dxi*(ri(i,j,k)-ri(i-1,j,k))
            hx(i,j,k)=coef*( pxa + g13*pza )
          end do
        end do
      end do
create boundary conditions at i=1
      call updated(hx,hx,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      if(leftdedge.eq.1) then                       !add d rrl
      do j=1,mp
        do k=1,L
          hx(1,j,k) = (ibcx-1)*hx(2,j,k) + ibcx*hx(0,j,k)

        end do
      end do
      endif
      call updated(hx,hx,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
c
compute y-flux at (i,j+-1/2,k)
      if (j3.eq.1) then
        do i=1,np
          do j=1+botdedge,mp                       !add d rrl
            do k=1,L
c            coef=0.66*0.09*0.5*(sb(i,j,k)+sb(i,j-1,k))*
c     .           sqrt(1.2*0.5*(rkb(i,j,k)+rkb(i,j-1,k)))
c     .      *0.5*(xv(i,j-1,k,nv)+xv(i,j,k,nv))
                coef=0.09*rturbprandtl*0.5*(sb(i,j,k)+sb(i,j-1,k))*
     .           sqrt(0.5*(rkb(i,j,k)+rkb(i,j-1,k)))
     .         *0.5*(xv(i,j-1,k,nv)+xv(i,j,k,nv))
            dcr1=3.*coef*dt*(dxi**2+j3*dyi**2+dzi**2)
            if(abs(dcr1).gt.0.20)
     .      write(6,*) 'diffusion too strong in theta',i,j,k
c     .coef=0.20/(3.*dt*(dxi**2+j3*dyi**2+dzi**2))
              g23=0.5*gmul(k)*(c23(i,j-j3)+c23(i,j))
              pza=0.25*(  pz(i,j-j3,k  ) + pz(i,j,k  )
     .                  + pz(i,j-j3,k+1) + pz(i,j,k+1) )
              pya=dyi*(ri(i,j,k)-ri(i,j-j3,k))
              hy(i,j,k)=coef*( pya + g23*pza )
            end do
          end do
        end do
create boundary conditions at j=1
        call updated(hy,hy,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
        if(botdedge.eq.1) then                       !add d rrl
        do k=1,l
          do i=1,np
            hy(i,1,k)= (ibcy-1)*hy(i,2,k) + ibcy*hy(i,0,k)
          end do
        end do
        endif
        call updated(hy,hy,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      endif

compute z-flux at (i,j,k+-1/2)
! i) include the dh/dx and dh/dz terms
      do k=2,L
        do j=1,mp
          do i=1,np-1*rightdedge                       !add d rrl
c             coef=0.66*0.09*0.5*(sb(i,j,k)+sb(i,j,k-1))*
c     .           sqrt(1.2*0.5*(rkb(i,j,k)+rkb(i,j,k-1)))
c     .      *0.5*(xv(i,j,k-1,nv)+xv(i,j,k,nv))
               coef=0.09*rturbprandtl*0.5*(sb(i,j,k)+sb(i,j,k-1))*
     .           sqrt(0.5*(rkb(i,j,k)+rkb(i,j,k-1)))
     .          *0.5*(xv(i,j,k-1,nv)+xv(i,j,k,nv))
            dcr1=3.*coef*dt*(dxi**2+j3*dyi**2+dzi**2)
            if(abs(dcr1).gt.0.20)
     .      write(6,*) 'diffusion too strong in theta',i,j,k
c     .coef=0.20/(3.*dt*(dxi**2+j3*dyi**2+dzi**2))
            g13=0.5*(gmul(k)+gmul(k-1))*c13(i,j)
            gii=0.5*(gi(i,j,k)+gi(i,j,k-1))
            hxa=0.25*(hx(i,j,k-1)+hx(i+1,j,k-1)+hx(i,j,k)+hx(i+1,j,k))
            hza=gii*pz(i,j,k)
            hz(i,j,k)= coef*gii*hza + g13*hxa
          end do
          if (rightdedge.eq.1) then                       !add d rrl
corporate b.c. for hx on i=n+1
c             coef=0.66*0.09*0.5*(sb(np,j,k)+sb(np,j,k-1))*
c     .           sqrt(1.2*0.5*(rkb(np,j,k)+rkb(np,j,k-1)))
c     .      *0.5*(xv(np,j,k-1,nv)+xv(np,j,k,nv))
                coef=0.09*rturbprandtl*0.5*(sb(np,j,k)+sb(np,j,k-1))*
     .           sqrt(0.5*(rkb(np,j,k)+rkb(np,j,k-1)))
     .         *0.5*(xv(np,j,k-1,nv)+xv(np,j,k,nv))
            dcr1=3.*coef*dt*(dxi**2+j3*dyi**2+dzi**2)
            if(abs(dcr1).gt.0.20)
c     .coef=0.20/(3.*dt*(dxi**2+j3*dyi**2+dzi**2))
     .      write(6,*) 'diffusion too strong in theta',i,j,k
             g13=0.5*(gmul(k)+gmul(k-1))*c13(np,j)
             gii=0.5*(gi(np,j,k)+gi(np,j,k-1))
             hx1 = (ibcx-1)*hx(np,j,k-1) + ibcx*hx(np+2,j,k-1)
             hx2 = (ibcx-1)*hx(np,j,k) + ibcx*hx(np+2,j,k)
             hxa=0.25*(hx(np,j,k-1)+hx1+hx(np,j,k)+hx2)
             hza=gii*pz(np,j,k)
             hz(np,j,k)= coef*gii*hza + g13*hxa
          end if
        end do
      end do
      k=1
      do j=1,mp
          do i=1,np-1*rightdedge                       !add d rrl
c             coef=0.66*0.09*0.5*(sb(i,j,k)+sb(i,j,k-1))*
c     .           sqrt(1.2*0.5*(rkb(i,j,k)+rkb(i,j,k-1)))
c     .      *0.5*(xv(i,j,k-1,nv)+xv(i,j,k,nv))
c             coef=0.66*0.09*sb(i,j,k)*
c     .           sqrt(1.2*rkb(i,j,k))
c     .          *xv(i,j,k,nv)
             coef=0.09*rturbprandtl*sb(i,j,k)*
     .           sqrt(rkb(i,j,k))
     .          *xv(i,j,k,nv)
c234567
            dcr1=3.*coef*dt*(dxi**2+j3*dyi**2+dzi**2)
            if(abs(dcr1).gt.0.20)
!     .coef=0.20/(3.*dt*(dxi**2+j3*dyi**2+dzi**2))
     .      write(6,*) 'diffusion too strong in theta',i,j,k
            g13=gmul(k)*c13(i,j)
            gii=gi(i,j,k)
            hxa=0.5*(hx(i,j,k)+hx(i+1,j,k))
            hza=gii*pz(i,j,k)
            hz(i,j,k)= coef*gii*hza + g13*hxa
          end do
          if (rightdedge.eq.1) then                       !add d rrl
corporate b.c. for hx on i=n+1
c             coef=0.66*0.09*0.5*(sb(np,j,k)+sb(np,j,k-1))*
c     .           sqrt(1.2*0.5*(rkb(np,j,k)+rkb(np,j,k-1)))
c     .      *0.5*(xv(np,j,k-1,nv)+xv(np,j,k,nv))
cc             coef=0.66*0.09*sb(i,j,k)*
c     .           sqrt(1.2*rkb(i,j,k))
c     .      *xv(i,j,k,nv)
             coef=0.09*rturbprandtl*sb(i,j,k)*
     .           sqrt(rkb(i,j,k))
     .      *xv(i,j,k,nv)

            dcr1=3.*coef*dt*(dxi**2+j3*dyi**2+dzi**2)
            if(abs(dcr1).gt.0.20)
c     .coef=0.20/(3.*dt*(dxi**2+j3*dyi**2+dzi**2))
     .      write(6,*) 'diffusion too strong in theta',i,j,k
             g13=gmul(k)*c13(np,j)
             gii=gi(np,j,k)
c             hx1 = (ibcx-1)*hx(np,j,k-1) + ibcx*hx(np+2,j,k-1)
             hx2 = (ibcx-1)*hx(np,j,k) + ibcx*hx(np+2,j,k)
             hxa=0.5*(hx(np,j,k)+hx2)
             hza=gii*pz(np,j,k)
             hz(np,j,k)= coef*gii*hza + g13*hxa
          end if
        end do
        k=l+1
        do j=1,mp
          do i=1,np-1*rightdedge                       !add d rrl
c             coef=0.66*0.09*0.5*(sb(i,j,k)+sb(i,j,k-1))*
c     .           sqrt(1.2*0.5*(rkb(i,j,k)+rkb(i,j,k-1)))
c     .      *0.5*(xv(i,j,k-1,nv)+xv(i,j,k,nv))
c             coef=0.66*0.09*sb(i,j,k-1)*
c     .           sqrt(1.2*rkb(i,j,k-1))
c     .      *xv(i,j,k-1,nv)
             coef=0.09*rturbprandtl*sb(i,j,k-1)*
     .           sqrt(rkb(i,j,k-1))
     .      *xv(i,j,k-1,nv)

            dcr1=3.*coef*dt*(dxi**2+j3*dyi**2+dzi**2)
            if(abs(dcr1).gt.0.20)
     .      write(6,*) 'diffusion too strong in theta',i,j,k
c     .coef=0.20/(3.*dt*(dxi**2+j3*dyi**2+dzi**2))
            g13=gmul(k-1)*c13(i,j)
            gii=gi(i,j,k-1)
            hxa=0.5*(hx(i,j,k-1)+hx(i+1,j,k-1))
            hza=gii*pz(i,j,k)
            hz(i,j,k)= coef*gii*hza + g13*hxa
          end do
          if (rightdedge.eq.1) then                       !add d rrl
corporate b.c. for hx on i=n+1
c             coef=0.66*0.09*0.5*(sb(np,j,k)+sb(np,j,k-1))*
c     .           sqrt(1.2*0.5*(rkb(np,j,k)+rkb(np,j,k-1)))
c     .      *0.5*(xv(np,j,k-1,nv)+xv(np,j,k,nv))
c             coef=0.66*0.09*sb(i,j,k-1)*
c     .           sqrt(1.2*rkb(i,j,k-1))
c     .      *xv(i,j,k-1,nv)
             coef=0.09*rturbprandtl*sb(i,j,k-1)*
     .           sqrt(rkb(i,j,k-1))
     .      *xv(i,j,k-1,nv)

            dcr1=3.*coef*dt*(dxi**2+j3*dyi**2+dzi**2)
            if(abs(dcr1).gt.0.20)
     .      write(6,*) 'diffusion too strong in theta',i,j,k
c     .coef=0.20/(3.*dt*(dxi**2+j3*dyi**2+dzi**2))

             g13=gmul(k-1)*c13(np,j)
             gii=gi(np,j,k-1)
             hx1 = (ibcx-1)*hx(np,j,k-1) + ibcx*hx(np+2,j,k-1)
c             hx2 = (ibcx-1)*hx(np,j,k) + ibcx*hx(np+2,j,k)
             hxa=.5*(hx(np,j,k-1)+hx1)
             hza=gii*pz(np,j,k)
             hz(np,j,k)= coef*gii*hza + g13*hxa
          end if
        end do

! ii) include the dh/dy term if 3D
      if(j3.eq.1) then
      call updated(hy,hy,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
        do k=2,L
          do i=1,np
            do j=1,mp-topdedge                       !add d rrl
              g23=0.5*(gmul(k)+gmul(k-1))*c23(i,j)
              hya=0.25*(hy(i,j,k-1)+hy(i,j+j3,k-1)+

     1                  hy(i,j,k)+hy(i,j+j3,k))
              hz(i,j,k)=hz(i,j,k) + g23*hya
            end do
          end do
        end do
      if (topdedge.eq.1) then                       !add d rrl
        do k=2,L
          do i=1,np
corporate b.c. for hy on j=m+1
            g23=0.5*(gmul(k)+gmul(k-1))*c23(i,mp)
            hy1= (ibcy-1)*hy(i,mp,k-1) + ibcy*hy(i,mp+2,k-1)
            hy2= (ibcy-1)*hy(i,mp,k) + ibcy*hy(i,mp+2,k)
            hya=0.25*(hy(i,mp,k-1)+hy1+hy(i,mp,k)+hy2)
            hz(i,mp,k)=hz(i,mp,k) + g23*hya
          end do
        end do
      end if
      k=1
         do i=1,np
            do j=1,mp-topdedge                       !add d rrl
              g23=gmul(k)*c23(i,j)
              hya=0.5*(
     1                  hy(i,j,k)+hy(i,j+j3,k))
              hz(i,j,k)=hz(i,j,k) + g23*hya
            end do
        end do
      if (topdedge.eq.1) then                       !add d rrl
          do i=1,np
corporate b.c. for hy on j=m+1
            g23=gmul(k)*c23(i,mp)
c            hy1= (ibcy-1)*hy(i,mp,k-1) + ibcy*hy(i,mp+2,k-1)
            hy2= (ibcy-1)*hy(i,mp,k) + ibcy*hy(i,mp+2,k)
            hya=0.5*(hy(i,mp,k)+hy2)
            hz(i,mp,k)=hz(i,mp,k) + g23*hya
          end do
      end if
      k=l+1
         do i=1,np
            do j=1,mp-topdedge                       !add d rrl
              g23=gmul(k-1)*c23(i,j)
              hya=0.5*(hy(i,j,k-1)+hy(i,j+j3,k-1))
              hz(i,j,k)=hz(i,j,k) + g23*hya
            end do
          end do
      if (topdedge.eq.1) then                       !add d rrl
          do i=1,np
corporate b.c. for hy on j=m+1
            g23=gmul(k-1)*c23(i,mp)
            hy1= (ibcy-1)*hy(i,mp,k-1) + ibcy*hy(i,mp+2,k-1)
c            hy2= (ibcy-1)*hy(i,mp,k) + ibcy*hy(i,mp+2,k)
            hya=0.5*(hy(i,mp,k-1)+hy1)
            hz(i,mp,k)=hz(i,mp,k) + g23*hya
          end do
      end if

      endif
create boundary conditions at k=1; for k=L see divergence below
c surface fluxes:
c      do j=1,mp
c      do i=1,np
c      hz(i,j,1)=-hz(i,j,2)
c      end do
c      end do
      call updated(hz,hz,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
  

compute Laplacian term by term
      do k=1,l
      do j=1,mp
      do i=1,np
      r(i,j,k)=0.
      end do
      end do
      end do
compute d/dx(dh/dx)
      do k=1,L
        do j=1,mp
          do i=1,np-rightdedge                       !add d rrl
            r(i,j,k) = dxi*(hx(i+1,j,k)-hx(i,j,k))
          end do
        end do
      end do
create boundary conditions on hx at i=n+1
      if (rightdedge.eq.1) then                       !add d rrl
         do k=1,L
            do j=1,mp
               hxn=hx(np,j,k)
               hxnp1 = (ibcx-1)*hxn + ibcx*hx(np+2,j,k)
               r(np,j,k) = dxi*(hxnp1-hxn)
            end do
         end do
      end if
c for 3D problem compute d/dy(dh/dy) and add to r
      if(j3.eq.1) then
        do k=1,L
          do i=1,np
            do j=1,mp-topdedge                       !add d rrl
              r(i,j,k) = r(i,j,k) + dyi*(hy(i,j+1,k)-hy(i,j,k))
            end do
          end do
        end do
       if (topdedge.eq.1) then                       !add d rrl
          do k=1,L
             do i=1,np
create boundary conditions at j=m+1
                hymp1= (ibcy-1)*hy(i,mp,k) + ibcy*hy(i,mp+2,k)
                r(i,mp,k) = r(i,mp,k) + dyi*(hymp1-hy(i,mp,k))
             end do
          end do
       end if
      endif
compute d/dz(dh/dz) and add to r
      do j=1,mp
        do i=1,np
          do k=1,L
            r(i,j,k) = r(i,j,k) + dzi*(hz(i,j,k+1)-hz(i,j,k))
          end do
corporate b.c. hz(i,j,L+1)=-hz(i,j,L) at k=L
          hzLp1 =-hz(i,j,L)
          r(i,j,L) = r(i,j,L) + dzi*(hzLp1-hz(i,j,L))
        end do
      end do

      call updated(fi,fi,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      do k=1,l
      do j=1,mp
      do i=1,np
       r(i,j,k)=r(i,j,k)*gi(i,j,k)
      fi(i,j,k)=fi(i,j,k)+2.*r(i,j,k)*dt
      enddo
      enddo
      enddo
      call updated(fi,fi,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
150   format(5x,e30.15)
151   format(e30.15,2x,e30.15)

      deallocate (ri)
      deallocate (r)
      deallocate (rkb)
      deallocate (hx)
      deallocate (hy)
      deallocate (hz)
      deallocate (pz)
      deallocate (srff)
      deallocate (u)
      deallocate (v)
      deallocate (w)
 
      return
      end subroutine lapdis
************************************************************************************
      subroutine lapdo(xv,ro,fox,d11,d22,d33,d12,d13,d23,il,iu,jl,ju,lls,nvp)
      use gridsetup
      use metryic
      use msga
      use turba

      Implicit None

      !JAS 3/7/06 added explicit declarations to comply with implicit
      !none
      integer,intent(in) :: il,iu,jl,ju,lls,nvp   
 
      real      xv(il:iu, jl:ju,lls,nvp),
     .          fox(il:iu, jl:ju,lls),
     .          ro(il:iu,jl:ju,lls)
      real            d13(il:iu+1, jl:ju,lls+1),
     .                d12(il:iu+1, jl:ju+1,lls),
     .                d23(il:iu, jl:ju+1,lls+1),
     .                d11(il:iu+1, jl:ju,lls),
     .                d22(il:iu, jl:ju+1,lls),
     .                d33(il:iu, jl:ju,lls+1) 
            real,allocatable::
     .         rkb(:,:,:),
     .          r(:, :,:),
     .          hx(:, :,:),
     .          hy(:, :,:),
     .          hz(:, :,:),
     .          pz(:, :,:),
     .        srff(:, :),
     .           u(:, :,:),
     .           v(:, :,:),
     .           w(:, :,:)
      
      !JAS 3/7/06 added explicit declarations to comply with implicit
      !none
      integer :: i,j,k,illim,iulim,jllim,julim,ip1,im1,jp1,jm1,kmz
      real :: hdxi,hdyi,hdzi,g33,g23,g13,gii,px,py,pxa,pya,pza
      real :: sba,sca,rkba
      real :: hxa,hya,hza,hx1,hy1,hx2,hy2,hxp1,hxp2,hyp1,hyp2,hzp1,hzp2
      real :: hxn,hxnp1,hymp1,hzLp1
      real :: d12ba,d13ba,d21ba,d23ba,d31ba,d32ba
      real :: sigmab11,sigmac11,sigmab12,sigmac12,sigmab22,sigmac22
      real :: sigmab23,sigmac23,sigmab13,sigmac13,sigmab31,sigmac31
      real :: sigmab32,sigmac32,sigmab33,sigmac33

      if(iturb.ge.1) then
      allocate (r(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (rkb(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (hx(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (hy(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (hz(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (pz(1-ih:np+ih+1, 1-ih:mp+ih+1,l+1))
      allocate (srff(1-ih:np+ih, 1-ih:mp+ih))
      allocate (u(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (v(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (w(1-ih:np+ih, 1-ih:mp+ih,l))
      endif

 
      hdxi=0.5*dxi
      hdyi=0.5*dyi
      hdzi=0.5*dzi

      do k=1,l
      do j=1,mp
      do i=1,np
      ro(i,j,k)=ro(i,j,k)/xv(i,j,k,nv)
      rkb(i,j,k)=xv(i,j,k,6)/xv(i,j,k,nv)
      enddo
      enddo
      enddo
      call updated(ro,ro,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      call updated(rkb,rkb,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)

compute z-derivatives at (i,j,k+-1/2)
      do k=2,L
        do j=1,mp
          do i=1,np
            pz(i,j,k)=dzi*(ro(i,j,k)-ro(i,j,k-1))
          end do
        end do
      end do
      do j=1,mp
        do i=1,np
          pz(i,j,L+1)=-pz(i,j,L)
        end do
      end do

      do j=1,mp
      do i=1,np
        srff(i,j)=0.
      enddo
      enddo

      do j=1,mp
      do i=1,np
      g33=(c13(i,j)*gmul(1))**2+(c23(i,j)*gmul(1))**2+gi(i,j,1)**2
      srff(i,j) = sqrt(g33)*srff(i,j)
      enddo
      enddo

       
      jllim = 1  + (j3-ibcy)*botdedge                       !add d rrl
      julim = mp + (ibcy-j3)*topdedge
      do 21 j=jllim,julim
         if (j3.eq.1) then
            if (topdedge.eq.1 .and. j.eq.mp) then
               jp1 = mp + 2
            else
               jp1 = j + 1
            end if
            if (botdedge.eq.1 .and. j.eq.1) then
               jm1 = -1
            else
               jm1 = j - 1
            end if
         else
            jp1=1
            jm1=1
         end if
c        jp1=j+j3-j/m*(m-1)
c        jm1=j-j3+(m-j)/(m-j3)*(m-j3)
         illim = 1  + (1-ibcx)*leftdedge                       !add d
         iulim = np + (ibcx-1)*rightdedge                       !add d
         do 21 i=illim,iulim
            if (rightdedge.eq.1 .and. i.eq.np) then
!add d rrl
               ip1 = np + 2
            else
               ip1 = i + 1
            end if
            if (leftdedge.eq.1 .and. i.eq.1) then
!add d rrl
               im1 = -1
            else
               im1 = i - 1
            end if
      g13=c13(i,j)*gmul(1)
      g23=c23(i,j)*gmul(1)
      g33=g13**2+g23**2+gi(i,j,1)**2
      px=hdxi*(ro(ip1,j,1)-ro(im1,j,1))
      py=hdyi*(ro(i,jp1,1)-ro(i,jm1,1))*j3
   21 pz(i,j,1)=-2.*(srff(i,j)+g13*px+g23*py)/g33-pz(i,j,2)
      if(ibcx.eq.0) then
         illim = 1*leftdedge   + np*(1-leftdedge)
!add d rrl
         iulim = np*rightdedge +  1*(1-rightdedge)
!add d rrl
         do 211 i=illim,iulim,np-1
            jllim = 1  + (j3-ibcy)*botdedge                       !add d
            julim = mp + (ibcy-j3)*topdedge                       !add d
            do 2111 j=jllim,julim
               if (j3.eq.1) then
                  if (topdedge.eq.1 .and. j.eq.mp) then
!add d rrl
                     jp1 = mp + 2
                  else
                     jp1 = j + 1
                  end if
                  if (botdedge.eq.1 .and. j.eq.1) then
!add d rrl
                     jm1 = -1
                  else
                     jm1 = j - 1
                  end if
               else
                  jp1=1
                  jm1=1
               end if
c        jp1=j+j3-j/m*(m-1)
c        jm1=j-j3+(m-j)/(m-j3)*(m-j3)
               g13=c13(i,j)*gmul(1)
               g23=c23(i,j)*gmul(1)
               g33=g13**2+g23**2+gi(i,j,1)**2
               py=hdyi*(ro(i,jp1,1)-ro(i,jm1,1))*j3
 2111      pz(i,j,1)=-2.*(srff(i,j)+g23*py)/(g33-g13*g13)-pz(i,j,2)
           if(ibcy.eq.0.and.j3.eq.1) then
              jllim = 1*botdedge  + mp*(1-botdedge)
!add d rrl
              julim = mp*topdedge +  1*(1-topdedge)
!add d rrl
              do 2112 j=jllim,julim,mp-j3
 2112         pz(i,j,1)=-2.*srff(i,j)/(gi(i,j,1)**2)-
     .                  pz(i,j,2)
           endif
 211    continue
      endif
      if(ibcy.eq.0.and.j3.eq.1) then
         jllim = 1*botdedge  + mp*(1-botdedge)
!add d rrl
         julim = mp*topdedge +  1*(1-topdedge)
!add d rrl
         do 212 j=jllim,julim,mp-j3
            illim = 1  + (1-ibcx)*leftdedge                       !add d
            iulim = np + (ibcx-1)*rightdedge                       !add
            do 2121 i=illim,iulim
               if (rightdedge.eq.1 .and. i.eq.np) then
!add d rrl
                  ip1 = np + 2
               else
                  ip1 = i + 1
               end if
               if (leftdedge.eq.1 .and. i.eq.1) then
!add d rrl
                  im1 = -1
               else
                  im1 = i - 1
               end if
c        ip1=i+1-i/n*(n-1)
c        im1=i-1+(n-i)/(n-1)*(n-1)
               g13=c13(i,j)*gmul(1)
               g23=c23(i,j)*gmul(1)
               g33=g13**2+g23**2+gi(i,j,1)**2
               px=hdxi*(ro(ip1,j,1)-ro(im1,j,1))
 2121      pz(i,j,1)=-2.*(srff(i,j)+g13*px)/(g33-g23*g23)-pz(i,j,2)
 212    continue
      endif
      if (rightdedge.eq.0 .and. topdedge.eq.0) then
!add d rrl
         call updated(pz,pz,np,mp,l+1,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      else if (rightdedge.eq.1 .and. topdedge.eq.0) then
!add d rrl
         call updated(pz,pz,np+1,mp,l+1,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      else if (rightdedge.eq.0 .and. topdedge.eq.1) then
!add d rrl
         call updated(pz,pz,np,mp+1,l+1,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      else
         call updated(pz,pz,np+1,mp+1,l+1,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      end if
compute x-flux at (i+-1/2,j,k)
      do j=1,mp
        do i=1+leftdedge,np                       !add d rrl
          do k=1,L
            g13=0.5*gmul(k)*(c13(i-1,j)+c13(i,j))
            pza=0.25*(pz(i-1,j,k)+pz(i,j,k)+pz(i-1,j,k+1)+pz(i,j,k+1))
            Pxa=dxi*(ro(i,j,k)-ro(i-1,j,k))
            hx(i,j,k)=( pxa + g13*pza )
          end do
        end do
      end do
create boundary conditions at i=1
      call updated(hx,hx,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      if(leftdedge.eq.1) then                       !add d rrl
      do j=1,mp
        do k=1,L
          hx(1,j,k) = (ibcx-1)*hx(2,j,k) + ibcx*hx(0,j,k)
        end do
      end do
      endif
      call updated(hx,hx,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
c
compute y-flux at (i,j+-1/2,k)
      if (j3.eq.1) then
        do i=1,np
          do j=1+botdedge,mp                       !add d rrl
            do k=1,L
              g23=0.5*gmul(k)*(c23(i,j-j3)+c23(i,j))
              pza=0.25*(  pz(i,j-j3,k  ) + pz(i,j,k  )
     .                  + pz(i,j-j3,k+1) + pz(i,j,k+1) )
              pya=dyi*(ro(i,j,k)-ro(i,j-j3,k))
              hy(i,j,k)=( pya + g23*pza )
            end do
          end do
        end do
create boundary conditions at j=1
        call updated(hy,hy,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
        if(botdedge.eq.1) then                       !add d rrl
        do k=1,l
          do i=1,np
            hy(i,1,k)= (ibcy-1)*hy(i,2,k) + ibcy*hy(i,0,k)
          end do
        end do
        endif
        call updated(hy,hy,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      endif

compute z-flux at (i,j,k+-1/2)
! i) include the dh/dx and dh/dz terms
      do k=2,L
        do j=1,mp
          do i=1,np-1*rightdedge                       !add d rrl
            g13=0.5*(gmul(k)+gmul(k-1))*c13(i,j)
            gii=0.5*(gi(i,j,k)+gi(i,j,k-1))
            hxa=0.25*(hx(i,j,k-1)+hx(i+1,j,k-1)+hx(i,j,k)+hx(i+1,j,k))
            hza=gii*pz(i,j,k)
            hz(i,j,k)= gii*hza + g13*hxa
          end do
          if (rightdedge.eq.1) then                       !add d rrl
corporate b.c. for hx on i=n+1
             g13=0.5*(gmul(k)+gmul(k-1))*c13(np,j)
             gii=0.5*(gi(np,j,k)+gi(np,j,k-1))
             hx1 = (ibcx-1)*hx(np,j,k-1) + ibcx*hx(np+2,j,k-1)
             hx2 = (ibcx-1)*hx(np,j,k) + ibcx*hx(np+2,j,k)
             hxa=0.25*(hx(np,j,k-1)+hx1+hx(np,j,k)+hx2)
             hza=gii*pz(np,j,k)
             hz(np,j,k)= gii*hza + g13*hxa
          end if
        end do
      end do
! ii) include the dh/dy term if 3D
      if(j3.eq.1) then
        do k=2,L
          do i=1,np
            do j=1,mp-topdedge                       !add d rrl
              g23=0.5*(gmul(k)+gmul(k-1))*c23(i,j)
              hya=0.25*(hy(i,j,k-1)+hy(i,j+j3,k-1)+
     1                  hy(i,j,k)+hy(i,j+j3,k))
              hz(i,j,k)=hz(i,j,k) + g23*hya
            end do
          end do
        end do
      if (topdedge.eq.1) then                       !add d rrl
        do k=2,L
          do i=1,np
corporate b.c. for hy on j=m+1
            g23=0.5*(gmul(k)+gmul(k-1))*c23(i,mp)
            hy1= (ibcy-1)*hy(i,mp,k-1) + ibcy*hy(i,mp+2,k-1)
            hy2= (ibcy-1)*hy(i,mp,k) + ibcy*hy(i,mp+2,k)
            hya=0.25*(hy(i,mp,k-1)+hy1+hy(i,mp,k)+hy2)
            hz(i,mp,k)=hz(i,mp,k) + g23*hya
          end do
        end do
      end if
      endif
create boundary conditions at k=1; for k=L see divergence below
c surface fluxes:
      do j=1,mp
      do i=1,np
      hz(i,j,1)=-hz(i,j,2)
      end do
      end do

compute u
      call updated(hz,hz,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      
      do k=1,l
      do j=1,mp
      do i=1+leftdedge,np                       !add d rrl
      d12ba=0.
      if(j3.eq.1) d12ba=0.5*(d12(i,j,k)+d12(i,j+1,k))
      d13ba=0.5*(d13(i,j,k)+d13(i,j,k+1))
      rkba=0.5*(rkb(i,j,k)+rkb(i-1,j,k))
      sba=0.5*(sb(i,j,k)+sb(i-1,j,k))
      sca=sc
      sigmab11=0.09*sba*d11(i,j,k)/sqrt(rkba)
      sigmac11=0.09*sca*0.2*d11(i,j,k)/sqrt(0.2*rkba)
      sigmab12=0.09*sba*(d12ba)/sqrt(rkba)
      sigmac12=0.09*sca*0.2*(d12ba)/sqrt(0.2*rkba)
      sigmab13=0.09*sba*d13ba/sqrt(rkba)
      sigmac13=0.09*sca*0.2*(d13ba)/sqrt(rkba)
      if(j3.eq.1) hyp1=hy(i,j+1,k)
      if(j3.eq.1) hyp2=hy(i-1,j+1,k)
      if((topdedge.eq.1.and.j.eq.mp).and.j3.eq.1)
!add d rrl
     .                             hyp1=hy(i,j,k)*(ibcy-1)
     .                                 +hy(i,mp+2,k)*ibcy
      if((topdedge.eq.1.and.j.eq.mp).and.j3.eq.1)
!add d rrl
     .                             hyp2=hy(i-1,j,k)*(ibcy-1)
     .                                 +hy(i-1,mp+2,k)*ibcy
      hya=0.
      if(j3.eq.1) hya=0.25*(hy(i,j,k)+hyp1+hy(i-1,j,k)+hyp2)
      if(k.ne.l) hzp1=hz(i,j,k+1)
      if(k.ne.l) hzp2=hz(i-1,j,k+1)
      if(k.eq.l) hzp1=-hz(i,j,k)
      if(k.eq.l) hzp2=-hz(i-1,j,k)
      hza=0.25*(hz(i,j,k)+hzp1+hz(i-1,j,k)+hzp2)
      u(i,j,k)=( (sigmab11+sigmac11)*hx(i,j,k)+
     .           (sigmab12+sigmac12)*hya+
     .           (sigmab13+sigmac13)*hza )
c     if(mpi_rank.eq.2) then
c     if(i.eq.1.and.k.eq.2) write(41,150) u(i,j,k)
c     if(i.eq.1.and.k.eq.2) write(41,151) hx(i,j,k),hza
c     if(i.eq.1.and.k.eq.2) write(41,151) sigmab11,sigmac11
c     if(i.eq.1.and.k.eq.2) write(41,151) sigmab13,sigmac13
c     endif
      enddo
      enddo
      enddo
      call updated(u,u,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)

      if(leftdedge.eq.1) then                       !add d rrl
      do j=1,mp
        do k=1,L
          u(1,j,k) = (ibcx-1)*u(2,j,k) + ibcx*u(0,j,k)
        end do
      end do
      endif
      call updated(u,u,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)

compute v
      if(j3.eq.1) then
      do k=1,l
      do j=1+botdedge,mp                       !add d rrl
      do i=1,np
      d21ba=0.5*(d12(i,j,k)+d12(i+1,j,k))
      d23ba=0.5*(d23(i,j,k)+d23(i,j,k+1))
      rkba=0.5*(rkb(i,j,k)+rkb(i,j-1,k))
      sba=0.5*(sb(i,j,k)+sb(i,j-1,k))
      sca=sc
      sigmab12=0.09*sba*d21ba/sqrt(rkba)
      sigmac12=0.09*sca*0.2*d21ba/sqrt(0.2*rkba)
      sigmab22=0.09*sba*d22(i,j,k)/sqrt(rkba)
      sigmac22=0.09*sca*0.2*d22(i,j,k)/sqrt(rkba)
      sigmab23=0.09*sba*d23ba/sqrt(rkba)
      sigmac23=0.09*sca*0.2*d23ba/sqrt(0.2*rkba)
      hxp1=hx(i+1,j,k)
      hxp2=hx(i+1,j-1,k)
      if(rightdedge.eq.1.and.i.eq.np) hxp1=hx(i,j,k)*(ibcx-1)
!add d rrl
     .                                   +hx(np+2,j,k)*ibcx
      if(rightdedge.eq.1.and.i.eq.np) hxp2=hx(i,j-1,k)*(ibcx-1)
!add d rrl
     .                                   +hx(np+2,j-1,k)*ibcx
      hxa=0.25*(hx(i,j,k)+hxp1+hx(i,j-1,k)+hxp2)
      if(k.ne.l) hzp1=hz(i,j,k+1)
      if(k.ne.l) hzp2=hz(i,j-1,k+1)
      if(k.eq.l) hzp1=-hz(i,j,k)
      if(k.eq.l) hzp2=-hz(i,j-1,k)
      hza=0.25*(hz(i,j,k)+hzp1+hz(i,j-1,k)+hzp2)
      v(i,j,k)=( (sigmab12+sigmac12)*hx(i,j,k)+
     .           (sigmab22+sigmac22)*hy(i,j,k)+
     .           (sigmab23+sigmac23)*hz(i,j,k) )
      enddo
      enddo
      enddo
      call updated(v,v,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      if(botdedge.eq.1) then                       !add d rrl
      do k=1,l
      do i=1,np
      v(i,1,k)= (ibcy-1)*v(i,2,k) + ibcy*v(i,0,k)
      end do
      end do
      endif
      call updated(v,v,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      endif
compute w
      do k=2,l
      do j=1,mp
      do i=1,np
      kmz=k-1
      if(k.eq.1) kmz=1
      d31ba=0.5*(d13(i,j,k)+d13(i+1,j,k))
      d32ba=0.
      if(j3.eq.1) d32ba=0.5*(d23(i,j,k)+d23(i,j+1,k))
      rkba=0.5*(rkb(i,j,k)+rkb(i,j,k-1))
      sba=0.5*(sb(i,j,k)+sb(i,j,kmz))
      sca=sc
      sigmab31=0.09*sba*d31ba/sqrt(rkba)
      sigmac31=0.09*sca*0.2*d31ba/sqrt(0.2*rkba)
      sigmab32=0.09*sba*d32ba/sqrt(rkba)
      sigmac32=0.09*sca*0.2*d32ba/sqrt(0.2*rkba)
      sigmab33=0.09*0.5*sba*d33(i,j,k)/sqrt(rkba)
      sigmac33=0.09*sca*0.2*d33(i,j,k)/sqrt(0.2*rkba)
      hxp1=hx(i+1,j,k)
      hxp2=hx(i+1,j,k-1)
      if(rightdedge.eq.1.and.i.eq.np) hxp1=hx(i,j,k)*(ibcx-1)
!add d rrl
     .                                   +hx(np+2,j,k)*ibcx
      if(rightdedge.eq.1.and.i.eq.np) hxp2=hx(i,j,k-1)*(ibcx-1)
!add d rrl
     .                                   +hx(np+2,j,k-1)*ibcx
      hxa=0.25*(hx(i,j,k)+hx(i,j,k-1)+hxp1+hxp2)
      if(j3.eq.1) hyp1=hy(i,j+1,k)
      if(j3.eq.1) hyp2=hy(i,j+1,k-1)
      if((topdedge.eq.1.and.j.eq.mp).and.j3.eq.1)
!add d rrl
     .                             hyp1=hy(i,j,k)*(ibcy-1)
     .                                 +hy(i,mp+2,k)*ibcy
      if((topdedge.eq.1.and.j.eq.mp).and.j3.eq.1)
!add d rrl
     .                             hyp2=hy(i,j,k-1)*(ibcy-1)
     .                                 +hy(i,mp+2,k-1)*ibcy
      hya=0.
      if(j3.eq.1) hya=0.25*(hy(i,j,k)+hyp1+hy(i,j,k-1)+hyp2)
      w(i,j,k)=( (sigmab31+sigmac31)*hxa+
     .           (sigmab32+sigmac32)*hya+
     .           (sigmab33+sigmac33)*hz(i,j,k) )
      enddo
      enddo
      enddo

      do j=1,mp
      do i=1,np
      w(i,j,1)=-w(i,j,2)
      enddo
      enddo

      do k=1,l
      do j=1,mp
      do i=1,np
      hx(i,j,k)=u(i,j,k)
      hy(i,j,k)=v(i,j,k)
      hz(i,j,k)=w(i,j,k)
      enddo
      enddo
      enddo
      call updated(hx,hx,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      call updated(hy,hy,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      call updated(hz,hz,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)

compute Laplacian term by term
      do k=1,l
      do j=1,mp
      do i=1,np
      r(i,j,k)=0.
      end do
      end do
      end do
compute d/dx(dh/dx)
      do k=1,L
        do j=1,mp
          do i=1,np-rightdedge                       !add d rrl
            r(i,j,k) = dxi*(hx(i+1,j,k)-hx(i,j,k))
          end do
        end do
      end do
create boundary conditions on hx at i=n+1
      if (rightdedge.eq.1) then                       !add d rrl
         do k=1,L
            do j=1,mp
               hxn=hx(np,j,k)
               hxnp1 = (ibcx-1)*hxn + ibcx*hx(np+2,j,k)
               r(np,j,k) = dxi*(hxnp1-hxn)
            end do
         end do
      end if
c for 3D problem compute d/dy(dh/dy) and add to r
      if(j3.eq.1) then
        do k=1,L
          do i=1,np
            do j=1,mp-topdedge                       !add d rrl
              r(i,j,k) = r(i,j,k) + dyi*(hy(i,j+1,k)-hy(i,j,k))
            end do
          end do
        end do
       if (topdedge.eq.1) then                       !add d rrl
          do k=1,L
             do i=1,np
create boundary conditions at j=m+1
                hymp1= (ibcy-1)*hy(i,mp,k) + ibcy*hy(i,mp+2,k)
                r(i,mp,k) = r(i,mp,k) + dyi*(hymp1-hy(i,mp,k))
             end do
          end do
       end if
      endif
compute d/dz(dh/dz) and add to r
      do j=1,mp
        do i=1,np
          do k=1,L-1
            r(i,j,k) = r(i,j,k) + dzi*(hz(i,j,k+1)-hz(i,j,k))
          end do
corporate b.c. hz(i,j,L+1)=-hz(i,j,L) at k=L
          hzLp1 =-hz(i,j,L)
          r(i,j,L) = r(i,j,L) + dzi*(hzLp1-hz(i,j,L))
        end do
      end do

      do k=1,l
      do j=1,mp
      do i=1,np
      fox(i,j,k)=fox(i,j,k)+2.*r(i,j,k)*gi(i,j,k)*dt
      ro(i,j,k)=ro(i,j,k)*xv(i,j,k,nv)
      enddo
      enddo
      enddo
      call updated(fox,fox,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      call updated(ro,ro,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
c     if(mpi_rank.eq.1) then
c     write(41,150) fox(25,1,2)
c     write(41,151) hz(25,1,2),hz(25,1,3)
c     write(41,151) hx(25,1,2),u(26,1,2)
c     if(mpi_rank.eq.2) write(41,151) u(1,1,2)
c     if(mpi_rank.eq.1) write(41,151) u(26,1,2)
150   format(e30.15)
151   format(e30.15,2x,e30.15)
c     endif

      deallocate (r)
      deallocate (rkb)
      deallocate (hx)
      deallocate (hy)
      deallocate (hz)
      deallocate (pz)
      deallocate (srff)
      deallocate (u)
      deallocate (v)
      deallocate (w)
 
      return
      end subroutine lapdo
!**************************end subroutine
!lapdo*****************************!
      !subroutine lapdos(xv,ro,fox,
      !.d11,d22,d33,d12,d13,d23,il,iu,jl,ju,lls,nvp)
      !JAS 3/7/06 changed routine signature as d-arrays are not used
      subroutine lapdos(xv,ro,fox,
     .il,iu,jl,ju,lls,nvp)
      use gridsetup
      use metryic
      use xve, only:xe
      use msga
      use turba

      Implicit None

      !JAS 3/7/06 added explicit declarations to comply with implicit
      !none
      integer,intent(in) :: il,iu,jl,ju,lls,nvp   

      real      xv(il:iu, jl:ju,lls,nvp),
     .          fox(il:iu, jl:ju,lls),
     .          ro(il:iu,jl:ju,lls)
      ! real            d13(il:iu+1, jl:ju,lls+1),
      !.                d12(il:iu+1, jl:ju+1,lls),
      !.                d23(il:iu, jl:ju+1,lls+1),
      !.                d11(il:iu+1, jl:ju,lls),
      !.                d22(il:iu, jl:ju+1,lls),
      !.                d33(il:iu, jl:ju,lls+1) 
            real,allocatable::
     .         rkb(:,:,:),
     .          r(:, :,:),
     .          hx(:, :,:),
     .          hy(:, :,:),
     .          hz(:, :,:),
     .          pz(:, :,:),
     .        srff(:, :),
     .           u(:, :,:),
     .           v(:, :,:),
     .           w(:, :,:)
      
      !JAS 3/7/06 added explicit declarations to comply with implicit
      !none
      integer :: i,j,k,illim,iulim,jllim,julim,ip1,im1,jp1,jm1
      real :: hdxi,hdyi,hdzi,g33,g23,g13,gii,px,py,pxa,pya,pza
      real :: hxa,hya,hza,hx1,hy1,hx2,hy2
      real :: hxn,hxnp1,hymp1,hzLp1
      real :: coef,dcr1
 
      if(iturb.ge.1) then
      allocate (r(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (rkb(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (hx(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (hy(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (hz(1-ih:np+ih, 1-ih:mp+ih,l+1))
      allocate (pz(1-ih:np+ih+1, 1-ih:mp+ih+1,l+1))
      allocate (srff(1-ih:np+ih, 1-ih:mp+ih))
      allocate (u(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (v(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (w(1-ih:np+ih, 1-ih:mp+ih,l))
      endif

 
      hdxi=0.5*dxi
      hdyi=0.5*dyi
      hdzi=0.5*dzi

      do k=1,l
      do j=1,mp
      do i=1,np
      ro(i,j,k)=ro(i,j,k)/xv(i,j,k,nv)
      rkb(i,j,k)=xv(i,j,k,6)/xv(i,j,k,nv)
      enddo
      enddo
      enddo
      call updated(ro,ro,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      call updated(rkb,rkb,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)

compute z-derivatives at (i,j,k+-1/2)
      do k=2,L
        do j=1,mp
          do i=1,np
            pz(i,j,k)=dzi*(ro(i,j,k)-ro(i,j,k-1))
          end do
        end do
      end do
      do j=1,mp
        do i=1,np
c          pz(i,j,L+1)=-pz(i,j,L)
          pz(i,j,L+1)=dzi*(xe(i,j,l,4)/xe(i,j,l,nv)-ro(i,j,l))   !rrl
        end do
      end do

      do j=1,mp
      do i=1,np
        srff(i,j)=0.
      enddo
      enddo

      do j=1,mp
      do i=1,np
      g33=(c13(i,j)*gmul(1))**2+(c23(i,j)*gmul(1))**2+gi(i,j,1)**2
      srff(i,j) = sqrt(g33)*srff(i,j)
      enddo
      enddo

       
      jllim = 1  + (j3-ibcy)*botdedge                       !add d rrl
      julim = mp + (ibcy-j3)*topdedge                       !add d rrl
      do 21 j=jllim,julim
         if (j3.eq.1) then
            if (topdedge.eq.1 .and. j.eq.mp) then
!add d rrl
               jp1 = mp + 2
            else
               jp1 = j + 1
            end if
            if (botdedge.eq.1 .and. j.eq.1) then
!add d rrl
               jm1 = -1
            else
               jm1 = j - 1
            end if
         else
            jp1=1
            jm1=1
         end if
c        jp1=j+j3-j/m*(m-1)
c        jm1=j-j3+(m-j)/(m-j3)*(m-j3)
         illim = 1  + (1-ibcx)*leftdedge                       !add d
         iulim = np + (ibcx-1)*rightdedge                       !add d
         do 21 i=illim,iulim
            if (rightdedge.eq.1 .and. i.eq.np) then
!add d rrl
               ip1 = np + 2
            else
               ip1 = i + 1
            end if
            if (leftdedge.eq.1 .and. i.eq.1) then
!add d rrl
               im1 = -1
            else
               im1 = i - 1
            end if
      g13=c13(i,j)*gmul(1)
      g23=c23(i,j)*gmul(1)
      g33=g13**2+g23**2+gi(i,j,1)**2
      px=hdxi*(ro(ip1,j,1)-ro(im1,j,1))
      py=hdyi*(ro(i,jp1,1)-ro(i,jm1,1))*j3
c   21 pz(i,j,1)=-2.*(srff(i,j)+g13*px+g23*py)/g33-pz(i,j,2)
   21 pz(i,j,1)=-1.*(srff(i,j)+g13*px+g23*py)/g33    !rrl (mirror
c ground)
      if(ibcx.eq.0) then
         illim = 1*leftdedge   + np*(1-leftdedge)
!add d rrl
         iulim = np*rightdedge +  1*(1-rightdedge)
!add d rrl
         do 211 i=illim,iulim,np-1
            jllim = 1  + (j3-ibcy)*botdedge                       !add d
            julim = mp + (ibcy-j3)*topdedge                       !add d
            do 2111 j=jllim,julim
               if (j3.eq.1) then
                  if (topdedge.eq.1 .and. j.eq.mp) then
!add d rrl
                     jp1 = mp + 2
                  else
                     jp1 = j + 1
                  end if
                  if (botdedge.eq.1 .and. j.eq.1) then
!add d rrl
                     jm1 = -1
                  else
                     jm1 = j - 1
                  end if
               else
                  jp1=1
                  jm1=1
               end if
c        jp1=j+j3-j/m*(m-1)
c        jm1=j-j3+(m-j)/(m-j3)*(m-j3)
               g13=c13(i,j)*gmul(1)
               g23=c23(i,j)*gmul(1)
               g33=g13**2+g23**2+gi(i,j,1)**2
               py=hdyi*(ro(i,jp1,1)-ro(i,jm1,1))*j3
c 2111      pz(i,j,1)=-2.*(srff(i,j)+g23*py)/(g33-g13*g13)-pz(i,j,2)
 2111      pz(i,j,1)=-1.*(srff(i,j)+g23*py)/(g33-g13*g13)  !rrl (mirror
c cround)
           if(ibcy.eq.0.and.j3.eq.1) then
              jllim = 1*botdedge  + mp*(1-botdedge)
!add d rrl
              julim = mp*topdedge +  1*(1-topdedge)
!add d rrl
              do 2112 j=jllim,julim,mp-j3
c 2112         pz(i,j,1)=-2.*srff(i,j)/(gi(i,j,1)**2)-
c     .                  pz(i,j,2)
 2112         pz(i,j,1)=-1.*srff(i,j)/(gi(i,j,1)**2)     !rrl (mirror
c ground)
           endif
 211    continue
      endif
      if(ibcy.eq.0.and.j3.eq.1) then
         jllim = 1*botdedge  + mp*(1-botdedge)
!add d rrl
         julim = mp*topdedge +  1*(1-topdedge)
!add d rrl
         do 212 j=jllim,julim,mp-j3
            illim = 1  + (1-ibcx)*leftdedge                       !add d
            iulim = np + (ibcx-1)*rightdedge                       !add
            do 2121 i=illim,iulim
               if (rightdedge.eq.1 .and. i.eq.np) then
!add d rrl
                  ip1 = np + 2
               else
                  ip1 = i + 1
               end if
               if (leftdedge.eq.1 .and. i.eq.1) then
!add d rrl
                  im1 = -1
               else
                  im1 = i - 1
               end if
c        ip1=i+1-i/n*(n-1)
c        im1=i-1+(n-i)/(n-1)*(n-1)
               g13=c13(i,j)*gmul(1)
               g23=c23(i,j)*gmul(1)
               g33=g13**2+g23**2+gi(i,j,1)**2
               px=hdxi*(ro(ip1,j,1)-ro(im1,j,1))
c 2121      pz(i,j,1)=-2.*(srff(i,j)+g13*px)/(g33-g23*g23)-pz(i,j,2)
 2121      pz(i,j,1)=-1.*(srff(i,j)+g13*px)/(g33-g23*g23) !rrl (mirror
c at ground)
 212    continue
      endif
      if (rightdedge.eq.0 .and. topdedge.eq.0) then
!add d rrl
         call updated(pz,pz,np,mp,l+1,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      else if (rightdedge.eq.1 .and. topdedge.eq.0) then
!add d rrl
         call updated(pz,pz,np+1,mp,l+1,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      else if (rightdedge.eq.0 .and. topdedge.eq.1) then
!add d rrl
         call updated(pz,pz,np,mp+1,l+1,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      else
         call updated(pz,pz,np+1,mp+1,l+1,1-ih,np+ih+1,1-ih,mp+ih+1,1)
      end if
compute x-flux at (i+-1/2,j,k)
      do j=1,mp
        do i=1+leftdedge,np                       !add d rrl
          do k=1,L
            !coef=.66*0.09*0.5*(sb(i,j,k)+sb(i-1,j,k))*
            coef=0.09*rturbprandtl*0.5*(sb(i,j,k)+sb(i-1,j,k))*
     .           sqrt(1.2*0.5*(rkb(i,j,k)+rkb(i-1,j,k)))
     .      *0.5*(xv(i-1,j,k,nv)+xv(i,j,k,nv))
            dcr1=3.*coef*dt*(dxi**2+j3*dyi**2+dzi**2)
            if(abs(dcr1).gt.0.20)
     .coef=0.20/(3.*dt*(dxi**2+j3*dyi**2+dzi**2))
            g13=0.5*gmul(k)*(c13(i-1,j)+c13(i,j))
            pza=0.25*(pz(i-1,j,k)+pz(i,j,k)+pz(i-1,j,k+1)+pz(i,j,k+1))
            Pxa=dxi*(ro(i,j,k)-ro(i-1,j,k))
            hx(i,j,k)=coef*( pxa + g13*pza )
          end do
        end do
      end do
create boundary conditions at i=1
      call updated(hx,hx,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      if(leftdedge.eq.1) then                       !add d rrl
      do j=1,mp
        do k=1,L
          hx(1,j,k) = (ibcx-1)*hx(2,j,k) + ibcx*hx(0,j,k)
        end do
      end do
      endif
      call updated(hx,hx,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
c
compute y-flux at (i,j+-1/2,k)
      if (j3.eq.1) then
        do i=1,np
          do j=1+botdedge,mp                       !add d rrl
            do k=1,L
             !coef=.66*0.09*0.5*(sb(i,j,k)+sb(i,j-1,k))*
             coef=0.09*rturbprandtl*0.5*(sb(i,j,k)+sb(i,j-1,k))*
     .           sqrt(1.2*0.5*(rkb(i,j,k)+rkb(i,j-1,k)))
     .      *0.5*(xv(i,j-1,k,nv)+xv(i,j,k,nv))
            dcr1=3.*coef*dt*(dxi**2+j3*dyi**2+dzi**2)
            if(abs(dcr1).gt.0.20)
     .coef=0.20/(3.*dt*(dxi**2+j3*dyi**2+dzi**2))
              g23=0.5*gmul(k)*(c23(i,j-j3)+c23(i,j))
              pza=0.25*(  pz(i,j-j3,k  ) + pz(i,j,k  )
     .                  + pz(i,j-j3,k+1) + pz(i,j,k+1) )
              pya=dyi*(ro(i,j,k)-ro(i,j-j3,k))
              hy(i,j,k)=coef*( pya + g23*pza )
            end do
          end do
        end do
create boundary conditions at j=1
        call updated(hy,hy,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
        if(botdedge.eq.1) then                       !add d rrl
        do k=1,l
          do i=1,np
            hy(i,1,k)= (ibcy-1)*hy(i,2,k) + ibcy*hy(i,0,k)
          end do
        end do
        endif
        call updated(hy,hy,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      endif

compute z-flux at (i,j,k+-1/2)
! i) include the dh/dx and dh/dz terms
      do k=2,L
        do j=1,mp
          do i=1,np-1*rightdedge                       !add d rrl
            !coef=0.66*0.09*0.5*(sb(i,j,k)+sb(i,j,k-1))*
            coef=0.09*rturbprandtl*0.5*(sb(i,j,k)+sb(i,j,k-1))*
     .           sqrt(1.2*0.5*(rkb(i,j,k)+rkb(i,j,k-1)))
     .      *0.5*(xv(i,j,k-1,nv)+xv(i,j,k,nv))
            dcr1=3.*coef*dt*(dxi**2+j3*dyi**2+dzi**2)
            if(abs(dcr1).gt.0.20)
     .coef=0.20/(3.*dt*(dxi**2+j3*dyi**2+dzi**2))
            g13=0.5*(gmul(k)+gmul(k-1))*c13(i,j)
            gii=0.5*(gi(i,j,k)+gi(i,j,k-1))
            hxa=0.25*(hx(i,j,k-1)+hx(i+1,j,k-1)+hx(i,j,k)+hx(i+1,j,k))
            hza=gii*pz(i,j,k)
            hz(i,j,k)= coef*gii*hza + g13*hxa
          end do
          if (rightdedge.eq.1) then                       !add d rrl
corporate b.c. for hx on i=n+1
             coef=.66*.09*.5*(sb(np,j,k)+sb(np,j,k-1))*
     .           sqrt(1.2*0.5*(rkb(np,j,k)+rkb(np,j,k-1)))
     .      *0.5*(xv(np,j,k-1,nv)+xv(np,j,k,nv))
            dcr1=coef*dt*(dxi**2+j3*dyi**2+dzi**2)
            if(abs(dcr1).gt.0.20)
     .coef=0.20/(dt*(dxi**2+j3*dyi**2+dzi**2))
             g13=0.5*(gmul(k)+gmul(k-1))*c13(np,j)
             gii=0.5*(gi(np,j,k)+gi(np,j,k-1))
             hx1 = (ibcx-1)*hx(np,j,k-1) + ibcx*hx(np+2,j,k-1)
             hx2 = (ibcx-1)*hx(np,j,k) + ibcx*hx(np+2,j,k)
             hxa=0.25*(hx(np,j,k-1)+hx1+hx(np,j,k)+hx2)
             hza=gii*pz(np,j,k)
             hz(np,j,k)= coef*gii*hza + g13*hxa
          end if
        end do
      end do
      k=1
      do j=1,mp
          do i=1,np-1*rightdedge                       !add d rrl
c            coef=0.66*0.09*0.5*(sb(i,j,k)+sb(i,j,k-1))*
c     .           sqrt(1.2*0.5*(rkb(i,j,k)+rkb(i,j,k-1)))
c     .      *0.5*(xv(i,j,k-1,nv)+xv(i,j,k,nv))
             !coef=0.66*0.09*sb(i,j,k)*
             coef=0.09*rturbprandtl*sb(i,j,k)*
     .           sqrt(1.2*rkb(i,j,k))
     .      *xv(i,j,k,nv)

            dcr1=3.*coef*dt*(dxi**2+j3*dyi**2+dzi**2)
            if(abs(dcr1).gt.0.20)
     .coef=0.20/(3.*dt*(dxi**2+j3*dyi**2+dzi**2))
            g13=gmul(k)*c13(i,j)
            gii=gi(i,j,k)
            hxa=0.5*(hx(i,j,k)+hx(i+1,j,k))
            hza=gii*pz(i,j,k)
            hz(i,j,k)= coef*gii*hza + g13*hxa
          end do
          if (rightdedge.eq.1) then                       !add d rrl
corporate b.c. for hx on i=n+1
c             coef=.66*0.09*.5*(sb(np,j,k)+sb(np,j,k-1))*
c     .           sqrt(1.2*0.5*(rkb(np,j,k)+rkb(np,j,k-1)))
c     .      *0.5*(xv(np,j,k-1,nv)+xv(np,j,k,nv))
             !coef=0.66*0.09*sb(i,j,k)*
             coef=0.09*rturbprandtl*sb(i,j,k)*
     .           sqrt(1.2*rkb(i,j,k))
     .      *xv(i,j,k,nv)
            dcr1=coef*dt*(dxi**2+j3*dyi**2+dzi**2)
            if(abs(dcr1).gt.0.20)
     .coef=0.20/(dt*(dxi**2+j3*dyi**2+dzi**2))
             g13=gmul(k)*c13(np,j)
             gii=0.5*gi(np,j,k)
c             hx1 = (ibcx-1)*hx(np,j,k-1) + ibcx*hx(np+2,j,k-1)
             hx2 = (ibcx-1)*hx(np,j,k) + ibcx*hx(np+2,j,k)
             hxa=0.5*(hx(np,j,k)+hx2)
             hza=gii*pz(np,j,k)
             hz(np,j,k)= coef*gii*hza + g13*hxa
          end if
        end do
        k=l+1
        do j=1,mp
          do i=1,np-1*rightdedge                       !add d rrl
c            coef=0.66*0.09*0.5*(sb(i,j,k)+sb(i,j,k-1))*
c     .           sqrt(1.2*0.5*(rkb(i,j,k)+rkb(i,j,k-1)))
c     .      *0.5*(xv(i,j,k-1,nv)+xv(i,j,k,nv))
             !coef=0.66*0.09*sb(i,j,k-1)*
             coef=0.09*rturbprandtl*sb(i,j,k-1)*
     .           sqrt(1.2*rkb(i,j,k-1))
     .      *xv(i,j,k-1,nv)
            dcr1=3.*coef*dt*(dxi**2+j3*dyi**2+dzi**2)
            if(abs(dcr1).gt.0.20)
     .coef=0.20/(3.*dt*(dxi**2+j3*dyi**2+dzi**2))
            g13=gmul(k-1)*c13(i,j)
            gii=gi(i,j,k-1)
            hxa=0.5*(hx(i,j,k-1)+hx(i+1,j,k-1))
            hza=gii*pz(i,j,k)
            hz(i,j,k)= coef*gii*hza + g13*hxa
          end do
          if (rightdedge.eq.1) then                       !add d rrl
corporate b.c. for hx on i=n+1
c             coef=.66*0.09*.5*(sb(np,j,k)+sb(np,j,k-1))*
c     .           sqrt(1.2*0.5*(rkb(np,j,k)+rkb(np,j,k-1)))
c     .      *0.5*(xv(np,j,k-1,nv)+xv(np,j,k,nv))
             !coef=0.66*0.09*sb(i,j,k-1)*
             coef=0.09*rturbprandtl*sb(i,j,k-1)*
     .           sqrt(1.2*rkb(i,j,k-1))
     .      *xv(i,j,k-1,nv)
            dcr1=coef*dt*(dxi**2+j3*dyi**2+dzi**2)
            if(abs(dcr1).gt.0.20)
     .coef=0.20/(dt*(dxi**2+j3*dyi**2+dzi**2))
             g13=gmul(k-1)*c13(np,j)
             gii=gi(np,j,k-1)
             hx1 = (ibcx-1)*hx(np,j,k-1) + ibcx*hx(np+2,j,k-1)
c             hx2 = (ibcx-1)*hx(np,j,k) + ibcx*hx(np+2,j,k)
             hxa=0.5*(hx1+hx(np,j,k-1))
             hza=gii*pz(np,j,k)
             hz(np,j,k)= coef*gii*hza + g13*hxa
          end if
        end do

! ii) include the dh/dy term if 3D
      if(j3.eq.1) then
        do k=2,L
          do i=1,np
            do j=1,mp-topdedge                       !add d rrl
              g23=0.5*(gmul(k)+gmul(k-1))*c23(i,j)
              hya=0.25*(hy(i,j,k-1)+hy(i,j+j3,k-1)+
     1                  hy(i,j,k)+hy(i,j+j3,k))
              hz(i,j,k)=hz(i,j,k) + g23*hya
            end do
          end do
        end do
      if (topdedge.eq.1) then                       !add d rrl
        do k=2,L
          do i=1,np
corporate b.c. for hy on j=m+1
            g23=0.5*(gmul(k)+gmul(k-1))*c23(i,mp)
            hy1= (ibcy-1)*hy(i,mp,k-1) + ibcy*hy(i,mp+2,k-1)
            hy2= (ibcy-1)*hy(i,mp,k) + ibcy*hy(i,mp+2,k)
            hya=0.25*(hy(i,mp,k-1)+hy1+hy(i,mp,k)+hy2)
            hz(i,mp,k)=hz(i,mp,k) + g23*hya
          end do
        end do
      end if
      k=1
      do i=1,np
            do j=1,mp-topdedge                       !add d rrl
              g23=gmul(k)*c23(i,j)
              hya=0.5*(
     1                  hy(i,j,k)+hy(i,j+j3,k))
              hz(i,j,k)=hz(i,j,k) + g23*hya
            end do
          end do
      if (topdedge.eq.1) then                       !add d rrl
          do i=1,np
corporate b.c. for hy on j=m+1
            g23=gmul(k)*c23(i,mp)
c            hy1= (ibcy-1)*hy(i,mp,k-1) + ibcy*hy(i,mp+2,k-1)
            hy2= (ibcy-1)*hy(i,mp,k) + ibcy*hy(i,mp+2,k)
            hya=0.5*(hy(i,mp,k)+hy2)
            hz(i,mp,k)=hz(i,mp,k) + g23*hya
          end do
        endif
        k=l+1
        do i=1,np
            do j=1,mp-topdedge                       !add d rrl
              g23=gmul(k-1)*c23(i,j)
              hya=0.5*(hy(i,j,k-1)+hy(i,j+j3,k-1))
              hz(i,j,k)=hz(i,j,k) + g23*hya
            end do
          end do
      if (topdedge.eq.1) then                       !add d rrl
          do i=1,np
corporate b.c. for hy on j=m+1
            g23=gmul(k-1)*c23(i,mp)
            hy1= (ibcy-1)*hy(i,mp,k-1) + ibcy*hy(i,mp+2,k-1)
c            hy2= (ibcy-1)*hy(i,mp,k) + ibcy*hy(i,mp+2,k)
            hya=0.5*(hy(i,mp,k-1)+hy1)
            hz(i,mp,k)=hz(i,mp,k) + g23*hya
          end do

         endif
      endif
create boundary conditions at k=1; for k=L see divergence below
c surface fluxes:
c      do j=1,mp
c      do i=1,np
c      hz(i,j,1)=-hz(i,j,2)
c      end do
c      end do

compute Laplacian term by term
      do k=1,l
      do j=1,mp
      do i=1,np
      r(i,j,k)=0.
      end do
      end do
      end do
compute d/dx(dh/dx)
      do k=1,L
        do j=1,mp
          do i=1,np-rightdedge                       !add d rrl
            r(i,j,k) = dxi*(hx(i+1,j,k)-hx(i,j,k))
          end do
        end do
      end do
create boundary conditions on hx at i=n+1
      if (rightdedge.eq.1) then                       !add d rrl
         do k=1,L
            do j=1,mp
               hxn=hx(np,j,k)
               hxnp1 = (ibcx-1)*hxn + ibcx*hx(np+2,j,k)
               r(np,j,k) = dxi*(hxnp1-hxn)
            end do
         end do
      end if
c for 3D problem compute d/dy(dh/dy) and add to r
      if(j3.eq.1) then
        do k=1,L
          do i=1,np
            do j=1,mp-topdedge                       !add d rrl
              r(i,j,k) = r(i,j,k) + dyi*(hy(i,j+1,k)-hy(i,j,k))
            end do
          end do
        end do
       if (topdedge.eq.1) then                       !add d rrl
          do k=1,L
             do i=1,np
create boundary conditions at j=m+1
                hymp1= (ibcy-1)*hy(i,mp,k) + ibcy*hy(i,mp+2,k)
                r(i,mp,k) = r(i,mp,k) + dyi*(hymp1-hy(i,mp,k))
             end do
          end do
       end if
      endif
compute d/dz(dh/dz) and add to r
      do j=1,mp
        do i=1,np
          do k=1,L
            r(i,j,k) = r(i,j,k) + dzi*(hz(i,j,k+1)-hz(i,j,k))
          end do
corporate b.c. hz(i,j,L+1)=-hz(i,j,L) at k=L
          hzLp1 =-hz(i,j,L)
          r(i,j,L) = r(i,j,L) + dzi*(hzLp1-hz(i,j,L))
        end do
      end do

      do k=1,l
      do j=1,mp
      do i=1,np
      fox(i,j,k)=fox(i,j,k)+2.*r(i,j,k)*gi(i,j,k)*dt
      ro(i,j,k)=ro(i,j,k)*xv(i,j,k,nv)
      enddo
      enddo
      enddo
      call updated(fox,fox,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      call updated(ro,ro,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
c     if(mpi_rank.eq.1) then
c     write(41,150) fox(25,1,2)
c     write(41,151) hz(25,1,2),hz(25,1,3)
c     write(41,151) hx(25,1,2),u(26,1,2)
c     if(mpi_rank.eq.2) write(41,151) u(1,1,2)
c     if(mpi_rank.eq.1) write(41,151) u(26,1,2)
150   format(e30.15)
151   format(e30.15,2x,e30.15)
c     endif

      deallocate (r)
      deallocate (rkb)
      deallocate (hx)
      deallocate (hy)
      deallocate (hz)
      deallocate (pz)
      deallocate (srff)
      deallocate (u)
      deallocate (v)
      deallocate (w)
 
      return
      end subroutine lapdos
      end module nonlocalfunction
