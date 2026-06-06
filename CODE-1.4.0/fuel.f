      !subroutine fuel(xv,xvfuel,il,iu,jl,ju,lls,nvp)
      !subroutine fuel(rkctemp,il,iu,jl,ju,lls,nvp)
      subroutine fuel(ift,i,j,k,rneteng,ff_sum,fw_sum,rhofold)
      use xve
      use xvo
      use metryic
      use fireteca
      use turba
      use turbb
      use constants
      use gridsetup
      use msga
      use pres

      Implicit None

      
      !JAS 3/7/06 added explicit declarationis to comply with implicit none
      integer,intent(in) :: i,j,k,ift
      real,intent(inout) ::rneteng ! defined for gas phase
      real,intent(inout)::ff_sum 
      real,intent(inout)::fw_sum
      real,intent(in)::rhofold
      real,external :: zcart
      real :: sigmac, th 
      real :: perchydroremaining,hydrofactor,tapsw !,taps,tapsw,xmap
      real :: slambdaof,psiwmaxold
      real :: rhowaterold,rhosold !,rhosnph
      real :: reactht,rmassloss,waterevpht,rwaterlossht,tmp
      real :: ff,fw
      
      real, external ::computePsi      
 
      real:: rkctemp !, rkctempref
      real::frhof,frhosies,frhowater,psif,psiw,thetasolid
      real::cf 

c it should be noted that xvfuel(i,j,k,5) corrspondes to xv(i,j,k,6)
c the .2 factor is included as an extrapolation for turbulence at the c
c scale.  The .2 factor should be replaced by a function of sc and sb.
      rkctemp=.2*xvfuel(ift,i,j,k,5)/xvb(i,j,k,nv)
c          rkctempref=.2*xvb(i,j,k,6)/xvb(i,j,k,nv)
c               if (tempg(i,j,k).gt.500)
c     &               write (*,*) rkctempref,rkctemp(i,j,k)

      !    zla=zcart(z(k),i,j)-zs(i,j)
      !    if(zla.le.fueldepth) then
       sigmac=0.
      if (rkctemp.gt.1.e-06) sigmac=sc*0.5*sqrt(rkctemp)
c   readdress the .09 s k**.5 issue             !(rrl)
c the traditional .09 has been left out of sigmac, however the effect of this 
c value is absorbed directly in cf and sigmac is only used in the reaction rate
c the .5 was inserted in the conversion to isotropic turb.  It 
c can be thought of as part of cf.


c computation of psif
      psif =computePsi (temps(ift,i,j,k),ifuel)

c computation of thetasolid and cf
      ! if (rhofinitial(i,j,k).gt.0.0) then
      perchydroremaining=max(0.,(rhof(ift,i,j,k)-rhohydrothresh*rhofinitial(ift,i,j,k))
     .  /(rhofinitial(ift,i,j,k)*(1.-rhohydrothresh)))
      cf=cfhydro*perchydroremaining+cfchar*(1-perchydroremaining)
       !else
       !  cf(i,j,k)=0.0
       !endif
      hydrofactor=exp(-1.*rno/rnfuel*psif*rhofold/xvb(i,j,k,7))
      thetasolid=perchydroremaining*(1.0-thetag)*hydrofactor+(1-perchydroremaining)*thetag

c if statements below are used to identify which part of the ramp
c the average temperature falls on.
c         taps=tramp
c         xmap=(temps(i,j,k)-taps)/(tcrit-taps)
c         if (xmap.lt.0) xmap=0.0
c         if (xmap.gt.2) xmap=2.0

c      if (temps(i,j,k).gt.taps.and.temps(i,j,k).lt.th)
c     .   psif(i,j,k)=(temps(i,j,k)-taps)/(th-taps)       !(rrl 6/1/99)

c      if (temps(i,j,k).ge.th) psif(i,j,k)=1.       

c similar process for the water evap. ramp

c      tapsw=tambient+tstep*(twvap-tambent)/(tcrit-tambient)
c      tapsw=tambient+.25*tstep*(twvap-tambient)/(tcrit-tambient)
      tapsw=tambientarray(i,j,k)+tstep
c      tapsw=tamb+80.
      th=2*twvap-tapsw
      if(temps(ift,i,j,k).le.tapsw) then
          psiw=0.
      else if(temps(ift,i,j,k).ge.th) then
          psiw=1.
      else
        psiw=(temps(ift,i,j,k)-tapsw)/(th-tapsw)        !(rrl 6/1/99)
      endif
      psiwmaxold=psiwmax(ift,i,j,k)
      psiwmax(ift,i,j,k)=max(psiwmax(ift,i,j,k),psiw)

c computation of ff
      slambdaof = rhofold*xvb(i,j,k,7)/(rhofold/rnfuel+xvb(i,j,k,7)/rno)**2
      ff= ffparam*cf*rhof(ift,i,j,k)*xvb(i,j,k,7)*sigmac*
     .          psif*slambdaof/(rhoref*sx*sx)     !(rrl 6/1/99)
      ff_sum = ff_sum+ff
c computation of fw
      if(psiwmax(ift,i,j,k).lt.1.and.rhowater(ift,i,j,k).gt.0.0) then
        fw=fw+min(rhowater(ift,i,j,k)*(psiwmax(ift,i,j,k)-psiwmaxold)/(1.0-psiwmax(ift,i,j,k))/dtp,rhowater(ift,i,j,k)/dtp)
        fw_sum = fw+fw_sum
      endif
c computation of frhof frhowater
        frhof=-ff*dtp*rnfuel
        frhowater=-fw*dtp
c computation of frhofsies 
        reactht=thetasolid*hf*ff
        rmassloss=-tcrit*rnfuel*cpwood*ff !rrl
        waterevpht=-fw*hwevap             !rrl
        rwaterlossht=-fw*cpwater*twvap    !rrl
      tmp=convht(ift,i,j,k)+reactht+rmassloss+waterevpht+rwaterlossht
        frhosies=tmp*dtp
c  update rhof, rhos, rhowater
      rhowaterold=rhowater(ift,i,j,k)
      rhosold=rhos(ift,i,j,k)
      rhof(ift,i,j,k)=rhof(ift,i,j,k)+frhof
      rhowater(ift,i,j,k)=max(0.0,rhowater(ift,i,j,k)+frhowater)
      rhos(ift,i,j,k)=rhof(ift,i,j,k)+rhowater(ift,i,j,k)
c update sies, temps
      sies(ift,i,j,k)=(sies(ift,i,j,k)*rhosold+frhosies+frhosiesrad(ift,i,j,k)*dtp)/rhos(ift,i,j,k)   !FP
      cpsolid(ift,i,j,k)=(rhof(ift,i,j,k)*cpwood+rhowater(ift,i,j,k)*cpwater)/rhos(ift,i,j,k)
      temps(ift,i,j,k)=sies(ift,i,j,k)/cpsolid(ift,i,j,k)
      ! variable for gas phase
      rneteng=rneteng+ff*hf*(1.-thetasolid)  !for gas phase
c update sies, temps
c      if (mpi_rank.eq.26.and.i.eq.16.and.j.eq.13.and.k.eq.5)then
c         write (*,*) 'inside fuel',temps(i,j,k),psiw(i,j,k),
c     &     psiwmax(i,j,k),th,twvap,tapsw
c         write (*,*) rhowater(i,j,k),
c     &        (psiwmax(i,j,k)-psiwmaxold), 
c     &        (1.00-psiwmax(i,j,k)),dtp 
c      endif
c      if (fw(i,j,k).gt.0.0) then
c       write (*,*) i,j,k,temps(i,j,k),rhos(i,j,k),fw(i,j,k),ff(i,j,k)
c     & psiwmax(i,j,k), psiwmaxold,tapsw,tamp,tstep
c23456
c       write (*,*) 'i,j,k,temp,rhowater(i,j,k),psiwmax(i,j,k), 
c     &psiwmaxold,tabsw'
c      pause
c       endif
      !endif
 
      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc computePsi(temps, ifuel) computes the pdf with several methods
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real function computePsi(temps,ifuel)
      use constants
      real ::temps
      integer :: ifuel
      real ::psif,frampmod, trampk,tfstepk,rlcrit,rlramp,rljoint1,tjoint1
      real ::scalefactor,rpsi,rltop,rlramp2,psijoint2,rljoint2,rpsi2,rl
      real er,xerf,c1,c2,c3
      real p,a1,a2,a3,terf,sum
      real erf
      data c1,c2,c3/0.5,0.0079,1.0/
      data p,a1,a2,a3/0.47047,0.3480242,-0.0958798,0.7478556/
c----------------------------------RRL's psif-------------------------------c
      if(ifuel==1)then            !Use RRL's method of calculating psif!
      frampmod=1.
      trampk=-(tcrit-tramp)*frampmod+tcrit
      trampk=tramp
      tfstepk=tfstep

      rlcrit=(tcrit-tfstepk)
      rlramp=(trampk-tfstepk)
c      write (*,*) 'rlcrit is',rlcrit,' rlramp is',rlramp

      rljoint1=psijoint*2.*(rlcrit-rlramp)+rlramp
      tjoint1=rljoint1+tfstepk

      scalefactor=sqrt(rlramp**2-(rljoint1-rlramp)**2)/psijoint
        rpsi=2.0/scalefactor*(rlcrit-rlramp)*rljoint1+scalefactor*psijoint

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
       ! temps is used here
       rl=temps-tfstepk
       if (rl.le.0.0 ) then
         psif=0.0
       elseif (rl.le.rljoint1) then
         psif=1./scalefactor*(rpsi-sqrt(rpsi**2-rl**2))
       elseif (rl.le.rljoint2) then
         psif=.5/(rlcrit-rlramp)*(rl-rlramp)
       elseif (rl.le.rltop) then
         psif=1./scalefactor*(rpsi2+sqrt(rpsi**2-(rl-rltop)**2))
       elseif (rl.gt.rltop) then
         psif=1.
       endif
c------------------------------- end RRL's psif-------------------------------c

      elseif(ifuel==2)then            !Use Michael Clark's method of calculating psif!
                  
c*******************************************************************
c Michael Clark 7/11/05
c 
c A subroutine that scales the error function and scales it
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
c  seems to slow things down maybe 15%
c******************************************************************* 
          
              xerf = c2*(temps-tcrit)
              terf = 1./(1.+p*abs(xerf))
              sum  = a1*terf+a2*terf*terf+a3*terf*terf*terf
              er   = 1-sum*exp(-xerf*xerf)
              
              if (xerf.lt.0.0) psif=-psif+c3
c
c-------------------------------------------------------------------i
      elseif(ifuel==3)then            !Use JAS's method of calculating psif!
c  reimplemented by JAS 2/2/06 to use intrinsic erf (much faster!)

              
              !look closely at the integral erf calculates to see why here!
c      psif(i,j,k) = c1*(er*((temps(i,j,k)-tcrit)/abs(temps(i,j,k)-tcrit))+c3)  

              if(temps.le.tfstep)then 
                psif=0.0
              elseif(temps.le.(2*(tcrit-tfstep)+tfstep))then
                er = erf(c2*(temps-tcrit))                         
                psif = c1*(c3+er)
              else 
                psif = 1.0
              endif  
           
c---------------------------------------------------------------------
      endif    !end if(ifuel==1,2,or 3)        
        computePsi=psif
        return
       end function computePsi
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

