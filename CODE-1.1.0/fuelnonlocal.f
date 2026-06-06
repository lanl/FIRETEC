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

      !JAS 3/7/06 added explicit declarationis to comply with implicit none
      integer,intent(in) :: il,iu,jl,ju,lls,nvp
      
      !JAS 3/7/06 added explicit declarationis to comply with implicit none
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
c the traditional .09 has been left out of sigmac, however the effect of this 
c value is absorbed directly in cf and sigmac is only used in the reaction rate
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

c----------------------------------RRL's psif-------------------------------c
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
c------------------------------- end RRL's psif-------------------------------c

      elseif(ifuel==2)then            !Use Michael Clark's method of calculating psif!
                  
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
      elseif(ifuel==3)then            !Use JAS's method of calculating psif!
c  reimplemented by JAS 2/2/06 to use intrinsic erf (much faster!)

              
              !look closely at the integral erf calculates to see why here!
c      psif(i,j,k) = c1*(er*((temps(i,j,k)-tcrit)/abs(temps(i,j,k)-tcrit))+c3)  

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

c      slambdaof=rhof(i,j,k)*xv(i,j,k,7)/(rhof(i,j,k)/rnfuel+ !(rrl 6/1/99)
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
c----------------------------------RRL's psif-------------------------------c
      if(ifuel==1)then            !Use RRL's method of calculating psif to calculate psig!

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
c------------------------------- end RRL's psif-------------------------------c

      elseif(ifuel==2)then            !Use Michael Clark's method of calculating psi!
                  
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
      elseif(ifuel==3)then            !Use JAS's method of calculating psi!
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

      slambdaog=xv(i,j,k,8)*xv(i,j,k,7)/(xv(i,j,k,8)/rng+     !(jjc 9/20/01)
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
      end
