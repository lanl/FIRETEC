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

      !JAS 3/7/06 added explicit declarations to comply with implicit none
      integer,intent(in) :: il,iu,jl,ju,lls,nvp

      real xv(il:iu,jl:ju,lls,nvp)
      real,allocatable:: u(:, :,:),
     .                   v(:, :,:), 
     .                   w(:, :,:)
c     .                  ,tbar(:,: ,:)
      
      !JAS 3/7/06 added explicit declarations to comply with implicit none
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
               if (rhodirt(i,j,k).ge.2.) sizescale(i,j,k)=.3 !FP

      rhogas=xv(i,j,k,nv)

c      usp=sqrt(u(i,j,k)*u(i,j,k)+v(i,j,k)*v(i,j,k)+w(i,j,k)*w(i,j,k))
c       rktemp=xv(i,j,k,6)/xv(i,j,k,nv)+
c     +            xv(i,j,k,5)/xv(i,j,k,nv)
c      re=ss*(usp+sqrt(rktemp))/2.e-05                                   !rrl
c      h=2*0.683*re**0.466*thermcondair/ss                          !rrl p335*2
                  ufuel=xvfuel(i,j,k,1)/xv(i,j,k,nv)
                  vfuel=xvfuel(i,j,k,2)/xv(i,j,k,nv)
                  wfuel=xvfuel(i,j,k,3)/xv(i,j,k,nv)
                  rktemp=xvfuel(i,j,k,4)/xv(i,j,k,nv)+
     +            1.2*(xvfuel(i,j,k,5)/xv(i,j,k,nv))

                  usp=sqrt(rktemp)+sqrt(ufuel*ufuel
     &                         +vfuel*vfuel+wfuel*wfuel)
                  re=sizescale(i,j,k)*usp/2.e-05 !FP                             !rrl
c                  h=0.683*re**0.466*thermcondair/ss                !rrl
c                  h=2*0.683*re**0.466*thermcondair/ss        !rrl*2 8/15/01
c                  h=2*0.683*re**0.466*thermcondair/ss  
                  h=.5*0.683*re**0.466*thermcondair/sizescale(i,j,k) !FP         !4/23/03
c      h=0.683*re**0.466*thermcondair/ss                          !rrl p335
c      av=2.*(rhof(i,j,k)*0.004+rhowater(i,j,k)*.0001)/ss    !rrl
      av=2.*(rhof(i,j,k)*rrhomicro)/sizescale(i,j,k) !FP
      sizescale(i,j,k)=sstemp  !FP               !added for rhodirt
      tambientarray(i,j,k)=xe(i,j,k,4)/xe(i,j,k,nv)
     &      *(pre(i,j,k)*1.0e-5)**(rg/cp)
      convhtb=h*av*(temps(i,j,k)-tempg(i,j,k))
c      if(mpi_rank.eq.0) write (*,*) 'energysource',j,k,temps(i,j,k),
c     .   tempg(i,j,k),xv(i,j,k,4),pr(i,j,k),rg,cp

      
      

c      qxt=rke*sigma/sb(i,j,k)*(tempg(i,j,k)**4-tambient**4)  !wss
      qxt=rke*sigma/sqrt(dx*dy)*(tempg(i,j,k)**4
     &                           -tambientarray(i,j,k)**4)  !wss
      if(irad.GE.1) qxt=0.                                  !rrl !KOO eq->GE for MC
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
      rnetmass=ff(i,j,k)*tcrit*cpwood+rwatergainht   !rwaterg.. term added jjc rrl 09/21/01
                                                           !cpwood added  jjc rrl 09/21/01
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
      end
