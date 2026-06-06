      subroutine fueltempnonlocal(xv,xvfuel,il,iu,jl,ju,lls,nvp)
      use turba
      use pres
      use fireteca
      use constants
      use gridsetup
      use msga
      use nonlocal

      Implicit None

      !JAS 3/7/06 added explicit declarations to comply with implicit none
      integer,intent(in) :: il,iu,jl,ju,lls,nvp

      real xv(il:iu,jl:ju,lls,nvp) 
      real xvfuel(il:iu,jl:ju,lls,5) 

      !JAS 3/7/06 added explicit declarations to comply with implicit none
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
               !if (rhodirt(i,j,k).gt.2.) ss=.3
               if (rhodirt(i,j,k).gt.2.) sizescale(i,j,k)=.3 !FP

               if (rhof(i,j,k).gt.1.e-04) then
                  uref=xv(i,j,k,1)/xv(i,j,k,nv)
                  vref=xv(i,j,k,2)/xv(i,j,k,nv)
                  wref=xv(i,j,k,3)/xv(i,j,k,nv)
                  rktempref=xv(i,j,k,6)/xv(i,j,k,nv)+
     +            xv(i,j,k,5)/xv(i,j,k,nv)
                  spref=sqrt(rktempref)
     &                 +sqrt(uref*uref+vref*vref+wref*wref) 
                  !reref=ss*spref/2.e-05                              !rrl
                  reref=sizescale(i,j,k)*spref/2.e-05                !FP              !rrl
                  !href=2*0.683*reref**0.466*thermcondair/ss                !rrl*2 8/15/01
                  href=2*0.683*reref**0.466*thermcondair/sizescale(i,j,k)  !FP               !rrl*2 8/15/01

                  u=xvfuel(i,j,k,1)/xv(i,j,k,nv)
                  v=xvfuel(i,j,k,2)/xv(i,j,k,nv)
                  w=xvfuel(i,j,k,3)/xv(i,j,k,nv)
                  rktemp=xvfuel(i,j,k,4)/xv(i,j,k,nv)+
     +            1.2*(xvfuel(i,j,k,5)/xv(i,j,k,nv))

                  sp=sqrt(rktemp)+sqrt(u*u+v*v+w*w) 
                  !re=ss*sp/2.e-05                              !rrl
                  re=sizescale(i,j,k)*sp/2.e-05                              !rrl
                  !h=.5*0.683*re**0.466*thermcondair/ss         !rrl
                  h=.5*0.683*re**0.466*thermcondair/sizescale(i,j,k) !FP         !rrl
c                  h=2*0.683*re**0.466*thermcondair/ss         !rrl*2 8/15/01
c                  h=2*0.683*re**0.466*thermcondair/ss         !rrl*2 8/15/01
c      av=2.*(rhof(i,j,k)*0.004+rhowater(i,j,k)*.0001)/ss      !rrl
c      av=2.*(rhof(i,j,k)*.004+rhodirt(i,j,k)*.0008)/ss      !rrl
      av=2.*(rhof(i,j,k)*rrhomicro)/sizescale(i,j,k)      !rrl
      sizescale(i,j,k)=sstemp
c      tg=xv(i,j,k,4)/xv(i,j,k,nv)*(pr(i,j,k)*1.0e5)**(rg/cp)  !rrl
      convht(i,j,k)=h*av*(tempg(i,j,k)-temps(i,j,k))    !wss
c      solidqxt=av*sigma*solidemisivity*(temps(i,j,k))
      if (xv(i,j,k,8).ne.0.d0) then
        reactht=(thetasolid(i,j,k)*ff(i,j,k)*dtp*rnhc/xv(i,j,k,8))
     &            *hfgas*fg(i,j,k)+hfsolid*ff(i,j,k)                  !jjc
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
      end
