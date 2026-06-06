c2345678***************************************************
       module higrad 
       Implicit None
       
       save

       contains
!*****************************************************************************************!
      subroutine rmoainner()
      use gridsetup
      use xvi
      use xve
      use forcinner
      use workavg
      use metryic
      use pres
      use constants
      use relax
      use weights
      use msga

      Implicit None

      !JAS 3/6/06 added explicit declarations to comply with implicit none
      integer :: its,i,j,k,kv,nvs !,ia,ja
!      real :: relt,tauv
! FIXME
      real*8 :: pr_double

      real,allocatable:: xvt(:,:,:,:)
      real,allocatable:: uab(:,:,:)
      real,allocatable:: vab(:,:,:)
      real,allocatable:: oab(:,:,:)
      real,allocatable:: ua(:,:,:)
      real,allocatable:: va(:,:,:)
      real,allocatable:: oa(:,:,:)
      real,allocatable:: u(:,:,:,:),
     .                   v(:,:,:,:),
     .                   o(:,:,:,:)
      integer,save::icounter=0
      real :: ug,vg ! geostrophic wind
      allocate (xvt(1-ih:np+ih,1-ih:mp+ih,l,3))
      allocate (uab(1-ih:np+ih+1,1-ih:mp+ih,l))
      allocate (vab(1-ih:np+ih,1-ih:mp+ih+1,l))
      allocate (oab(1-ih:np+ih,1-ih:mp+ih,l+1))
      allocate (ua(1-ih:np+ih+1,1-ih:mp+ih,l))
      allocate (va(1-ih:np+ih,1-ih:mp+ih+1,l))
      allocate (oa(1-ih:np+ih,1-ih:mp+ih,l+1))
      allocate (u(1-ih:np+ih, 1-ih:mp+ih,l,0:2))
      allocate (v(1-ih:np+ih, 1-ih:mp+ih,l,0:2))
      allocate (o(1-ih:np+ih, 1-ih:mp+ih,l,0:2))

      do its=1,nts
       time=time+dts
       icounter=icounter+1
       if(its.gt.1)then
        do k=1,l
         do j=1,mp
          do i=1,np
           u(i,j,k,0)=u(i,j,k,1)
           v(i,j,k,0)=v(i,j,k,1)
           o(i,j,k,0)=o(i,j,k,1)
          enddo
         enddo
        enddo
       else
        f1=0.0
        f2=0.0
        f3=0.0
        f1avg=0.0
        f2avg=0.0
        f3avg=0.0
        do k=1,l
         do j=1,mp
          do i=1,np
           u(i,j,k,0)=xv(i,j,k,1)/xv(i,j,k,5)
           v(i,j,k,0)=xv(i,j,k,2)/xv(i,j,k,5)
           o(i,j,k,0)=(xv(i,j,k,3)*gi(i,j,k)
     &               +c13(i,j)*gmul(k)*xv(i,j,k,1)
     &               +c23(i,j)*gmul(k)*xv(i,j,k,2))/xv(i,j,k,5)


          enddo
         enddo
        enddo
        if(isplit.eq.1)then
         call velprdmoa(u,v,o,f1,f2,f3,uab,vab,oab,gc1,gc2,gc3,1-ih,np+ih,1-ih,mp+ih,l,0,2)
        endif !if(isplit.eq.1)
       endif !if(its.gt.1)
       do k=1,l
        do j=1,mp
         do i=1,np
          u(i,j,k,1)=xv(i,j,k,1)/xv(i,j,k,5)
          v(i,j,k,1)=xv(i,j,k,2)/xv(i,j,k,5)
          o(i,j,k,1)=(xv(i,j,k,3)*gi(i,j,k)
     &               +c13(i,j)*gmul(k)*xv(i,j,k,1)
     &               +c23(i,j)*gmul(k)*xv(i,j,k,2))/xv(i,j,k,5)
         enddo
        enddo
       enddo
       call velprdf(u,v,o,ua,va,oa,gc1s,gc2s,gc3s,1-ih,np+ih,1-ih,mp+ih,l,0,2)

       if(islip.eq.1) then
        do j=1,mp
         do i=1,np
          oa(i,j,2)=0.
          oab(i,j,2)=0.
         enddo
        enddo
       endif !(islip.eq.1)

       if(isplit.eq.1)then
        if(its.eq.1)then
         do kv=1,3
          do k=1,l
           do j=1,mp
            do i=1,np
             xvt(i,j,k,kv)=xv(i,j,k,kv)
            enddo
           enddo
          enddo
         enddo
         call donorcell(uab,vab,oab,xvt(1-ih,1-ih,1,1),
     .                  xe(1-ih,1-ih,1,1),1-ih,np+ih,1-ih,mp+ih,l)
         call donorcell(uab,vab,oab,xvt(1-ih,1-ih,1,2),
     .                  xe(1-ih,1-ih,1,2),1-ih,np+ih,1-ih,mp+ih,l)
         call donorcell(uab,vab,oab,xvt(1-ih,1-ih,1,3),
     .                  xe(1-ih,1-ih,1,3),1-ih,np+ih,1-ih,mp+ih,l)
        endif !if(its.eq.1)
       else !else (isplit.ne.1)
        call donorcell(ua,va,oa,xv(1-ih,1-ih,1,1),
     .                 xe(1-ih,1-ih,1,1),1-ih,np+ih,1-ih,mp+ih,l)
        call donorcell(ua,va,oa,xv(1-ih,1-ih,1,2),
     .                 xe(1-ih,1-ih,1,2),1-ih,np+ih,1-ih,mp+ih,l)
        call donorcell(ua,va,oa,xv(1-ih,1-ih,1,3),
     .                 xe(1-ih,1-ih,1,3),1-ih,np+ih,1-ih,mp+ih,l)
       endif !if(isplit.eq.1)

       if(isplit.eq.1) then
        if(its.eq.1) then
         do kv=1,3
          do k=1,l
           do j=1,mp
            do i=1,np
             xvt(i,j,k,kv)=(xvt(i,j,k,kv)-xv(i,j,k,kv))/real(nts)
            enddo
           enddo
          enddo
         enddo
        endif !if(is.eq.1)
        do kv=1,3
         do k=1,l
          do j=1,mp
           do i=1,np
            xv(i,j,k,kv)=xv(i,j,k,kv)+xvt(i,j,k,kv)
           enddo
          enddo
         enddo
        enddo
       endif !if(isplit.eq.1)

       call donorcell(ua,va,oa,xv(1-ih,1-ih,1,5),
     .                xe(1-ih,1-ih,1,nv),1-ih,np+ih,1-ih,mp+ih,l)
       call donorcell(ua,va,oa,xv(1-ih,1-ih,1,4),
     .                xe(1-ih,1-ih,1,4),1-ih,np+ih,1-ih,mp+ih,l)

        do k=1,l
         do j=1,mp
          do i=1,np
!           pr(i,j,k)=(xv(i,j,k,4)*rg/prrcp)**(cp/cv)-pre(i,j,k)
           pr_double=(xv(i,j,k,4)*rg/prrcp)**(cp/cv)-pre(i,j,k)
           pr(i,j,k)=real(pr_double)

          enddo
         enddo
        enddo

!       if(its.eq.nts)then
!        call rmaxmin1(pr,'pr',1-ih,np+ih,1-ih,mp+ih,l)
!       endif

       call gradmoa(gc1s,gc2s,gc3s)


       if(its.eq.nts)then
        call rmaxmin1(pr,'pr',1-ih,np+ih,1-ih,mp+ih,l)
       endif 

       do k=1,l
        do j=1,mp
         do i=1,np
c gravity effects 
          f1(i,j,k)=f1(i,j,k)-dts*(xv(i,j,k,5)-xe(i,j,k,nv))*g
     +               *sin(slopeangle*3.14159/180.) 
     +               *cos(slopeazimuth*3.14159/180.) 
          f2(i,j,k)=f2(i,j,k)-dts*(xv(i,j,k,5)-xe(i,j,k,nv))*g
     +               *sin(slopeangle*3.14159/180.) 
     +               *sin(slopeazimuth*3.14159/180.) 
          f3(i,j,k)=f3(i,j,k)-dts*(xv(i,j,k,5)-xe(i,j,k,nv))*g
     +               *cos(slopeangle*3.14159/180.) 
         enddo
        enddo
       enddo

       if(icorio.ge.1)then
        do k=1,l
         do j=1,mp
          do i=1,np
          if (icorio.eq.1) then
           ug=xe(i,j,k,1)
           vg=xe(i,j,k,2)
          else if (icorio.eq.2) then!icorio.eq.2
           ug=xe(i,j,l,1)
           vg=xe(i,j,l,2)
          end if
          ! FP switches dt to dts
           f1(i,j,k)= dts*(fcor3*(xv(i,j,k,2)-vg)
     &                    -fcor2*(xv(i,j,k,3)) )+f1(i,j,k)
           f2(i,j,k)=-dts*(fcor3*(xv(i,j,k,1)-ug) )+f2(i,j,k)
           f3(i,j,k)= dts*(fcor2*(xv(i,j,k,1)-ug))+f3(i,j,k)
          enddo
         enddo
        enddo
       end if !icorio.ge.1
       if (ilspgf.ge.1) then !ilspgf
        do k=1,l
         do j=1,mp
          do i=1,np
           ug=xe(i,j,l,1)
           vg=xe(i,j,l,2)
          ! FP switches dt to dts
           f1(i,j,k)= dts*fcor3*ug*sintheta(i,j,k)+f1(i,j,k)
           f2(i,j,k)= dts*fcor3*vg*sintheta(i,j,k)+f2(i,j,k)
          !write(6,*) 'corio ',i,j,k,sintheta(i,j,k),ug,vg,fcor3
          enddo
         enddo
        enddo
       endif !ilspgf.ge.1 
       do k=1,l
         do j=1,mp
          do i=1,np
           xv(i,j,k,1)=xv(i,j,k,1)+f1(i,j,k)
           xv(i,j,k,2)=xv(i,j,k,2)+f2(i,j,k)
           xv(i,j,k,3)=xv(i,j,k,3)+f3(i,j,k)
          enddo
         enddo
        enddo

       ! HERE FP DIVIDED RELAXXV BY nts BECAUSE RELAXATION IS DONE MORE OFTEN
       do nvs=1,5
        if(nvs.ne.5) then
         do k=1,l
          do j=1,mp
           do i=1,np
            xv(i,j,k,nvs)=xv(i,j,k,nvs)*(1-relaxxv(i,j,k,min(nvs,4))/nts)+
     &                    xe(i,j,k,nvs)*relaxxv(i,j,k,min(nvs,4))/nts
           enddo
          enddo
         enddo
        else ! case of nv stuff FP
         do k=1,l
          do j=1,mp
           do i=1,np
            xv(i,j,k,5)=xv(i,j,k,5)*(1-relaxxv(i,j,k,min(nvs,4))/nts)+
     &      xe(i,j,k,nv)*relaxxv(i,j,k,min(nvs,4))/nts
           enddo
          enddo
         enddo
        endif
       enddo !do nvs=1,5

       if(ibctopbot.eq.0)then
        do j=1,mp
         do i=1,np
          if(islip.eq.1) xv(i,j,1,1)=0.
          if(islip.eq.1) xv(i,j,1,2)=0.
          xv(i,j,1,3)=-(xv(i,j,1,1)*c13(i,j)*gmul(1)+
     &                  xv(i,j,1,2)*c23(i,j)*gmul(1))/gi(i,j,1)
          xv(i,j,l,3)=0.
         enddo
        enddo
       endif !if(ibctopbot.eq.0)

       do k=1,l
        do j=1,mp
         do i=1,np
          uavg(i,j,k)=uavg(i,j,k)+wt(its+1)*xv(i,j,k,1)/xv(i,j,k,5)
          vavg(i,j,k)=vavg(i,j,k)+wt(its+1)*xv(i,j,k,2)/xv(i,j,k,5)
          wavg(i,j,k)=wavg(i,j,k)+wt(its+1)*xv(i,j,k,3)/xv(i,j,k,5)
          f1avg(i,j,k)=f1avg(i,j,k)+wtf(its+1)*f1(i,j,k)
          f2avg(i,j,k)=f2avg(i,j,k)+wtf(its+1)*f2(i,j,k)
          f3avg(i,j,k)=f3avg(i,j,k)+wtf(its+1)*f3(i,j,k)
         enddo
        enddo
       enddo

      enddo !do its=1,nts

      do k=1,l
       do j=1,mp
        do i=1,np
         uavg(i,j,k)=(uavg(i,j,k)/real(nts))/gi(i,j,k)
         vavg(i,j,k)=(vavg(i,j,k)/real(nts))/gi(i,j,k)
         wavg(i,j,k)=(wavg(i,j,k)/real(nts))/gi(i,j,k)
         oavg(i,j,k)=uavg(i,j,k)*c13(i,j)*gmul(k)+
     &               vavg(i,j,k)*c23(i,j)*gmul(k)+
     &               wavg(i,j,k)*gi(i,j,k)
        enddo
       enddo
      enddo

      deallocate   (xvt)
      deallocate   (uab)
      deallocate   (vab)
      deallocate   (oab)
      deallocate   (ua)
      deallocate   (va)
      deallocate   (oa)
      deallocate   (u)
      deallocate   (v)
      deallocate   (o)

      return
      end subroutine rmoainner
!*****************************************************************************************!
      subroutine advvel
      use gridsetup
      use advo
      use workavg
      use msga
      use xve   !this use staement allows access to the environmental array
      use metryic !this use statement allows access to the gi array

      Implicit None

      !JAS 3/7/06 added explicit declarations to comply with implicit none
      integer :: i,j,k

      call updated(uavg,uavg,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      call updated(vavg,vavg,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)

      !compute covariant cell face velocites in x direction
      do k=1,l
       do j=1,mp
        do i=1,np+1                         !added d rrl
         uab(i,j,k)=0.5*(uavg(i,j,k)+uavg(i-1,j,k))*gc1
        enddo
       enddo
      enddo
      if(leftdedge.eq.1)then !if noncyclic bcs set domain edge to environmental value
       do k=1,l
        do j=1,mp
      ! FP below it is not correct because xe is cell centered,
      ! but ze don't know (xe(0) so a constant gradient is assumed
         uab(1,j,k)=xe(1,j,k,1)/xe(1,j,k,nv)*gc1/gi(i,j,k)
c         uab(1,j,k)=0.5*gc1*(
c     +    (xe(1,j,k,1)/xe(1,j,k,nv))/gi(1,j,k)
c     +    +(xe(0,j,k,1)/xe(0,j,k,nv))/gi(0,j,k))
         !uab(1,j,k)=uab(2,j,k)
        enddo
       enddo
      endif
      if(rightdedge.eq.1)then !if noncyclic bcs set domain edge to environmental value
       do k=1,l
        do j=1,mp
         uab(np+1,j,k)=xe(np,j,k,1)/xe(np,j,k,nv)*gc1/gi(i,j,k)
!            uab(np+1,j,k)=0.5*gc1*(
!     +     (xe(np+1,j,k,1)/xe(np+1,j,k,nv))/gi(np+1,j,k)
!     +     +(xe(np,j,k,1)/xe(np,j,k,nv))/gi(np,j,k))
         if (ibclatopen.eq.1.) uab(np+1,j,k)=max(0.0,uab(np,j,k))
         if (ibclatopen.eq.2.) uab(np+1,j,k)=uab(np,j,k)
        enddo
       enddo
      endif

      !compute covariant cell face velocites in y direction
      if(j3.ne.0) then
       do k=1,l
        do j=1,mp+1
         do i=1,np
          vab(i,j,k)=0.5*(vavg(i,j,k)+vavg(i,j-1,k))*gc2
         enddo
        enddo
       enddo
       if(botdedge.eq.1)then !if noncyclic bcs set domain edge to environmental value
        do k=1,l
         do i=1,np
          vab(i,1,k)=xe(i,1,k,2)/xe(i,1,k,nv)*gc2/gi(i,j,k)
!          vab(i,1,k)=0.5*gc2*(
!     +         (xe(i,1,k,2)/xe(i,1,k,nv))/gi(i,1,k)
!     +         +(xe(i,0,k,2)/xe(i,0,k,nv))/gi(i,0,k))
         if (ibclatopen.eq.1.) vab(i,1,k)=min(0.0,vab(i,2,k))
         if (ibclatopen.eq.2.) vab(i,1,k)=vab(i,2,k)
          !vab(i,1,k)=vab(i,2,k)
         enddo
        enddo
       endif
       if(topdedge.eq.1)then !if noncyclic bcs set domain edge to environmental value
        do k=1,l
         do i=1,np
          vab(i,mp+1,k)=xe(i,mp,k,2)/xe(i,mp,k,nv)*gc2/gi(i,j,k)
!          vab(i,mp+1,k)=0.5*gc2*(
!     +            (xe(i,mp+1,k,2)/xe(i,mp+1,k,nv))/gi(i,mp+1,k)
!     +            +(xe(i,mp,k,2)/xe(i,mp,k,nv))/gi(i,mp,k))
         if (ibclatopen.eq.1.) vab(i,mp+1,k)=max(0.0,vab(i,mp,k))
         if (ibclatopen.eq.2.) vab(i,mp+1,k)=vab(i,mp,k)
          !vab(i,mp+1,k)=vab(i,mp,k)
         enddo
        enddo
       endif
      else
      vab=0.0
      endif

      !compute covariant cell face velocites in z direction
      do k=2,l
       do j=1,mp
        do i=1,np
          oab(i,j,k)=0.5*(oavg(i,j,k)+oavg(i,j,k-1))*gc3
        enddo
       enddo
      enddo
      if(ibctopbot.eq.0) then
       do j=1,mp
        do i=1,np
         !oab(i,j,1)=xe(i,j,1,3)/xe(i,mp,1,5)*gc3
         !oab(i,j,l+1)=xe(i,j,l,3)/xe(i,mp,l,5)*gc3
         oab(i,j,1)=oab(i,j,2)
         oab(i,j,l+1)=oab(i,j,l)
        enddo
       enddo
      else
       oab(:,:,1)=0.0
       oab(:,:,l+1)=0.0
      endif

      return
      end subroutine advvel
!*****************************************************************************************!
      subroutine advec(xv,il,iu,jl,ju,lls,nvp)
      use gridsetup
      use advo
      use xve
      use metryic
      use msga

      Implicit None


      !JAS 3/7/06 added explicit declarations to comply with implicit none
      integer,intent(in) :: il,iu,jl,ju,lls,nvp

      real xv(il:iu,jl:ju,lls,nvp)
      real, allocatable::xv2d(:,:,:,:),
     .                   xe2d(:,:,:,:)

      !JAS 3/7/06 added explicit declarations to comply with implicit none
      integer :: i,k,ila,iua,jla,jua,lla,nvpa,kv!,ierr

      if(j3.eq.0) then
        allocate (xv2d(1-ih:np+ih,1-ih:mp+ih,l,nv2d))
        allocate (xe2d(1-ih:np+ih,1-ih:mp+ih,l,nv2d))
      endif ! for 2d

      ila=1
      iua=1
      jla=1
      jua=1
      lla=1
      nvpa=1
      if(j3.eq.1) then
       ila=il
       iua=iu
       jla=jl
       jua=ju
       lla=lls
       nvpa=nvp
c walltimer !KOO
       if(iwallclock.EQ.1) wtime2=MPI_Wtime()

       call mpdatanew3d(uab,vab,oab,xv,h,xe,
     .ila,iua,jla,jua,lla,nvpa)
  
      if(iwallclock.EQ.1) call computeWallTime(wtime2,11,'mpdata3dnew')

      else ! 3d
       xv2d(1:np,1,1:l,1)=xv(1:np,1,1:l,1)
       xe2d(1:np,1,1:l,1)=xe(1:np,1,1:l,1)
       do kv=2,nv2d
        do k=1,l
         do i=1,np
          xv2d(i,1,k,kv)=xv(i,1,k,kv+1)
          xe2d(i,1,k,kv)=xe(i,1,k,kv+1)
         enddo
        enddo
       enddo
       ila=il
       iua=iu
       jla=jl
       jua=ju
       lla=lls
 
       call mpdatanew2d(uab,oab,xv2d,h,xe2d,
     .ila,iua,jla,jua,lla,nv2d)

       xv(1:np,1,1:l,1)=xv2d(1:np,1,1:l,1)
       do kv=2,nv-1
        do k=1,l
         do i=1,np
          xv(i,1,k,kv+1)=xv2d(i,1,k,kv)
         enddo
        enddo
       enddo
      endif ! 2d

      if(j3.eq.0) deallocate (xv2d) ! 2d
      if(j3.eq.0) deallocate (xe2d) ! 2d

      return
      end subroutine advec
!*****************************************************************************************!
      subroutine gradmoa(gcx,gcy,gcz)
      use gridsetup
      use pres
      use forcinner
      use metryic
      use msga
      use xvi

      Implicit None

      !JAS 3/7/06 added explict declarations to comply with implicit none
      real :: gcx,gcy,gcz
      real :: ghx,ghy,ghz,g13,g23,g33
      integer :: i,j,k

      real,allocatable:: px(:, :,:),
     .                   py(:, :,:),
     .                   pz(:, :,:)

      !JAS not used --real,pointer::xjac0(:,:)
      ghx=0.5*gcx
      ghy=0.5*gcy
      ghz=0.5*gcz
      allocate (px(1-ih:np+ih,1-ih:mp+ih,l))
      allocate (py(1-ih:np+ih,1-ih:mp+ih,l))
      allocate (pz(1-ih:np+ih,1-ih:mp+ih,l))


      call updated(pr,pr,np,mp,l,1-ih,np+ih,1-ih,mp+ih,0)
c in order to get rid of leftedge and rightedge we need to make sure pressure is updated correctly
      do i=1+leftdedge,np-rightdedge                       !add d rrl
      do k=1,l
      do j=1,mp
       px(i,j,k)=    ghx*(pr(i+1,j,k)-pr(i-1,j,k))
      enddo
      enddo
      enddo 

      if(leftdedge.eq.1) then                       !add d rrl
      do k=1,l
      do j=1,mp
      px(1,j,k)=(1-ibcx)*gcx*(pr(2,j,k)-pr(1,j,k))
     1            +ibcx*ghx*(pr(2,j,k)-pr(-1,j,k))
      enddo
      enddo
      endif
      if(rightdedge.eq.1) then                       !add d rrl
      do k=1,l
      do j=1,mp
      px(np,j,k)=(1-ibcx)*gcx*(pr(np,j,k)-pr(np-1,j,k))
     1              +ibcx*ghx*(pr(np+2,j,k)-pr(np-1,j,k))
      enddo
      enddo
      endif

      if(j3.ne.0) then
      do k=1,l
      do i=1,np
      do j=1+botdedge,mp-topdedge                       !add d rrl
   28 py(i,j,k)=     ghy*(pr(i,j+1,k)-pr(i,j-1,k))
      enddo
      enddo
      enddo
      if(botdedge.eq.1) then                       !add d rrl
      do k=1,l
      do i=1,np
      py(i,1,k)=(1-ibcy)*gcy*(pr(i,2,k)-pr(i,1,k))
     1            +ibcy*ghy*(pr(i,2,k)-pr(i,-1,k))
      enddo
      enddo
      endif
      if(topdedge.eq.1) then                       !add d rrl
      do k=1,l
      do i=1,np
      py(i,mp,k)=(1-ibcy)*gcy*(pr(i,mp,k)-pr(i,mp-1,k))
     1            +ibcy*ghy*(pr(i,mp+2,k)-pr(i,mp-1,k))
      enddo
      enddo
      endif
      else
      do k=1,l
      do i=1,np
      py(i,1,k)=0.
      enddo
      enddo
      endif
      do 3 k=2,l-1
      do 3 j=1,mp
      do 3 i=1,np
    3 pz(i,j,k)=ghz*(pr(i,j,k+1)-pr(i,j,k-1))
      do 38 j=1,mp
      do 38 i=1,np
c     
c      pz(i,j,1)= gcz*(pr(i,j,2)-pr(i,j,1))     
      pz(i,j,1)= ghz*(pr(i,j,2)-pr(i,j,1))     !modified by rrl 10/02/01
   38 pz(i,j,l)= ghz*(pr(i,j,l)-pr(i,j,l-1))
c top and bottom b.c.
      if(islip.eq.1) xv(:,:,1,1)=0.
      if(islip.eq.1) xv(:,:,1,2)=0.
      if(islip.eq.1) xv(:,:,1,3)=0.
      if(ibctopbot.eq.0) then
      do k=1,l,l-1
      do j=1,mp
      do i=1,np
        g13=c13(i,j)*gmul(k)
        g23=c23(i,j)*gmul(k)
        g33= gi(i,j,k)
        pz(i,j,k)=(g13*xv(i,j,k,1)+g23*xv(i,j,k,2)+g33*xv(i,j,k,3)-
     .             g13*px(i,j,k)-g23*py(i,j,k))
     .           /(g33*g33+ g13*g13+g23*g23)
      enddo
      enddo
      enddo
      endif

      do 10 k=1,l
      do 10 j=1,mp
      do 10 i=1,np
        g13=c13(i,j)*gmul(k)
        g23=c23(i,j)*gmul(k)
        g33= gi(i,j,k)

        f1(i,j,k)=-(    px(i,j,k)               +g13*pz(i,j,k))
        f2(i,j,k)=-(                   py(i,j,k)+g23*pz(i,j,k))
        f3(i,j,k)=-(g33*pz(i,j,k))

   10 continue

      deallocate (px)
      deallocate (py)
      deallocate (pz)
      return
      end subroutine gradmoa

!*****************************************************************************************!
      subroutine donorcell(u1,u2,u3,xf,xe,il,iu,jl,ju,lls)
      use metryic
      use gridsetup
      use msga

      Implicit None

      !JAS 3/7/06 added explicit declarations to comply with implicit none
      integer,intent(in) :: il,iu,jl,ju,lls

      real u1(il:iu+1,jl:ju,lls),
     .     u2(il:iu,jl:ju+1,lls),
     .     u3(il:iu,jl:ju,lls+1)
      real xf(il:iu,jl:ju,lls),
     .     xe(il:iu,jl:ju,lls)
      real,allocatable:: f1(:,:,:),
     .                   f2(:,:,:),
     .                   f3(:,:,:)
      !JAS 3/7/06 added explicit declarations to comply with implicit none
      integer :: i,j,k

      allocate (f1(1-ih:np+ih+1,1-ih:mp+ih,l))
      allocate (f2(1-ih:np+ih,1-ih:mp+ih+1,l))
      allocate (f3(1-ih:np+ih,1-ih:mp+ih,l+1))
      f1=0.0
      f2=0.0
      f3=0.0

c      write (*,*) 'inside donorcell'
      call updated(xf,xe,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      !call updated(u1,u1,np+1,mp,l,1-ih,np+ih+1,1-ih,mp+ih,1)
      !call updated(u2,u2,np,mp+1,l,1-ih,np+ih,1-ih,mp+ih+1,1)
      !call updated(u3,u3,np,mp,l+1,1-ih,np+ih,1-ih,mp+ih,1)

      do k=1,l
       do j=1,mp
        do i=1,np+1
         f1(i,j,k)=donor(xf(i-1,j,k),xf(i,j,k),u1(i,j,k))
        enddo
       enddo
      enddo
c      write (*,*) 'after left and right edge'

      if(j3.eq.1) then
       do k=1,l
        do j=1,mp+1
         do i=1,np
          f2(i,j,k)=donor(xf(i,j-1,k),xf(i,j,k),u2(i,j,k))
         enddo
        enddo
       enddo
      else
       do k=1,l
        do i=1,np
         f2(i,1,k)=0.0
         f2(i,2,k)=0.0
        enddo
       enddo
      endif
c      write (*,*) 'after top and bottom edge'

      do k=2,l
       do j=1,mp
        do i=1,np
         f3(i,j,k)=donor(xf(i,j,k-1),xf(i,j,k),u3(i,j,k))
        enddo
       enddo
      enddo
c      write (*,*) 'before ibctopbot'
      if(ibctopbot.eq.0) then
       do j=1,mp
        do i=1,np
         f3(i,j, 1)=-f3(i,j,2)
         f3(i,j,l+1)=-f3(i,j,l)
        enddo
       enddo
      else
       do j=1,mp
        do i=1,np
         f3(i,j, 1)=0.
         f3(i,j,l+1)=0.
        enddo
       enddo
      endif
c      write (*,*) 'after ibctopbot'
c      TODO : check if it should be updated or update with rrl
!
         call updated(f1,f1,np,mp,l,1-ih,np+ih+1,1-ih,mp+ih,2)
         call updated(f2,f2,np,mp,l,1-ih,np+ih,1-ih,mp+ih+1,3)

c      write (*,*) 'before xf'
      do k=1,l
       do j=1,mp
        do i=1,np
         xf(i,j,k)=xf(i,j,k)-( f1(i+1,j,k)-f1(i,j,k)
     &                        +f2(i,j+1,k)-f2(i,j,k)
     &                        +f3(i,j,k+1)-f3(i,j,k) )/h(i,j,k)
        enddo
       enddo
      enddo
c      write (*,*) 'after xf'

c      FP : this should be update"d" because cell centered 

      call updated(xf,xf,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)

      deallocate (f1)
      deallocate (f2)
      deallocate (f3)

      return
      end subroutine donorcell
!*****************************************************************************************!
      subroutine mpdatanew3d(u1,u2,u3,x,h,xe,il,iu,jl,ju,lls,nvp)
      use gridsetup
      use msga

      Implicit None

      !JAS 3/7/06 added explicit declarations to comply with implicit none
      integer,intent(in) :: il,iu,jl,ju,lls,nvp

      real,dimension(il:iu+1,jl:ju,lls) :: u1
      real,dimension(il:iu,jl:ju+1,lls) :: u2
      real,dimension(il:iu,jl:ju,lls+1) :: u3
      real,dimension(il:iu,jl:ju,lls,nvp) :: x,xe
      real,dimension(il:iu,jl:ju,lls) :: h
      real,allocatable:: v1(:,:,:,:)
      real,allocatable:: v2(:,:,:,:)
      real,allocatable:: v3(:,:,:,:)
      real,allocatable:: f1(:,:,:,:)
      real,allocatable:: f2(:,:,:,:)
      real,allocatable:: f3(:,:,:,:)
      real,allocatable:: f1o(:,:,:,:)
      real,allocatable:: f2o(:,:,:,:)
      real,allocatable:: f3o(:,:,:,:)
      real,allocatable:: cp(:,:,:,:)
      real,allocatable:: cn(:,:,:,:)
      real,allocatable:: mxo(:,:,:,:)
      real,allocatable:: mno(:,:,:,:)
      real,allocatable::  a(:,:,:,:)
      real,allocatable:: mx(:,:,:,:)
      real,allocatable:: mn(:,:,:,:)

      !JAS 3/7/06 added explicit declarations to comply with implicit none
      integer :: n3,n3m,kv,i,j,k,ip,im,jp,jm,kp,km,itrfct
      integer :: itr,ilft,illim,iulim,jllim,julim,ia,ja
      real :: ep,rhoin,rhoout,ain,aout,rmxuse,rmnuse,a1p,a2p,a1n,a2n
      real :: tmpp,tmpn,tmp,v1d,v2d,c1,c2

      n3=l+1
      n3m=n3-1
      ep=1.e-10

      if(j3.eq.1) then
      allocate (v1(1-ih:np+ih+1,1-ih:mp+ih, l,nv))
      allocate (v2(1-ih:np+ih,1-ih:mp+ih+1, l,nv))
      allocate (v3(1-ih:np+ih,1-ih:mp+ih, l+1,nv))
      allocate (f1(1-ih:np+ih+1,1-ih:mp+ih, l,nv))
      allocate (f2(1-ih:np+ih,1-ih:mp+ih+1, l,nv))
      allocate (f3(1-ih:np+ih,1-ih:mp+ih, l+1,nv))
      allocate (f1o(1-ih:np+ih+1,1-ih:mp+ih, l,nv))
      allocate (f2o(1-ih:np+ih,1-ih:mp+ih+1, l,nv))
      allocate (f3o(1-ih:np+ih,1-ih:mp+ih, l+1,nv))
      allocate (cp(1-ih:np+ih, 1-ih:mp+ih, l,nv))
      allocate (cn(1-ih:np+ih, 1-ih:mp+ih, l,nv))
      allocate (mxo(1-ih:np+ih, 1-ih:mp+ih, l,nv))
      allocate (mno(1-ih:np+ih, 1-ih:mp+ih, l,nv))
      allocate (a(1-ih:np+ih, 1-ih:mp+ih, l,nv-1))
      allocate (mx(1-ih:np+ih, 1-ih:mp+ih, l,nv-1))
      allocate (mn(1-ih:np+ih, 1-ih:mp+ih, l,nv-1))
      endif
c
c
c     transfer data from shared arrays to msg arrays
c
c updated (array,dim in x, dim in y, dim in z, ngp in x, ngp in y,
c0-east west and north-south, 1 all directions, 2 east-west only,
c3-north-south only

      do kv=1,nv
      call updated(x(1-ih,1-ih,1,kv),x(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih,1)

      do k=1,n3m+1
       do j=1,mp
        do i=1,np
            v3(i,j,k,kv) = u3(i,j,k)
          enddo
        enddo
      end do

      do k=1,n3m
        do j=1,mp
          do i=1,np+1
            v1(i,j,k,kv) = u1(i,j,k)
          end do
        end do

        do i=1,np
          do j=1,mp+1
            v2(i,j,k,kv) = u2(i,j,k)
          end do
        end do
      end do
      enddo

      if(nonos.eq.1) then
      do kv=1,nv-1
      do k=1,l
      do j=1,mp
      do i=1,np
      a(i,j,k,kv)=x(i,j,k,kv)/x(i,j,k,nv)
      enddo
      enddo
      enddo
      call updated(a(1-ih,1-ih,1,kv),a(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih,1)

         do k=1,n3m
            km=max0(k-1,1  )
            kp=min0(k+1,n3m)
            do j=1,mp
               if (botdedge.eq.1 .and. j.eq.1) then                       !add d rrl
                  jm = ibcy*(-1) + ibyo*1
               else
                  jm = j - 1
               end if
               if (topdedge.eq.1 .and. j.eq.mp) then                       !add d rrl
                  jp = ibcy*(mp+2) + ibyo*mp
               else
                  jp = j + 1
               end if
               do i=1,np
                  if (leftdedge.eq.1 .and. i.eq.1) then                       !add d rrl
                     im = ibcx*(-1) + ibxo*1
                  else
                     im = i - 1
                  end if
                  if (rightdedge.eq.1 .and. i.eq.np) then                       !add d rrl
                     ip = ibcx*(np+2) + ibxo*np
                  else
                     ip = i + 1
                  end if
           mx(i,j,k,kv)=amax1(a(im,j,k,kv),a(i,j,k,kv),a(ip,j,k,kv),
     .      a(i,jm,k,kv),a(i,jp,k,kv),a(i,j,kp,kv),a(i,j,km,kv))
           mn(i,j,k,kv)=amin1(a(im,j,k,kv),a(i,j,k,kv),a(ip,j,k,kv),
     .      a(i,jm,k,kv),a(i,jp,k,kv),a(i,j,kp,kv),a(i,j,km,kv))
               end do
            end do
         end do
      enddo
      end if


      if(nonosold.eq.1) then
         do kv=1,nv
         do k=1,n3m
            km=max0(k-1,1  )
            kp=min0(k+1,n3m)
            do j=1,mp
               if (botdedge.eq.1 .and. j.eq.1) then                       !add d rrl
                  jm = ibcy*(-1) + ibyo*1
               else
                  jm = j - 1
               end if
               if (topdedge.eq.1 .and. j.eq.mp) then                       !add d rrl
                  jp = ibcy*(mp+2) + ibyo*mp
               else
                  jp = j + 1
               end if
               do i=1,np
                  if (leftdedge.eq.1 .and. i.eq.1) then                       !add d rrl
                     im = ibcx*(-1) + ibxo*1
                  else
                     im = i - 1
                  end if
                  if (rightdedge.eq.1 .and. i.eq.np) then                       !add d rrl
                     ip = ibcx*(np+2) + ibxo*np
                  else
                     ip = i + 1
                  end if
           mxo(i,j,k,kv)=amax1(x(im,j,k,kv),x(i,j,k,kv),x(ip,j,k,kv),
     .             x(i,jm,k,kv),x(i,jp,k,kv),x(i,j,kp,kv),x(i,j,km,kv))
           mno(i,j,k,kv)=amin1(x(im,j,k,kv),x(i,j,k,kv),x(ip,j,k,kv),
     .             x(i,jm,k,kv),x(i,jp,k,kv),x(i,j,kp,kv),x(i,j,km,kv))
               end do
            end do
         end do
      enddo
      end if


      c1=1.
      c2=0.
                         do 30 itr=1,iord
      do kv=1,nv
c     if(itr.eq.2) then
c     do k=1,n3m
c     do j=1,mp
c     do i=1,np+1
c        v1(i,j,k,kv)=0.
c     end do
c     end do
c     end do
c     endif

      call updated(x(1-ih,1-ih,1,kv),x(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih,1)

       ilft=1+leftdedge*1                       !add d rrl
         do k=1,n3m
            do j=1,mp
               do i=ilft,np
                  f1(i,j,k,kv)=donor(c1*x(i-1,j,k,kv)+c2,
     .                               c1*x(i,j,k,kv)+c2,
     .                 v1(i,j,k,kv))
               end do
            end do
         end do


       if (rightdedge.eq.0) then                       !add d rrl
          call updated(f1(1-ih,1-ih,1,kv),f1(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,2)
       else
          call updated(f1(1-ih,1-ih,1,kv),f1(1-ih,1-ih,1,kv),np+1,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,2)
       end if


      if (leftdedge.eq.1) then                       !add d rrl
         do k=1,n3m
            do j=1,mp
               f1(1 ,j,k,kv)= ibcx*f1(-1,j,k,kv)
     .              +ibxo*donor(c1*xe(1,j,k,kv)+c2,
     .                          c1*x(1,j,k,kv)+c2,
     .              v1(1,j,k,kv))
            end do
         end do
      end if

      if (rightdedge.eq.1) then                       !add d rrl
         do k=1,n3m
            do j=1,mp
               f1(np+1,j,k,kv)= ibcx*f1(np+3,j,k,kv)
     .  +ibxo*donor(c1*x(np,j,k,kv)+c2,
     .              c1*xe(np,j,k,kv)+c2,
     .         v1(np+1,j,k,kv))
            end do
         end do
      end if


      if (topdedge.eq.0) then                       !add d rrl
          call updated(v2(1-ih,1-ih,1,kv),v2(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih+1,3)
      else
          call updated(v2(1-ih,1-ih,1,kv),v2(1-ih,1-ih,1,kv),np,mp+1,l,
     .1-ih,np+ih,1-ih,mp+ih+1,3)
      end if

      if (botdedge.eq.0) then                       !add d rrl
         do k=1,n3m
            do j=1,mp
               do i=1,np
                  f2(i,j,k,kv)=donor(c1*x(i,j-1,k,kv)+c2,
     .                               c1*x(i,j,k,kv)+c2,
     .                 v2(i,j,k,kv))
               end do
            end do
         end do
      end if
      if (botdedge.eq.1) then                       !add d rrl
         do k=1,n3m
            do j=2,mp
               do i=1,np
                  f2(i,j,k,kv)=donor(c1*x(i,j-1,k,kv)+c2,
     .                              c1*x(i,j,k,kv)+c2,
     .                 v2(i,j,k,kv))
               end do
            end do
         end do
      end if

      if (topdedge.eq.0) then                       !add d rrl
         call updated(f2(1-ih,1-ih,1,kv),f2(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih+1,3)
      else
         call updated(f2(1-ih,1-ih,1,kv),f2(1-ih,1-ih,1,kv),np,mp+1,l,
     .1-ih,np+ih,1-ih,mp+ih+1,3)
      end if


      if (botdedge.eq.1) then                       !add d rrl
         do k=1,n3m
            do i=1,np
               f2(i,1 ,k,kv)= ibcy*f2(i,-1,k,kv)
     .              +ibyo*donor(c1*xe(i,1,k,kv)+c2,
     .                          c1*x(i,1,k,kv)+c2,
     .              v2(i,1,k,kv))

            end do
         end do
      end if
      if (topdedge.eq.1) then                       !add d rrl
         do k=1,n3m
            do i=1,np
               f2(i,mp+1,k,kv)= ibcy*f2(i,mp+3,k,kv)
     .              +ibyo*donor(c1*x(i,mp,k,kv)+c2,
     .                          c1*xe(i,mp,k,kv)+c2,
     .              v2(i,mp+1,k,kv))
            end do
         end do
      end if


      do 333 k=2,n3m
      do 333 j=1,mp
      do 333 i=1,np
  333 f3(i,j,k,kv)=donor(c1*x(i,j,k-1,kv)+c2,c1*x(i,j,k,kv)+c2,
     .v3(i,j,k,kv))

         if(ibctopbot.eq.0) then
         do j=1,mp
            do i=1,np
               f3(i,j, 1,kv)=-f3(i,j,2,kv)
               f3(i,j,n3,kv)=-f3(i,j,n3m,kv)
            enddo
         enddo
         else
         do j=1,mp
            do i=1,np
               f3(i,j, 1,kv)=0.
               f3(i,j,n3,kv)=0.
            enddo
         enddo
         endif


      if (rightdedge.eq.0) then                       !add d rrl
         call updated(f1(1-ih,1-ih,1,kv),f1(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,2)
      else
         call updated(f1(1-ih,1-ih,1,kv),f1(1-ih,1-ih,1,kv),np+1,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,2)
      endif
      if(topdedge.eq.0) then                       !add d rrl
         call updated(f2(1-ih,1-ih,1,kv),f2(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih+1,3)
      else
         call updated(f2(1-ih,1-ih,1,kv),f2(1-ih,1-ih,1,kv),np,mp+1,l,
     .1-ih,np+ih,1-ih,mp+ih+1,3)
      end if

      do k=1,n3m
      do j=1,mp
      do i=1,np
         x(i,j,k,kv)=x(i,j,k,kv)-( f1(i+1,j,k,kv)-f1(i,j,k,kv)
     .                            +f2(i,j+1,k,kv)-f2(i,j,k,kv)
     .                            +f3(i,j,k+1,kv)-f3(i,j,k,kv) )
     .  /h(i,j,k)
      end do
      end do
      end do

      enddo


      if(itr.eq.iord) go to 6
      c1=0.
      c2=1.
      do kv=1,nv
      do k=1,n3m
      do j=1,mp
      do i=1,np+1
         f1(i,j,k,kv)=v1(i,j,k,kv)
         v1(i,j,k,kv)=0.
      end do
      end do
      end do


      do k=1,n3m
      do j=1,mp+1
      do i=1,np
         f2(i,j,k,kv)=v2(i,j,k,kv)
         v2(i,j,k,kv)=0.
      end do
      end do
      end do


      do k=1,n3
      do j=1,mp
      do i=1,np
         f3(i,j,k,kv)=v3(i,j,k,kv)
         v3(i,j,k,kv)=0.
      end do
      end do
      end do


compute antidiffusive velocities in x direction
      call updated(x(1-ih,1-ih,1,kv),x(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih,1)

      if (topdedge.eq.0) then                       !add d rrl
         call updated(f2(1-ih,1-ih,1,kv),f2(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih+1,1)
      else
         call updated(f2(1-ih,1-ih,1,kv),f2(1-ih,1-ih,1,kv),np,mp+1,l,
     .1-ih,np+ih,1-ih,mp+ih+1,1)
      end if
      call updated(f3(1-ih,1-ih,1,kv),f3(1-ih,1-ih,1,kv),np,mp,l+1,
     .1-ih,np+ih,1-ih,mp+ih,1)

      illim = 1 + 1*leftdedge                       !add d rrl
      iulim = np
      jllim = 1  + 1*botdedge                       !add d rrl
      julim = mp - 1*topdedge                       !add d rrl

      do k=2,n3-2
      do j=jllim,julim
      do i=illim,iulim
      v1(i,j,k,kv)=vdyf(x(i-1,j,k,kv),x(i,j,k,kv),f1(i,j,k,kv),
     *              .5*(h(i-1,j,k)+h(i,j,k)))
     * +vcorr(f1(i,j,k,kv),
     *      f2(i-1,j,k,kv)+f2(i-1,j+1,k,kv)
     *     +f2(i,j+1,k,kv)+f2(i,j,k,kv),
     *      x(i-1,j-1,k,kv),x(i,j-1,k,kv),
     *      x(i-1,j+1,k,kv),x(i,j+1,k,kv),
     *  .5*(h(i-1,j,k)+h(i,j,k)))
     * +vcorr(f1(i,j,k,kv),
     *      f3(i-1,j,k,kv)+f3(i-1,j,k+1,kv)
     *     +f3(i,j,k+1,kv)+f3(i,j,k,kv),
     *       x(i-1,j,k-1,kv),x(i,j,k-1,kv),
     *       x(i-1,j,k+1,kv),x(i,j,k+1,kv),
     *   .5*(h(i-1,j,k)+h(i,j,k)))
      end do
      end do
      end do
      if (rightdedge.eq.0) then                       !add d rrl
         call updated(f1(1-ih,1-ih,1,kv),f1(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,1)
      else
         call updated(f1(1-ih,1-ih,1,kv),f1(1-ih,1-ih,1,kv),np+1,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,1)
      end if

      if(ibcy.eq.1) then

         if (botdedge.eq.1) then                       !add d rrl
            illim = 1  + 1*leftdedge                       !add d rrl
            iulim = np
            do k=2,n3-2
               do i=illim,iulim
                  v1(i,1,k,kv)=vdyf(x(i-1,1,k,kv),x(i,1,k,kv),
     *            f1(i,1,k,kv),
     *         .5*(h(i-1,1,k)+h(i,1,k)))
     *     +vcorr(f1(i,1,k,kv),
     *            f2(i-1,1,k,kv)+f2(i-1,2,k,kv)
     *           +f2(i,2,k,kv)+f2(i,1,k,kv),
     *             x(i-1,-1,k,kv),x(i,-1,k,kv),
     *             x(i-1,2,k,kv),x(i,2,k,kv),
     *         .5*(h(i-1,1,k)+h(i,1,k)))
     *     +vcorr(f1(i,1,k,kv),
     *            f3(i-1,1,k,kv)+f3(i-1,1,k+1,kv)
     *           +f3(i,1,k+1,kv)+f3(i,1,k,kv),
     *             x(i-1,1,k-1,kv),x(i,1,k-1,kv),
     *             x(i-1,1,k+1,kv),x(i,1,k+1,kv),
     *             .5*(h(i-1,1,k)+h(i,1,k)))
               enddo
            enddo
         end if

         if (rightdedge.eq.0) then                       !add d rrl
            call updated(v1(1-ih,1-ih,1,kv),v1(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,1)
         else
            call updated(v1(1-ih,1-ih,1,kv),v1(1-ih,1-ih,1,kv),np+1,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,1)
         end if

         if (topdedge.eq.1) then                       !add d rrl
            illim = 1 + 1*leftdedge                       !add d rrl
            iulim = np
            do k=2,n3-2
               do i=illim,iulim
                  v1(i,mp,k,kv)=v1(i,mp+1,k,kv)
               enddo
            enddo
         end if
      end if


      if(idiv.eq.1) then
         illim = 1 + 1*leftdedge                       !add d rrl
         iulim = np
         jllim = 1  + (1-ibcy)*botdedge                       !add d rrl
         julim = mp + (-1+ibcy)*topdedge                       !add d rrl
      do 511 k=2,n3-2
      do 511 j=jllim,julim
      do 511 i=illim,iulim
      v1d=-vdiv1(f1(i-1,j,k,kv),f1(i,j,k,kv),f1(i+1,j,k,kv),
     *                 .5*(h(i-1,j,k)+h(i,j,k)))
     *  -vdiv2(f1(i,j,k,kv),f2(i-1,j+1,k,kv),
     *         f2(i,j+1,k,kv),f2(i-1,j,k,kv),
     *   f2(i,j,k,kv),   .5*(h(i-1,j,k)+h(i,j,k)))
     *  -vdiv2(f1(i,j,k,kv),f3(i-1,j,k+1,kv),
     *         f3(i,j,k+1,kv),f3(i-1,j,k,kv),
     *   f3(i,j,k,kv),   .5*(h(i-1,j,k)+h(i,j,k)))
  511 v1(i,j,k,kv)=v1(i,j,k,kv)+
     *     (pp(v1d)*x(i-1,j,k,kv)-pn(v1d)*x(i,j,k,kv))
      endif

compute antidiffusive velocities in y direction
         illim = 1  + 1*leftdedge                       !add d rrl
         iulim = np - 1*rightdedge                       !add d rrl
         jllim = 1 + 1*botdedge                       !add d rrl
         julim = mp

      do k=2,n3-2
      do j=jllim,julim
      do i=illim,iulim
      v2(i,j,k,kv)=vdyf(x(i,j-1,k,kv),x(i,j,k,kv),f2(i,j,k,kv),
     *              .5*(h(i,j-1,k)+h(i,j,k)))
     * +vcorr(f2(i,j,k,kv),
     *      f1(i,j-1,k,kv)+f1(i,j,k,kv)
     *     +f1(i+1,j,k,kv)+f1(i+1,j-1,k,kv),
     *       x(i-1,j-1,k,kv),x(i-1,j,k,kv),
     *       x(i+1,j-1,k,kv),x(i+1,j,k,kv),
     *               .5*(h(i,j-1,k)+h(i,j,k)))
     * +vcorr(f2(i,j,k,kv),
     *      f3(i,j-1,k,kv)+f3(i,j,k,kv)
     *     +f3(i,j,k+1,kv)+f3(i,j-1,k+1,kv),
     *       x(i,j-1,k-1,kv),x(i,j,k-1,kv),
     *       x(i,j-1,k+1,kv),x(i,j,k+1,kv),
     *               .5*(h(i,j-1,k)+h(i,j,k)))
      end do
      end do
      end do

      if(ibcx.eq.1) then
         jllim = 1 + 1*botdedge                       !add d rrl
         julim = mp
         if (leftdedge.eq.1) then                       !add d rrl
            do k=2,n3-2
               do j=jllim,julim
                  v2(1,j,k,kv)=vdyf(x(1,j-1,k,kv),
     *                     x(1,j,k,kv),f2(1,j,k,kv),
     *                 .5*(h(1,j-1,k)+h(1,j,k)))
     *             +vcorr(f2(1,j,k,kv),
     *                    f1(1,j-1,k,kv)+f1(1,j,k,kv)
     *                   +f1(2,j,k,kv)+f1(2,j-1,k,kv),
     *                     x(-1,j-1,k,kv),x(-1,j,k,kv),
     *                     x(2,j-1,k,kv),x(2,j,k,kv),
     *             .5*(h(1,j-1,k)+h(1,j,k)))
     *             +vcorr(f2(1,j,k,kv),
     *                    f3(1,j-1,k,kv)+f3(1,j,k,kv)
     *                   +f3(1,j,k+1,kv)+f3(1,j-1,k+1,kv),
     *                     x(1,j-1,k-1,kv),x(1,j,k-1,kv),
     *                     x(1,j-1,k+1,kv),x(1,j,k+1,kv),
     *             .5*(h(1,j-1,k)+h(1,j,k)))
               end do
            end do
         end if

         if (topdedge.eq.0) then                       !add d rrl
            call updated(v2(1-ih,1-ih,1,kv),v2(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih+1,2)
         else
            call updated(v2(1-ih,1-ih,1,kv),v2(1-ih,1-ih,1,kv),np,mp+1,l,
     .1-ih,np+ih,1-ih,mp+ih+1,2)
         end if

         if (rightdedge.eq.1) then                       !add d rrl
            do k=2,n3-2
               do j=jllim,julim
                  v2(np,j,k,kv)=v2(np+1,j,k,kv)
               end do
            end do
         end if
      end if


      if(idiv.eq.1) then
         illim = 1  + (1-ibcx)*leftdedge                       !add d rrl
         iulim = np + (-1+ibcx)*rightdedge                       !add d rrl
         jllim = 1 + 1*botdedge                       !add d rrl
         julim = mp
         do k=2,n3-2
            do j=jllim,julim
               do i=illim,iulim
                  v2d=-vdiv1(f2(i,j-1,k,kv),
     *                       f2(i,j,k,kv),f2(i,j+1,k,kv),
     *             .5*(h(i,j-1,k)+h(i,j,k)))
     *                -vdiv2(f2(i,j,k,kv),
     *                       f1(i+1,j-1,k,kv),f1(i+1,j,k,kv),
     *                       f1(i,j-1,k,kv),f1(i,j,k,kv),
     *                 .5*(h(i,j-1,k)+h(i,j,k)))
     *                -vdiv2(f2(i,j,k,kv),f3(i,j-1,k+1,kv),
     *                 f3(i,j,k+1,kv),f3(i,j-1,k,kv),
     *                 f3(i,j,k,kv),    .5*(h(i,j-1,k)+h(i,j,k)))
                  v2(i,j,k,kv)=v2(i,j,k,kv)+(pp(v2d)*x(i,j-1,k,kv)-
     *                                       pn(v2d)*x(i,j,k,kv))
               end do
            end do
         end do
      endif

compute antidiffusive velocities in z direction

      illim = 1  + 1*leftdedge                       !add d rrl
      iulim = np - 1*rightdedge                       !add d rrl
      jllim = 1  + 1*botdedge                       !add d rrl
      julim = mp - 1*topdedge                       !add d rrl


      do 53 k=2,n3m
      do 53 j=jllim,julim
      do 53 i=illim,iulim

   53 v3(i,j,k,kv)=vdyf(x(i,j,k-1,kv),x(i,j,k,kv),f3(i,j,k,kv),
     *                 .5*(h(i,j,k-1)+h(i,j,k)))
     * +vcorr(f3(i,j,k,kv),
     *        f1(i,j,k-1,kv)+f1(i,j,k,kv)
     *       +f1(i+1,j,k,kv)+f1(i+1,j,k-1,kv),
     *         x(i-1,j,k-1,kv),x(i-1,j,k,kv),
     *         x(i+1,j,k-1,kv),x(i+1,j,k,kv),
     *               .5*(h(i,j,k-1)+h(i,j,k)))
     * +vcorr(f3(i,j,k,kv),
     *        f2(i,j,k-1,kv)+f2(i,j+1,k-1,kv)
     *       +f2(i,j+1,k,kv)+f2(i,j,k,kv),
     *         x(i,j-1,k-1,kv),x(i,j-1,k,kv),
     *         x(i,j+1,k-1,kv),x(i,j+1,k,kv),
     *               .5*(h(i,j,k-1)+h(i,j,k)))
c     if(mpi_rank.eq.1.and.kv.eq.4) then 
c     print*,f3(1,5,4,4),'v3'
c     print*,v3(1,5,4,4),'v3'
c     endif


      if(ibcx.eq.1) then
         jllim = 1  + 1*botdedge                       !add d rrl
         julim = mp - 1*topdedge                       !add d rrl
      if (leftdedge.eq.1) then                       !add d rrl
         do k=2,n3m
            do j=jllim,julim
               v3(1,j,k,kv)=vdyf(x(1,j,k-1,kv),
     *                           x(1,j,k,kv),f3(1,j,k,kv),
     *             .5*(h(1,j,k-1)+h(1,j,k)))
     *             +vcorr(f3(1,j,k,kv),
     *                    f1(1,j,k-1,kv)+f1(1,j,k,kv)
     *                   +f1(2,j,k,kv)+f1(2,j,k-1,kv),
     *                     x(-1,j,k-1,kv),x(-1,j,k,kv),
     *                     x(2,j,k-1,kv),x(2,j,k,kv),
     *             .5*(h(1,j,k-1)+h(1,j,k)))
     *             +vcorr(f3(1,j,k,kv),
     *                    f2(1,j,k-1,kv)+f2(1,j+1,k-1,kv)
     *                   +f2(1,j+1,k,kv)+f2(1,j,k,kv),
     *                     x(1,j-1,k-1,kv),x(1,j-1,k,kv),
     *                     x(1,j+1,k-1,kv),x(1,j+1,k,kv),
     *             .5*(h(1,j,k-1)+h(1,j,k)))
            enddo
         enddo
      endif
      call updated(v3(1-ih,1-ih,1,kv),v3(1-ih,1-ih,1,kv),np,mp,l+1,
     .1-ih,np+ih,1-ih,mp+ih,2)

      if (rightdedge.eq.1) then                       !add d rrl
         do k=2,n3m
            do j=jllim,julim
               v3(np,j,k,kv)=v3(np+1,j,k,kv)
c               ja=(mpos-1)*mp + j
c               v3(np,j,k)=vdyf(x(np+1,j,k-1),x(np+1,j,k),f3(np+1,j,k),
c     *             .5*(h(np+1,j,k-1)+h(np+1,j,k)))
c     *             +vcorr(f3(np+1,j,k),
c     *             f1(np+1,j,k-1)+f1(np+1,j,k)+f1(np+2,j,k)+
c     *             f1(np+2,j,k-1),
c     *             x(np-2,j,k-1),x(np-2,j,k),x(np+2,j,k-1),x(np+2,j,k),
c     *             .5*(h(np+1,j,k-1)+h(np+1,j,k)))
c     *             +vcorr(f3(np+1,j,k),
c     *             f2(np+1,j,k-1)+f2(np+1,j+1,k-1)+
c     *             f2(np+1,j+1,k)+f2(np+1,j,k),
c     *             x(np+1,j-1,k-1),x(np+1,j-1,k),x(np+1,j+1,k-1),
c     *             x(np+1,j+1,k),
c     *             .5*(h(np+1,j,k-1)+h(np+1,j,k)))
            enddo
         enddo
      endif
      end if

      if(ibcy.eq.1) then
         illim = 1  + 1*leftdedge                       !add d rrl
         iulim = np - 1*rightdedge                       !add d rrl
         if (botdedge.eq.1) then                       !add d rrl
            do k=2,n3m
               do i=illim,iulim

                  v3(i,1,k,kv)=vdyf(x(i,1,k-1,kv),
     *             x(i,1,k,kv),f3(i,1,k,kv),
     *             .5*(h(i,1,k-1)+h(i,1,k)))
     *             +vcorr(f3(i,1,k,kv),
     *                    f1(i,1,k-1,kv)+f1(i,1,k,kv)
     *                   +f1(i+1,1,k,kv)+f1(i+1,1,k-1,kv),
     *                     x(i-1,1,k-1,kv),x(i-1,1,k,kv),
     *                     x(i+1,1,k-1,kv),x(i+1,1,k,kv),
     *             .5*(h(i,1,k-1)+h(i,1,k)))
     *             +vcorr(f3(i,1,k,kv),
     *                    f2(i,1,k-1,kv)+f2(i,2,k-1,kv)
     *                   +f2(i,2,k,kv)+f2(i,1,k,kv),
     *                     x(i,-1,k-1,kv),x(i,-1,k,kv),
     *                     x(i,2,k-1,kv),x(i,2,k,kv),
     *             .5*(h(i,1,k-1)+h(i,1,k)))
               enddo
            enddo
         end if

      call updated(v3(1-ih,1-ih,1,kv),v3(1-ih,1-ih,1,kv),np,mp,l+1,
     .1-ih,np+ih,1-ih,mp+ih,3)

         if (topdedge.eq.1) then                       !add d rrl
            do k=2,n3m
               do i=illim,iulim
                  v3(i,mp,k,kv)=v3(i,mp+1,k,kv)

c                  v3(i,mp,k)=vdyf(x(i,mp+1,k-1),x(i,mp+1,k),
c     *             f3(i,mp+1,k),
c     *             .5*(h(i,mp+1,k-1)+h(i,mp+1,k)))
c     *             +vcorr(f3(i,mp+1,k),
c     *             f1(i,mp+1,k-1)+f1(i,mp+1,k)+f1(i+1,mp+1,k)+
c     *             f1(i+1,mp+1,k-1),
c     *             x(i-1,mp+1,k-1),x(i-1,mp+1,k),x(i+1,mp+1,k-1),
c     *             x(i+1,mp+1,k),
c     *             .5*(h(i,mp+1,k-1)+h(i,mp+1,k)))
c     *             +vcorr(f3(i,mp+1,k),
c     *             f2(i,mp+2,k-1)+f2(i,mp+3,k-1)+f2(i,mp+3,k)+
c     *             f2(i,mp+2,k),
c     *             x(i,mp-1,k-1),x(i,mp-1,k),x(i,mp+2,k-1),x(i,mp+2,k),
c     *             .5*(h(i,mp+1,k-1)+h(i,mp+1,k)))
               enddo
            enddo
         end if



      if(ibcx.eq.1) then
         if (leftdedge.eq.1 .and. botdedge.eq.1) then                       !add d rrl
            do k=2,n3m
               v3(1,1,k,kv)=vdyf(x(1,1,k-1,kv),
     *          x(1,1,k,kv),f3(1,1,k,kv),
     *              .5*(h(1,1,k-1)+h(1,1,k)))
     *              +vcorr(f3(1,1,k,kv),
     *                     f1(1,1,k-1,kv)+f1(1,1,k,kv)
     *                    +f1(2,1,k,kv)+f1(2,1,k-1,kv),
     *                      x(-1,1,k-1,kv),x(-1,1,k,kv),
     *                      x(2,1,k-1,kv),x(2,1,k,kv),
     *              .5*(h(1,1,k-1)+h(1,1,k)))
     *              +vcorr(f3(1,1,k,kv),
     *                     f2(1,1,k-1,kv)+f2(1,2,k-1,kv)
     *                    +f2(1,2,k,kv)+f2(1,1,k,kv),
     *                      x(1,-1,k-1,kv),x(1,-1,k,kv),
     *                      x(1,2,k-1,kv),x(1,2,k,kv),
     *              .5*(h(1,1,k-1)+h(1,1,k)))
            end do
         end if
         call updated(v3(1-ih,1-ih,1,kv),v3(1-ih,1-ih,1,kv),np,mp,l+1,
     .1-ih,np+ih,1-ih,mp+ih,1)
         if (rightdedge.eq.1 .and. topdedge.eq.1) then                       !add d rrl
            do k=2,n3m
               v3(np,mp,k,kv)=v3(np+1,mp+1,k,kv)
            enddo
         endif
         if (rightdedge.eq.1 .and. botdedge.eq.1) then                       !add d rrl
            do k=2,n3m
               v3(np,1,k,kv)=v3(np+1,1,k,kv)
            enddo
         endif
         if (leftdedge.eq.1 .and. topdedge.eq.1) then                       !add d rrl
            do k=2,n3m
               v3(1,mp,k,kv)=v3(1,mp+1,k,kv)
            enddo
         endif
      end if
      end if


      if(idiv.eq.1) then
         illim = 1  + (1-ibcx)*leftdedge                       !add d rrl
         iulim = np + (-1+ibcx)*rightdedge                       !add d rrl
       jllim = 1  + (1-ibcy)*botdedge                       !add d rrl
       julim = mp + (-1+ibcy)*topdedge                       !add d rrl

      do 531 k=2,n3m
      do 531 j=jllim,julim
      do 531 i=illim,iulim
         ia=(npos-1)*np + i
         ja=(mpos-1)*mp + j



         v2d=-vdiv1(f3(i,j,k-1,kv),
     *              f3(i,j,k,kv),f3(i,j,k+1,kv),
     *        .5*(h(i,j,k-1)+h(i,j,k)))
     *        -vdiv2(f3(i,j,k,kv),
     *               f1(i+1,j,k-1,kv),f1(i+1,j,k,kv),
     *               f1(i,j,k-1,kv),
     *               f1(i,j,k,kv),  .5*(h(i,j,k-1)+h(i,j,k)))
     *        -vdiv2(f3(i,j,k,kv),
     *               f2(i,j+1,k-1,kv),f2(i,j+1,k,kv),
     *        f2(i,j,k-1,kv),
     *     f2(i,j,k,kv),    .5*(h(i,j,k-1)+h(i,j,k)))
  531 v3(i,j,k,kv)=v3(i,j,k,kv)
     *    +(pp(v2d)*x(i,j,k-1,kv)-pn(v2d)*x(i,j,k,kv))

      endif

      if (rightdedge.eq.0) then                       !add d rrl
         call updated(v1(1-ih,1-ih,1,kv),v1(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,1)
      else
         call updated(v1(1-ih,1-ih,1,kv),v1(1-ih,1-ih,1,kv),np+1,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,1)
      end if

       if(ibcx.eq.1) then
          if (leftdedge.eq.1) then                       !add d rrl
             do k=1,n3m
                do j=1,mp
                   v1(1,j,k,kv)=v1(-1,j,k,kv)
                end do
             end do
          end if
          if (rightdedge.eq.1) then                       !add d rrl
             do k=1,n3m
                do j=1,mp
                   v1(np+1,j,k,kv)=v1(np+3,j,k,kv)
                enddo
             enddo
          end if
       end if

      if (topdedge.eq.0) then                       !add d rrl
         call updated(v2(1-ih,1-ih,1,kv),v2(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih+1,3)
      else
         call updated(v2(1-ih,1-ih,1,kv),v2(1-ih,1-ih,1,kv),np,mp+1,l,
     .1-ih,np+ih,1-ih,mp+ih+1,3)
      end if

       if(ibcy.eq.1) then
          if (botdedge.eq.1) then                       !add d rrl
             do k=1,n3m
                do i=1,np
                   v2(i,1,k,kv)=v2(i,-1,k,kv)
                end do
             end do
          end if
          if (topdedge.eq.1) then                       !add d rrl
             do k=1,n3m
                do i=1,np
                   v2(i,mp+1,k,kv)=v2(i,mp+3,k,kv)
                enddo
             enddo
          end if
       end if
       enddo

c     if(mpi_rank.eq.1) print*,v3(1,5,4,4),'v3',f3(1,5,4,4)
                  if(nonosold.eq.1) then
c                 non-osscilatory option
      do kv=1,nv
      do 401 k=1,n3m
      km=max0(k-1,1  )
      kp=min0(k+1,n3m)
      do 401 j=1,mp
         if (botdedge.eq.1 .and. j.eq.1) then                       !add d rrl
            jm = ibcy*(-1) + ibyo*1
         else
            jm = j - 1
         end if
         if (topdedge.eq.1 .and. j.eq.mp) then                       !add d rrl
            jp = ibcy*(mp+2) + ibyo*mp
         else
            jp = j + 1
         end if
      do 401 i=1,np
         if (leftdedge.eq.1 .and. i.eq.1) then                       !add d rrl
            im = ibcx*(-1) + ibxo*1
         else
            im = i - 1
         end if
         if (rightdedge.eq.1 .and. i.eq.np) then                       !add d rrl
            ip = ibcx*(np+2) + ibxo*np
         else
            ip = i + 1
         end if
         mxo(i,j,k,kv)=amax1(x(im,j,k,kv),x(i,j,k,kv),
     .                       x(ip,j,k,kv),mxo(i,j,k,kv),
     .                       x(i,jm,k,kv),x(i,jp,k,kv),
     .                       x(i,j,kp,kv),x(i,j,km,kv))
 401     mno(i,j,k,kv)=amin1(x(im,j,k,kv),x(i,j,k,kv),
     .                       x(ip,j,k,kv),mno(i,j,k,kv),
     .                       x(i,jm,k,kv),x(i,jp,k,kv),
     .                       x(i,j,kp,kv),x(i,j,km,kv))


      illim = 1
      iulim = np + 1*rightdedge                       !add d rrl
      jllim = 1
      julim = mp


      do 402 k=1,n3m
      do 402 j=jllim,julim
      do 402 i=illim,iulim
  402 f1(i,j,k,kv)=donor(c2,c2,v1(i,j,k,kv))


      illim = 1
      iulim = np
      jllim = 1
      julim = mp + 1*topdedge                       !add d rrl


      do 403 k=1,n3m
      do 403 j=jllim,julim
      do 403 i=illim,iulim
  403 f2(i,j,k,kv)=donor(c2,c2,v2(i,j,k,kv))



      do 4033 k=1,n3
      do 4033 j=1,mp
      do 4033 i=1,np
 4033 f3(i,j,k,kv)=donor(c2,c2,v3(i,j,k,kv))

      if (rightdedge.eq.0) then                       !add d rrl
         call updated(f1(1-ih,1-ih,1,kv),f1(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,2)
      else
         call updated(f1(1-ih,1-ih,1,kv),f1(1-ih,1-ih,1,kv),np+1,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,2)
      end if
      if (topdedge.eq.0) then                       !add d rrl
         call updated(f2(1-ih,1-ih,1,kv),f2(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih+1,3)
      else
         call updated(f2(1-ih,1-ih,1,kv),f2(1-ih,1-ih,1,kv),np,mp+1,l,
     .1-ih,np+ih,1-ih,mp+ih+1,3)
      end if


      do 444 k=1,n3m
      do 444 j=1,mp
      do 444 i=1,np
      cp(i,j,k,kv)=(mxo(i,j,k,kv)-x(i,j,k,kv))*h(i,j,k)/
     1( pn(f1(i+1,j,k,kv))+pp(f1(i,j,k,kv))
     2 +pn(f2(i,j+1,k,kv))+pp(f2(i,j,k,kv))
     3 +pn(f3(i,j,k+1,kv))+pp(f3(i,j,k,kv))+ep)

      cn(i,j,k,kv)=(x(i,j,k,kv)-mno(i,j,k,kv))*h(i,j,k)/
     1( pp(f1(i+1,j,k,kv))+pn(f1(i,j,k,kv))
     2 +pp(f2(i,j+1,k,kv))+pn(f2(i,j,k,kv))
     3 +pp(f3(i,j,k+1,kv))+pn(f3(i,j,k,kv))+ep)
 444  continue


      call updated(cp(1-ih,1-ih,1,kv),cp(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih,0)
      call updated(cn(1-ih,1-ih,1,kv),cn(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih,0)


      illim = 1 + 1*leftdedge                       !add d rrl
      do k=1,n3m
        do j=1,mp
          do i=illim,np
            v1(i,j,k,kv)= pp(v1(i,j,k,kv))
     *                *amin1(1.,cp(i,j,k,kv),cn(i-1,j,k,kv))
     *                 -pn(v1(i,j,k,kv))
     *                *amin1(1.,cp(i-1,j,k,kv),cn(i,j,k,kv))
          end do
        end do
      end do
      if (ibcx.eq.1) then
         if (rightdedge.eq.0) then                       !add d rrl
            call updated(v1(1-ih,1-ih,1,kv),v1(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,1)
         else
            call updated(v1(1-ih,1-ih,1,kv),v1(1-ih,1-ih,1,kv),np+1,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,1)
         end if
         do k=1,n3m
            do j=1,mp
               if (leftdedge.eq.1) then                       !add d rrl
                  v1(1 ,j,k,kv)=v1(-1,j,k,kv)
               end if
               if (rightdedge.eq.1) then                       !add d rrl
                  v1(np+1,j,k,kv)=v1(np+3  ,j,k,kv)
               end if
            end do
         end do
      end if

      jllim = 1 + 1*botdedge                       !add d rrl
      do k=1,n3m
        do j=jllim,mp
          do i=1,np
            v2(i,j,k,kv)= pp(v2(i,j,k,kv))
     .                *amin1(1.,cp(i,j,k,kv),cn(i,j-1,k,kv))
     .                -pn(v2(i,j,k,kv))
     .                *amin1(1.,cp(i,j-1,k,kv),cn(i,j,k,kv))
          end do
        end do
      end do

      if (ibcy.eq.1) then
         if (topdedge.eq.0) then                       !add d rrl
            call updated(v2(1-ih,1-ih,1,kv),v2(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih+1,1)
         else
            call updated(v2(1-ih,1-ih,1,kv),v2(1-ih,1-ih,1,kv),np,mp+1,l,
     .1-ih,np+ih,1-ih,mp+ih+1,1)
         end if

        do k=1,n3m
          do i=1,np
             if (botdedge.eq.1) then                       !add d rrl
                v2(i, 1,k,kv)=v2(i,-1,k,kv)
             end if
             if (topdedge.eq.1) then                       !add d rrl
                v2(i,mp+1,k,kv)=v2(i,mp+3,k,kv)
             end if
          end do
        end do
      end if

      do k=2,n3m
        do j=1,mp
          do i=1,np
            v3(i,j,k,kv)= pp(v3(i,j,k,kv))
     .             *amin1(1.,cp(i,j,k,kv),cn(i,j,k-1,kv))
     *                -pn(v3(i,j,k,kv))
     .             *amin1(1.,cp(i,j,k-1,kv),cn(i,j,k,kv))
          end do
        end do
      end do

      enddo
                  endif

      if(nonos.eq.1) then
      do 1000 itrfct=1,nfct

      if(itrfct.eq.1) then
      do kv=1,nv

      illim = 1
      iulim = np + 1*rightdedge                       !add d rrl
      jllim = 1
      julim = mp
      do 502 k=1,l
      do 502 j=jllim,julim
      do 502 i=illim,iulim
  502 f1(i,j,k,kv)=donor(c2,c2,v1(i,j,k,kv))

      illim = 1
      iulim = np
      jllim = 1
      julim = mp + 1*topdedge                       !add d rrl
      do 503 k=1,n3m
      do 503 j=jllim,julim
      do 503 i=illim,iulim
  503 f2(i,j,k,kv)=donor(c2,c2,v2(i,j,k,kv))

      do 5033 k=1,n3
      do 5033 j=1,mp
      do 5033 i=1,np
 5033 f3(i,j,k,kv)=donor(c2,c2,v3(i,j,k,kv))
      if (rightdedge.eq.0) then                       !add d rrl
         call updated(f1(1-ih,1-ih,1,kv),f1(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,2)
      else
         call updated(f1(1-ih,1-ih,1,kv),f1(1-ih,1-ih,1,kv),np+1,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,2)
      end if
      if (topdedge.eq.0) then                       !add d rrl
         call updated(f2(1-ih,1-ih,1,kv),f2(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih+1,3)
      else
         call updated(f2(1-ih,1-ih,1,kv),f2(1-ih,1-ih,1,kv),np,mp+1,l,
     .1-ih,np+ih,1-ih,mp+ih+1,3)
      end if
      enddo
      endif
correction coefficients for variables nv-1
      do kv=1,nv-1
      do k=1,l
      do j=1,mp
      do i=1,np
      if(abs(mx(i,j,k,kv)).lt.ep) mx(i,j,k,kv)=0.
      if(abs(mn(i,j,k,kv)).lt.ep) mn(i,j,k,kv)=0.
      if(k.eq.1) f3(i,j,k,kv)=-f3(i,j,k+1,kv)
      if(k.eq.1) f3(i,j,k,nv)=-f3(i,j,k+1,nv)
      if(k.eq.l) f3(i,j,k+1,kv)=-f3(i,j,k,kv)
      if(k.eq.l) f3(i,j,k+1,nv)=-f3(i,j,k,nv)
      rhoin=  (pn(f1(i+1,j,k,nv))+pp(f1(i,j,k,nv))
     .        +pn(f2(i,j+1,k,nv))+pp(f2(i,j,k,nv))
     .        +pn(f3(i,j,k+1,nv))+pp(f3(i,j,k,nv)))
      rhoout=-(pp(f1(i+1,j,k,nv))+pn(f1(i,j,k,nv))
     .        +pp(f2(i,j+1,k,nv))+pn(f2(i,j,k,nv))
     .        +pp(f3(i,j,k+1,nv))+pn(f3(i,j,k,nv)))
      ain=    (pn(f1(i+1,j,k,kv))+pp(f1(i,j,k,kv))
     .        +pn(f2(i,j+1,k,kv))+pp(f2(i,j,k,kv))
     .        +pn(f3(i,j,k+1,kv))+pp(f3(i,j,k,kv)))
      aout=  -(pp(f1(i+1,j,k,kv))+pn(f1(i,j,k,kv))
     .        +pp(f2(i,j+1,k,kv))+pn(f2(i,j,k,kv))
     .        +pp(f3(i,j,k+1,kv))+pn(f3(i,j,k,kv)))
          cp(i,j,k,kv)=
     . pp(mx(i,j,k,kv)*x(i,j,k,nv)-x(i,j,k,kv))*h(i,j,k)
     .       /(ain-pp(mx(i,j,k,kv))*rhoout
     .            +pn(mx(i,j,k,kv))*rhoin+ep)
         cn(i,j,k,kv)=
     . pp(x(i,j,k,kv)-mn(i,j,k,kv)*x(i,j,k,nv))*h(i,j,k)
     .     /(-aout+pp(mn(i,j,k,kv))*rhoin
     .            -pn(mn(i,j,k,kv))*rhoout+ep)
      enddo
      enddo
      enddo
      call updated(cp(1-ih,1-ih,1,kv),cp(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih,0)
      call updated(cn(1-ih,1-ih,1,kv),cn(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih,0)
      enddo



correction coeffiecients for rho
      do kv=1,nv-1
      do k=1,l
      do j=1,mp
      do i=1,np
c  JLW - Fix from Jon Reisner 9/2005
c     a1p=cp(i,j,k,kv)+amax1(0.,sign(1., mx(i,j,k,kv)))
c     a2p=cn(i,j,k,kv)+amax1(0.,sign(1.,-mn(i,j,k,kv)))
c     a1n=cn(i,j,k,kv)+amax1(0.,sign(1., mn(i,j,k,kv)))
c     a2n=cp(i,j,k,kv)+amax1(0.,sign(1.,-mx(i,j,k,kv)))
      rmxuse=-1.*mx(i,j,k,kv)
      if(mx(i,j,k,kv).eq.0.) rmxuse=0.
      rmnuse=-1.*mn(i,j,k,kv)
      if(mn(i,j,k,kv).eq.0.) rmnuse=0.
      a1p=cp(i,j,k,kv)+amax1(0.,sign(1., mx(i,j,k,kv)))
      a2p=cn(i,j,k,kv)+amax1(0.,sign(1.,rmnuse))
      a1n=cn(i,j,k,kv)+amax1(0.,sign(1., mn(i,j,k,kv)))
      a2n=cp(i,j,k,kv)+amax1(0.,sign(1.,rmxuse))
      tmpp=amin1(a1p,a2p)
      tmpn=amin1(a1n,a2n)
      if(kv.eq.1) cp(i,j,k,nv)=tmpp
      if(kv.eq.1) cn(i,j,k,nv)=tmpn
      cp(i,j,k,nv)=amin1(tmpp,cp(i,j,k,nv))
      cn(i,j,k,nv)=amin1(tmpn,cn(i,j,k,nv))
      enddo
      enddo
      enddo
      enddo
      call updated(cp(1-ih,1-ih,1,nv),cp(1-ih,1-ih,1,nv),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih,0)
      call updated(cn(1-ih,1-ih,1,nv),cn(1-ih,1-ih,1,nv),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih,0)

c     if(mpi_rank.eq.1) then
c     print*,v1(1,5,3,4),'v1'
c     print*,v1(2,5,3,4),'v1'
c     print*,v2(1,5,3,4),'v2'
c     print*,v2(1,6,3,4),'v2'
c     print*,v3(1,5,3,4),'v3'
c     print*,v3(1,5,4,4),'v3'
c     print*, pp(v1(1,5,3,4))
c    *        *amin1(1.,cp(1,5,3,4),cn(0,5,3,4))
c    *                 -pn(v1(1,5,3,4))
c    *        *amin1(1.,cp(0,5,3,4),cn(1,5,3,4))
c     endif

      do kv=1,nv
      do k=1,l
        do j=1,mp
          do i=1+leftdedge,np                        !add d rrl
          v1(i,j,k,kv)= pp(v1(i,j,k,kv))
     *        *amin1(1.,cp(i,j,k,kv),cn(i-1,j,k,kv))
     *                 -pn(v1(i,j,k,kv))
     *        *amin1(1.,cp(i-1,j,k,kv),cn(i,j,k,kv))
          end do
        end do
      end do
      if (ibcx.eq.1) then
         if (rightdedge.eq.0) then                       !add d rrl
            call updated(v1(1-ih,1-ih,1,kv),v1(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,1)
         else
            call updated(v1(1-ih,1-ih,1,kv),v1(1-ih,1-ih,1,kv),np+1,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,1)
         end if
         do k=1,n3m
            do j=1,mp
               if (leftdedge.eq.1) then                       !add d rrl
                  v1(1 ,j,k,kv)=v1(-1,j,k,kv)
               end if                       !add d rrl
               if (rightdedge.eq.1) then                       !add d rrl
                  v1(np+1,j,k,kv)=v1(np+3  ,j,k,kv)
               end if
            end do
         end do
      end if

      do k=1,l
        do j=1+botdedge,mp                        !add d rrl
          do i=1,np
            v2(i,j,k,kv)= pp(v2(i,j,k,kv))
     *             *amin1(1.,cp(i,j,k,kv),cn(i,j-1,k,kv))
     *                   -pn(v2(i,j,k,kv))
     *             *amin1(1.,cp(i,j-1,k,kv),cn(i,j,k,kv))
          end do
        end do
      end do

      if (ibcy.eq.1) then
         if (topdedge.eq.0) then                       !add d rrl
            call updated(v2(1-ih,1-ih,1,kv),v2(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih+1,1)
         else
            call updated(v2(1-ih,1-ih,1,kv),v2(1-ih,1-ih,1,kv),np,mp+1,l,
     .1-ih,np+ih,1-ih,mp+ih+1,1)
         end if

        do k=1,n3m
          do i=1,np
             if (botdedge.eq.1) then                       !add d rrl
                v2(i, 1,k,kv)=v2(i,-1,k,kv)
             end if
             if (topdedge.eq.1) then                       !add d rrl
                v2(i,mp+1,k,kv)=v2(i,mp+3,k,kv)
             end if
          end do
        end do
      end if

      do k=2,l
        do j=1,mp
          do i=1,np
            v3(i,j,k,kv)= pp(v3(i,j,k,kv))
     *             *amin1(1.,cp(i,j,k,kv),cn(i,j,k-1,kv))
     *                   -pn(v3(i,j,k,kv))
     *             *amin1(1.,cp(i,j,k-1,kv),cn(i,j,k,kv))
          end do
        end do
      end do
      enddo
c     if(mpi_rank.eq.1) then
c     print*,v1(1,5,3,4),'v1'
c     print*,v1(2,5,3,4),'v1'
c     print*,cn(1,5,3,4),'cn'
c     print*,cn(0,5,3,4),'cn'
c     print*,cp(1,5,3,4),'cp'
c     print*,cp(0,5,3,4),'cp'
c     print*,v2(1,5,3,4),'v2'
c     print*,v2(1,6,3,4),'v2'
c     print*,v3(1,5,3,4),'v3'
c     print*,v3(1,5,4,4),'v3'
c     endif

      if(itrfct.lt.nfct) then
      do kv=1,nv
      do 602 k=1,l
      do 602 j=1,mp
      do 602 i=1,np+1
      tmp=f1(i,j,k,kv)
      f1o(i,j,k,kv)=donor(c2,c2,v1(i,j,k,kv))
  602 f1(i,j,k,kv)=tmp-f1o(i,j,k,kv)
      do 603 k=1,l
      do 603 j=1,mp+1
      do 603 i=1,np
      tmp=f2(i,j,k,kv)
      f2o(i,j,k,kv)=donor(c2,c2,v2(i,j,k,kv))
  603 f2(i,j,k,kv)=tmp-f2o(i,j,k,kv)
      do 6033 k=1,n3
      do 6033 j=1,mp
      do 6033 i=1,np
      tmp=f3(i,j,k,kv)
      f3o(i,j,k,kv)=donor(c2,c2,v3(i,j,k,kv))
 6033 f3(i,j,k,kv)=tmp-f3o(i,j,k,kv)

      if(rightdedge.eq.0) then                       !add d rrl
       call updated(f1o(1-ih,1-ih,1,kv),f1o(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,2)
      else
       call updated(f1o(1-ih,1-ih,1,kv),f1o(1-ih,1-ih,1,kv),np+1,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,2)
      endif
      if(botdedge.eq.0) then                       !add d rrl
       call updated(f2o(1-ih,1-ih,1,kv),f2o(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih+1,3)
      else
       call updated(f2o(1-ih,1-ih,1,kv),f2o(1-ih,1-ih,1,kv),np,mp+1,l,
     .1-ih,np+ih,1-ih,mp+ih+1,3)
      endif

      do k=1,l
      do j=1,mp
      do i=1,np
      x(i,j,k,kv)=x(i,j,k,kv)-( f1o(i+1,j,k,kv)-f1o(i,j,k,kv)
     .                         +f2o(i,j+1,k,kv)-f2o(i,j,k,kv)
     .                         +f3o(i,j,k+1,kv)-f3o(i,j,k,kv) )
     .                          /h(i,j,k)
      enddo
      enddo
      enddo
      enddo
      endif
1000  continue
      endif

   30                      continue
    6 continue

      do kv=1,nv
      call updated(x(1-ih,1-ih,1,kv),x(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih,1)
      enddo

      deallocate (v1)
      deallocate (v2)
      deallocate (v3)
      deallocate (f1)
      deallocate (f2)
      deallocate (f3)
      deallocate (f1o)
      deallocate (f2o)
      deallocate (f3o)
      deallocate (cp)
      deallocate (cn)
      deallocate (mxo)
      deallocate (mno)
      deallocate (mx)
      deallocate (mn)
      deallocate (a)
2500  continue

      return
      end subroutine mpdatanew3d
!*****************************************************************************************!
      subroutine velprdf(u,v,o,ua,va,oa,gcx,gcy,gcz,il,iu,jl,ju,lls,kkl,kkh)
      use metryic
      use gridsetup
      use msga
      use xve   !this use staement allows access to the environmental array
      Implicit None

      !JAS 3/6/06 added explicit declarations to comply with implicit none
      integer,intent(in) :: il,iu,jl,ju,lls,kkl,kkh
      real :: gcx,gcy,gcz
      integer :: i,j,k

      real u(il:iu,jl:ju,lls,kkl:kkh),
     .     v(il:iu,jl:ju,lls,kkl:kkh),
     .     o(il:iu,jl:ju,lls,kkl:kkh)
      real ua(il:iu+1,jl:ju,lls),
     .     va(il:iu,jl:ju+1,lls),
     .     oa(il:iu,jl:ju,lls+1)

      do k=1,l
       do j=1,mp
        do i=1,np
         u(i,j,k,2)=0.5*(3.*u(i,j,k,1)-u(i,j,k,0))
         v(i,j,k,2)=0.5*(3.*v(i,j,k,1)-v(i,j,k,0))
         o(i,j,k,2)=0.5*(3.*o(i,j,k,1)-o(i,j,k,0))
        enddo
       enddo
      enddo

      do k=1,l
       do j=1,mp
        do i=1,np
         u(i,j,k,2)=u(i,j,k,2)*gcx/gi(i,j,k)
         v(i,j,k,2)=v(i,j,k,2)*gcy/gi(i,j,k)
         o(i,j,k,2)=o(i,j,k,2)*gcz/gi(i,j,k)
        enddo
       enddo
      enddo

      !compute covariant cell face velocites in z direction
      do k=2,l
       do j=1,mp
        do i=1,np
         oa(i,j,k)=0.5*(o(i,j,k,2)+o(i,j,k-1,2))
        enddo
       enddo
      enddo
      if(ibctopbot.eq.0) then
       do j=1,mp
        do i=1,np
         !oa(i,j,1)=xe(i,j,1,3)/xe(i,mp,1,5)*gcz/gi(i,j,k)
         !oa(i,j,l+1)=xe(i,j,l,3)/xe(i,mp,l,5)*gcz/gi(i,j,k)
         oa(i,j,1)=-oa(i,j,2)
         oa(i,j,l+1)=-oa(i,j,l)
        enddo
       enddo
      else
       do j=1,mp
        do i=1,np
         oa(i,j,1)=0.0
         oa(i,j,l+1)=0.0
        enddo
       enddo
      endif !if(ibctopbot.eq.0)

      call updated(u(1-ih,1-ih,1,2),u(1-ih,1-ih,1,2),np,mp,l,1-ih,np+ih,1-ih,mp+ih,2)
      call updated(v(1-ih,1-ih,1,2),v(1-ih,1-ih,1,2),np,mp,l,1-ih,np+ih,1-ih,mp+ih,0)

      !compute covariant cell face velocites in y direction
      if(j3.ne.0)then
       do k=1,l
        do j=1,mp+1
         do i=1,np
          va(i,j,k)=0.5*(v(i,j,k,2)+v(i,j-1,k,2))
         enddo
        enddo
       enddo
       if(botdedge.eq.1)then
        do k=1,l
         do i=1,np
           ! here we assume 0 gradient between xe
           va(i,1,k)=xe(i,1,k,2)/xe(i,1,k,nv)*gcy/gi(i,j,k)
!          va(i,1,k)=0.5*gcy*(
!     +          (xe(i,1,k,2)/xe(i,1,k,nv))/gi(i,1,k)
!     +          +(xe(i,0,k,2)/xe(i,0,k,nv))/gi(i,0,k))
          !va(i,1,k)=va(i,2,k)
         enddo
        enddo
       endif
       if(topdedge.eq.1)then
        do k=1,l
         do i=1,np
           ! here we assume 0 gradient between xe
          va(i,mp+1,k)=xe(i,mp,k,2)/xe(i,mp,k,nv)*gcy/gi(i,j,k)
!          va(i,mp+1,k)=0.5*gcy*(
!     +          (xe(i,mp+1,k,2)/xe(i,mp+1,k,nv))/gi(i,mp+1,k)
!     +          +(xe(i,mp,k,2)/xe(i,mp,k,nv))/gi(i,mp,k))
          !va(i,mp+1,k)=va(i,mp,k)
         enddo
        enddo
       endif
      else
       va=0.0
      endif

      !compute covariant cell face velocites in y direction
      do k=1,l
       do j=1,mp
        do i=1,np+1
         ua(i,j,k)=0.5*(u(i,j,k,2)+u(i-1,j,k,2))
        enddo
       enddo
      enddo
      if(leftdedge.eq.1) then
       do k=1,l
        do j=1,mp
          ua(1,j,k)=xe(1,j,k,1)/xe(1,j,k,nv)*gcx/gi(i,j,k)
            ! here we assume 0 gradient between xe
!         ua(1,j,k)=0.5*gcx*(
!     +          (xe(1,j,k,1)/xe(1,j,k,nv))/gi(1,j,k)
!     +          +(xe(0,j,k,1)/xe(0,j,k,nv))/gi(0,j,k))
         !ua(1,j,k)=ua(2,j,k)
        enddo
       enddo
      endif
      if(rightdedge.eq.1) then
       do k=1,l
        do j=1,mp
           ua(np+1,j,k)=xe(np,j,k,1)/xe(np,j,k,nv)*gcx/gi(i,j,k)
           ! here we assume 0 gradient between xe
!         ua(np+1,j,k)=0.5*gcx*(
!     +            (xe(np+1,j,k,1)/xe(np+1,j,k,nv))/gi(np+1,j,k)
!     +            +(xe(np,j,k,1)/xe(np,j,k,nv))/gi(np,j,k))
         !ua(np+1,j,k)=ua(np,j,k)
        enddo
       enddo
      endif

       return
       end subroutine velprdf
!*****************************************************************************************!
      subroutine velprdmoa(u,v,o,fu,fv,fo,ua,va,oa,gcx,gcy,gcz,il,iu,jl,ju,lls,kkl,kkh)
      use gridsetup
      use metryic
      use msga

      Implicit None

      !JAS 3/6/06 added explicit declarations to comply with implicit none
      integer,intent(in) :: il,iu,jl,ju,lls,kkl,kkh
      real :: gcx,gcy,gcz
      integer :: i,j,k,lc
      real :: cp,cn,c1p,c1n,cnp,cnn,cmp,cmn

      real u(il:iu, jl:ju,lls,kkl:kkh),
     .     v(il:iu, jl:ju,lls,kkl:kkh),
     .     o(il:iu, jl:ju,lls,kkl:kkh)
      real fu(il:iu, jl:ju,lls),
     .     fv(il:iu, jl:ju,lls),
     .     fo(il:iu, jl:ju,lls)
      real ua(il+1:iu, jl:ju,lls),
     .     va(il:iu, jl+1:ju,lls),
     .     oa(il:iu, jl:ju,lls+1)
      real,allocatable:: ux(:, :,:),
     .                   vx(:, :,:),
     .                   ox(:, :,:)
      real,allocatable:: uy(:,:,:),
     .                   vy(:,:,:),
     .                   oy(:, :,:)
      real,allocatable:: uz(:, :,:),
     .                   vz(:, :,:),
     .                   oz(:, :,:)
      allocate (ux(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (vx(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (ox(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (uy(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (vy(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (oy(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (uz(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (vz(1-ih:np+ih, 1-ih:mp+ih,l))
      allocate (oz(1-ih:np+ih, 1-ih:mp+ih,l))

      do k=1,l
      do j=1,mp
      do i=1,np
      u(i,j,k,1)=u(i,j,k,0)+0.5*fu(i,j,k)
      v(i,j,k,1)=v(i,j,k,0)+0.5*fv(i,j,k)
      o(i,j,k,1)=o(i,j,k,0)+0.5*fo(i,j,k)
      enddo
      enddo
      enddo

      do j=1,mp
      do i=1,np
      o(i,j,1,0)=0.
      o(i,j,l,0)=0.
      o(i,j,1,1)=0.
      o(i,j,l,1)=0.
      enddo
      enddo

      lc=1
      call updated(u(1-ih,1-ih,1,lc),u(1-ih,1-ih,1,lc),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih,0)
      call updated(v(1-ih,1-ih,1,lc),v(1-ih,1-ih,1,lc),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih,0)
      call updated(o(1-ih,1-ih,1,lc),o(1-ih,1-ih,1,lc),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih,0)

      do 1 i=1+leftdedge,np-rightdedge                         !added d rrl
      do 1 j=1,mp
      do 1 k=1,l
      cp=gcx*amax1(0., u(i,j,k,0))
      cn=gcx*amin1(0., u(i,j,k,0))
      ux(i,j,k)= cp*(u(i,j,k,lc)-u(i-1,j,k,lc))
     .           +cn*(u(i+1,j,k,lc)-u(i,j,k,lc))
      vx(i,j,k)= cp*(v(i,j,k,lc)-v(i-1,j,k,lc))
     .           +cn*(v(i+1,j,k,lc)-v(i,j,k,lc))
      ox(i,j,k)= cp*(o(i,j,k,lc)-o(i-1,j,k,lc))
     .           +cn*(o(i+1,j,k,lc)-o(i,j,k,lc))
    1 continue

      if(leftdedge.eq.1) then                         !added d rrl
      do k=1,l
      do j=1,mp
      c1p=gcx*amax1(0., u(1,j,k,0))
      c1n=gcx*amin1(0., u(1,j,k,0))
      ux(1,j,k)=         c1n*(u(2,j,k,lc)-u(1,j,k,lc))
     1              +ibcx*c1p*(u(1,j,k,lc)-u(-1,j,k,lc))
      vx(1,j,k)=         c1n*(v(2,j,k,lc)-v(1,j,k,lc))
     1              +ibcx*c1p*(v(1,j,k,lc)-v(-1,j,k,lc))
      ox(1,j,k)=         c1n*(o(2,j,k,lc)-o(1,j,k,lc))
     1              +ibcx*c1p*(o(1,j,k,lc)-o(-1,j,k,lc))
      enddo
      enddo
      endif

      if(rightdedge.eq.1) then                         !added d rrl
      do k=1,l
      do j=1,mp
      cnp=gcx*amax1(0., u(np,j,k,0))
      cnn=gcx*amin1(0., u(np,j,k,0))
      ux(np,j,k)=         cnp*(u(np,j,k,lc)-u(np-1,j,k,lc))
     1              +ibcx*cnn*(u(np+2,j,k,lc)-u(np,j,k,lc))
      vx(np,j,k)=         cnp*(v(np,j,k,lc)-v(np-1,j,k,lc))
     1              +ibcx*cnn*(v(np+2,j,k,lc)-v(np,j,k,lc))
      ox(np,j,k)=         cnp*(o(np,j,k,lc)-o(np-1,j,k,lc))
     1              +ibcx*cnn*(o(np+2,j,k,lc)-o(np,j,k,lc))
      enddo
      enddo
      endif

      if(j3.ne.0) then
      do 2  k=1,l
      do 2  i=1,np
      do j=1+botdedge,mp-topdedge                         !added d rrl
      cp=gcy*amax1(0., v(i,j,k,0))
      cn=gcy*amin1(0., v(i,j,k,0))
      uy(i,j,k)= cp*(u(i,j,k,lc)-u(i,j-1,k,lc))
     .          +cn*(u(i,j+1,k,lc)-u(i,j,k,lc))
      vy(i,j,k)= cp*(v(i,j,k,lc)-v(i,j-1,k,lc))
     .          +cn*(v(i,j+1,k,lc)-v(i,j,k,lc))
      oy(i,j,k)= cp*(o(i,j,k,lc)-o(i,j-1,k,lc))
     .          +cn*(o(i,j+1,k,lc)-o(i,j,k,lc))
      enddo
2     continue
      if(botdedge.eq.1) then                         !added d rrl
      do k=1,l
      do i=1,np
      c1p=gcy*amax1(0., v(i,1,k,0))
      c1n=gcy*amin1(0., v(i,1,k,0))
      uy(i,1,k)=         c1n*(u(i,2,k,lc)-u(i,1,k,lc))
     1             +ibcy*c1p*(u(i,1,k,lc)-u(i,-1,k,lc))
      vy(i,1,k)=         c1n*(v(i,2,k,lc)-v(i,1,k,lc))
     1             +ibcy*c1p*(v(i,1,k,lc)-v(i,-1,k,lc))
      oy(i,1,k)=         c1n*(o(i,2,k,lc)-o(i,1,k,lc ))
     1             +ibcy*c1p*(o(i,1,k,lc)-o(i,-1,k,lc))
      enddo
      enddo
      endif
      if(topdedge.eq.1) then                         !added d rrl
      do k=1,l
      do i=1,np
      cmp=gcy*amax1(0., v(i,mp,k,0))
      cmn=gcy*amin1(0., v(i,mp,k,0))
      uy(i,mp,k)=         cmp*(u(i,mp,k,lc)-u(i,mp-1,k,lc))
     1             +ibcy*cmn*(u(i,mp+2,k,lc)-u(i,mp,k,lc))
      vy(i,mp,k)=         cmp*(v(i,mp,k,lc)-v(i,mp-1,k,lc))
     1             +ibcy*cmn*(v(i,mp+2,k,lc)-v(i,mp,k,lc))
      oy(i,mp,k)=         cmp*(o(i,mp,k,lc)-o(i,mp-1,k,lc))
     1             +ibcy*cmn*(o(i,mp+2,k,lc)-o(i,mp,k,lc))
      enddo
      enddo
      endif
      else
      do i=1,np
      do k=1,l
      uy(i,1,k)=0.
      vy(i,1,k)=0.
      oy(i,1,k)=0.
      enddo
      enddo
      endif

      do 3 k=2,l-1
      do 3 j=1,mp
      do 3 i=1,np
      cp=gcz*amax1(0., o(i,j,k,0))
      cn=gcz*amin1(0., o(i,j,k,0))
      uz(i,j,k)= cp*(u(i,j,k,lc)-u(i,j,k-1,lc))
     .          +cn*(u(i,j,k+1,lc)-u(i,j,k,lc))
      vz(i,j,k)= cp*(v(i,j,k,lc)-v(i,j,k-1,lc))
     .          +cn*(v(i,j,k+1,lc)-v(i,j,k,lc))
      oz(i,j,k)= cp*(o(i,j,k,lc)-o(i,j,k-1,lc))
     .          +cn*(o(i,j,k+1,lc)-o(i,j,k,lc))
    3 continue
      do j=1,mp
      do i=1,np
      c1n=gcz*amin1(0., o(i,j,1,0))
      c1p=gcz*amax1(0., o(i,j,l,0))
      uz(i,j,1)=        c1n*(u(i,j,2,lc)-u(i,j,1,lc))
      uz(i,j,l)=        c1p*(u(i,j,l,lc)-u(i,j,l-1,lc))
      vz(i,j,1)=        c1n*(v(i,j,2,lc)-v(i,j,1,lc))
      vz(i,j,l)=        c1p*(v(i,j,l,lc)-v(i,j,l-1,lc))
      oz(i,j,1)=        c1n*(o(i,j,2,lc)-o(i,j,1,lc))
      oz(i,j,l)=        c1p*(o(i,j,l,lc)-o(i,j,l-1,lc))
      enddo
      enddo
      do k=1,l
      do j=1,mp
      do i=1,np
      u(i,j,k,1)=u(i,j,k,1)-ux(i,j,k)-uy(i,j,k)-uz(i,j,k)
      v(i,j,k,1)=v(i,j,k,1)-vx(i,j,k)-vy(i,j,k)-vz(i,j,k)
      o(i,j,k,1)=o(i,j,k,1)-ox(i,j,k)-oy(i,j,k)-oz(i,j,k)
      enddo
      enddo
      enddo

      do k=1,l
      do j=1,mp
      do i=1,np
      u(i,j,k,1)=u(i,j,k,1)*gcx/gi(i,j,k)
      v(i,j,k,1)=v(i,j,k,1)*gcy/gi(i,j,k)
      o(i,j,k,1)=o(i,j,k,1)*gcz/gi(i,j,k)
      enddo
      enddo
      enddo


        do k=2,l
         do j=1,mp
          do i=1,np
           oa(i,j,k)=0.5*(o(i,j,k,1)+o(i,j,k-1,1))
          enddo
         enddo
        enddo

       do j=1,mp
       do i=1,np
       oa(i,j,1)=oa(i,j,2)
       oa(i,j,l+1)=oa(i,j,l)
       enddo
       enddo

       if(j3.ne.0) then
       call updated(v(1-ih,1-ih,1,1),v(1-ih,1-ih,1,1),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih,0)

        do k=1,l
         do j=1+botdedge,mp                         !added d rrl
          do i=1,np
           va(i,j,k)=0.5*(v(i,j,k,1)+v(i,j-1,k,1))
          enddo
         enddo
        enddo
        if(topdedge.eq.1) then                         !added d rrl
        call updated(va,va,np,mp+1,l,
     .1-ih,np+ih,1-ih,mp+ih+1,0)
        else
        call updated(va,va,np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih+1,0)
        endif
        if(botdedge.eq.1) then                         !added d rrl
        do k=1,l
        do i=1,np
        va(i,1,k)=va(i,-1,k)*ibcy+(1-ibcy)*va(i,2,k)
        enddo
        enddo
        endif
        if(topdedge.eq.1) then                         !added d rrl
        do k=1,l
        do i=1,np
        va(i,mp+1,k)=va(i,mp+3,k)*ibcy+(1-ibcy)*va(i,mp,k)
        enddo
        enddo
        endif
        else
        do i=1,np
        do j=1,mp+1
        do k=1,l
        va(i,j,k)=0.
        enddo
        enddo
        enddo
        endif

       call updated(u(1-ih,1-ih,1,1),u(1-ih,1-ih,1,1),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih,0)
        do k=1,l
         do j=1,mp
          do i=1+leftdedge,np                         !added d rrl
           ua(i,j,k)=0.5*(u(i,j,k,1)+u(i-1,j,k,1))
c          ua(i,j,k)=0.
          enddo
         enddo
        enddo

        if(rightdedge.eq.1) then                         !added d rrl
         call updated(ua,ua,np+1,mp,l,1-ih,np+ih+1,1-ih,mp+ih,0)
        else
         call updated(ua,ua,np,mp,l,1-ih,np+ih+1,1-ih,mp+ih,0)
        endif

        if(leftdedge.eq.1) then                         !added d rrl
        do k=1,l
        do j=1,mp
        ua(1,j,k)=ua(-1,j,k)*ibcx+(1-ibcx)*ua(2,j,k)
        enddo
        enddo
        endif

        if(rightdedge.eq.1) then                         !added d rrl
        do k=1,l
        do j=1,mp
        ua(np+1,j,k)=ua(np+3,j,k)*ibcx+(1-ibcx)*ua(np,j,k)
        enddo
        enddo
        endif
       deallocate (ux)
       deallocate (vx)
       deallocate (ox)
       deallocate (uy)
       deallocate (vy)
       deallocate (oy)
       deallocate (uz)
       deallocate (vz)
       deallocate (oz)

       return
       end subroutine velprdmoa
!*****************************************************************************************!
      subroutine mpdatanew2d(u1,u2,x,h,xe,il,iu,jl,ju,lls,nv2dp)
      use gridsetup
      use msga

      Implicit None

      !JAS 3/7/06 added explicit declarations to comply with implicit none
      integer,intent(in) :: il,iu,jl,ju,lls,nv2dp

      real,dimension(il:iu+1,jl:ju,lls) :: u1
      real,dimension(il:iu,jl:ju,lls+1) :: u2
      real,dimension(il:iu,jl:ju,lls,nv2dp) :: x,xe
      real,dimension(il:iu,jl:ju,lls) :: h
      real,allocatable:: v1(:,:, :,:)
      real,allocatable:: v2(:,  :, :,:)
      real,allocatable:: f1(:,:, :,:)
      real,allocatable:: f2(:,  :, :,:),
     .                   cp(:,:,:,:),
     .                   cn(:,:,:,:),
     .                   mx(:,:,:,:),
     .                   mn(:,:,:,:)
      real,allocatable::   f1o(:,:,:,:),
     .                     f2o(:,:,:,:)
      real,allocatable:: mxo(:, :,:,:),
     .                   mno(:, :,:,:)
      real,allocatable::   a(:, :,:,:)

      !JAS 3/7/06 added explicit declarations to comply with implicit none
      integer :: n1,n1m,n2,n2m,kv,i,j,k,ip,im,jp,jm,itrfct
      integer :: illim,iulim,ibc,ibo
      real :: ep,rhoin,rhoout,ain,aout,rmxuse,rmnuse,a1p,a2p,a1n,a2n
      real :: tmpp,tmpn,tmp,v1d,v2d,c1,c2

      n1=n+1
      n2=l+1
      n1m=n1-1
      n2m=n2-1
      ep=1.e-10

      if(j3.eq.0) then
      allocate (v1(1-ih:np+ih+1,1-ih:mp+ih, l,nv2d))
      allocate (v2(1-ih:np+ih,  1-ih:mp+ih, l+1,nv2d))
      allocate (f1(1-ih:np+ih+1,1-ih:mp+ih, l,nv2d))
      allocate (f2(1-ih:np+ih,  1-ih:mp+ih, l+1,nv2d))
      allocate (cp(1-ih:np+ih, 1-ih:mp+ih, l,nv2d))
      allocate (cn(1-ih:np+ih, 1-ih:mp+ih, l,nv2d))
      allocate (mx(1-ih:np+ih, 1-ih:mp+ih, l,nv2d-1))
      allocate (mn(1-ih:np+ih, 1-ih:mp+ih, l,nv2d-1))
      allocate (mxo(1-ih:np+ih, 1-ih:mp+ih, l,nv2d))
      allocate (mno(1-ih:np+ih, 1-ih:mp+ih, l,nv2d))
      allocate (f1o(1-ih:np+ih+1, 1-ih:mp+ih,l,nv2d))
      allocate (f2o(1-ih:np+ih, 1-ih:mp+ih,l+1,nv2d))
      allocate (a(1-ih:np+ih, 1-ih:mp+ih,l,nv2d-1))
      endif
      f1=0.
      f2=0.

      do kv=1,nv2d
      call update(x(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih,1)
      enddo

      ibc=ibcx
      ibo=1-ibc

      do kv=1,nv2d
      illim = 1
      iulim = np+1
      do j=1,l
        do i=illim,iulim
          v1(i,1,j,kv) = u1(i,1,j)
        end do
      end do

      do i=1,np
        do j=1,n2
          v2(i,1,j,kv) = u2(i,1,j)
        end do
      enddo
      enddo

      if(nonosold.eq.1) then
      do kv=1,nv2d
      do j=1,n2m
      jm=max0(j-1,1  )
      jp=min0(j+1,n2m)
      do i=1,np
         if (leftedge.eq.1 .and. i.eq.1) then
            im = ibc*(-1) + ibo*1
         else
            im = i - 1
         end if
         if (rightedge.eq.1 .and. i.eq.np) then
            ip = ibc*(np+2) + ibo*np
         else
            ip = i + 1
         end if
c      im=ibc*(i-1+(n1-i)/n1m*(n1-2))+ibo*max0(i-1,1  )
c      ip=ibc*(i+1    -i /n1m*(n1-2))+ibo*min0(i+1,n1m)
      mxo(i,1,j,kv)=amax1(x(im,1,j,kv),x(i,1,j,kv),
     .  x(ip,1,j,kv),x(i,1,jm,kv),x(i,1,jp,kv))
      mno(i,1,j,kv)=amin1(x(im,1,j,kv),x(i,1,j,kv),
     .  x(ip,1,j,kv),x(i,1,jm,kv),x(i,1,jp,kv))
      end do
      end do
      enddo
      endif

      if(nonos.eq.1) then
      do kv=1,nv2d-1
      do j=1,n2m
      do i=1,np
      a(i,1,j,kv)=x(i,1,j,kv)/x(i,1,j,nv2d)
      enddo
      enddo
      call update(a(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih,1)
      enddo

      do kv=1,nv2d-1
      do j=1,n2m
      jm=max0(j-1,1  )
      jp=min0(j+1,n2m)
      do i=1,np
c     im=ibcx*(i-1+(n1-i)/n1m*(n1-2))+ibxo*max0(i-1,1  )
c     ip=ibcx*(i+1    -i /n1m*(n1-2))+ibxo*min0(i+1,n1m)
         if (leftedge.eq.1 .and. i.eq.1) then
            im = ibc*(-1) + ibo*1
         else
            im = i - 1
         end if
         if (rightedge.eq.1 .and. i.eq.np) then
            ip = ibc*(np+2) + ibo*np
         else
            ip = i + 1
         end if
      mx(i,1,j,kv)=amax1(a(im,1,j,kv),a(i,1,j,kv),
     .a(ip,1,j,kv),a(i,1,jm,kv),a(i,1,jp,kv))
      mn(i,1,j,kv)=amin1(a(im,1,j,kv),a(i,1,j,kv),
     .a(ip,1,j,kv),a(i,1,jm,kv),a(i,1,jp,kv))
      end do
      end do

      enddo
      endif



      c1=1.
      c2=0.
                         do 3 k=1,iord

      do kv=1,nv2d

      illim = 1 + 1*leftedge
      iulim = np
      do 331 j=1,n2-1
      do 331 i=illim,iulim
  331 f1(i,1,j,kv)=donor(c1*x(i-1,1,j,kv)+c2,c1*x(i,1,j,kv)+c2,
     .                v1(i,1,j,kv))
      call update(x(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih,1)

      if (rightedge.eq.0) then
         call update(f1(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,1)
      else
         call update(f1(1-ih,1-ih,1,kv),np+1,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,1)
      end if

      if (leftedge.eq.1) then
         do j=1,n2-1
            f1(1,1,j,kv)=ibo*donor(c1*xe(1,1,j,kv)+c2,
     .                             c1*x(1,1,j,kv)+c2,v1(1,1,j,kv))
     .           +ibc*f1(-1,1,j,kv)
         end do
      end if
      if (rightedge.eq.1) then
         do j=1,n2-1
            f1(np+1,1,j,kv)=ibo*donor(c1*x(np,1,j,kv)+c2,
     .                                c1*xe(np,1,j,kv)+c2,
     .           v1(np+1,1,j,kv))+ibc*f1(np+3,1,j,kv)
         enddo
      end if

      do 332 j=2,n2-1
      do 332 i=1,np
  332 f2(i,1,j,kv)=donor(c1*x(i,1,j-1,kv)+c2,c1*x(i,1,j,kv)+c2,
     .v2(i,1,j,kv))

        if(ibctopbot.eq.0) then
        do i=1,np
          f2(i,1, 1,kv)=-f2(i,1,  2,kv)
          f2(i,1,n2,kv)=-f2(i,1,n2m,kv)
        end do
        else
        do i=1,np
          f2(i,1, 1,kv)=0.
          f2(i,1,n2,kv)=0.
        end do
        endif

      if (rightedge.eq.0) then
         call update(f1(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,1)
      else
         call update(f1(1-ih,1-ih,1,kv),np+1,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,1)
      end if

      do 333 j=1,n2-1
      do 333 i=1,np
      x(i,1,j,kv)=x(i,1,j,kv)-(f1(i+1,1,j,kv)-f1(i,1,j,kv)
     .                        +f2(i,1,j+1,kv)-f2(i,1,j,kv))/
     .         h(i,1,j)
c     if(mpi_rank.eq.0.and.(j.eq.1.and.kv.eq.1)) 
c    .print*,x(i,1,j,kv),u1(i,1,j)
 333  continue

      call update(x(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih,1)
      enddo
c     do i=1,np
c      print*,x(i,1,1,nv2d),u1(i,1,1),u1(i+1,1,1)
c     enddo
c     call combine(x(1-ih,1-ih,1,1),100,tmpa)
c     if (mpi_rank.eq.0) then      
c     nh=25
c     do j=1,l
c     is=0
c     do i=1,n
c     if(tmpa(i,j).le.0.) diff=tmpa(n-is,j)+tmpa(i,j)
c     if(tmpa(i,j).ge.0.) diff=tmpa(n-is,j)-tmpa(i,j)
c     is=is+1
c     enddo
c     enddo
c     endif
c     print*,'yes'


      if(k.eq.iord) go to 6
      c1=0.
      c2=1.
      do kv=1,nv2d

      illim = 1
      iulim = np + 1*rightedge
      do 49 j=1,n2-1
      do 49 i=illim,iulim
      f1(i,1,j,kv)=v1(i,1,j,kv)
   49 v1(i,1,j,kv)=0.


      do 50 j=1,n2
      do 50 i=1,np
      f2(i,1,j,kv)=v2(i,1,j,kv)
   50 v2(i,1,j,kv)=0.

      if (rightedge.eq.0) then
         call update(f1(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,1)
      else
         call update(f1(1-ih,1-ih,1,kv),np+1,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,1)
      end if
      call update(f2(1-ih,1-ih,1,kv),np,mp,l+1,
     .1-ih,np+ih,1-ih,mp+ih,1)
      call update(x(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih,1)
      illim = 1 + 1*leftedge
      iulim = np
      do 51 j=2,n2-2
      do 51 i=illim,iulim
   51 v1(i,1,j,kv)=vdyf(x(i-1,1,j,kv),x(i,1,j,kv),f1(i,1,j,kv),.5*
     .                 (h(i-1,1,j)+h(i,1,j)))
     * +vcorr(f1(i,1,j,kv), f2(i-1,1,j,kv)
     .       +f2(i-1,1,j+1,kv)+f2(i,1,j+1,kv)+
     .        f2(i,1,j,kv),
     *        x(i-1,1,j-1,kv),x(i,1,j-1,kv),
     .        x(i-1,1,j+1,kv),x(i,1,j+1,kv),
     *        .5*(h(i-1,1,j)+h(i,1,j)))
c     if(kv.eq.1) then
c     if(mpi_rank.eq.0) print*,v1(7,9,4,1),'v1'
c     if(mpi_rank.eq.0) print*,v1(8,9,4,1),'v1'
c     endif
      if(idiv.eq.1) then
      illim = 1 + 1*leftedge
      iulim = np
      do 511 j=2,n2-2
      do 511 i=illim,iulim
      v1d=-vdiv1(f1(i-1,1,j,kv),f1(i,1,j,kv),f1(i+1,1,j,kv),.5*
     .           (h(i-1,1,j)+h(i,1,j)))
     *    -vdiv2(f1(i,1,j,kv),f2(i-1,1,j+1,kv),
     .           f2(i,1,j+1,kv),f2(i-1,1,j,kv),
     *           f2(i,1,j,kv), .5*(h(i-1,1,j)+h(i,1,j)))
  511 v1(i,1,j,kv)=v1(i,1,j,kv)+(pp(v1d)*x(i-1,1,j,kv)
     .                          -pn(v1d)*x(i,1,j,kv))
      endif

      illim = 1  + 1*leftedge
      iulim = np - 1*rightedge
      do 52 j=2,n2-1
      do 52 i=illim,iulim
   52 v2(i,1,j,kv)=vdyf(x(i,1,j-1,kv),x(i,1,j,kv),f2(i,1,j,kv),.5*
     .              (h(i,1,j-1)+h(i,1,j)))
     * +vcorr(f2(i,1,j,kv), f1(i,1,j-1,kv)
     .       +f1(i,1,j,kv)+f1(i+1,1,j,kv)+
     .        f1(i+1,1,j-1,kv),
     *         x(i-1,1,j-1,kv),x(i-1,1,j,kv),
     *         x(i+1,1,j-1,kv),x(i+1,1,j,kv),
     *        .5*(h(i,1,j-1)+h(i,1,j)))

      if(ibc.eq.1) then
         if (leftedge.eq.1) then
            do j=2,n2-1
               v2(1,1,j,kv)=vdyf(x(1,1,j-1,kv),x(1,1,j,kv),
     .         f2(1,1,j,kv),.5*
     .                       (h(1,1,j-1)+h(1,1,j)))
     *              +vcorr(f2(1,1,j,kv),
     .                     f1(1,1,j-1,kv)+f1(1,1,j,kv)+f1(2,1,j,kv)+
     .                     f1(2,1,j-1,kv),x(-1,1,j-1,kv),x(-1,1,j,kv),
     .                      x(2,1,j-1,kv),x(2,1,j,kv),
     *                  .5*(h(1,1,j-1)+h(1,1,j)))
            end do
         end if
         call update(v2(1-ih,1-ih,1,kv),np,mp,l+1,
     .1-ih,np+ih,1-ih,mp+ih,1)
         if (rightedge.eq.1) then
            do j=2,n2-1
               v2(np,1,j,kv)=v2(np+1,1,j,kv)
            enddo
         end if
      endif

      if(idiv.eq.1) then
         illim = 1  + (1-ibc)*leftedge
         iulim = np + (ibc-1)*rightedge
      do 521 j=2,n2-1
      do 521 i=illim,iulim
      v2d=-vdiv1(f2(i,1,j-1,kv),f2(i,1,j,kv),f2(i,1,j+1,kv),
     .           .5*(h(i,1,j-1)+h(i,1,j)))
     *    -vdiv2(f2(i,1,j,kv),f1(i+1,1,j-1,kv),
     .           f1(i+1,1,j,kv),f1(i,1,j-1,kv),
     .           f1(i,1,j,kv),.5*(h(i,1,j-1)+h(i,1,j)))
  521 v2(i,1,j,kv)=v2(i,1,j,kv)
     .    +(pp(v2d)*x(i,1,j-1,kv)-pn(v2d)*x(i,1,j,kv))
      endif

      if(isor.eq.3) then
         illim = 1  + 2*leftedge
         iulim = np - 1*rightedge
      do 61 j=2,n2-2
      do 61 i=illim,iulim
   61 v1(i,1,j,kv)=v1(i,1,j,kv)     +vcor31(f1(i,1,j,kv),
     1        x(i-2,1,j,kv),x(i-1,1,j,kv),x(i,1,j,kv),x(i+1,1,j,kv),
     .        .5*(h(i-1,1,j)+h(i,1,j)))
c
      if(ibc.eq.1) then
         if (leftedge.eq.1) then
            do j=2,n2-2
               v1(2,1,j,kv)=v1(2,1,j,kv) +vcor31(f1(2,1,j,kv),
     1              x(-1,1,j,kv),x(1,1,j,kv),x(2,1,j,kv),x(3,1,j,kv),
     .              .5*(h(1,1,j)+h(2,1,j)))
            end do
         end if
         if (rightedge.eq.1) then
            do j=2,n2-2
               v1(np,1,j,kv)=v1(np,1,j,kv)
     .    +vcor31(f1(np,1,j,kv),x(np-2,1,j,kv),
     .          x(np-1,1,j,kv),
     1          x(np,1,j,kv),x(np+2,1,j,kv),.5*(h(np-1,1,j)+h(np,1,j)))
            enddo
         end if
      endif

      illim = 1  + (2-ibc)*leftedge
      iulim = np + (ibc-1)*rightedge
      do 62 j=2,n2-2
      do 62 i=illim,iulim
   62 v1(i,1,j,kv)=v1(i,1,j,kv)
     1 +vcor32(f1(i,1,j,kv),f2(i-1,1,j,kv)
     1                     +f2(i-1,1,j+1,kv)+f2(i,1,j+1,kv)+
     .                      f2(i,1,j,kv),
     *                       x(i,1,j-1,kv),x(i-1,1,j+1,kv),
     *                       x(i-1,1,j-1,kv),x(i,1,j+1,kv),
     *                   .5*(h(i-1,1,j)+h(i,1,j)))


      illim = 1  + (1-ibc)*leftedge
      iulim = np + (ibc-1)*rightedge
      do 63 j=3,n2-2
      do 63 i=illim,iulim
   63 v2(i,1,j,kv)=v2(i,1,j,kv)     +vcor31(f2(i,1,j,kv),
     1              x(i,1,j-2,kv),x(i,1,j-1,kv),
     .              x(i,1,j,kv),x(i,1,j+1,kv),
     .        .5*(h(i,1,j-1)+h(i,1,j)))


      illim = 1  + 1*leftedge
      iulim = np - 1*rightedge
      do 64 j=3,n2-2
      do 64 i=illim,iulim
   64 v2(i,1,j,kv)=v2(i,1,j,kv)
     1 +vcor32(f2(i,1,j,kv),f1(i,1,j-1,kv)
     .        +f1(i+1,1,j-1,kv)+f1(i+1,1,j,kv)+
     .         f1(i,1,j,kv),
     *          x(i+1,1,j-1,kv),x(i-1,1,j,kv),
     .          x(i-1,1,j-1,kv),x(i+1,1,j,kv),
     *                  .5*(h(i,1,j-1)+h(i,1,j)))
      if(ibc.eq.1) then
         if (leftedge.eq.1) then
            do j=3,n2-2
               v2(1,1,j,kv)=v2(1,1,j,kv)
     1              +vcor32(f2(1,1,j,kv),
     .                      f1(1,1,j-1,kv)+f1(2,1,j-1,kv)+
     .                      f1(2,1,j,kv)+f1(1,1,j,kv),
     *                       x(2,1,j-1,kv),x(-1,1,j,kv),
     *                       x(-1,1,j-1,kv),x(2,1,j,kv),
     *              .5*(h(1,1,j-1)+h(1,1,j)))
            end do
         end if
         call update(v2(1-ih,1-ih,1,kv),np,mp,l+1,
     .1-ih,np+ih,1-ih,mp+ih,1)
         if (rightedge.eq.1) then
            do j=3,n2-2
               v2(np,1,j,kv)=v2(np+1,1,j,kv)
            end do
         end if
      endif
      endif

      if (ibc.eq.1) then
         if (rightedge.eq.0) then
            call update(v1(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,1)
         else
            call update(v1(1-ih,1-ih,1,kv),np+1,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,1)
         end if
         if (leftedge.eq.1) then
            do j=1,n2m
               v1( 1,1,j,kv)=v1(-1,1,j,kv)
            end do
         end if
c     the following update is probably unnecessary ; when I have time
c     I'll check to make sure and remove it if that's the case.
         if (rightedge.eq.0) then
            call update(v1(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,1)
         else
            call update(v1(1-ih,1-ih,1,kv),np+1,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,1)
         end if
         if (rightedge.eq.1) then
            do j=1,n2m
               v1(np+1,1,j,kv)=v1(np+3,1,j,kv)
            end do
         end if
      end if
      enddo
c     if(mpi_rank.eq.0) print*,v1(25,1,3,4),v1(24,1,3,4),'v1'
c     if(mpi_rank.eq.0) print*,f1(25,1,3,4),f1(24,1,3,4),'v1'

                  if(nonosold.eq.1) then
c                 non-osscilatory option
      do kv=1,nv2d
      do 401 j=1,n2m
      jm=max0(j-1,1  )
      jp=min0(j+1,n2m)
      do 401 i=1,np
         if (leftedge.eq.1 .and. i.eq.1) then
            im = ibc*(-1) + ibo*1
         else
            im = i - 1
         end if
         if (rightedge.eq.1 .and. i.eq.np) then
            ip = ibc*(np+2) + ibo*np
         else
            ip = i + 1
         end if
c      im=ibc*(i-1+(n1-i)/n1m*(n1-2))+ibo*max0(i-1,1  )
c      ip=ibc*(i+1    -i /n1m*(n1-2))+ibo*min0(i+1,n1m)
      mxo(i,1,j,kv)=amax1(x(im,1,j,kv),x(i,1,j,kv),x(ip,1,j,kv),
     .         x(i,1,jm,kv),x(i,1,jp,kv),
     .        mxo(i,1,j,kv))
  401 mno(i,1,j,kv)=amin1(x(im,1,j,kv),x(i,1,j,kv),x(ip,1,j,kv),
     .         x(i,1,jm,kv),x(i,1,jp,kv),
     .        mno(i,1,j,kv))

      illim = 1
      iulim = np + 1*rightedge
      do 402 j=1,n2m
      do 402 i=illim,iulim
  402 f1(i,1,j,kv)=donor(c2,c2,v1(i,1,j,kv))
      do 403 j=1,n2
      do 403 i=1,np
  403 f2(i,1,j,kv)=donor(c2,c2,v2(i,1,j,kv))

      if (rightedge.eq.0) then
         call update(f1(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,1)
      else
         call update(f1(1-ih,1-ih,1,kv),np+1,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,1)
      end if
      call update(f2(1-ih,1-ih,1,kv),np,mp,l+1,
     .1-ih,np+ih,1-ih,mp+ih,1)

      do 404 j=1,n2m
      do 404 i=1,np
      if(j.eq.1) f2(i,1,j,kv)=-f2(i,1,j+1,kv)
      if(j.eq.n2m) f2(i,1,j+1,kv)=-f2(i,1,j,kv)
      cp(i,1,j,kv)=(mxo(i,1,j,kv)-x(i,1,j,kv))*h(i,1,j)/
     1(pn(f1(i+1,1,j,kv))+pp(f1(i,1,j,kv))
     1+pn(f2(i,1,j+1,kv))+pp(f2(i,1,j,kv))+ep)
      cn(i,1,j,kv)=(x(i,1,j,kv)-mno(i,1,j,kv))*h(i,1,j)/
     1(pp(f1(i+1,1,j,kv))+pn(f1(i,1,j,kv))
     1+pp(f2(i,1,j+1,kv))+pn(f2(i,1,j,kv))+ep)
  404 continue

      call update(cp(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih,1)
      call update(cn(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih,1)

      illim = 1 + 1*leftedge
      iulim = np
      do j=1,n2m
        do i=illim,iulim
          v1(i,1,j,kv)=pp(v1(i,1,j,kv))
     *             *amin1(1.,cp(i,1,j,kv),cn(i-1,1,j,kv))
     *                -pn(v1(i,1,j,kv))
     *             *amin1(1.,cp(i-1,1,j,kv),cn(i,1,j,kv))
        end do
      end do
      if (ibc.eq.1) then
          if (rightedge.eq.0) then
             call update(v1(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,1)
          else
             call update(v1(1-ih,1-ih,1,kv),np+1,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,1)
          end if
         if (leftedge.eq.1) then
            do j=1,n2m
               v1( 1,1,j,kv)=v1(-1,1,j,kv)
            end do
         end if
c     following update is probably unnecessary
          if (rightedge.eq.0) then
             call update(v1(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,1)
          else
             call update(v1(1-ih,1-ih,1,kv),np+1,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,1)
          end if
         if (rightedge.eq.1) then
            do j=1,n2m
               v1(np+1,1,j,kv)=v1(np+3,1,j,kv)
            end do
         end if
      end if
      do j=2,n2m
        do i=1,np
          v2(i,1,j,kv)=pp(v2(i,1,j,kv))
     *             *amin1(1.,cp(i,1,j,kv),cn(i,1,j-1,kv))
     *                -pn(v2(i,1,j,kv))
     *             *amin1(1.,cp(i,1,j-1,kv),cn(i,1,j,kv))
        end do
      end do
      enddo
                  endif

      if(nonos.eq.1) then
      do 1000 itrfct=1,nfct

      if(itrfct.eq.1) then
      do kv=1,nv2d
      if (rightedge.eq.0) then
         call update(v1(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,1)
      else
         call update(v1(1-ih,1-ih,1,kv),np+1,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,1)
      end if

      illim = 1
      iulim = np + 1*rightedge
c     iulim = np + 1
      do 502 j=1,n2m
      do 502 i=illim,iulim
  502 f1(i,1,j,kv)=donor(c2,c2,v1(i,1,j,kv))

      do 5033 j=1,n2
      do 5033 i=1,np
 5033 f2(i,1,j,kv)=donor(c2,c2,v2(i,1,j,kv))
      if (rightedge.eq.0) then
         call update(f1(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,1)
      else
         call update(f1(1-ih,1-ih,1,kv),np+1,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,1)
      end if
      call update(f2(1-ih,1-ih,1,kv),np,mp,l+1,
     .1-ih,np+ih,1-ih,mp+ih,1)
      enddo
      endif
correction coefficients for variables nv-1
      do kv=1,nv2d-1
      do 504 j=1,n2m
      do 504 i=1,np
      if(abs(mx(i,1,j,kv)).lt.ep) mx(i,1,j,kv)=0.
      if(abs(mn(i,1,j,kv)).lt.ep) mn(i,1,j,kv)=0.
      if(j.eq.1) f2(i,1,j,nv2d)=-f2(i,1,j+1,nv2d)
      if(j.eq.1) f2(i,1,j,kv)=-f2(i,1,j+1,kv)
      if(j.eq.n2m) f2(i,1,j+1,nv2d)=-f2(i,1,j,nv2d)
      if(j.eq.n2m) f2(i,1,j+1,kv)=-f2(i,1,j,kv)
      rhoin=
     1pn(f1(i+1,1,j,nv2d))+pp(f1(i,1,j,nv2d))+
     1pn(f2(i,1,j+1,nv2d))+pp(f2(i,1,j,nv2d))
      rhoout=-
     1(pp(f1(i+1,1,j,nv2d))+pn(f1(i,1,j,nv2d))+
     1 pp(f2(i,1,j+1,nv2d))+pn(f2(i,1,j,nv2d)))
      aout=-
     1(pp(f1(i+1,1,j,kv))+pn(f1(i,1,j,kv))
     1+pp(f2(i,1,j+1,kv))+pn(f2(i,1,j,kv)))
       ain=
     1 pn(f1(i+1,1,j,kv))+pp(f1(i,1,j,kv))
     1+pn(f2(i,1,j+1,kv))+pp(f2(i,1,j,kv))
c     if(i.eq.24.and.j.eq.3.and.(mpi_rank.eq.0))
c    .print*,ain,aout,rhoout,rhoin,'1'
c     if(i.eq.24.and.j.eq.3.and.(mpi_rank.eq.0))
c    .print*,f1(i+1,1,j,kv),f1(i,1,j,kv)
c     if(i.eq.24.and.j.eq.3.and.(mpi_rank.eq.0))
c    .print*,v1(i+1,1,j,kv),v1(i,1,j,kv)

      cp(i,1,j,kv)=pp(mx(i,1,j,kv)*x(i,1,j,nv2d)-x(i,1,j,kv))*h(i,1,j)/
     1(ain-pp(mx(i,1,j,kv))*rhoout+pn(mx(i,1,j,kv))*rhoin+ep)

504   cn(i,1,j,kv)=pp(x(i,1,j,kv)-mn(i,1,j,kv)*x(i,1,j,nv2d))*h(i,1,j)/
     1(-aout+pp(mn(i,1,j,kv))*rhoin-pn(mn(i,1,j,kv))*rhoout+ep)

      call update(cp(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih,1)
      call update(cn(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih,1)
      enddo

correction coeffiecients for rho
      do kv=1,nv2d-1
      do j=1,n2m
      do i=1,np
      if(abs(mx(i,1,j,kv)).lt.ep) mx(i,1,j,kv)=0.
      if(abs(mn(i,1,j,kv)).lt.ep) mn(i,1,j,kv)=0.
c  JLW - Fix from Jon Reisner 9/2005
c     a1p=cp(i,1,j,kv)+amax1(0.,sign(1.,mx(i,1,j,kv)))
c     a2p=cn(i,1,j,kv)+amax1(0.,sign(1.,-mn(i,1,j,kv)))
c     a1n=cn(i,1,j,kv)+amax1(0.,sign(1.,mn(i,1,j,kv)))
c     a2n=cp(i,1,j,kv)+amax1(0.,sign(1.,-mx(i,1,j,kv)))
      rmxuse=-1.*mx(i,1,j,kv)
      if(mx(i,1,j,kv).eq.0.) rmxuse=0.
      rmnuse=-1.*mn(i,1,j,kv)
      if(mn(i,1,j,kv).eq.0.) rmnuse=0.
      a1p=cp(i,1,j,kv)+amax1(0.,sign(1., mx(i,1,j,kv)))
      a2p=cn(i,1,j,kv)+amax1(0.,sign(1.,rmnuse))
      a1n=cn(i,1,j,kv)+amax1(0.,sign(1., mn(i,1,j,kv)))
      a2n=cp(i,1,j,kv)+amax1(0.,sign(1.,rmxuse))
      tmpp=amin1(a1p,a2p)
      tmpn=amin1(a1n,a2n)
      if(kv.eq.1) cp(i,1,j,nv2d)=tmpp
      if(kv.eq.1) cn(i,1,j,nv2d)=tmpn
      cp(i,1,j,nv2d)=amin1(tmpp,cp(i,1,j,nv2d))
      cn(i,1,j,nv2d)=amin1(tmpn,cn(i,1,j,nv2d))
c     if(i.eq.24.and.j.eq.3.and.(mpi_rank.eq.0))
c    .print*,cp(24,1,3,nv2d),cn(24,1,3,nv2d)
c     if(i.eq.24.and.j.eq.3.and.(mpi_rank.eq.0))
c    .print*,cp(24,1,3,kv),cn(24,1,3,kv)
c     cp(i,j,nv)=amax1(0.,cp(i,j,nv))
c     cn(i,j,nv)=amax1(0.,cn(i,j,nv))
      enddo
      enddo
      enddo
      call update(cp(1-ih,1-ih,1,nv2d),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih,1)
      call update(cn(1-ih,1-ih,1,nv2d),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih,1)


      do kv=1,nv2d
c     if (mpi_rank.eq.0.and.kv.eq.nv2d) 
c    .print*,v1(24,1,3,kv),v1(25,1,3,kv)
c     if (mpi_rank.eq.0.and.kv.eq.nv2d) 
c    .print*,v2(24,1,3,kv),v2(24,1,4,kv)
c     if (mpi_rank.eq.0.and.kv.eq.nv2d) 
c    .print*,cp(24,1,3,kv),cn(24,1,3,kv)
c     if (mpi_rank.eq.0.and.kv.eq.nv2d) 
c    .print*,cp(23,1,3,kv),cn(23,1,3,kv)
      do j=1,n2m
          do i=1+leftedge,np
          v1(i,1,j,kv)=pp(v1(i,1,j,kv))
     .*amin1(1.,cp(i,1,j,kv),cn(i-1,1,j,kv))
     *                -pn(v1(i,1,j,kv))
     .*amin1(1.,cp(i-1,1,j,kv),cn(i,1,j,kv))
          end do
        end do
c     if (mpi_rank.eq.0.and.kv.eq.nv2d) 
c    .print*,v1(24,1,3,kv),v1(25,1,3,kv)
c     if (mpi_rank.eq.0.and.kv.eq.nv2d) 
c    .print*,v2(24,1,3,kv),v2(24,1,4,kv)
c     if (mpi_rank.eq.0.and.kv.eq.nv2d) 
c    .print*,x(24,1,3,kv)
      if (ibcx.eq.1) then
         if (rightedge.eq.0) then
            call update(v1(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,1)
         else
            call update(v1(1-ih,1-ih,1,kv),np+1,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,1)
         end if
         do j=1,n2m
               if (leftedge.eq.1) then
                  v1(1 ,1,j,kv)=v1(-1,1,j,kv)
               end if
               if (rightedge.eq.1) then
                  v1(np+1,1,j,kv)=v1(np+3  ,1,j,kv)
               end if
         end do
      end if
c     if (rightedge.eq.0) then
c        call update(v1(1-ih,1-ih,1,kv),np+1,mp,l,np,mp,1)
c     else
c        call update(v1(1-ih,1-ih,1,kv),np+1,mp,l,np+1,mp,1)
c     end if

      do j=2,n2m
          do i=1,np
            v2(i,1,j,kv)= pp(v2(i,1,j,kv))
     *             *amin1(1.,cp(i,1,j,kv),cn(i,1,j-1,kv))
     *                   -pn(v2(i,1,j,kv))
     *             *amin1(1.,cp(i,1,j-1,kv),cn(i,1,j,kv))
          end do
        end do
c     do j=1,n2
c     do i=1,np
c     f2(i,1,j,kv)=0.
c     enddo
c     enddo
c     do j=1,n2m
c     do i=1,np+1
c     f1(i,1,j,kv)=0.
c     enddo
c     enddo
      end do

      if(itrfct.lt.nfct) then
      print*,'yes'
      do kv=1,nv2d
      do 602 j=1,n2m
      do 602 i=1,np+1*rightedge
      tmp=f1(i,1,j,kv)
      f1o(i,1,j,kv)=donor(c2,c2,v1(i,1,j,kv))
  602 f1(i,1,j,kv)=tmp-f1o(i,1,j,kv)
      do 6033 j=1,n2
      do 6033 i=1,np
      tmp=f2(i,1,j,kv)
      f2o(i,1,j,kv)=donor(c2,c2,v2(i,1,j,kv))
 6033 f2(i,1,j,kv)=tmp-f2o(i,1,j,kv)

      if(rightedge.eq.0) then
       call update(f1(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,1)
       call update(f1o(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,1)
      else
       call update(f1(1-ih,1-ih,1,kv),np+1,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,1)
       call update(f1o(1-ih,1-ih,1,kv),np+1,mp,l,
     .1-ih,np+ih+1,1-ih,mp+ih,1)
      endif

      do j=1,n2m
      do i=1,np
      x(i,1,j,kv)=x(i,1,j,kv)-( f1o(i+1,1,j,kv)-f1o(i,1,j,kv)
     .                         +f2o(i,1,j+1,kv)-f2o(i,1,j,kv) )
     .                          /h(i,1,j)
      enddo
      enddo
      enddo
      endif
1000  continue
      endif
    3                      continue
    6 continue

      do kv=1,nv2d
      call update(x(1-ih,1-ih,1,kv),np,mp,l,
     .1-ih,np+ih,1-ih,mp+ih,1)
      enddo

      deallocate (v1)
      deallocate (v2)
      deallocate (f1)
      deallocate (f2)
      deallocate (cp)
      deallocate (cn)
      deallocate (mx)
      deallocate (mn)
      deallocate (mxo)
      deallocate (mno)
      deallocate (f1o)
      deallocate (f2o)
      deallocate (a)

      return
      end subroutine mpdatanew2d

!*****************************************************************************************!
      real function pp(y)
      Implicit None
      real :: y

      pp = amax1(0.0,y)

      end function pp
      !*************************************
      real function pn(y)
      Implicit None
      real :: y

      pn = -amin1(0.0,y)

      end function pn
      !*************************************
      real function donor(y1,y2,a10)
      Implicit None
      real :: y1,y2,a10

      donor = pp(a10)*y1-pn(a10)*y2

      end function donor
      !*************************************
      real function vdiv1(a1,a2,a3,r)
      Implicit None
      real :: a1,a2,a3,r

      vdiv1 = 0.25*a2*(a3-a1)/r

      end function vdiv1
      !*************************************
      real function vdiv2(aa,b1,b2,b3,b4,r)
      Implicit None
      real :: aa,b1,b2,b3,b4,r

      vdiv2 = 0.25*aa*(b1+b2-b3-b4)/r

      end function vdiv2
      !*************************************
      real function rat2(z1,z2)
      Implicit None
      real :: z1,z2

      rat2 = (z2-z1)*0.5

      end function rat2
      !*************************************
      real function rat4(z0,z1,z2,z3)
      Implicit None
      real :: z0,z1,z2,z3

      rat4 = (z3+z2-z1-z0)*0.25

      end function rat4
      !*************************************
      real function vdyf(x1,x2,aa,r)
      Implicit None
      real :: x1,x2,aa,r

      vdyf = (abs(aa)-aa**2/r)*rat2(x1,x2)

      end function vdyf
      !*************************************
      real function vcorr(aa,b,y0,y1,y2,y3,r)
      Implicit None
      real :: aa,b,y0,y1,y2,y3,r

      vcorr = -0.125*aa*b/r*rat4(y0,y1,y2,y3)

      end function vcorr
      !*************************************
      real function vcor31(aa,x0,x1,x2,x3,r)
      Implicit None
      real :: aa,x0,x1,x2,x3,r

      vcor31 = -(aa -3.*abs(aa)*aa/r+2.*aa**3/r**2)/3.
     & *rat4(x1,x2,x0,x3)

      end function vcor31
      !*************************************
      real function vcor32(aa,b,y0,y1,y2,y3,r)
      Implicit None
      real :: aa,b,y0,y1,y2,y3,r

      vcor32 = 0.25*b/r*(abs(aa)-2.*aa**2/r)*rat4(y0,y1,y2,y3)

      end function vcor32
      !*************************************
!*****************************************************************************************!
      end module higrad
!*****************************************************************************************!
