!***************************************************
       module higrad
       Implicit None

       save

       contains
       
       !-rmoainner: higrad main routine, see just below
       !-donorcell: donorcell advection scheme for individual variable
       !-gradmoa: pressure force for moa
       
       
!*****************************************************************************************!
!  RMOAINNER is the MOA inner loop on small time steps
!    The method of average does the time integration on the xv field
!    using a donor cell scheme for advection, and computing pressure, gravity,
!    coriolis, lspgf, and relaxation to ambiant variable
!    . the sum of the forcing (except the advection) is added in f1avg, f2avg,
!    f3avg (using weights wtf)
!    . the mean advective velocities uavg, vavg and oavg are averaged centered
!    contravariant velocities / gi for advection scheme on large time
!    step
!***********************************************************************

      subroutine rmoainner()
      use gridsetup
      use xvi
      use xve
      use advo
      use forcinner
      use workavg
      use constants
      use metryic
      use pres
      use lspgf
      use relax
      use weights
      use msga
      use advectionScheme
      use xvo  ! to get rhovapor fraction
      Implicit None

      integer :: its,i,j,k,kv,kve, kref 
      real ::rhovaporFrac

      real,allocatable:: u(:,:,:,:),  ! centered contravariant coordinates/gi 
     .                   v(:,:,:,:),
     .                   o(:,:,:,:)
               ! last index 0 is t-dts, 1 is current time t, 2 is t + 0.5 *  dts

      real :: uc,vc ! reference wind for coriolis
      real :: piOver180, g1,g2, g3 ! projected gravity coeff

      piOver180 = pi/180.0
      g1 = g*sin(slopeangle*piOver180)*cos(slopeazimuth*piOver180)
      g2 = g*sin(slopeangle*piOver180)*sin(slopeazimuth*piOver180)
      g3 = g*cos(slopeangle*piOver180)
       
      allocate (u(1-ih:np+ih, 1-ih:mp+ih,l,0:2)) ! ghost cell required here to compute u1,u2,u3
      allocate (v(1-ih:np+ih, 1-ih:mp+ih,l,0:2))
      allocate (o(1-ih:np+ih, 1-ih:mp+ih,l,0:2))
      !FP added this section for computing of gas thermal properties,
      !which are assumed to be constant within the higrad loop
      do k=1,l
        do j=1,mp
          do i=1,np
            rhovaporFrac=real(irhovapor)*xvb(i,j,k,8)/xvb(i,j,k,nv)
            call updateGasThermalProperties(rhovaporFrac)
            rg_over_prrcp(i,j,k)=rg_over_prrcp_gas
            cp_over_cv(i,j,k)=cp_over_cv_gas
          enddo
        enddo
      enddo
      ! initialization of u,v,w (for last index 1), f1avg,uavg...
      f1avg=0.0
      f2avg=0.0
      f3avg=0.0
      do k=1,l
        do j=1,mp
         do i=1,np
           u(i,j,k,1)=xv(i,j,k,1)/xv(i,j,k,5)/gi(i,j,k)
           v(i,j,k,1)=xv(i,j,k,2)/xv(i,j,k,5)/gi(i,j,k)
           o(i,j,k,1)=(xv(i,j,k,3)*gi(i,j,k)
     &               +c13(i,j)*gmul(k)*xv(i,j,k,1)
     &       +c23(i,j)*gmul(k)*xv(i,j,k,2))/xv(i,j,k,5)/gi(i,j,k)
           uavg(i,j,k)=wt(1)*u(i,j,k,1)
           vavg(i,j,k)=wt(1)*v(i,j,k,1)
           wavg(i,j,k)=wt(1)*xv(i,j,k,3)/xv(i,j,k,5)/gi(i,j,k)
         enddo
       enddo
      enddo
      !if (mpi_rank.eq.0) write(6,*) 'begin inner loop'

      ! beginning of the inner loop
      do its=1,nts
       time=time+dts
       !if (mpi_rank.eq.0) write(6,*) 'inner loop : iteration ',its
       do k=1,l
         do j=1,mp
           do i=1,np
             ! past u,v,o at t-dts: NB for its=1, it is for t...
             u(i,j,k,0)=u(i,j,k,1)
             v(i,j,k,0)=v(i,j,k,1)
             o(i,j,k,0)=o(i,j,k,1)
             ! current u,v, o
             u(i,j,k,1)=xv(i,j,k,1)/xv(i,j,k,5)/gi(i,j,k)
             v(i,j,k,1)=xv(i,j,k,2)/xv(i,j,k,5)/gi(i,j,k)
             o(i,j,k,1)=(xv(i,j,k,3)*gi(i,j,k)
     &               +c13(i,j)*gmul(k)*xv(i,j,k,1)
     &         +c23(i,j)*gmul(k)*xv(i,j,k,2))/xv(i,j,k,5)/gi(i,j,k)
             ! time integrated u,v,o at t+0.5*dts
             u(i,j,k,2)=0.5*(3.*u(i,j,k,1)-u(i,j,k,0))
             v(i,j,k,2)=0.5*(3.*v(i,j,k,1)-v(i,j,k,0))
             o(i,j,k,2)=0.5*(3.*o(i,j,k,1)-o(i,j,k,0))
           enddo
         enddo
       enddo
      ! computation of advective velocities at t+0.5dts: u1,u2,u3
      call setAdvectiveVelocities(u(1-ih,1-ih,1,2),v(1-ih,1-ih,1,2),o(1-ih,1-ih,1,2)
     +               ,gc1s,gc2s,gc3s,1-ih,np+ih,1-ih,mp+ih,l,0)

      ! computation of advective tendencies (donorcell scheme)

      do kv=1,5
        kve=kv
        if (kv==5) kve=nv
        call updated(xv(1-ih,1-ih,1,kv),xe(1-ih,1-ih,1,kve),np,mp,l,1-ih,np+ih,1-ih,mp+ih,0,0)
        call donorcell(xv(1-ih,1-ih,1,kv),1-ih,np+ih,1-ih,mp+ih,l)
      enddo 

      ! compute pressure force and put it in f1,f2,f3
      call gradmoa()
      if(its.eq.nts) call rmaxmin1(pr,'pr',1-ih,np+ih,1-ih,mp+ih,l)

      ! forces : gravity
      do k=1,l
        do j=1,mp
          do i=1,np
            f1(i,j,k)=f1(i,j,k)-dts*(xv(i,j,k,5)-xe(i,j,k,nv))*g1
            f2(i,j,k)=f2(i,j,k)-dts*(xv(i,j,k,5)-xe(i,j,k,nv))*g2
            f3(i,j,k)=f3(i,j,k)-dts*(xv(i,j,k,5)-xe(i,j,k,nv))*g3
          enddo
        enddo
      enddo

      ! forces : coriolis
      if(icorio.ge.1)then
        do k=1,l
         do j=1,mp
          do i=1,np
           kref = k+(l-k)*(icorio-1) ! k or l depending icorio=1 or 2
           uc=xe(i,j,kref,1)
           vc=xe(i,j,kref,2)
           f1(i,j,k)= dts*(fcor3*(xv(i,j,k,2)-vc)
     &                    -fcor2*(xv(i,j,k,3)) )+f1(i,j,k)
           f2(i,j,k)=-dts*(fcor3*(xv(i,j,k,1)-uc) )+f2(i,j,k)
           f3(i,j,k)= dts*(fcor2*(xv(i,j,k,1)-uc))+f3(i,j,k)
          enddo
         enddo
        enddo
      end if !icorio.ge.1
      ! forces : lspgf
      if (ilspgf.ge.1) then !ilspgf
        do k=1,l
         do j=1,mp
          do i=1,np
           f1(i,j,k)= dts*rhoug*sinthetaf(i,j,k)+f1(i,j,k)
           f2(i,j,k)= dts*rhovg*sinthetaf(i,j,k)+f2(i,j,k)
          enddo
         enddo
        enddo
       endif !ilspgf.ge.1 

      ! update xv (1,2,3)
       do k=1,l
         do j=1,mp
          do i=1,np
           xv(i,j,k,1)=xv(i,j,k,1)+f1(i,j,k)
           xv(i,j,k,2)=xv(i,j,k,2)+f2(i,j,k)
           xv(i,j,k,3)=xv(i,j,k,3)+f3(i,j,k)
          enddo
         enddo
        enddo
      !if (mpi_rank.eq.0) write(6,*) 'end add forces '

      ! Relaxation
       ! HERE FP DIVIDED RELAXXV BY nts BECAUSE RELAXATION IS DONE MORE OFTEN
      do kv=1,5
        kve = kv
        if(kv.eq.5) kve = nv
         do k=1,l
          do j=1,mp
           do i=1,np
            xv(i,j,k,kv)=xv(i,j,k,kv)*(1-relaxxv(i,j,k,min(kv,4))/nts)+
     &                    xe(i,j,k,kve)*relaxxv(i,j,k,min(kv,4))/nts
           enddo
          enddo
         enddo
       enddo !do kv=1,5

       if(ibctopbot.eq.0)then
        write(6,*) 'ibctopbot=0 is not supported'
        stop
!        do j=1,mp
!         do i=1,np
!          if(islip.eq.1) xv(i,j,1,1)=0.
!          if(islip.eq.1) xv(i,j,1,2)=0.
!          xv(i,j,1,3)=-(xv(i,j,1,1)*c13(i,j)*gmul(1)+
!     &                  xv(i,j,1,2)*c23(i,j)*gmul(1))/gi(i,j,1)
!          xv(i,j,l,3)=0.
!         enddo
!        enddo
       endif !if(ibctopbot.eq.0)

      ! computation of uavg, f1avg...
      do k=1,l
        do j=1,mp
         do i=1,np
          ! nb the first term in uavg is added in compress (wt(1))
          uavg(i,j,k)=uavg(i,j,k)+wt(its+1)*xv(i,j,k,1)/xv(i,j,k,5)/gi(i,j,k)
          vavg(i,j,k)=vavg(i,j,k)+wt(its+1)*xv(i,j,k,2)/xv(i,j,k,5)/gi(i,j,k)
          wavg(i,j,k)=wavg(i,j,k)+wt(its+1)*xv(i,j,k,3)/xv(i,j,k,5)/gi(i,j,k)
          f1avg(i,j,k)=f1avg(i,j,k)+wtf(its+1)*f1(i,j,k)
          f2avg(i,j,k)=f2avg(i,j,k)+wtf(its+1)*f2(i,j,k)
          f3avg(i,j,k)=f3avg(i,j,k)+wtf(its+1)*f3(i,j,k)
         enddo
        enddo
       enddo
      !if (mpi_rank.eq.0) write(6,*) 'end inner loop : iteration ',its
      enddo !do its=1,nts

      ! compute oavg
      do k=1,l
       do j=1,mp
        do i=1,np
         oavg(i,j,k)=uavg(i,j,k)*c13(i,j)*gmul(k)+
     &               vavg(i,j,k)*c23(i,j)*gmul(k)+
     &               wavg(i,j,k)*gi(i,j,k)
        enddo
       enddo
      enddo

      deallocate   (u)
      deallocate   (v)
      deallocate   (o)
      return
      end subroutine rmoainner

!*********************************************************************
! compute the donorcell advection scheme on field xf base on adective
! velocities u1,u2,u3  (face contravariant *dt/dx*1/gi)
! this scheme is used in the inner loop of moa and in force.f
! (forcing advection in compress)
!*********************************************************************

      subroutine donorcell(xf,il,iu,jl,ju,lls)
      use metryic
      use gridsetup
      use advo
      use advectionScheme
      use msga

      Implicit None

      integer,intent(in) :: il,iu,jl,ju,lls
      real xf(il:iu,jl:ju,lls)
      integer :: i,j,k

      do k=1,l
       do j=1,mp
        do i=1,np+1
         fd1(i,j,k)=donor(xf(i-1,j,k),xf(i,j,k),u1(i,j,k))
        enddo
       enddo
      enddo

      do k=1,l
       do j=1,mp+1
        do i=1,np
         fd2(i,j,k)=donor(xf(i,j-1,k),xf(i,j,k),u2(i,j,k))
        enddo
       enddo
      enddo

      do k=2,l
       do j=1,mp
        do i=1,np
         fd3(i,j,k)=donor(xf(i,j,k-1),xf(i,j,k),u3(i,j,k))
        enddo
       enddo
      enddo

      if(ibctopbot.eq.0) then
       stop
!       do j=1,mp
!        do i=1,np
!         fd3(i,j, 1)=-fd3(i,j,2)
!         fd3(i,j,l+1)=-fd3(i,j,l)
!        enddo
!       enddo
      else
       do j=1,mp
        do i=1,np
         fd3(i,j, 1)=0.
         fd3(i,j,l+1)=0.
        enddo
       enddo
      endif

! KOO to make CODE2 same as CODE1 by adding these two update (and adding
! halo cell in defineArray), but 
! FP09/2009: this "fix" was incorrect because CODE1 was incorrect 
!       call updated(fd1,fd1,np,mp,l,1-ih,np+ih+1,1-ih,mp+ih,2,0) 
!       call updated(fd2,fd2,np,mp,l,1-ih,np+ih,1-ih,mp+ih+1,3,0) 
 
      do k=1,l
       do j=1,mp
        do i=1,np
         xf(i,j,k)=xf(i,j,k)-( fd1(i+1,j,k)-fd1(i,j,k)
     &                        +fd2(i,j+1,k)-fd2(i,j,k)
     &                        +fd3(i,j,k+1)-fd3(i,j,k) )/h(i,j,k)
        enddo
       enddo
      enddo
      
      return
      end subroutine donorcell


!*****************************************************************************************!
!     gradmoa computes the pressure force on the small time steps (see
!     gc1s, gc2s and gc3s below) and add those terms to f1, f2, f3
!     NB : the lower bc on pr assumes that current is pr(0)=pr(1)
!     and not that pr(1/2)=pr(1)....
!*****************************************************************************************!


      subroutine gradmoa()
      use constants
      use gridsetup
      use pres
      use forcinner
      use metryic
      use msga
      use xvi
      Implicit None

      real :: ghx,ghy,ghz,g13,g23,g33
      integer :: i,j,k
      real*8 :: pr_double 
      real,allocatable:: px(:, :,:),
     .                   py(:, :,:),
     .                   pz(:, :,:)

      ghx=0.5*gc1s
      ghy=0.5*gc2s
      ghz=0.5*gc3s
      allocate (px(1-ih:np+ih,1-ih:mp+ih,l))
      allocate (py(1-ih:np+ih,1-ih:mp+ih,l))
      allocate (pz(1-ih:np+ih,1-ih:mp+ih,l))

      do k=1,l
        do j=1,mp
          do i=1,np
! FP: NB this double precision pressure was added, but it does not seem 
! to affect the results at least on early time steps (09/2019)
!            pr(i,j,k)=(xv(i,j,k,4)*rg/prrcp)**(cp/cv)-pre(i,j,k)
            !pr_double=(xv(i,j,k,4)*rg/prrcp)**(cp/cv)-pre(i,j,k)
            pr_double=(xv(i,j,k,4)
     +    *rg_over_prrcp(i,j,k))**(cp_over_cv(i,j,k))-pre(i,j,k) 
            pr(i,j,k)=real(pr_double)

          enddo
        enddo
      enddo
!      call rmaxmin1(pr,'pr',1-ih,np+ih,1-ih,mp+ih,l)
!      if(mpi_rank.eq.0) then
!       write (6,*) 'cp/cv,rg/cp,cp0,rg_over_prrcp0'
!     + ,cp_over_cv_gas,rg_over_cp_gas,cp_gas,rg_over_prrcp_gas
!      endif
      call updated(pr,pr,np,mp,l,1-ih,np+ih,1-ih,mp+ih,0,0)
      
      do i=1+leftdedge,np-rightdedge                       
      do k=1,l
      do j=1,mp
        px(i,j,k)=    ghx*(pr(i+1,j,k)-pr(i-1,j,k))

      enddo
      enddo
      enddo
      if(leftdedge.eq.1) then
        do k=1,l
        do j=1,mp
         px(1,j,k)=gc1s*(pr(2,j,k)-pr(1,j,k))
        enddo
        enddo
      endif
      if(rightdedge.eq.1) then                     
        do k=1,l
        do j=1,mp
          px(np,j,k)=gc1s*(pr(np,j,k)-pr(np-1,j,k))
        enddo
        enddo
      endif

      do k=1,l
      do i=1,np
      do j=1+botdedge,mp-topdedge                  
        py(i,j,k)=     ghy*(pr(i,j+1,k)-pr(i,j-1,k))
      enddo
      enddo
      enddo
      if(botdedge.eq.1) then                     
        do k=1,l
        do i=1,np
          py(i,1,k)=gc2s*(pr(i,2,k)-pr(i,1,k))
        enddo
        enddo
      endif
      if(topdedge.eq.1) then
        do k=1,l
        do i=1,np
          py(i,mp,k)=gc2s*(pr(i,mp,k)-pr(i,mp-1,k))
        enddo
        enddo
      endif
      do k=2,l-1
      do j=1,mp
      do i=1,np
        pz(i,j,k)=ghz*(pr(i,j,k+1)-pr(i,j,k-1))
      enddo
      enddo
      enddo
      do j=1,mp
      do i=1,np
!      pz(i,j,1)= gcz*(pr(i,j,2)-pr(i,j,1))     
        pz(i,j,1)= ghz*(pr(i,j,2)-pr(i,j,1))     !modified by rrl 10/02/01
        pz(i,j,l)= ghz*(pr(i,j,l)-pr(i,j,l-1))
      enddo
      enddo
! top and bottom b.c.
      if(islip.eq.1) then
         stop
         ! should not be here!
!         xv(:,:,1,1)=0.
!         xv(:,:,1,2)=0.
!         xv(:,:,1,3)=0.
      endif
      if(ibctopbot.eq.0) then
      stop
!      do k=1,l,l-1
!      do j=1,mp
!      do i=1,np
!        g13=c13(i,j)*gmul(k)
!        g23=c23(i,j)*gmul(k)
!        g33= gi(i,j,k)
!        pz(i,j,k)=(g13*xv(i,j,k,1)+g23*xv(i,j,k,2)+g33*xv(i,j,k,3)-
!     .             g13*px(i,j,k)-g23*py(i,j,k))
!     .           /(g33*g33+ g13*g13+g23*g23)
!      enddo
!      enddo
!      enddo
      endif

      do k=1,l
      do j=1,mp
      do i=1,np
        g13=c13(i,j)*gmul(k)
        g23=c23(i,j)*gmul(k)
        g33= gi(i,j,k)

        f1(i,j,k)=-(    px(i,j,k)               +g13*pz(i,j,k))
        f2(i,j,k)=-(                   py(i,j,k)+g23*pz(i,j,k))
        f3(i,j,k)=-(g33*pz(i,j,k))
      enddo
      enddo
      enddo

      deallocate (px)
      deallocate (py)
      deallocate (pz)
      return
      end subroutine gradmoa

!*****************************************************************************************!




!*****************************************************************************************!
      end module higrad
!*****************************************************************************************!
   
