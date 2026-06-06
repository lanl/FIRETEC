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
!    . the sum of the forcing (except the advection) is added in force(:,:,:,1), force(:,:,:,2),
!    force(:,:,:,3) (using weights wtf)
!    . the mean advective velocities uavg, vavg and oavg are averaged centered
!    contravariant velocities / gi for advection scheme on large time
!    step
!***********************************************************************
  subroutine rmoainner
    use gridlist_variables, only : nts,dts,l,ih,irhovapor,icorio,ilspgf, &
      slopeangle,slopeazimuth,ifire
    use forcings, only : force,uavg,vavg,oavg,wt,wtf
    use xvall, only : xv,xe,xvtmp,relaxxv,iuvel,ivvel,iwvel,itemp,irho,nv, &
      ivapor,iO2
    use thermo_variables, only : cp_over_cv_gas,rg_over_prrcp_gas,pr
    use metric_variables, only : gi,gmul,c13,c23
    use lspgf_variables, only : sinthetaf
    use gridsetup, only : np,mp,time
    use constants, only : g,pi,fcor2,fcor3
    use forcinner
    use workavg
    use metric_variables_old
    Implicit None

    ! Local Variables
    integer :: its,i,j,k,kv,kve, kref 
    real :: uc,vc ! reference wind for coriolis
    real :: piOver180, g1,g2, g3 ! projected gravity coeff
    real,allocatable :: u(:,:,:,:),&  ! centered contravariant coordinates/gi 
                        v(:,:,:,:),&
                        o(:,:,:,:)

    ! Executable Code
    piOver180 = pi/180.0
    g1 = g*sin(slopeangle*piOver180)*cos(slopeazimuth*piOver180)
    g2 = g*sin(slopeangle*piOver180)*sin(slopeazimuth*piOver180)
    g3 = g*cos(slopeangle*piOver180)
     
    allocate(u(1-ih:np+ih, 1-ih:mp+ih,l,0:2))
    allocate(v(1-ih:np+ih, 1-ih:mp+ih,l,0:2))
    allocate(o(1-ih:np+ih, 1-ih:mp+ih,l,0:2))
    !FP added this section for computing of gas thermal properties,
    !which are assumed to be constant within the higrad loop
    do k=1,l
      do j=1,mp
        do i=1,np
          call updateGasThermProps(ifire*xv(i,j,k,iO2)/xv(i,j,k,irho), &
            irhovapor*xv(i,j,k,ivapor)/xv(i,j,k,irho))
          rg_over_prrcp(i,j,k)=rg_over_prrcp_gas
          cp_over_cv(i,j,k)=cp_over_cv_gas
        enddo
      enddo
    enddo
    ! initialization of u,v,w (for last index 1), force(:,:,:,1),uavg...
    force(:,:,:,1)=0.0
    force(:,:,:,2)=0.0
    force(:,:,:,3)=0.0
    do j=1,mp
      do i=1,np
        u(i,j,:,1)=xv(i,j,:,iuvel)/xv(i,j,:,irho)/gi(i,j,:)
        v(i,j,:,1)=xv(i,j,:,ivvel)/xv(i,j,:,irho)/gi(i,j,:)
        o(i,j,:,1)=(xv(i,j,:,iwvel)*gi(i,j,:) &
          +c13(i,j)*gmul(:)*xv(i,j,:,iuvel) &
          +c23(i,j)*gmul(:)*xv(i,j,:,ivvel))/xv(i,j,:,irho)/gi(i,j,:)
        uavg(i,j,:)=wt(1)*u(i,j,:,1)
        vavg(i,j,:)=wt(1)*v(i,j,:,1)
        wavg(i,j,:)=wt(1)*xv(i,j,:,iwvel)/xv(i,j,:,irho)/gi(i,j,:)
        xvtmp(i,j,:,1)=xv(i,j,:,iuvel)
        xvtmp(i,j,:,2)=xv(i,j,:,ivvel)
        xvtmp(i,j,:,3)=xv(i,j,:,iwvel)
        xvtmp(i,j,:,4)=xv(i,j,:,itemp)
        xvtmp(i,j,:,5)=xv(i,j,:,irho)
      enddo
    enddo

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
            u(i,j,k,1)=xvtmp(i,j,k,1)/xvtmp(i,j,k,5)/gi(i,j,k)
            v(i,j,k,1)=xvtmp(i,j,k,2)/xvtmp(i,j,k,5)/gi(i,j,k)
            o(i,j,k,1)=(xvtmp(i,j,k,3)*gi(i,j,k)&
              +c13(i,j)*gmul(k)*xvtmp(i,j,k,1)&
              +c23(i,j)*gmul(k)*xvtmp(i,j,k,2))/xvtmp(i,j,k,5)/gi(i,j,k)
            ! time integrated u,v,o at t+0.5*dts
            u(i,j,k,2)=0.5*(3.*u(i,j,k,1)-u(i,j,k,0))
            v(i,j,k,2)=0.5*(3.*v(i,j,k,1)-v(i,j,k,0))
            o(i,j,k,2)=0.5*(3.*o(i,j,k,1)-o(i,j,k,0))
          enddo
        enddo
      enddo
    
      ! computation of advective velocities at t+0.5dts: u1,u2,u3
      call setAdvectiveVelocities(u(:,:,:,2),v(:,:,:,2),o(:,:,:,2) &
                  ,1-ih,np+ih,1-ih,mp+ih,l,dts,0)

      ! computation of advective tendencies (donorcell scheme)

      do kv=1,5
        kve=kv
        if (kv==5) kve=nv
        call update(xvtmp(1-ih,1-ih,1,kv),xe(1-ih,1-ih,1,kve),np,mp,l,1-ih,np+ih,1-ih,mp+ih,0,0)
        call donorcell(xvtmp(1-ih,1-ih,1,kv),1-ih,np+ih,1-ih,mp+ih,l)
      enddo 

      ! compute pressure force and put it in f1,f2,f3
      call gradmoa()
      if(its.eq.nts) call rmaxmin(pr,'pr',1-ih,np+ih,1-ih,mp+ih,l,1)

      ! forces : gravity
      do k=1,l
        do j=1,mp
          do i=1,np
            f1(i,j,k)=f1(i,j,k)-dts*(xvtmp(i,j,k,5)-xe(i,j,k,nv))*g1
            f2(i,j,k)=f2(i,j,k)-dts*(xvtmp(i,j,k,5)-xe(i,j,k,nv))*g2
            f3(i,j,k)=f3(i,j,k)-dts*(xvtmp(i,j,k,5)-xe(i,j,k,nv))*g3
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
              f1(i,j,k)= dts*(fcor3*(xvtmp(i,j,k,2)-vc)  &
                         -fcor2*(xvtmp(i,j,k,3)) )+f1(i,j,k)
              f2(i,j,k)=-dts*(fcor3*(xvtmp(i,j,k,1)-uc) )+f2(i,j,k)
              f3(i,j,k)= dts*(fcor2*(xvtmp(i,j,k,1)-uc))+f3(i,j,k)
            enddo
          enddo
        enddo
      endif !icorio.ge.1

      ! forces : lspgf
      if (ilspgf.ge.1) then !ilspgf
        do k=1,l
          do j=1,mp
            do i=1,np
              f1(i,j,k)= dts*xe(1,1,1,iuvel)*sinthetaf(i,j,k)+f1(i,j,k)
              f2(i,j,k)= dts*xe(1,1,1,ivvel)*sinthetaf(i,j,k)+f2(i,j,k)
            enddo
          enddo
        enddo
      endif !ilspgf.ge.1 

      ! update xv (1,2,3)
      do k=1,l
        do j=1,mp
          do i=1,np
            xvtmp(i,j,k,1)=xvtmp(i,j,k,1)+f1(i,j,k)
            xvtmp(i,j,k,2)=xvtmp(i,j,k,2)+f2(i,j,k)
            xvtmp(i,j,k,3)=xvtmp(i,j,k,3)+f3(i,j,k)
          enddo
        enddo
      enddo

      ! Relaxation
       ! HERE FP DIVIDED RELAXXV BY nts BECAUSE RELAXATION IS DONE MORE OFTEN
      do kv=1,5
        kve = kv
        if(kv.eq.5) kve = nv
        do k=1,l
          do j=1,mp
            do i=1,np
              xvtmp(i,j,k,kv)=xvtmp(i,j,k,kv)*(1-relaxxv(i,j,k,min(kv,4))/nts)+&
                         xe(i,j,k,kve)*relaxxv(i,j,k,min(kv,4))/nts
            enddo
          enddo
        enddo
      enddo !do kv=1,5

      ! computation of uavg, force(:,:,:,1)...
      do k=1,l
        do j=1,mp
          do i=1,np
            ! nb the first term in uavg is added in compress (wt(1))
            uavg(i,j,k)=uavg(i,j,k)+wt(its+1)*xvtmp(i,j,k,1)/xvtmp(i,j,k,5)/gi(i,j,k)
            vavg(i,j,k)=vavg(i,j,k)+wt(its+1)*xvtmp(i,j,k,2)/xvtmp(i,j,k,5)/gi(i,j,k)
            wavg(i,j,k)=wavg(i,j,k)+wt(its+1)*xvtmp(i,j,k,3)/xvtmp(i,j,k,5)/gi(i,j,k)
            force(i,j,k,1)=force(i,j,k,1)+wtf(its+1)*f1(i,j,k)
            force(i,j,k,2)=force(i,j,k,2)+wtf(its+1)*f2(i,j,k)
            force(i,j,k,3)=force(i,j,k,3)+wtf(its+1)*f3(i,j,k)
          enddo
        enddo
      enddo
    enddo !do its=1,nts

    ! compute oavg
    do k=1,l
      do j=1,mp
        do i=1,np
          oavg(i,j,k)=uavg(i,j,k)*c13(i,j)*gmul(k)+&
                  vavg(i,j,k)*c23(i,j)*gmul(k)+&
                  wavg(i,j,k)*gi(i,j,k)
        enddo
      enddo
    enddo
    deallocate(u,v,o)
  
  end subroutine rmoainner

  !----------------------------------------------------------------
  ! compute the donorcell advection scheme on field xf base on
  ! advective velocities u1,u2,u3 (face contravariant *dt/dx*1/gi)
  ! this scheme is used in the inner loop of moa and in force.f
  ! (forcing advection in compress)
  !----------------------------------------------------------------
  subroutine donorcell(xf,il,iu,jl,ju,lls)
    use gridlist_variables, only : l
    use forcings, only : u1,u2,u3
    use metric_variables, only : gi
    use gridsetup, only : np,mp
    Implicit None
  
    ! Local Variables
    integer,intent(in) :: il,iu,jl,ju,lls
    real,intent(inout) :: xf(il:iu,jl:ju,lls)
  
    integer :: i,j,k
    real,external :: donor
    real,allocatable :: fd(:,:,:,:)
  
    ! Executable Code
    allocate(fd(3,np+1,mp+1,l+1)); fd=0.
    do k=1,l
      do j=1,mp
        do i=1,np+1
          fd(1,i,j,k)=donor(xf(i-1,j,k),xf(i,j,k),u1(i,j,k))
        enddo
      enddo
    enddo
    do k=1,l
      do j=1,mp+1
        do i=1,np
          fd(2,i,j,k)=donor(xf(i,j-1,k),xf(i,j,k),u2(i,j,k))
        enddo
      enddo
    enddo
    do k=2,l
      do j=1,mp
        do i=1,np
          fd(3,i,j,k)=donor(xf(i,j,k-1),xf(i,j,k),u3(i,j,k))
        enddo
      enddo
    enddo
  
    do k=1,l
      do j=1,mp
        do i=1,np
          xf(i,j,k)=xf(i,j,k)-gi(i,j,k)*(fd(1,i+1,j,k)-fd(1,i,j,k) &
            +fd(2,i,j+1,k)-fd(2,i,j,k)+fd(3,i,j,k+1)-fd(3,i,j,k))
        enddo
      enddo
    enddo
    deallocate(fd)
      
  end subroutine donorcell

  !----------------------------------------------------------------
  ! gradmoa computes the pressure force on the small time steps and 
  ! add those terms to the overall forcing terms
  ! NB : the lower bc on pr assumes that current is pr(0)=pr(1)
  ! and not that pr(1/2)=pr(1)....
  !----------------------------------------------------------------
  subroutine gradmoa()
    use gridlist_variables, only : l,ih,dts
    use xvall, only : xvtmp,iO2,ivapor,iuvel,ivvel,iwvel
    use thermo_variables, only : rg_over_prrcp_gas,cp_over_cv_gas, &
      pr,pre
    use metric_variables, only : gmul,gi,c13,c23
    use gridsetup, only : np,mp,dxi,dyi,dzi,dti
    use msga_variables, only : leftdedge,rightdedge,botdedge, &
      topdedge
    use forcinner
    Implicit None

    ! Local Variables
    real :: ghx,ghy,ghz,g13,g23,g33
    integer :: i,j,k
    real*8 :: pr_double 
    real,allocatable:: px(:, :,:),&
                      py(:, :,:),&
                      pz(:, :,:)

    ! Executable Code
    allocate (px(1-ih:np+ih,1-ih:mp+ih,l))
    allocate (py(1-ih:np+ih,1-ih:mp+ih,l))
    allocate (pz(1-ih:np+ih,1-ih:mp+ih,l))

    do k=1,l
      do j=1,mp
        do i=1,np
          pr_double=(xvtmp(i,j,k,4)&
           *rg_over_prrcp(i,j,k))**(cp_over_cv(i,j,k))-pre(i,j,k) 
          pr(i,j,k)=real(pr_double)
        enddo
      enddo
    enddo
    call update(pr,pr,np,mp,l,1-ih,np+ih,1-ih,mp+ih,0,0)
    
    do i=1+leftdedge,np-rightdedge                       
      do k=1,l
        do j=1,mp
          px(i,j,k)=0.5*(pr(i+1,j,k)-pr(i-1,j,k))*dts*dxi
        enddo
      enddo
    enddo
    if(leftdedge.eq.1) then
      do k=1,l
        do j=1,mp
          px(1,j,k)=(pr(2,j,k)-pr(1,j,k))*dts*dxi
        enddo
      enddo
    endif
    if(rightdedge.eq.1) then                     
      do k=1,l
        do j=1,mp
          px(np,j,k)=(pr(np,j,k)-pr(np-1,j,k))*dts*dxi
        enddo
      enddo
    endif

    do k=1,l
      do i=1,np
        do j=1+botdedge,mp-topdedge                  
          py(i,j,k)=0.5*(pr(i,j+1,k)-pr(i,j-1,k))*dts*dyi
        enddo
      enddo
    enddo
    if(botdedge.eq.1) then                     
      do k=1,l
        do i=1,np
          py(i,1,k)=(pr(i,2,k)-pr(i,1,k))*dts*dyi
        enddo
      enddo
    endif
    if(topdedge.eq.1) then
      do k=1,l
        do i=1,np
          py(i,mp,k)=(pr(i,mp,k)-pr(i,mp-1,k))*dts*dyi
        enddo
      enddo
    endif

    do k=2,l-1
      do j=1,mp
        do i=1,np
          pz(i,j,k)=0.5*(pr(i,j,k+1)-pr(i,j,k-1))*dts*dzi
        enddo
      enddo
    enddo
    
    do j=1,mp
      do i=1,np
        pz(i,j,1)= 0.5*(pr(i,j,2)-pr(i,j,1))*dts*dzi
        pz(i,j,l)= 0.5*(pr(i,j,l)-pr(i,j,l-1))*dts*dzi
      enddo
    enddo
    
    do k=1,l
      do j=1,mp
        do i=1,np
          g13=c13(i,j)*gmul(k)
          g23=c23(i,j)*gmul(k)
          g33= gi(i,j,k)

          f1(i,j,k)=-(px(i,j,k)+g13*pz(i,j,k))
          f2(i,j,k)=-(py(i,j,k)+g23*pz(i,j,k))
          f3(i,j,k)=-(g33*pz(i,j,k))
        enddo
      enddo
    enddo

    deallocate (px,py,pz)
  end subroutine gradmoa

!*****************************************************************************************!
end module higrad
!*****************************************************************************************!
   
