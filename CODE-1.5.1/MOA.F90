!-----------------------------------------------------------------------
! The MOA inner loop on small time steps
!    The method of average does the time integration on the xv fields
!    using a donor cell scheme for advection, and computing pressure, 
!    gravity, coriolis, lspgf, and relaxation to ambiant variable
!    the sum of the forcing (except the advection) is added in 
!    force(:,:,:,1), force(:,:,:,2),force(:,:,:,3) (using weights wtf)
!    the mean advective velocities uavg, vavg and oavg are averaged 
!    centered contravariant velocities / gi for advection scheme on 
!    large time step
!-----------------------------------------------------------------------
subroutine MOA
  use gridlist_variables, only : nts,l,ih,wp
  use vegetation_gridlist_variables, only : nfuel
  use forcings, only : forceSE_xv,forceSI_xv,forceLE_xv, &
    forceSE_xvfuel,forceSI_xvfuel,forceLE_xvfuel,xvLim,xvfuelLim, &
    uavg,vavg,wavg,oavg,wt,wtf
  use xvall, only : xv,xvtmp,xvfuel,xvfueltmp,nv,nvfuel,iuvel,ivvel, &
    iwvel,irho
  use fuel_variables, only : lfuel
  use metric_variables, only : gi,gmul,c13,c23
  use gridsetup, only : np,mp,dts
  Implicit None

  ! Local Variables
  integer :: its,i,j,k,kv,iuv,ift
  real(wp),allocatable :: u(:,:,:,:),v(:,:,:,:),o(:,:,:,:)  ! centered contravariant coordinates/gi 

  ! Executable Code
  allocate(u(1-ih:np+ih, 1-ih:mp+ih,l,0:2))
  allocate(v(1-ih:np+ih, 1-ih:mp+ih,l,0:2))
  allocate(o(1-ih:np+ih, 1-ih:mp+ih,l,0:2))

  ! initialization of u,v,w (for last index 1), force(:,:,:,1),uavg...
  xvtmp=xv
  xvfueltmp=xvfuel
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
    enddo
  enddo

  ! beginning of the inner loop
  do its=1,nts
    do k=1,l
      do j=1,mp
        do i=1,np
         ! past u,v,o at t-dts: NB for its=1, it is for t...
          u(i,j,k,0)=u(i,j,k,1)
          v(i,j,k,0)=v(i,j,k,1)
          o(i,j,k,0)=o(i,j,k,1)
          ! current u,v, o
          u(i,j,k,1)=xvtmp(i,j,k,iuvel)/xvtmp(i,j,k,irho)/gi(i,j,k)
          v(i,j,k,1)=xvtmp(i,j,k,ivvel)/xvtmp(i,j,k,irho)/gi(i,j,k)
          o(i,j,k,1)=(xvtmp(i,j,k,iwvel)*gi(i,j,k) &
            +c13(i,j)*gmul(k)*xvtmp(i,j,k,iuvel) &
            +c23(i,j)*gmul(k)*xvtmp(i,j,k,ivvel))/xvtmp(i,j,k,irho) &
            /gi(i,j,k)
          ! time integrated u,v,o at t+0.5*dts
          u(i,j,k,2)=0.5_wp*(3._wp*u(i,j,k,1)-u(i,j,k,0))
          v(i,j,k,2)=0.5_wp*(3._wp*v(i,j,k,1)-v(i,j,k,0))
          o(i,j,k,2)=0.5_wp*(3._wp*o(i,j,k,1)-o(i,j,k,0))
        enddo
      enddo
    enddo
  
    ! computation of advective velocities at t+0.5dts: u1,u2,u3
    call setAdvectiveVelocities(u(:,:,:,2),v(:,:,:,2),o(:,:,:,2) &
                ,1-ih,np+ih,1-ih,mp+ih,l,dts)

    ! computation of advective tendencies (donorcell scheme)
    do kv=1,nv
      call donorcell(xvtmp(:,:,:,kv),1-ih,np+ih,1-ih,mp+ih,l)
      call loop_update_flags(kv,iuv)
      call update(xvtmp(:,:,:,kv),np,mp,l,1-ih,np+ih,1-ih,mp+ih,0,0,iuv)
    enddo 

    call implicit_domain_updates
    call small_implicit_forcings

    ! Apply small step forcings
    do kv=1,nv
      do k=1,l
        do j=1,mp
          do i=1,np
            xvtmp(i,j,k,kv)=xvtmp(i,j,k,kv)+forceSE_xv(i,j,k,kv)*dts
            xvtmp(i,j,k,kv)=xvtmp(i,j,k,kv)+forceSI_xv(i,j,k,kv)*dts
          enddo
        enddo
      enddo
      call loop_update_flags(kv,iuv)
      call update(xvtmp(:,:,:,kv),np,mp,l,1-ih,np+ih,1-ih,mp+ih,1,0,iuv) 
    enddo
    call forcingLimits(xvtmp,xvLim,1-ih,np+ih,1-ih,mp+ih,l,nv)
    xvfueltmp=xvfueltmp+forceSE_xvfuel*dts
    xvfueltmp=xvfueltmp+forceSI_xvfuel*dts
    do ift=1,nfuel
      call forcingLimits(xvfueltmp(ift,:,:,:,:),xvfuelLim, &
        1,np,1,mp,lfuel,nvfuel)
    enddo

    ! Incorporate into large step forcings
    forceLE_xv    =forceLE_xv    +wtf(its)*forceSI_xv
    forceLE_xvfuel=forceLE_xvfuel+wtf(its)*forceSI_xvfuel

    ! computation of uavg, force(:,:,:,1)...
    do k=1,l
      do j=1,mp
        do i=1,np
          ! nb the first term in uavg is added in compress (wt(1))
          uavg(i,j,k)=uavg(i,j,k)+wt(its+1)*xvtmp(i,j,k,iuvel)/ &
            xvtmp(i,j,k,irho)/gi(i,j,k)
          vavg(i,j,k)=vavg(i,j,k)+wt(its+1)*xvtmp(i,j,k,ivvel)/ &
            xvtmp(i,j,k,irho)/gi(i,j,k)
          wavg(i,j,k)=wavg(i,j,k)+wt(its+1)*xvtmp(i,j,k,iwvel)/ &
            xvtmp(i,j,k,irho)/gi(i,j,k)
        enddo
      enddo
    enddo
  enddo !do its=1,nts

  ! compute oavg
  do k=1,l
    do j=1,mp
      do i=1,np
        oavg(i,j,k)=uavg(i,j,k)*c13(i,j)*gmul(k)+ &
                vavg(i,j,k)*c23(i,j)*gmul(k)+ &
                wavg(i,j,k)*gi(i,j,k)
      enddo
    enddo
  enddo
  deallocate(u,v,o)

end subroutine MOA
