!-----------------------------------------------------------------------
! forcing applies the large step explicit forcings 
!-----------------------------------------------------------------------
subroutine forcing
  use gridlist_variables, only : l,ih,dt
  use vegetation_gridlist_variables, only : nfuel
  use forcings, only : forceLE_xv,forceLE_xvfuel,xvLim,xvfuelLim
  use fuel_variables, only : lfuel
  use gridsetup, only : np,mp
  use forcings, only : uavg,vavg,oavg
  use xvall, only : xv,xvrho,nv,xvfuel,nvfuel
  use gridsetup, only : np,mp
  use metric_variables, only : h
  Implicit None

  ! Local Variables 
  integer :: i,j,k,kv,iuv,ift
 
  ! Executable Code
  call setAdvectiveVelocities(uavg,vavg,oavg,1-ih,np+ih,1-ih,mp+ih,l, &
    dt) ! this is advvel in code. 
  call mpdata(xv,xvrho,h,1-ih,np+ih,1-ih,mp+ih,l,nv) ! LES forcings
  do kv=1,nv
    !!!!! Transported Forcings !!!!!!
    call loop_update_flags(kv,iuv)
    call update(forceLE_xv(:,:,:,kv),np,mp,l, &
      1-ih,np+ih,1-ih,mp+ih,0,0,iuv)
    call donorcell(forceLE_xv(:,:,:,kv),1-ih,np+ih,1-ih,mp+ih,l)
    do k=1,l
      do j=1,mp
        do i=1,np
          xv(i,j,k,kv)=xv(i,j,k,kv)+forceLE_xv(i,j,k,kv)*dt
        enddo
      enddo
    enddo
    call forcingLimits(xv(:,:,:,kv),xvLim(:,kv), &
      1-ih,np+ih,1-ih,mp+ih,l,1)
    call update(xv(:,:,:,kv),np,mp,l,1-ih,np+ih,1-ih,mp+ih,1,0,iuv) 
  enddo

  !!!!! Fuel Forcings !!!!!
  xvfuel=xvfuel+forceLE_xvfuel*dt
  do ift=1,nfuel
    call forcingLimits(xvfuel(ift,:,:,:,:),xvfuelLim, &
      1,np,1,mp,lfuel,nvfuel)
  enddo

end subroutine forcing

!-----------------------------------------------------------------------
! Small clippings to xv arrays in case beyond the physically possible
!-----------------------------------------------------------------------
subroutine forcingLimits(xv,xvLimits,il,iu,jl,ju,lls,nvp)
  use gridlist_variables, only : wp
  Implicit None

  ! Local Variables
  integer,intent(in) :: il,iu,jl,ju,lls,nvp
  real(wp),intent(inout) :: xv(il:iu,jl:ju,lls,nvp)
  real(wp),intent(in) :: xvLimits(2,nvp)
  integer :: i,j,k,kv

  ! Executable Code
  do kv=1,nvp
    do k=1,lls
      do j=jl,ju
        do i=il,iu
          xv(i,j,k,kv)=max(xvLimits(1,kv),min(xvLimits(2,kv), &
            xv(i,j,k,kv)))
        enddo
      enddo
    enddo
  enddo

end subroutine forcingLimits

!-----------------------------------------------------------------------
! small_explicit_forcings calculates all the forcing terms that will be 
! applied in the implicit stepping method but will only be calculated 
! once on the large time step
!-----------------------------------------------------------------------
subroutine small_explicit_forcings
  use gridlist_variables, only : iturb,ih,l
  use vegetation_gridlist_variables, only : nfuel,idiffsies
  use radiation_gridlist_variables, only : irad
  use thermo_variables, only : tempg
  use gridsetup, only : np,mp
  use forcings, only : forceSE_xv,forceSE_xvfuel
  use fuel_variables, only : lfuel,temps
  use xvall, only : itheta,nv,xv,isies,irhof,irhow,xvfuel
  use msga_variables, only : mpi_rank,ierror
  Implicit None

  ! Local Variables

  ! Executable Code
  !!!!! Transported Forcings !!!!!!
  ! compute macro drag contribution
  call dragm(forceSE_xv,xv,1-ih,np+ih,1-ih,mp+ih,l,nv)
     
  ! compute forcing terms associated with subgrid turbulence
  if(iturb.eq.2) call calcForceTurb(xv,forceSE_xv, &
    1-ih,np+ih,1-ih,mp+ih,l,nv)

  !!!!! Fuel Forcings !!!!!
  ! Radiation
  if(irad.eq.1)then ! Diffusive Radiation Scheme
    if(mpi_rank.eq.0) print*,'Diffusive Radiation Not implemented'
    call mpi_finalize(ierror)
    STOP
  elseif(irad.eq.2)then ! Monte-Carlo Radiation Scheme
    call firerad_MC(forceSE_xv(:,:,:,itheta), &
      forceSE_xvfuel(:,:,:,:,isies), &
      xvfuel(:,:,:,:,irhof)+xvfuel(:,:,:,:,irhow), &
      1-ih,np+ih,1-ih,mp+ih,1,np,1,mp,l,nfuel)
  elseif(irad.eq.3)then ! Radiation Sink
    call radiationSink(forceSE_xv(:,:,:,itheta), &
      forceSE_xvfuel(:,:,:,:,isies),xvfuel(:,:,:,:,irhof), &
      xvfuel(:,:,:,:,irhow),tempg,temps, &
      1-ih,np+ih,1-ih,mp+ih,1,np,1,mp,l,lfuel,nfuel)
  endif
      
  if(idiffsies.eq.1) call diffusesies(xvfuel(:,:,:,:,irhof), &
    xvfuel(:,:,:,:,irhow),forceSE_xvfuel(:,:,:,:,irhof), &
    forceSE_xvfuel(:,:,:,:,isies),1,np,1,mp,lfuel,nfuel,ih)

end subroutine small_explicit_forcings

!-----------------------------------------------------------------------
! large_explicit_forcings calculates all the forcing terms that will be
! applied only once explicitly on the large time steps
!-----------------------------------------------------------------------
subroutine large_explicit_forcings
  use gridlist_variables, only : iheatsource,ifbrand,ih,l
  use heatsource_gridlist_variables, only : starttime_hs
  use vegetation_gridlist_variables, only : nfuel,ifire
  use radiation_gridlist_variables, only : isolar
  use forcings, only : forceSE_xv,forceLE_xv,forceSE_xvfuel, &
    forceLE_xvfuel
  use gridsetup, only : np,mp,time
  use xvall, only : xv,xvrho,xvfuel,nv,iuvel,ivvel,iwvel,irho, &
    nvfuel,irhof,isies,iH2O,irhow,irhof
  use fuel_variables, only : temps,lfuel
  use thermo_variables, only : tempg
  Implicit None

  ! Local Variables

  ! Executable Code
  !!!!! Transported Forcings !!!!!!
  forceLE_xv=forceLE_xv+forceSE_xv

  !!!!! Fuel Forcings !!!!!!
  forceLE_xvfuel=forceLE_xvfuel+forceSE_xvfuel

  ! Ignite new cells
  if(ifire.eq.1) call ignite(forceLE_xvfuel(:,:,:,:,isies),xvfuel, &
    1,np,1,mp,lfuel,nfuel,nvfuel)
  if(iheatsource.ge.1) then
    if(any(time.ge.starttime_hs))then
      call heatSource(iheatsource, &
        forceLE_xv,1,np,1,mp,l,nv)
    endif
  endif

  !!!!! Firebrand Forcings !!!!!
  if(ifbrand.eq.1) call firebrand(xvrho(:,:,:,iuvel), &
    xvrho(:,:,:,ivvel),xvrho(:,:,:,iwvel),1-ih,np+ih,1-ih,mp+ih,0,l, &
    xv(:,:,:,irho),tempg(:,:,:),xvfuel(:,:,:,:,irhof), &
    temps(:,:,:,:),1,np,1,mp,lfuel,nfuel)
  
  !!!!! Solar heating and evaporation Forcings !!!!!
  if(isies.ne.1)then
    call evapBulk(forceLE_xv,forceLE_xvfuel,xv,xvrho,xvfuel, &
      1-ih,np+ih,1-ih,mp+ih,0,l,1,np,1,mp,lfuel,nv,nvfuel,nfuel)
    if(isolar.eq.1)then
      call solarHeating(forceLE_xvfuel(:,:,:,:,isies), &
        xvfuel(:,:,:,:,irhof),xvfuel(:,:,:,:,irhow), &
        xvrho(:,:,:,iH2O),1-ih,np+ih,1-ih,mp+ih,l, &
        1,np,1,mp,lfuel,nfuel)  
    endif
  endif

end subroutine large_explicit_forcings

!-----------------------------------------------------------------------
! small_implicit_forcings calculates all the forcing terms that will be
! applied implicitly on the small time steps
!-----------------------------------------------------------------------
subroutine small_implicit_forcings
  use gridlist_variables, only : ih,l,icorio,ilspgf,iemissions, &
    slopeangle,slopeazimuth,wp
  use vegetation_gridlist_variables, only : nfuel,ifire
  use forcings, only : forceSI_xv,forceSI_xvfuel
  use xvall, only : xe,xvtmp,xvfueltmp,xvrhotmp,nv,nvfuel,relaxxv, &
    itheta,isies,iuvel,ivvel,iwvel,irho,irhof,iEmitStart,iEmitStop
  use lspgf_variables, only : sinthetaf
  use gridsetup, only : np,mp,dts
  use constants, only : g,pi,fcor2,fcor3
  use emission_gridlist_variables, only : nEmit,nAero,nMAero
  use fuel_variables, only : temps,lfuel
  use gridlist_variables, only : wp
  use thermo_variables, only : tempg
  Implicit None

  ! Local Variables
  integer :: i,j,k,kref 
  real(wp) :: uc,vc ! reference wind for coriolis
  real(wp) :: piOver180,g1,g2,g3 ! projected gravity coeff

  ! Executable Code
  
  !!!!! Transported Forcings !!!!!!
  ! forces : pressure
  call gradmoa()

  ! forces : gravity
  piOver180 = pi/180.0_wp
  g1 = g*sin(slopeangle*piOver180)*cos(slopeazimuth*piOver180)
  g2 = g*sin(slopeangle*piOver180)*sin(slopeazimuth*piOver180)
  g3 = g*cos(slopeangle*piOver180)
  do k=1,l
    do j=1,mp
      do i=1,np
        forceSI_xv(i,j,k,iuvel)=forceSI_xv(i,j,k,iuvel) &
          -(xvtmp(i,j,k,irho)-xe(i,j,k,irho))*g1
        forceSI_xv(i,j,k,ivvel)=forceSI_xv(i,j,k,ivvel) &
          -(xvtmp(i,j,k,irho)-xe(i,j,k,irho))*g2
        forceSI_xv(i,j,k,iwvel)=forceSI_xv(i,j,k,iwvel) &
          -(xvtmp(i,j,k,irho)-xe(i,j,k,irho))*g3
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
          forceSI_xv(i,j,k,iuvel)=forceSI_xv(i,j,k,iuvel) &
            +fcor3*(xvtmp(i,j,k,ivvel)-vc)  &
            -fcor2*(xvtmp(i,j,k,iwvel))
          forceSI_xv(i,j,k,ivvel)=forceSI_xv(i,j,k,ivvel) &
            -fcor3*(xvtmp(i,j,k,iuvel)-uc)
          forceSI_xv(i,j,k,iwvel)=forceSI_xv(i,j,k,iwvel) &
            +fcor2*(xvtmp(i,j,k,iuvel)-uc)
        enddo
      enddo
    enddo
  endif !icorio.ge.1

  ! forces : lspgf
  if (ilspgf.ge.1) then !ilspgf
    forceSI_xv(:,:,:,iuvel)=forceSI_xv(:,:,:,iuvel) &
      +xe(1,1,1,iuvel)*sinthetaf
    forceSI_xv(:,:,:,ivvel)=forceSI_xv(:,:,:,ivvel) &
      +xe(1,1,1,ivvel)*sinthetaf
  endif !ilspgf.ge.1 
    
  ! Relaxation
  call boundary(forceSI_xv,xvtmp,xe,relaxxv,dts,1-ih,np+ih,1-ih,mp+ih, &
    1,np,1,mp,l,nv)

  !!!!! Vegetation Forcings !!!!!!
  if(isies.ne.1)then
    call convection(forceSI_xv(:,:,:,itheta), &
      forceSI_xvfuel(:,:,:,:,isies),xvrhotmp,xvfueltmp, &
      1-ih,np+ih,1-ih,mp+ih,1,np,1,mp,0,l,lfuel,nv,nvfuel,nfuel)
    if(ifire.eq.1)then
      call evapPDF(forceSI_xv,forceSI_xvfuel,xvfueltmp, &
        1-ih,np+ih,1-ih,mp+ih,l,1,np,1,mp,lfuel,nv,nvfuel,nfuel)
      call burnFuel(forceSI_xv,forceSI_xvfuel,xvtmp,xvfueltmp, &
        1-ih,np+ih,1-ih,mp+ih,l,1,np,1,mp,lfuel,nv,nvfuel,nfuel)
    endif
  endif

  ! Emission Forcings
  select case (iemissions)
    case (1)
      call pointSource(forceSI_xv(:,:,:,iEmitStart:iEmitStop), &
        1-ih,np+ih,1-ih,mp+ih,l,nEmit,nAero,nMAero)
    case (2)
      call efSource(forceSI_xv(:,:,:,iEmitStart:iEmitStop), &
        forceSI_xvfuel(:,:,:,:,irhof),1-ih,np+ih,1-ih,mp+ih,l, &
        1,np,1,mp,lfuel,nfuel,nEmit)
    case (3)
      call ZBEST(forceSI_xv,xvtmp,-forceSI_xvfuel(:,:,:,:,irhof), &
        temps(:,1:np,1:mp,:),xvfueltmp(:,:,:,:,irhof), &
        1-ih,np+ih,1-ih,mp+ih,l,1,np,1,mp,lfuel,nv,nfuel,dts)
  end select

  if((iemissions.gt.0).and.(nAero.gt.0))then
    call stokesDrag(forceSI_xv(:,:,:,iEmitStart+nEmit:iEmitStop), &
      xvtmp(:,:,:,iEmitStart+nEmit:iEmitStop),tempg,xvtmp(:,:,:,irho), &
      1-ih,np+ih,1-ih,mp+ih,l)
    call deposit(forceSI_xv(:,:,:,iEmitStart+nEmit:iEmitStop), &
      xvtmp(:,:,:,iEmitStart+nEmit:iEmitStop),tempg(:,:,1:lfuel), &
      xvfueltmp(:,:,:,:,irhof),1-ih,np+ih,1-ih,mp+ih,1,np,1,mp, &
      l,lfuel,nfuel)
  endif

end subroutine small_implicit_forcings
