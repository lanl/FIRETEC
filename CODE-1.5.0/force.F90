!-----------------------------------------------------------------------
! forcing applies the large step explicit forcings 
!-----------------------------------------------------------------------
subroutine forcing
  use gridlist_variables, only : l,ih
  use forcings, only : forceLE_xv,forceLE_xvfuel,xvLim,xvfuelLim
  use fuel_variables, only : lfuel
  use gridsetup, only : np,mp
  use gridlist_variables, only : iord
  use forcings, only : uavg,vavg,oavg
  use xvall, only : xv,nv,xvfuel,nvfuel,xe
  use gridsetup, only : np,mp,dt
  use metric_variables_old
  Implicit None

  ! Local Variables 
  integer :: kv
 
  ! Executable Code
  call setAdvectiveVelocities(uavg,vavg,oavg,1-ih,np+ih,1-ih,mp+ih,l, &
    dt,iord) ! this is advvel in code. 
  call mpdata(xv,h,xe,1-ih,np+ih,1-ih,mp+ih,l,nv) ! LES forcings
  
  do kv=1,nv
    !!!!! Transported Forcings !!!!!!
    call update(forceLE_xv(:,:,:,kv),forceLE_xv(:,:,:,kv),np,mp,l, &
      1-ih,np+ih,1-ih,mp+ih,0,0)
    call donorcell(forceLE_xv(:,:,:,kv),1-ih,np+ih,1-ih,mp+ih,l)
    xv(1:np,1:mp,:,kv)=xv(1:np,1:mp,:,kv)+forceLE_xv(1:np,1:mp,:,kv)*dt
    call forcingLimits(xv(1:np,1:mp,:,kv),xvLim(kv,:),1,np,1,mp,l,1)
    call update(xv(:,:,:,kv),xv(:,:,:,kv),np,mp,l, &
      1-ih,np+ih,1-ih,mp+ih,0,0) 
  enddo

  !!!!! Fuel Forcings !!!!!
  xvfuel=xvfuel+forceLE_xvfuel*dt
  call forcingLimits(xvfuel,xvfuelLim,1,np,1,mp,lfuel,nvfuel)

end subroutine forcing

!-----------------------------------------------------------------------
! Small clippings to xv arrays in case beyond the physically possible
!-----------------------------------------------------------------------
subroutine forcingLimits(xv,xvLimits,il,iu,jl,ju,lls,nvp)
  use gridlist_variables, only : prec
  Implicit None

  ! Local Variables
  integer,intent(in) :: il,iu,jl,ju,lls,nvp
  real(prec),intent(inout) :: xv(il:iu,jl:ju,lls,nvp)
  real(prec),intent(in) :: xvLimits(nvp,2)
  integer :: i,j,k,kv

  ! Executable Code
  do kv=1,nvp
    do i=il,iu
      do j=jl,ju
        do k=1,lls
          xv(i,j,k,kv)=max(xvLimits(kv,1),min(xvLimits(kv,2), &
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
  use gridlist_variables, only : iturb,ih,l,idiffsies,irad,nfuel
  use gridsetup, only : np,mp
  use forcings, only : forceSE_xv,forceSE_xvfuel
  use fuel_variables, only : lfuel
  use thermo_variables, only : tempg
  use xvall, only : itemp,nv,xv,isies,irhof,irhow,xvfuel
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
    call firerad_MC(forceSE_xv(1:np,1:mp,:,itemp),1,np,1,mp,l)
  elseif(irad.eq.3)then ! Radiation Sink
    call radiationSink(forceSE_xv(1:np,1:mp,:,itemp),1,np,1,mp,l)
  endif
      
  if(idiffsies.eq.1) call diffusesies(xvfuel(:,:,:,:,isies), &
    xvfuel(:,:,:,:,irhof),xvfuel(:,:,:,:,irhow), &
    forceSE_xvfuel(:,:,:,:,irhof),forceSE_xvfuel(:,:,:,:,isies), &
    tempg(:,:,1:lfuel),1,np,1,mp,lfuel,nfuel,ih)

end subroutine small_explicit_forcings

!-----------------------------------------------------------------------
! large_explicit_forcings calculates all the forcing terms that will be
! applied only once explicitly on the large time steps
!-----------------------------------------------------------------------
subroutine large_explicit_forcings
  use gridlist_variables, only : ifire,iheatsource,ifbrand,ih,l,nfuel
  use forcings, only : forceSE_xv,forceLE_xv,forceSE_xvfuel, &
    forceLE_xvfuel
  use gridsetup, only : np,mp
  use xvall, only : xv,xvrho,xvfuel,nv,iuvel,ivvel,iwvel,irho, &
    nvfuel,irhof,isies
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
  if(iheatsource.ge.1) call heatSource(iheatsource, &
    forceLE_xv(1:np,1:mp,:,:),1,np,1,mp,l,nv)

  !!!!! Firebrand Forcings !!!!!
  if(ifbrand.eq.1) call firebrand(xvrho(:,:,:,iuvel), &
    xvrho(:,:,:,ivvel),xvrho(:,:,:,iwvel),1-ih,np+ih,1-ih,mp+ih,0,l, &
    xv(1:np,1:mp,:,irho),tempg(1:np,1:mp,:),xvfuel(:,:,:,:,irhof), &
    temps(:,1:np,1:mp,:),1,np,1,mp,lfuel,nfuel)

end subroutine large_explicit_forcings

!-----------------------------------------------------------------------
! small_implicit_forcings calculates all the forcing terms that will be
! applied implicitly on the small time steps
!-----------------------------------------------------------------------
subroutine small_implicit_forcings
  use gridlist_variables, only : l,icorio,ilspgf,nfuel,iemissions, &
    slopeangle,slopeazimuth,dts,ifire,prec
  use forcings, only : forceSI_xv,forceSI_xvfuel
  use xvall, only : xe,xvtmp,xvfueltmp,xvrho,nv,nvfuel,relaxxv,itemp, &
    isies,iuvel,ivvel,iwvel,irho,irhof,irhow,iEmitStart,iEmitStop,iH2O
  use lspgf_variables, only : sinthetaf
  use gridsetup, only : np,mp
  use constants, only : g,pi,fcor2,fcor3
  use emission_gridlist_variables, only : nEmit,nAero,nMAero
  use fuel_variables, only : temps,lfuel
  Implicit None

  ! Local Variables
  integer :: i,j,k,kref 
  real(prec) :: uc,vc ! reference wind for coriolis
  real(prec) :: piOver180,g1,g2,g3 ! projected gravity coeff

  ! Executable Code
  
  !!!!! Transported Forcings !!!!!!
  ! forces : pressure
  call gradmoa()

  ! forces : gravity
  piOver180 = pi/180.0
  g1 = g*sin(slopeangle*piOver180)*cos(slopeazimuth*piOver180)
  forceSI_xv(:,:,:,iuvel)=forceSI_xv(:,:,:,iuvel) &
    -(xvtmp(:,:,:,irho)-xe(:,:,:,irho))*g1
  g2 = g*sin(slopeangle*piOver180)*sin(slopeazimuth*piOver180)
  forceSI_xv(:,:,:,ivvel)=forceSI_xv(:,:,:,ivvel) &
    -(xvtmp(:,:,:,irho)-xe(:,:,:,irho))*g2
  g3 = g*cos(slopeangle*piOver180)
  forceSI_xv(:,:,:,iwvel)=forceSI_xv(:,:,:,iwvel) &
    -(xvtmp(:,:,:,irho)-xe(:,:,:,irho))*g3

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
  call boundary(forceSI_xv(1:np,1:mp,:,:),xvtmp(1:np,1:mp,:,:), &
    xe(1:np,1:mp,:,:),relaxxv,dts,1,np,1,mp,l,nv)

  !!!!! Fuel Forcings !!!!!!
  if (ifire.eq.1) then
    call convection(forceSI_xv(1:np,1:mp,1:lfuel,itemp), &
      forceSI_xvfuel(:,:,:,:,isies),xvrho(1:np,1:mp,1:lfuel,:), &
      xvfueltmp,1,np,1,mp,lfuel,nv,nvfuel,nfuel)
    call burnFuel(forceSI_xv(1:np,1:mp,1:lfuel,:),forceSI_xvfuel, &
      xvtmp(1:np,1:mp,1:lfuel,:),xvfueltmp, &
      1,np,1,mp,lfuel,nv,nvfuel,nfuel)
  endif

  ! Emission Forcings
  select case (iemissions)
    case (1)
      call pointSource(forceSI_xv(1:np,1:mp,:,iEmitStart:iEmitStop), &
        1,np,1,mp,l,nEmit,nAero,nMAero)
    case (2)
      call efSource(forceSI_xv(1:np,1:mp,1:lfuel,iEmitStart:iEmitStop),&
        -forceSI_xvfuel(:,:,:,:,irhof),-forceSI_xvfuel(:,:,:,:,irhow), &
        1,np,1,mp,lfuel,nfuel,nEmit,iH2O-iEmitStart+1)
    case (3)
      call ZBEST(forceSI_xv(1:np,1:mp,1:lfuel,:), &
        xvtmp(1:np,1:mp,1:lfuel,:),-forceSI_xvfuel(:,:,:,:,irhof), &
        temps(:,1:np,1:mp,:),-forceSI_xvfuel(:,:,:,:,irhow), &
        xvfueltmp(:,:,:,:,irhof),1,np,1,mp,lfuel,nv,nfuel,dts)
  end select

end subroutine small_implicit_forcings
