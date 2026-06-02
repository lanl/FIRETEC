!-----------------------------------------------------------------------
! Allocate global arrays
!-----------------------------------------------------------------------
subroutine definearray
  use gridlist_variables, only : nts,ih,iturb,ilspgf,l,n,m, &
    uswitch,vswitch
  use forcings, only : force,forceSE_xv,forceSI_xv,forceLE_xv,xvLim, &
    uavg,vavg,oavg,wt,wtf,u1,u2,u3
  use xvall, only : xv,xe,xvtmp,relaxxv,nv,xvrho,uprofile
  use metric_variables, only : x,y,z,zedge,zs,gmul,c13,c23,gi
  use thermo_variables, only : tempg,pr,pre,cp_gas,cv_gas,mw_gas
  use lspgf_variables, only : sintheta,sinthetaf,flspgf
  use gridsetup, only : np,mp
  use workavg
  use forcinner
  use metric_variables_old
  Implicit None

  ! Executable Code
  call xvpointers
  
  allocate(xe(1-ih:np+ih,1-ih:mp+ih,l,nv)); xe=0.
  allocate(xv(1-ih:np+ih,1-ih:mp+ih,l,nv)); xv=0.
  allocate(xvLim(nv,2)); xvLim(:,1)=-1.e20; xvLim(:,2)=1.e20
  allocate(wt(nts+1)); wt=1./nts
  allocate(wtf(nts)); wtf=1./nts
  allocate(uavg(1-ih:np+ih,1-ih:mp+ih,l))
  allocate(vavg(1-ih:np+ih,1-ih:mp+ih,l))
  allocate(wavg(1-ih:np+ih,1-ih:mp+ih,l))
  allocate(oavg(1-ih:np+ih,1-ih:mp+ih,l))
  allocate(force(1-ih:np+ih,1-ih:mp+ih,l,nv))
  allocate(forceSE_xv(1-ih:np+ih,1-ih:mp+ih,l,nv))
  allocate(forceSI_xv(1-ih:np+ih,1-ih:mp+ih,l,nv))
  allocate(forceLE_xv(1-ih:np+ih,1-ih:mp+ih,l,nv))
  allocate(f1(np,mp,l))
  allocate(f2(np,mp,l))
  allocate(f3(np,mp,l))
  allocate(rg_over_prrcp(np,mp,l)) !array for higrad
  allocate(cp_over_cv(np,mp,l)) !arrays for higrad
  allocate(c13(1-ih:np+ih,1-ih:mp+ih))
  allocate(c23(1-ih:np+ih,1-ih:mp+ih))
  allocate(gi(1-ih:np+ih,1-ih:mp+ih,l))
  allocate(gmul(l))
  allocate(h(1-ih:np+ih,1-ih:mp+ih,l))
  allocate(zs(1-ih:np+ih,1-ih:mp+ih)); zs=0.
  allocate(zsio(n,m))
  allocate(x(1-ih:n+ih))
  allocate(y(1-ih:m+ih))
  allocate(z(l))
  allocate(zedge(l+1))
  allocate(relaxxv(np,mp,l,4)); relaxxv=0.
  allocate(xvrho(1-ih:np+ih,1-ih:mp+ih,0:l,nv))
  allocate(xvtmp(1-ih:np+ih,1-ih:mp+ih,l,nv))
  if(iturb.gt.0) call defineTurbArray
  if (uswitch.eq.2.or.vswitch.eq.2) & 
    allocate(uprofile(1-ih:np+ih,1-ih:mp+ih,l)) !FP
  allocate(pr(1-ih:np+ih,1-ih:mp+ih,l))
  allocate(mw_gas(np,mp,l))
  allocate(cp_gas(np,mp,l))
  allocate(cv_gas(np,mp,l))
  allocate(tempg(1-ih:np+ih,1-ih:mp+ih,l))
  if(ilspgf.ge.1) then
    allocate(sintheta(1-ih:np+ih,1-ih:mp+ih,l))
    allocate(sinthetaf(1-ih:np+ih,1-ih:mp+ih,l))
    allocate(flspgf(np,mp))
  endif
  allocate(pre(1-ih:np+ih,1-ih:mp+ih,l))
  allocate(u1(1-ih:np+ih+1,1-ih:mp+ih,l))
  allocate(u2(1-ih:np+ih,1-ih:mp+ih+1,l))
  allocate(u3(1-ih:np+ih,1-ih:mp+ih,l+1))

end subroutine definearray

!-----------------------------------------------------------------------
! defineFuelArray is the primary allocation function for fuel 
! arrays used throughout the simulation
!-----------------------------------------------------------------------
subroutine defineFuelArray
  use gridlist_variables, only : ih,nfuel,ifire
  use gridsetup, only : np,mp
  use fuel_variables, only : lfuel,sizeScale,actualFuelDepth,temps, &
    rhoFuelInitial ! fcorr,rhoFuelInitial
  use xvall, only : xvfuel,xvfueltmp,nvfuel
  use forcings, only : forceSE_xvfuel,forceSI_xvfuel,forceLE_xvfuel, &
    xvfuelLim
  Implicit None

  allocate(xvfuel(nfuel,np,mp,lfuel,nvfuel))
  allocate(xvfueltmp(nfuel,np,mp,lfuel,nvfuel))
  allocate(xvfuelLim(nvfuel,2))
  xvfuelLim(:,1)=-1e20; xvfuelLim(:,2)=1e20
  allocate(forceSE_xvfuel(nfuel,np,mp,lfuel,nvfuel))
  allocate(forceSI_xvfuel(nfuel,np,mp,lfuel,nvfuel))
  allocate(forceLE_xvfuel(nfuel,np,mp,lfuel,nvfuel))
  allocate(rhoFuelInitial(nfuel,np,mp,lfuel));
  allocate(sizeScale(nfuel,np,mp,lfuel));
  allocate(actualFuelDepth(nfuel,np,mp));
  if(ifire.eq.1)then
    allocate(temps(nfuel,1-ih:np+ih,1-ih:mp+ih,lfuel))
    ! allocate(fcorr(nfuel,2)); fcorr=1. ! Fuel correction factor for fuel depth in bottom cell
  endif

end subroutine defineFuelArray

!-----------------------------------------------------------------------
! defineRadArray is the primary allocation function for radiation 
! arrays used throughout the simulation
!-----------------------------------------------------------------------
subroutine defineRadArray
  use gridlist_variables, only : n,m,l,nfuel,irad,iradeastflux
  use gridsetup, only : np,mp
  use fuel_variables, only : lfuel
  use radiation_variables, only : sourcesol,sourcegas,kfindex, &
    papvift,papvgas,papvtotAZ,volumeAZ,zCartEdgeAZGlobal,EAZ, &
    nphotemisAZsplit,nphotonAZgather,nphotonAZ, &
    nphotonEastFlux,eastFlux,firad,fsiesrad
  Implicit None

  allocate(firad(np,mp,l)); firad=0.
  allocate(fsiesrad(nfuel,np,mp,lfuel)); fsiesrad=0.
  if(irad.eq.0) return
  allocate(sourcesol(nfuel,np,mp,lfuel+nfuel)); sourcesol=0.
  allocate(sourcegas(np,mp,l)); sourcegas=0.
  
  allocate(kfindex(np,mp,nfuel))
  allocate(papvift(nfuel,np,mp,lfuel+nfuel)); papvift = 0.0  !projected area per unit volume of the solid
  allocate(papvgas(np,mp,l+nfuel)); papvgas = 0.0 !projected area per unit volume of the gas
  allocate(papvtotAZ(n,m,l+nfuel)); papvtotAZ = 0.0 !FPAZ
  allocate(volumeAZ(np,mp,l+nfuel)); volumeAZ = 0.0 !FPAZ volume of cell
  allocate(zCartEdgeAZGlobal(4,n,m,0:l+nfuel))
  allocate(nphotemisAZsplit(n,m,l+nfuel))
  allocate(nphotonAZgather(n,m,l+nfuel))
  allocate(nphotonAZ(n,m,l+nfuel))
  if(iradeastflux.eq.1)then
    allocate(nphotonEastFlux(n,m,l+nfuel)); nphotonEastFlux=0.
    allocate(eastFlux(np,mp,l)); eastFlux=0.
  endif
  allocate(EAZ(n,m,l+nfuel))

end subroutine defineRadArray

!-----------------------------------------------------------------------
! xvpointers sets up an array of indexes in the xv array for a more
! dynamic capability in the transporting array
!-----------------------------------------------------------------------
subroutine xvpointers
  use gridlist_variables, only : iturb,iemissions,inonlocal,ifire
  use emission_gridlist_variables, only : nEmit,nAero,nMAero
  use emission_general_variables, only : spEmit,spAero
  use xvall
  Implicit None

  ! Local Variables
  integer :: i,j
  integer :: indx

  ! Executable Code
  ! Transported array pointers
  iuvel = 1
  ivvel = 2
  iwvel = 3
  itemp = 4
  nv = 4
  if(iturb.gt.0)then
    nv = nv+1
    ika = nv
    if(iturb.gt.1)then
      nv = nv+1
      ikb = nv
    endif
  endif
  if(iemissions.gt.0)then
    iEmitStart = nv+1
    if(nEmit.gt.0)then
      iEmitStop  = nv+nEmit
      nv = nv+nEmit
      do i=1,nEmit
        if(spEmit(i).eq.'O2')  iO2  = iEmitStart-1+i
        if(spEmit(i).eq.'H2O') iH2O = iEmitStart-1+i
      enddo
    endif
    if(nAero.gt.0)then
      iEmitStop  = nv+nAero*nMAero
      nv = nv + nAero*nMAero
    endif
  endif
  if(inonlocal.eq.1)then
    nv = nv+1
    imixfrac = nv
  endif
  nv = nv+1
  irho = nv
  
  allocate(xv_list(nv))
  xv_list(iuvel) = 'u'
  xv_list(ivvel) = 'v'
  xv_list(iwvel) = 'w'
  xv_list(itemp) = 'theta'
  if(iturb.gt.0)then
    xv_list(ika) = 'ka'
    if(iturb.gt.1) xv_list(ikb) = 'kb'
  endif
  if(iemissions.gt.0)then
    if(nEmit.gt.0)then
      do i=1,nEmit
        xv_list(iEmitStart-1+i) = spEmit(i)
      enddo
    endif
    if(nAero.gt.0)then
      do i=1,nAero
        do j=1,nMAero
          indx=index(spAero(i),' ')
          write(xv_list(iEmitStart+nEmit+nMAero*(i-1)+j-1),"(A,I0)") &
            spAero(i)(1:indx),j
        enddo
      enddo
    endif
  endif
  if(inonlocal.eq.1) xv_list(imixfrac) = 'mixfrac'
  xv_list(irho) = 'rho'

  ! Non-advected solid array pointers
  irhof = 1
  irhow = 2
  nvfuel = 2
  if(ifire.eq.1)then
    isies = nvfuel+1
    ipsiw = nvfuel+2
    nvfuel = nvfuel+2
  endif

  allocate(xvfuel_list(nvfuel))
  xvfuel_list(irhof) = 'rhoFuel'
  xvfuel_list(irhow) = 'rhoWater'
  if(ifire.eq.1)then
    xvfuel_list(isies) = 'sies'
    xvfuel_list(ipsiw) = 'psiwmax'
  endif

end subroutine xvpointers
