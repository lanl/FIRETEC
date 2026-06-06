subroutine definearray
  use gridlist_variables, only : nts,ih,ifire,iturb,inonlocal,idiffsies,irhovapor, &
    nfuel,iemissions,ilspgf,iord,iwindfieldin,l,n,m,uswitch,vswitch
  use forcings, only : force,tkewght,uavg,vavg,oavg,wt,wtf,u1,u2,u3
  use xvall, only : xv,xe,xvtmp,relaxxv,xvdataold,xvdatanew,nv, &
    xvrho,uprofile
  use metric_variables, only : x,y,z,zedge,zs,gmul,c13,c23,gi
  use thermo_variables, only : tempg,pr,pre
  use lspgf_variables, only : sintheta,sinthetaf,flspgf
  use gridsetup, only : np,mp
  use turb_variables, only : sa,sb,saxy,saz,zonehts,zonedzs,zonerhos, &
    iftwght,fuelinds,zoneus,zonevs,deltaus,deltavs,sqrtG_Kxy, &
    sqrtG_Kz,sqrtG_KG33,K_axy,K_az,K_b,rtke_abc
  use workavg
  use forcinner
  use metric_variables_old
  use turba
  use xvbin
  use msga_variables, only : mpi_rank
  Implicit None

  ! Executable Code
  allocate(xe(1-ih:np+ih,1-ih:mp+ih,l,nv)); xe=0.
  allocate(xv(1-ih:np+ih,1-ih:mp+ih,l,nv)); xv=0.
  allocate(wt(nts+1)); wt=1./nts
  allocate(wtf(nts+1)); wtf=1.
  allocate(uavg(1-ih:np+ih,1-ih:mp+ih,l))
  allocate(vavg(1-ih:np+ih,1-ih:mp+ih,l))
  allocate(wavg(1-ih:np+ih,1-ih:mp+ih,l))
  allocate(oavg(1-ih:np+ih,1-ih:mp+ih,l))
  allocate(force(1-ih:np+ih,1-ih:mp+ih,l,np))
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
  allocate(relaxxv(1-ih:np+ih,1-ih:mp+ih,l,4)); relaxxv=0.
  if(iturb.ge.1) then
    allocate(sa(1-ih:np+ih,1-ih:mp+ih,l))
    allocate(saxy(1-ih:np+ih,1-ih:mp+ih,l))
    allocate(saz(1-ih:np+ih,1-ih:mp+ih,l))
    allocate(sb(1-ih:np+ih,1-ih:mp+ih,l))
    !allocate(sc(1-ih:np+ih,1-ih:mp+ih,l))
  endif
  allocate(xvrho(1-ih:np+ih,1-ih:mp+ih,0:l,nv))
  allocate(rtke_abc(1-ih:np+ih,1-ih:mp+ih,0:l))
  allocate(K_axy(1-ih:np+ih,1-ih:mp+ih,l))
  allocate(K_az(1-ih:np+ih,1-ih:mp+ih,l))
  allocate(K_b(1-ih:np+ih,1-ih:mp+ih,l))
  allocate(sqrtG_Kxy(1-ih:np+ih,1-ih:mp+ih,l))
  allocate(sqrtG_Kz(1-ih:np+ih,1-ih:mp+ih,l))
  allocate(sqrtG_KG33(1-ih:np+ih,1-ih:mp+ih,l))
  allocate(AB11c(np+1,mp,l),AB12c(np,mp+1,l),AB13c(np,mp,l+1))
  allocate(AB21c(np+1,mp,l),AB22c(np,mp+1,l),AB23c(np,mp,l+1))
  allocate(AB31c(np+1,mp,l),AB32c(np,mp+1,l),AB33c(np,mp,l+1))
  allocate(Ci(np+1,mp,l),Cj(np,mp+1,l),Ck(np,mp,l+1))
  allocate (xvtmp(1-ih:np+ih,1-ih:mp+ih,l,5))
  if (uswitch.eq.2.or.vswitch.eq.2) & 
    allocate (uprofile(1-ih:np+ih,1-ih:mp+ih,l)) !FP
  allocate (pr(1-ih:np+ih,1-ih:mp+ih,l))
  if(ifire.eq.1) then
    allocate(tempg(1-ih:np+ih,1-ih:mp+ih,l)); tempg=0.
  endif
  if(ilspgf.ge.1) then
    allocate (sintheta(1-ih:np+ih,1-ih:mp+ih,l))
    allocate (sinthetaf(1-ih:np+ih,1-ih:mp+ih,l))
    allocate (flspgf(np,mp))
  endif
  allocate (pre(1-ih:np+ih,1-ih:mp+ih,l))
  allocate (u1(1-ih:np+ih+1,1-ih:mp+ih,l))
  allocate (u2(1-ih:np+ih,1-ih:mp+ih+1,l))
  allocate (u3(1-ih:np+ih,1-ih:mp+ih,l+1))
  if(iemissions.eq.1) call setUpChemistryTables
  if(iwindfieldin.eq.1)then
    allocate (xvbdataold(1-ih:np+ih,1-ih:mp+ih,l,nv))
    allocate (xvbdatanew(1-ih:np+ih,1-ih:mp+ih,l,nv))
  endif

end subroutine definearray

!----------------------------------------------------------------
! defineFuelArray is the primary allocation function for fuel 
! arrays used throughout the simulation
!----------------------------------------------------------------
subroutine defineFuelArray
  use gridlist_variables, only : ih,nfuel,idiffsies
  use gridsetup, only : np,mp
  use fuel_variables, only : lfuel,rhoFuel,rhoWater,sizeScale,actualFuelDepth,psiwmax, &
    temps,sies,fcorr,convht,rhoFuelInitial,siesdiff
  Implicit None

  allocate(rhoFuel(nfuel,np,mp,lfuel));
  allocate(rhoFuelInitial(nfuel,np,mp,lfuel));
  allocate(rhoWater(nfuel,np,mp,lfuel));
  allocate(sizeScale(nfuel,np,mp,lfuel));
  allocate(actualFuelDepth(nfuel,np,mp));
  allocate(sies(nfuel,1-ih:np+ih,1-ih:mp+ih,lfuel))
  allocate(psiwmax(nfuel,np,mp,lfuel)); psiwmax=0.
  allocate(temps(nfuel,1-ih:np+ih,1-ih:mp+ih,lfuel))
  allocate(fcorr(nfuel,2)); fcorr=1. ! Fuel correction factor for fuel depth in bottom cell
  allocate(convht(nfuel,np,mp,lfuel))
  if(idiffsies.eq.1) allocate(siesdiff(nfuel,np,mp,lfuel))

end subroutine defineFuelArray

!----------------------------------------------------------------
! defineRadArray is the primary allocation function for radiation 
! arrays used throughout the simulation
!----------------------------------------------------------------
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

!----------------------------------------------------------------
! xvpointers sets up an array of indexes in the xv array for a more
! dynamic capability in the transporting array
!----------------------------------------------------------------
subroutine xvpointers
  use gridlist_variables, only : ifire,iturb,iemissions,inonlocal,irhovapor
  use xvall
  Implicit None

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
  if(ifire.eq.1)then
    nv = nv+1
    iO2 = nv
  endif
  if(irhovapor.eq.1)then
    nv = nv+1
    ivapor = nv
  endif
  if(iemissions.eq.1)then
    nv = nv+2
    iM0 = nv-1
    iM1 = nv
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
  xv_list(itemp) = 'temp'
  if(iturb.gt.0)then
    xv_list(ika) = 'ka'
    if(iturb.gt.1) xv_list(ikb) = 'kb'
  endif
  if(ifire.eq.1) xv_list(iO2) = 'O2'
  if(irhovapor.eq.1) xv_list(ivapor) = 'rhovapor'
  if(iemissions.eq.1)then
    xv_list(iM0) = 'M0'
    xv_list(iM1) = 'M1'
  endif
  if(inonlocal.eq.1) xv_list(imixfrac) = 'mixfrac'
  xv_list(irho) = 'rho'

end subroutine xvpointers
