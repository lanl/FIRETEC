!-----------------------------------------------------------------------
! defineFuelArray is the primary allocation function for fuel 
! arrays used throughout the simulation
!-----------------------------------------------------------------------
subroutine defineFuelArray
  use gridlist_variables, only : ih
  use vegetation_gridlist_variables, only : nfuel
  use gridsetup, only : np,mp
  use fuel_variables, only : lfuel,sizeScale,actualFuelDepth,temps, &
    rhoFuelInitial,kfindex,kfrank
  use xvall, only : xvfuel,xvfueltmp,nvfuel,isies
  use forcings, only : forceSE_xvfuel,forceSI_xvfuel,forceLE_xvfuel, &
    xvfuelLim
  Implicit None

  allocate(xvfuel(nfuel,np,mp,lfuel,nvfuel))
  allocate(xvfueltmp(nfuel,np,mp,lfuel,nvfuel))
  allocate(xvfuelLim(2,nvfuel))
  xvfuelLim(1,:)=-1e20; xvfuelLim(2,:)=1e20
  allocate(forceSE_xvfuel(nfuel,np,mp,lfuel,nvfuel))
  allocate(forceSI_xvfuel(nfuel,np,mp,lfuel,nvfuel))
  allocate(forceLE_xvfuel(nfuel,np,mp,lfuel,nvfuel))
  allocate(rhoFuelInitial(nfuel,np,mp,lfuel));
  allocate(sizeScale(nfuel,np,mp,lfuel));
  allocate(actualFuelDepth(nfuel,np,mp));
  allocate(kfindex(np,mp,nfuel))
  allocate(kfrank(np,mp,nfuel))
  if(isies.ne.1)then
    allocate(temps(nfuel,1-ih:np+ih,1-ih:mp+ih,lfuel))
  endif

end subroutine defineFuelArray
