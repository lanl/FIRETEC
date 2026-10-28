!-----------------------------------------------------------------------
! defineRadArray is the primary allocation function for radiation 
! arrays used throughout the simulation
!-----------------------------------------------------------------------
subroutine defineRadArray
  use gridlist_variables, only : n,m,l,wp
  use radiation_gridlist_variables, only : irad
  use vegetation_gridlist_variables, only : nfuel
  use gridsetup, only : np,mp
  use fuel_variables, only : lfuel
  use radiation_montecarlo_variables, only : sourcesol,sourcegas, &
    papvift,papvgas,papvtotAZ,volumeAZ,zCartEdgeAZGlobal, &
    EAZ,nphotemisAZsplit,nphotonAZgather,nphotonAZ
  use radiation_general_variables, only : firad,fsiesrad
  Implicit None

  allocate(firad(np,mp,l)); firad=0._wp
  allocate(fsiesrad(nfuel,np,mp,lfuel)); fsiesrad=0._wp
  if(irad.eq.0) return
  allocate(sourcesol(nfuel,np,mp,lfuel+nfuel)); sourcesol=0._wp
  allocate(sourcegas(np,mp,l)); sourcegas=0._wp
  
  allocate(papvift(nfuel,np,mp,lfuel+nfuel)); papvift = 0._wp  !projected area per unit volume of the solid
  allocate(papvgas(np,mp,l+nfuel)); papvgas = 0._wp !projected area per unit volume of the gas
  allocate(papvtotAZ(n,m,l+nfuel)); papvtotAZ = 0._wp !FPAZ
  allocate(volumeAZ(np,mp,l+nfuel)); volumeAZ = 0._wp !FPAZ volume of cell
  allocate(zCartEdgeAZGlobal(4,n,m,0:l+nfuel))
  allocate(nphotemisAZsplit(n,m,l+nfuel))
  allocate(nphotonAZgather(n,m,l+nfuel))
  allocate(nphotonAZ(n,m,l+nfuel))
  allocate(EAZ(n,m,l+nfuel))

end subroutine defineRadArray
