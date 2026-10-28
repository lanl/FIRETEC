!-----------------------------------------------------------------------
! rinitfire initializes all arrays and variables related to 
! FIRETEC for use throughout the simulation
!-----------------------------------------------------------------------
subroutine rinitVeg
  use, intrinsic :: iso_fortran_env, only : real32
  use gridlist_variables, only : n,m,wp
  use vegetation_gridlist_variables, only : nfuel,ivegread,cpwood, &
    ignfile,ifire
  use gridsetup, only : np,mp
  use msga_variables, only : mpi_rank
  use xvall, only : xvfuel,irhof,irhow,isies,ipsiw
  use fuel_variables, only : lfuel,sizeScale,actualFuelDepth,temps, &
    rhoFuelInitial,min_rhoFuel,kfindex,kfrank
  use ign_variables, only : ignLocation,ignTime,nIgn,targetTemp, &
    rampRate
  use forcings, only : xvfuelLim
  use thermo_variables, only : cpwater
  use metric_variables, only : zedge,zs
  use thermo_variables, only : cpwater,tempg
  use rmaxmin_mod, only : rmaxmin
  Implicit None
  
  ! Local Variables
  integer :: i,j,k,ift
  integer :: fsize
  real(wp) :: cpsolid
  real(wp),external :: zcart
  real(wp),allocatable :: afdtmp(:)
  real(wp),allocatable :: rmoist(:,:,:,:),fuelDepth(:,:,:,:)
  namelist/ignitelist/ nIgn,targetTemp,rampRate

  ! Executable Code
  if(ivegread.eq.1)then
    inquire(file='treesrhof.dat',size=fsize)
    ! Legacy vegetation files use 4-byte real payloads.
    lfuel=(fsize-8)/(wp*n*m*nfuel)
    if(mpi_rank.eq.0) print*,'lfuel for this simulation is ',lfuel
    call defineFuelArray
    if(mpi_rank.eq.0)then
      open(unit=92,file='treesrhof.dat',form='unformatted')
      open(unit=93,file='treesss.dat',form='unformatted')
      open(unit=94,file='treesmoist.dat',form='unformatted')
      open(unit=95,file='treesfueldepth.dat',form='unformatted')
    endif
    allocate(rmoist(nfuel,np,mp,lfuel),fuelDepth(nfuel,np,mp,lfuel))
    rmoist = 0.0_wp
    fuelDepth = 0.0_wp

    do ift=1,nfuel
      call readio(xvfuel(ift,:,:,:,irhof),92,lfuel)
      call readio(sizeScale(ift,:,:,:),93,lfuel)
      call readio(rmoist(ift,:,:,:),94,lfuel)
      call readio(fuelDepth(ift,:,:,:),95,1)
     
      if(mpi_rank.eq.0) print*,'Fuel Type ',ift 
      call rmaxmin(xvfuel(ift,:,:,:,irhof),'rhoFuel')
      call rmaxmin(sizeScale(ift,:,:,:),'szsc')
      call rmaxmin(rmoist(ift,:,:,:),'moist')
      call rmaxmin(fuelDepth(ift,:,:,1),'afd')
    enddo
    xvfuel(:,:,:,:,irhof)=max(min_rhoFuel,xvfuel(:,:,:,:,irhof))
    xvfuel(:,:,:,:,irhow)=rmoist*xvfuel(:,:,:,:,irhof)
    actualFuelDepth=max(0.000001_wp,fuelDepth(:,:,:,1))
    sizeScale=max(0.0005_wp,sizeScale)
    if(mpi_rank.eq.0)then
      close(92)
      close(93)
      close(94)
      close(95)
    endif
    deallocate(rmoist,fuelDepth)
    xvfuel(:,:,:,:,irhof)=max(min_rhoFuel,xvfuel(:,:,:,:,irhof))
  else ! Default vegetation  if (ivegread.eq.0) then
    lfuel=1
    call defineFuelArray
    sizeScale=0.0005_wp
    actualFuelDepth=0.7_wp
    do j=1,mp
      do i=1,np
        xvfuel(:,i,j,1,irhof)=1._wp*0.7_wp/(zcart(zedge(2),i,j)-zs(i,j))
        xvfuel(:,i,j,1,irhow)=0.05_wp*xvfuel(:,i,j,1,irhof)
      enddo
    enddo
    do ift=1,nfuel
      if(mpi_rank.eq.0) print*,'Fuel Type ',ift 
      call rmaxmin(xvfuel(ift,:,:,:,irhof),'rhoFuel')
      call rmaxmin(sizeScale(ift,:,:,:),'sizeScale')
      call rmaxmin(xvfuel(ift,:,:,:,irhow),'rhoWater')
    enddo
  endif ! end of ivegread loop

  xvfuelLim(1,irhof)=min_rhoFuel
  xvfuelLim(1,irhow)=0.0_wp
  
  ! sort individual fuel zones by afd
  allocate(afdtmp(nfuel)); afdtmp = 0.0_wp
  do j=1,mp
    do i=1,np
      afdtmp=actualFuelDepth(:,i,j)
      do k=1,nfuel
        ift=minloc(afdtmp,dim=1)
        kfindex(i,j,k)=ift
        kfrank(i,j,ift)=k
        afdtmp(ift)=maxval(actualFuelDepth(:,i,j))+1
      enddo
    enddo
  enddo 
  deallocate(afdtmp)

  if(isies.ne.1)then
    !Initialize temperatures
    rhoFuelInitial=xvfuel(:,:,:,:,irhof)
    do ift=1,nfuel
      do k=1,lfuel
        do j=1,mp
          do i=1,np
            temps(ift,i,j,k)=tempg(i,j,k)
            if(xvfuel(ift,i,j,k,irhof)+xvfuel(ift,i,j,k,irhow).gt.0)then
              cpsolid=(xvfuel(ift,i,j,k,irhof)*cpwood+ &
                xvfuel(ift,i,j,k,irhow)*cpwater)/ &
                (xvfuel(ift,i,j,k,irhof)+xvfuel(ift,i,j,k,irhow))
              xvfuel(ift,i,j,k,isies)=temps(ift,i,j,k)*cpsolid
            else
              xvfuel(ift,i,j,k,isies)=0._wp
            endif
          enddo
        enddo
      enddo
    enddo
    xvfuelLim(1,isies)=0.0_wp
  endif

  if(ifire.eq.1)then
    ! Initialize moisture PDF 
    xvfuel(:,:,:,:,ipsiw)=0._wp
    xvfuelLim(1,ipsiw)=0._wp
    xvfuelLim(2,ipsiw)=1._wp

    ! Import ignitions
    open(unit=1000,file=ignfile,form='formatted',status='old')
    read(1000,nml=ignitelist)
    allocate(ignLocation(nIgn,3))
    allocate(ignTime(nIgn))
    if(mpi_rank.eq.0) print*,'RampRate',rampRate, &
      'TempTarget',targetTemp,'nIgn',nIgn
    do i=1,nIgn
      read(1000,*) ignLocation(i,1),ignLocation(i,2),ignLocation(i,3), &
        ignTime(i)
    enddo
    close(1000)
  endif
end subroutine rinitVeg
