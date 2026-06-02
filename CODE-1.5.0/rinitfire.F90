!-----------------------------------------------------------------------
! rinitfire initializes all arrays and variables related to 
! FIRETEC for use throughout the simulation
!-----------------------------------------------------------------------
subroutine rinitfire
  use gridlist_variables, only : n,m,nfuel,ivegread,ignfile, &
    tambient,cpwood,ifire,prec
  use gridsetup, only : np,mp
  use msga_variables, only : mpi_rank
  use xvall, only : xvfuel,irhof,irhow,isies,ipsiw
  use fuel_variables, only : lfuel,sizeScale,actualFuelDepth,temps, &
    rhoFuelInitial,min_rhoFuel
  use ign_variables, only : ignLocation,ignTime,nIgn,targetTemp, &
    rampRate
  use forcings, only : xvfuelLim
  use thermo_variables, only : cpwater
  use metric_variables, only : zedge,zs
  use gridlist_variables, only : cpwood
  use thermo_variables, only : cpwater
  use zcart_function
  Implicit None
  
  ! Local Variables
  integer :: i,j,k,ift
  integer :: fsize
  real(prec) :: cpsolid
  real(prec),allocatable :: rmoist(:,:,:,:),fuelDepth(:,:,:,:)
  !real(prec),external :: zcart
  namelist/ignitelist/ nIgn,targetTemp,rampRate

  ! Executable Code
  if(ivegread.eq.1)then
    inquire(file='treesrhof.dat',size=fsize)
    ! Divide by 4 bytes
    lfuel=(fsize-8)/(4*n*m*nfuel)
#ifdef DBL_PREC
    lfuel = int(lfuel/2)
#endif
    if(mpi_rank.eq.0) print*,'lfuel for this simulation is ',lfuel
    call defineFuelArray
    if(mpi_rank.eq.0)then
      open(unit=92,file='treesrhof.dat',form='unformatted')
      open(unit=93,file='treesss.dat',form='unformatted')
      open(unit=94,file='treesmoist.dat',form='unformatted')
      open(unit=95,file='treesfueldepth.dat',form='unformatted')
    endif
    allocate(rmoist(nfuel,np,mp,lfuel),fuelDepth(nfuel,np,mp,lfuel))
    rmoist = 0.0
    fuelDepth = 0.0

    do ift=1,nfuel
      call readio(xvfuel(ift,:,:,:,irhof),92,lfuel)
      call readio(sizeScale(ift,:,:,:),93,lfuel)
      call readio(rmoist(ift,:,:,:),94,lfuel)
      call readio(fuelDepth(ift,:,:,:),95,lfuel)
     
      if(mpi_rank.eq.0) print*,'Fuel Type ',ift 
      call rmaxmin_3D(xvfuel(ift,:,:,:,irhof),'rhoFuel',1,np,1,mp,lfuel,1)
      call rmaxmin_3D(sizeScale(ift,:,:,:),'szsc',1,np,1,mp,lfuel,1)
      call rmaxmin_3D(rmoist(ift,:,:,:),'moist',1,np,1,mp,lfuel,1)
      call rmaxmin_3D(fuelDepth(ift,:,:,1),'afd',1,np,1,mp,1,1)
    enddo
    xvfuel(:,:,:,:,irhof)=max(min_rhoFuel,xvfuel(:,:,:,:,irhof))
    xvfuel(:,:,:,:,irhow)=rmoist*xvfuel(:,:,:,:,irhof)
    actualFuelDepth=max(0.000001,fuelDepth(:,:,:,1))
    sizeScale=max(0.0005,sizeScale)
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
    sizeScale=0.0005
    actualFuelDepth=0.7
    do j=1,mp
      do i=1,np
        xvfuel(:,i,j,1,irhof)=1.*0.7/(zcart(zedge(2),i,j)-zs(i,j))
        xvfuel(:,i,j,1,irhow)=0.05*xvfuel(:,i,j,1,irhof)
      enddo
    enddo
    do ift=1,nfuel
      if(mpi_rank.eq.0) print*,'Fuel Type ',ift 
      call rmaxmin_3D(xvfuel(ift,:,:,:,irhof),'rhoFuel',1,np,1,mp,lfuel,1)
      call rmaxmin_3D(sizeScale(ift,:,:,:),'sizeScale',1,np,1,mp,lfuel,1)
      call rmaxmin_3D(xvfuel(ift,:,:,:,irhow),'rhoWater',1,np,1,mp,lfuel,1)
    enddo
  endif ! end of ivegread loop

  xvfuelLim(irhof,1)=min_rhoFuel
  xvfuelLim(irhow,1)=0.

  if(ifire.eq.1)then
    ! Initialize temperatures  
    rhoFuelInitial=xvfuel(:,:,:,:,irhof)
    xvfuel(:,:,:,:,ipsiw)=0.
    xvfuelLim(ipsiw,1)=0.
    xvfuelLim(ipsiw,2)=1.
    do ift=1,nfuel
      do k=1,lfuel
        do j=1,mp
          do i=1,np
            temps(ift,i,j,k)=tambient
            if(xvfuel(ift,i,j,k,irhof)+xvfuel(ift,i,j,k,irhow).gt.0)then
              cpsolid=(xvfuel(ift,i,j,k,irhof)*cpwood+ &
                xvfuel(ift,i,j,k,irhow)*cpwater)/ &
                (xvfuel(ift,i,j,k,irhof)+xvfuel(ift,i,j,k,irhow))
              xvfuel(ift,i,j,k,isies)=temps(ift,i,j,k)*cpsolid
            else
              xvfuel(ift,i,j,k,isies)=0.
            endif
          enddo
        enddo
      enddo
    enddo
    xvfuelLim(isies,1)=0.

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

end subroutine rinitfire
