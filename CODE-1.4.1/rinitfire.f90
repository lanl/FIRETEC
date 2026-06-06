!----------------------------------------------------------------
! rinitfire initializes all arrays and variables related to 
! FIRETEC for use throughout the simulation
!----------------------------------------------------------------
subroutine rinitfire
  use gridlist_variables, only : n,m,l,ih,nfuel,ivegread,ignfile, &
    irhovapor,tambient,cpwood,ifire
  use gridsetup, only : np,mp
  use msga_variables, only : mpi_rank
  use fuel_variables, only : lfuel,rhoFuel,rhoWater,sizeScale,actualFuelDepth,sies, &
    temps,rhoFuelInitial
  use ign_variables, only : ignLocation,ignTime,nIgn,targetTemp, &
    rampRate
  use thermo_variables, only : specificHumidity,tempg,cpwater
  use metric_variables, only : zedge,zs
  use xvall, only : xe,iO2,irho,ivapor
  use gridlist_variables, only : cpwood
  use thermo_variables, only : cpwater
  Implicit None
  
  ! Local Variables
  integer :: i,j,k,ift
  integer :: fsize
  real :: cpsolid
  real,allocatable :: rmoist(:,:,:,:),fuelDepth(:,:,:,:)
  real,external :: zcart
  namelist/ignitelist/ nIgn,targetTemp,rampRate

  ! Executable Code
  ! Set environmental higrad arrays
  xe(:,:,:,iO2)=0.233*xe(:,:,:,irho)*(1.-specificHumidity)
  if(irhovapor.eq.1) xe(:,:,:,ivapor)=specificHumidity*xe(:,:,:,irho)
  
  if(ivegread.eq.1)then
    inquire(file='treesrhof.dat',size=fsize)
    ! Divide by 4 bytes
    lfuel=(fsize-8)/(4*n*m*nfuel)
    if(mpi_rank.eq.0) print*,'lfuel for this simulations is ',lfuel
    call defineFuelArray
    if(mpi_rank.eq.0)then
      open(unit=92,file='treesrhof.dat',form='unformatted')
      open(unit=93,file='treesss.dat',form='unformatted')
      open(unit=94,file='treesmoist.dat',form='unformatted')
      open(unit=95,file='treesfueldepth.dat',form='unformatted')
    endif
    allocate(rmoist(nfuel,np,mp,lfuel),fuelDepth(nfuel,np,mp,lfuel))
    do ift=1,nfuel
      call readio(rhoFuel(ift,:,:,:),92,lfuel)
      call readio(sizeScale(ift,:,:,:),93,lfuel)
      call readio(rmoist(ift,:,:,:),94,lfuel)
      call readio(fuelDepth(ift,:,:,:),95,lfuel)
     
      if(mpi_rank.eq.0) print*,'Fuel Type ',ift 
      call rmaxmin(rhoFuel(ift,:,:,:),'rhoFuel',1,np,1,mp,lfuel,1)
      call rmaxmin(sizeScale(ift,:,:,:),'szsc',1,np,1,mp,lfuel,1)
      call rmaxmin(rmoist(ift,:,:,:),'moist',1,np,1,mp,lfuel,1)
      call rmaxmin(fuelDepth(ift,:,:,1),'afd',1,np,1,mp,1,1)
    enddo
    rhoWater=rmoist*rhoFuel
    actualFuelDepth=max(1.e-6,fuelDepth(:,:,:,1))
    sizeScale=max(5.e-4,sizeScale)
    if(mpi_rank.eq.0)then
      close(92)
      close(93)
      close(94)
      close(95)
    endif
    deallocate(rmoist,fuelDepth)
  else ! Default vegetation  if (ivegread.eq.0) then
    lfuel=1
    call defineFuelArray
    sizeScale=0.0005
    actualFuelDepth=0.7
    do j=1,mp
      do i=1,np
        rhoFuel(:,i,j,1)=1.*0.7/(zcart(zedge(2),i,j)-zs(i,j))
        rhoWater(:,i,j,1)=0.05*rhoFuel(:,i,j,1)
      enddo
    enddo
    do ift=1,nfuel
      if(mpi_rank.eq.0) print*,'Fuel Type ',ift 
      call rmaxmin(rhoFuel(ift,:,:,:),'rhoFuel',1,np,1,mp,lfuel,1)
      call rmaxmin(sizeScale(ift,:,:,:),'sizeScale',1,np,1,mp,lfuel,1)
      call rmaxmin(rhoWater(ift,:,:,:),'rhoWater',1,np,1,mp,lfuel,1)
    enddo
  endif ! end of ivegread loop

  if(ifire.eq.1)then
    
    ! Initialize temperatures  
    rhoFuelInitial=rhoFuel
    call rmaxmin(tempg,'tempg',1-ih,np+ih,1-ih,mp+ih,l,1)
    do ift=1,nfuel
      do k=1,lfuel
        do j=1,mp
          do i=1,np
            cpsolid=(rhoFuel(ift,i,j,k)*cpwood+rhoWater(ift,i,j,k)*cpwater)/(rhoFuel(ift,i,j,k)+rhoWater(ift,i,j,k))
            temps(ift,i,j,k)=tambient
            sies(ift,i,j,k)=temps(ift,i,j,k)*cpsolid
          enddo
        enddo
      enddo
    enddo

    ! Import ignitions
    open(unit=1000,file=ignfile,form='formatted',status='old')
    read(1000,nml=ignitelist)
    allocate(ignLocation(nIgn,3))
    allocate(ignTime(nIgn))
    if(mpi_rank.eq.0) print*,'RampRate',rampRate,'TempTarget',targetTemp,'nIgn',nIgn
    do i=1,nIgn
      read(1000,*) ignLocation(i,1),ignLocation(i,2),ignLocation(i,3),igntime(i)
    enddo
    close(1000)
  endif

end subroutine rinitfire
