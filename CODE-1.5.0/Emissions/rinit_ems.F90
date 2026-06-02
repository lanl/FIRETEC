!-----------------------------------------------------------------------
! rinitEmission initializes all arrays and variables related to 
! smoke, emissions, and plume evolution throughout the simulation
!-----------------------------------------------------------------------
subroutine rinitEmission(iemissions,xe,il,iu,jl,ju,lls,nv)
  use gridlist_variables, only : prec
  use emission_gridlist_variables, only : nMAero
  use emission_pointSource_variables, only : nEmitPoints, &
    emitPointLocation,emitPointRate,emitPointSpecies,nAeroPoints, &
    aeroPointLocation,aeroPointRate,aeroPointSpecies
  use emission_general_variables, only : spEmit
  use msga_variables, only : mpi_rank
  use thermo_variables, only : specificHumidity
  use xvall, only : iO2,iH2O,irho,iEmitStart,iEmitStop
  use forcings, only : xvLim
  Implicit None

  ! Local Variables
  integer,intent(in) :: iemissions
  integer,intent(in) :: il,iu,jl,ju,lls,nv
  real(prec),intent(inout) :: xe(il:iu,jl:ju,lls,nv)

  integer :: i
  logical :: check
  character(len=10) :: text
  character(len=1)  :: equals

  ! Executable Code
  ! Initialize ambient conditions
  if(ANY(spEmit.eq."O2")) xe(:,:,:,iO2)=0.233*xe(:,:,:,irho)* &
    (1.-specificHumidity)
  if(ANY(spEmit.eq."H2O")) xe(:,:,:,iH2O)=specificHumidity* &
    xe(:,:,:,irho)
  
  ! Initialize source terms
  if(iemissions.eq.1) then ! Read in emission point source data
    ! Gaseous Emission Point Source
    inquire(file='emission_points.dat',exist=check)
    if(check)then
      open(unit=48,file='emission_points.dat',form='formatted', &
        status='old')
      read(48,*) text,equals,nEmitPoints
      read(48,*) ! Space for header
      allocate(emitPointLocation(nEmitPoints,3))
      allocate(emitPointRate(nEmitPoints))
      allocate(emitPointSpecies(nEmitPoints))
      if(mpi_rank.eq.0) print*, &
        'Source Point Emissions, # of Emission Points',nEmitPoints
      do i=1,nEmitPoints
        read(48,*) emitPointSpecies(i),emitPointRate(i), &
          emitPointLocation(i,:)
      enddo
      close(48)
    endif

    ! Aerosol Emission Point Source
    inquire(file='aerosol_points.dat',exist=check)
    if(check)then
      open(unit=48,file='aerosol_points.dat',form='formatted', &
        status='old')
      read(48,*) text,equals,nAeroPoints
      read(48,*) ! Space for header
      allocate(aeroPointLocation(nAeroPoints,3))
      allocate(aeroPointRate(nAeroPoints,nMAero))
      allocate(aeroPointSpecies(nAeroPoints))
      if(mpi_rank.eq.0) print*,'Source Point Aerosol Emissions,', &
        '# of Emission Points',nAeroPoints
      do i=1,nAeroPoints
        read(48,*) aeroPointSpecies(i),aeroPointRate(i,:), &
          aeroPointLocation(i,:)
        aeroPointRate(i,:)=10**aeroPointRate(i,:)
      enddo
      close(48)
    endif
  elseif(iemissions.eq.2)then

  elseif(iemissions.eq.3)then

  endif
  
  xvLim(iEmitStart:iEmitStop,1)=0.

end subroutine rinitEmission
