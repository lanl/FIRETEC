!-----------------------------------------------------------------------
! Reads basic user input for emissions in HIGRAD/FIRETEC
!-----------------------------------------------------------------------
subroutine input_emissions(iemissions,ifire)
  use gridlist_variables, only : wp,pressground,Tambient
  use emission_gridlist_variables
  use emission_general_variables, only : spEmit,cpEmit,mwEmit,spAero, &
    rhoAero,specificHumidity
  use emission_factor_variables, only : efEmit,efAero
  use QueryGridlist_mod, only : QueryGridlist_variable, &
    QueryGridlist_array
  use constants, only : tolerance
  Implicit None
  
  ! Local Variables
  integer,intent(in) :: ifire
  integer,intent(inout) :: iemissions
  integer :: isp
  
  real(wp) :: pw
  real(wp),allocatable :: cpEmit_tmp(:),mwEmit_tmp(:)
  real(wp),allocatable :: efEmit_tmp(:)
  character(len=6),allocatable :: spEmit_tmp(:)

  ! Executable Code
  if((ifire.eq.1).and.(iemissions.eq.0)) iemissions=2
  
  !-----Gaseous Species-----!
  call QueryGridlist_variable('nEmit',nEmit,48)
  allocate(spEmit(1:nEmit))
  call QueryGridlist_array('Species',spEmit,nEmit,48)
  allocate(cpEmit(1:nEmit)); cpEmit(:)=1040.0_wp  ! J/kg*K (Initializing all values as N2)
  allocate(mwEmit(1:nEmit)); mwEmit(:)=28.014_wp ! kg/kmol (Initializing all values as N2)
  call QueryGridlist_array('cpEmit',cpEmit,nEmit,48)
  call QueryGridlist_array('mwEmit',mwEmit,nEmit,48)
  
  ! Fixed Emissions Source Factor
  if(iemissions.eq.2) then
    allocate(efEmit(nEmit)); efEmit(:)=1.6_wp ! kg/kg (Initializing all values as CO2)
    do isp=1,nEmit
      if(spEmit(isp).eq."O2") efEmit(isp)=0._wp
    enddo
    call QueryGridlist_array('efEmit',efEmit,nEmit,48)
  endif
 
  !-----Water as a gaseous species-----!
  call QueryGridlist_variable('relativehumidity',relativeHumidity,48)
  ! Set specific humidity
  ! water pressure saturation in Pa (from Humidity Conversion Formulas, Vaisala)
  ! OK at 0.083% within the range 50-100°C
  pw=relativeHumidity*611.64_wp* &
    exp(17.62_wp*(Tambient-273.15_wp)/(Tambient-30.13_wp)) ! water vapor pressure (Pa)
  specificHumidity=0.62199_wp*pw/(pressground-0.378*pw)
  if((relativeHumidity.gt.tolerance).and.ALL(spEmit.ne."H2O")) then
    if(iemissions.eq.0) iemissions=2
    
    allocate(spEmit_tmp(nEmit),cpEmit_tmp(nEmit),mwEmit_tmp(nEmit))
    spEmit_tmp=spEmit
    cpEmit_tmp=cpEmit
    mwEmit_tmp=mwEmit
    deallocate(spEmit,cpEmit,mwEmit)

    nEmit=nEmit+1
    allocate(spEmit(nEmit),cpEmit(nEmit),mwEmit(nEmit))
    spEmit(1)="H2O";      spEmit(2:nEmit)=spEmit_tmp(:)
    cpEmit(1)=1990.0_wp;  cpEmit(2:nEmit)=cpEmit_tmp(:)  ! J/kg*K
    mwEmit(1)=18.02_wp;   mwEmit(2:nEmit)=mwEmit_tmp(:)  ! kg/kmol
    deallocate(spEmit_tmp,cpEmit_tmp,mwEmit_tmp)
   
    if((iemissions.eq.2).and.(ifire.eq.1))then
      allocate(efEmit_tmp(nEmit-1))
      efEmit_tmp=efEmit
      deallocate(efEmit); allocate(efEmit(nEmit))
      efEmit(1)=0.56_wp
      efEmit(2:nEmit)=efEmit_tmp(:)
      deallocate(efEmit_tmp)
    endif
  endif
      
  ! Fire needs oxygen to resolve it 
  if(ifire.eq.1.and.ALL(spEmit.ne."O2"))then ! Oxygen needed for fire
    allocate(spEmit_tmp(nEmit),cpEmit_tmp(nEmit),mwEmit_tmp(nEmit))
    spEmit_tmp=spEmit
    cpEmit_tmp=cpEmit
    mwEmit_tmp=mwEmit
    deallocate(spEmit,cpEmit,mwEmit)

    nEmit=nEmit+1
    allocate(spEmit(nEmit),cpEmit(nEmit),mwEmit(nEmit))
    spEmit(1)="O2";   spEmit(2:nEmit)=spEmit_tmp(:)
    cpEmit(1)=918.0_wp;   cpEmit(2:nEmit)=cpEmit_tmp(:)  ! J/kg*K
    mwEmit(1)=31.998_wp; mwEmit(2:nEmit)=mwEmit_tmp(:)  ! kg/kmol
    deallocate(spEmit_tmp,cpEmit_tmp,mwEmit_tmp)
   
    if(iemissions.eq.2)then
      allocate(efEmit_tmp(nEmit-1))
      efEmit_tmp=efEmit
      deallocate(efEmit); allocate(efEmit(nEmit))
      efEmit(1)=-0.5448_wp/0.4552_wp
      efEmit(2:nEmit)=efEmit_tmp(:)
      deallocate(efEmit_tmp)
    endif
  endif

  !------Aerosol Species-------!
  call QueryGridlist_variable('nAero',nAero,48)    ! Number of aerosol species
  call QueryGridlist_variable('nMAero',nMAero,48)  ! Number of moments transport for the aerosol size distribution (must be at least 2)
  allocate(spAero(nAero))
  call QueryGridlist_array('SpeciesAero',spAero,nAero,48)
  allocate(rhoAero(nAero)); rhoAero(:)=1850.0_wp ! kg/m3 (Initializing as incipient soot)
  call QueryGridlist_array('rhoAero',rhoAero,nAero,48)
  
  ! Fixed Emissions Source Factor
  if(iemissions.eq.2) then
    allocate(efAero(nAero)); efAero(:)=1.e-2_wp ! kg/kg (Initializing all values)
    call QueryGridlist_array('efAero',efAero,nAero,48)
  endif

end subroutine input_emissions
