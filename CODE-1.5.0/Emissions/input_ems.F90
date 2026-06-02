!-----------------------------------------------------------------------
! Reads basic user input for emissions in HIGRAD/FIRETEC
!-----------------------------------------------------------------------
subroutine input_emissions(iemissions,ifire)
  use gridlist_variables, only : prec
  use emission_gridlist_variables
  use emission_general_variables, only : spEmit,cpEmit,mwEmit,spAero, &
    rhoAero
  use emission_factor_variables, only : efEmit,efAero
  Implicit None
  
  ! Local Variables
  integer,intent(in) :: iemissions,ifire
  
  character(len=6),allocatable :: spEmit_tmp(:)
  real(prec),allocatable :: cpEmit_tmp(:),mwEmit_tmp(:)
  real(prec),allocatable :: efEmit_tmp(:)

  ! Executable Code
  !-----Gaseous Species-----!
  call QueryGridlist_integer('nEmit',nEmit,48)
  allocate(spEmit(nEmit))
  call QueryGridlist_string_array('Species',spEmit,nEmit,48)
  allocate(cpEmit(nEmit)); cpEmit(:)=1040.  ! J/kg*K (Initializing all values as N2)
  allocate(mwEmit(nEmit)); mwEmit(:)=28.014 ! kg/kmol (Initializing all values as N2)
  call QueryGridlist_real_array('cpEmit',cpEmit,nEmit,48)
  call QueryGridlist_real_array('mwEmit',mwEmit,nEmit,48)
  
  ! Fixed Emissions Source Factor
  if(iemissions.eq.2) then
    allocate(efEmit(nEmit)); efEmit(:)=1.6 ! kg/kg (Initializing all values as CO2)
    call QueryGridlist_real_array('efEmit',efEmit,nEmit,48)
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
    cpEmit(1)=918.;   cpEmit(2:nEmit)=cpEmit_tmp(:)  ! J/kg*K
    mwEmit(1)=31.998; mwEmit(2:nEmit)=mwEmit_tmp(:)  ! kg/kmol
    deallocate(spEmit_tmp,cpEmit_tmp,mwEmit_tmp)
    
    if(iemissions.eq.2)then
      allocate(efEmit_tmp(nEmit-1))
      efEmit_tmp=efEmit
      deallocate(efEmit); allocate(efEmit(nEmit))
      efEmit(1)=0.
    endif
  endif

  ! print*,'TEST',spEmit,cpEmit,mwEmit ! debug print statement for assigned vals 

  !------Aerosol Species-------!
  call QueryGridlist_integer('nAero',nAero,48)    ! Number of aerosol species
  call QueryGridlist_integer('nMAero',nMAero,48)  ! Number of moments transport for the aerosol size distribution (must be at least 2)
  allocate(spAero(nAero))
  call QueryGridlist_string_array('SpeciesAero',spAero,nAero,48)
  allocate(rhoAero(nAero)); rhoAero(:)=1850. ! kg/m3 (Initializing as incipient soot)
  call QueryGridlist_real_array('rhoAero',rhoAero,nAero,48)
  
  ! Fixed Emissions Source Factor
  if(iemissions.eq.2) then
    allocate(efAero(nAero)); efAero(:)=1.e-2 ! kg/kg (Initializing all values)
    call QueryGridlist_real_array('efAero',efAero,nAero,48)
  endif

end subroutine input_emissions
