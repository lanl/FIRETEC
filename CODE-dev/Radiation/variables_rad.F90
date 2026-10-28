!-----------------------------------------------------------------------
! Gridlist variables related to vegetation set to defaults here 
! and changed only in the gridlist itself
!-----------------------------------------------------------------------
module radiation_gridlist_variables
  use gridlist_variables, only : wp
  use, intrinsic :: iso_fortran_env, only : int64
  Implicit None
  
  ! Executable Code
  integer :: irad       = 2           ! Radiation options
  integer :: icallrad   = 5           ! Frequency of radiation call
  integer :: isootmodel = 0           ! Options for radiating soot model
  integer :: isolar     = 0           ! Heat fuels with solar radiation
  real(wp) :: plotLat   = 35.8800_wp  ! Latitude of the plot (default is Los Alamos)                 
  real(wp) :: plotLong  = -106.3031_wp! Longitude of the plot (default is Los Alamos)
  real(wp) :: crad      = 50._wp      ! Soot radiation constant
  integer(int64) :: unixTimeStart= 1491483600 ! Unix time (default is 6 April 2017 at 7:00 AM in Los Alamos)

end module radiation_gridlist_variables

!-----------------------------------------------------------------------
! General radiation arrays and associated variables
!-----------------------------------------------------------------------
module radiation_general_variables
  use gridlist_variables, only : wp

  Implicit None
  
  real(wp) :: papvthresh
  real(wp),allocatable :: firad (:,:,:),fsiesrad(:,:,:,:)

end module radiation_general_variables

!-----------------------------------------------------------------------
! Radiation arrays and associated variables
!-----------------------------------------------------------------------
module radiation_montecarlo_variables
  use gridlist_variables, only : wp

  Implicit None
  
  integer :: lmc
  integer :: nphotbatch0=1e6
  integer :: nphotbatch
  real(wp) :: Emin=1000.0_wp    ! Minimal energy for a photon
  real(wp) :: Etot
  real(wp),allocatable :: sourcesol(:,:,:,:),sourcegas(:,:,:)
  real(wp),allocatable :: papvift(:,:,:,:),papvgas(:,:,:)
  real(wp),allocatable :: papvtotAZ(:,:,:),volumeAZ(:,:,:)
  real(wp),allocatable :: zCartEdgeAZGlobal(:,:,:,:)
  real(wp),allocatable :: EAZ(:,:,:)
  integer,allocatable :: nphotemisAZsplit(:,:,:)
  integer,allocatable :: nphotonAZgather(:,:,:),nphotonAZ(:,:,:)

end module radiation_montecarlo_variables

!-----------------------------------------------------------------------
! Solar heating arrays and associated variables
!-----------------------------------------------------------------------
module solar_variables
  use gridlist_variables, only : wp
  Implicit None

  ! Executable Code
  real(wp) :: aod500      = 0.1_wp   ! Aerosol optical depth at 500-550 nm
  real(wp) :: col_d_w     = 0.015_wp ! Total depth of column wpipitable water [m]
  real(wp) :: col_d_O3    = 0.003_wp ! Total depth of column wpipitable ozone [m-atm]
  real(wp) :: cloud_frac  = 0.05_wp  ! Total cloud cover fraction 
  real(wp) :: leaf_albedo = 0.10_wp  ! Estimated albedo for vegetation
  real(wp) :: leaf_emis   = 0.97_wp  ! Estimated emissivity for vegetation

end module solar_variables
