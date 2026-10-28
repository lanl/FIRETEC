!-----------------------------------------------------------------------
! Gridlist variables related to vegetation set to defaults here 
! and changed only in the gridlist itself
!-----------------------------------------------------------------------
module vegetation_gridlist_variables
  use gridlist_variables, only : wp
  Implicit None
  
  ! Executable Code
  integer :: ifire        = 1          ! Firetec capability
  integer :: ifuel        = 1          ! Option for burning fuel pdf choice
  integer :: nfuel        = 1          ! Number of fuel types for multiple-fuels
  integer :: ivegread     = 0          ! Option for fuel bed read-in
  integer :: inonlocal    = 0          ! Option for non-local burning
  integer :: idiffsies    = 0          ! Option for diffusing ignitions for creeping fires
  real(wp) :: rhoMicro  = 500._wp       ! Fuel's true density (kg/m3)
  real(wp) :: cpwood    = 2500._wp      ! Fuel's heat capacity (J/kg K)
  character(len=60) :: ignfile = 'ignite.dat' ! Name of ignition file

end module vegetation_gridlist_variables

!-----------------------------------------------------------------------
! Fuel arrays and associated variables
!-----------------------------------------------------------------------
module fuel_variables
  use gridlist_variables, only : wp
  Implicit None

  ! Executable Code
  integer :: lfuel 
  real(wp) :: min_rhoFuel=0.000001_wp
  real(wp) :: rnfuel=0.4552_wp   ! stoichiometric fuel constant (Drysdale pp. 179)
  real(wp) :: rno=0.5448_wp      ! stoichiometric oxygen constant (Drysdale pp. 179)
  real(wp) :: tcrit=600.0_wp      ! temperature of pyrolysis
  real(wp) :: twvap=373.0_wp      ! temperature of water vaporization
  real(wp) :: tfstep=310.0_wp     ! tail of the psi ramp curve
  real(wp) :: hf=8913480      ! heat of reaction for simple wood  (J/Kg of products)
  real(wp) :: hwevap=2257000  ! energy needed for water evap at 373 K (Moran)(J/Kg)
  integer,allocatable :: kfindex(:,:,:)
  integer,allocatable :: kfrank(:,:,:)
  real(wp),allocatable :: sizeScale(:,:,:,:),actualFuelDepth(:,:,:)
  real(wp),allocatable :: temps(:,:,:,:)
  real(wp),allocatable :: rhoFuelInitial(:,:,:,:)

end module fuel_variables
