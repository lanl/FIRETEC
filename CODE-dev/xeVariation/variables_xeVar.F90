!-----------------------------------------------------------------------
! U/V Ramp up variables and arrays
!-----------------------------------------------------------------------
module uvRamp_variables
  use gridlist_variables, only : wp

  Implicit None

  real(wp)    :: uramp            = 0._wp    ! Ramp increase of ambient x direction wind magnitude (m/s)
  real(wp)    :: utarget          = 0._wp    ! Target increase of ambient x direction wind magnitude (m/s)
  real(wp)    :: uramptime        = 30.0_wp  ! Time to reach target increase ambient x direction wind (s)
  real(wp)    :: vramp            = 0._wp    ! Ramp increase of ambient y direction wind magnitude (m/s)
  real(wp)    :: vtarget          = 0._wp    ! Target increase of ambient y direction wind magnitude (m/s)
  real(wp)    :: vramptime        = 30.0_wp  ! Time to reach target increase ambient y direction wind (s)

end module uvRamp_variables

!-----------------------------------------------------------------------
! Windfield input/output variables and arrays
!-----------------------------------------------------------------------
module windfieldio
  use gridlist_variables, only : wp

  Implicit None

  !!! Gridlist Variables !!!
  integer :: windspeedupfactor= 1       ! Factor of speed up between wind run and fire run
  integer :: itwindfield      = 80000   ! Timestep to start saving wind data
  integer :: itinterp         = 10      ! # of timesteps between windfield interpolation
  integer :: ibcells          = 5       ! # of x cells saved on boundary
  integer :: jbcells          = 5       ! # of y cells saved on boundary
  integer :: is               = 0       ! Lower x indice for production xvdata file on subdomain (0 = whole domain is used)
  integer :: ie               = 0       ! Upper x indice for production xvdata file on subdomain (0 = whole domain is used)
  integer :: js               = 0       ! Lower y indice for production xvdata file on subdomain (0 = whole domain is used)
  integer :: je               = 0       ! Upper y indice for production xvdata file on subdomain (0 = whole domain is used)
  character(len=257) :: windfieldstartfile= 'windfieldstart'! Name of windfield start file
  character(len=257) :: xvdataname        = 'xvdata'        ! Name of xvdata files
  
  !!! Modular Variables !!!
  integer :: nvwind
  integer :: itabsold,itabsnew
  integer,allocatable :: w2f_index(:)
  real(wp),allocatable :: xvdataold(:,:,:,:),xvdatanew(:,:,:,:)
  character(len=12),allocatable :: xvwind_list(:)
  
end module windfieldio

!-----------------------------------------------------------------------
! Sensor arrays and variables
!-----------------------------------------------------------------------
module sensor_variables
  use gridlist_variables, only : wp

  Implicit None

  !! Gridlist Variables !!
  integer :: nSensors
  real(wp),allocatable :: SensorTimes(:)
  real(wp),allocatable :: SensorX(:),SensorY(:),SensorZ(:)
  real(wp),allocatable :: SensorUvel(:),SensorVvel(:)

  !! Module Variables !!
  real(wp),allocatable :: senw(:)
  real(wp),allocatable :: senDataOld(:,:),senDataNew(:,:)
  character(len=257) :: senDataName = 'sensor.dat'  ! Name of sensor data files

end module sensor_variables
