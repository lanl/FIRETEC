!-----------------------------------------------------------------------
! U/V Ramp up variables and arrays
!-----------------------------------------------------------------------
module uvRamp_variables
  use gridlist_variables, only : prec

  Implicit None

  real(prec)    :: uramp            = 10.     ! Target of ambient x direction wind magnitude (m/s)
  real(prec)    :: uramptime        = 30.     ! Time to reach target ambient x direction wind (s)
  real(prec)    :: vramp            = 10.     ! Target of ambient y direction wind magnitude (m/s)
  real(prec)    :: vramptime        = 30.     ! Time to reach target ambient y direction wind (s)

end module uvRamp_variables

!-----------------------------------------------------------------------
! Windfield input/output variables and arrays
!-----------------------------------------------------------------------
module windfieldio
  use gridlist_variables, only : prec

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
  real(prec),allocatable :: xvdataold(:,:,:,:),xvdatanew(:,:,:,:)
  character(len=12),allocatable :: xvwind_list(:)
  
end module windfieldio

!-----------------------------------------------------------------------
! Sensor arrays and variables
!-----------------------------------------------------------------------
module sensor_variables
  use gridlist_variables, only : prec

  Implicit None

  !! Gridlist Variables !!
  integer :: nSensors,freqSensor
  real(prec),allocatable :: SensorTimes(:)
  real(prec),allocatable :: SensorX(:),SensorY(:),SensorZ(:)
  real(prec),allocatable :: SensorUvel(:),SensorVvel(:)

  !! Module Variables !!
  real(prec),allocatable :: senw(:)
  real(prec),allocatable :: senDataOld(:,:),senDataNew(:,:)
  character(len=257) :: senDataName = 'sensor.dat'  ! Name of sensor data files

end module sensor_variables
