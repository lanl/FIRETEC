!-----------------------------------------------------------------------
! Gridlist variables related to sensors set to defaults here 
! and changed only in the gridlist itself
!-----------------------------------------------------------------------
module sensor_gridlist_variables
  use gridlist_variables, only : prec
   
  Implicit None
  
  ! Module Variables
  integer :: se_frq_write     = 1       ! Frequency in time cycles to write sensor data
  integer :: numsensors       = 5       ! default number of sensors
  integer :: ncycle           = 0       ! time step to begin sensors
  integer :: locfile          = 0       ! if providing text file for sensor locations
  character(len=260) :: sensorfile  = 'sensor' ! Name of the sensor output files
  character(len=60) :: sensorlocs   = 'sensorlocs.txt' !File that contains three columns - X, Y, and Z locations for all sensors - must be tab-delimited

  real(prec),allocatable :: se_xloc(:),se_yloc(:),se_zloc(:)  

end module sensor_gridlist_variables

!-----------------------------------------------------------------------
! General emission variables used throughout the emissions module
!-----------------------------------------------------------------------
module sensor_general_variables
  use gridlist_variables, only : prec
  
  Implicit None

  ! Module Variables
  integer,allocatable    :: se_proc(:)
  integer,allocatable    :: se_i(:),se_j(:),se_k(:)
  real(prec),allocatable :: se_x(:),se_y(:),se_z(:)
  real(prec),allocatable :: xv_values(:,:)

end module sensor_general_variables
