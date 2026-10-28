!-----------------------------------------------------------------------
! Gridlist variables set to defaults here and changed only in the 
! gridlist itself
!-----------------------------------------------------------------------
module gridlist_variables
  use, intrinsic :: iso_fortran_env, only : real32,real64
  
  Implicit None

#ifdef DBL_PREC
  integer,parameter :: wp = 8
#else
  integer,parameter :: wp = 4
#endif
  
  integer :: irst             = 0         ! Restart flag
  integer :: nt               = 10000     ! Total number of large timesteps
  integer :: nts              = 10        ! Number of small timesteps (higrad) in a large timestep
  real(wp):: dt               = 0.01_wp   ! Temporal size of a large timestep
  integer :: n                = 200       ! Number of cells in the x dimension
  integer :: m                = 200       ! Number of cells in the y dimension
  integer :: l                = 40        ! Number of cells in the z dimension
  real(wp):: dx               = 2._wp     ! Size of cell in x dimension
  real(wp):: dy               = 2._wp     ! Size of cell in y dimension
  real(wp):: dz               = 15._wp    ! Size of cell in z dimension before stretching
  real(wp):: aa1              = 0.1_wp    ! Stretching parameter for z direction
  integer :: nprocx           = 8         ! Number of processors in the x dimension
  integer :: nprocy           = 8         ! Number of processors in the y dimension
  integer :: ih               = 1         ! Number of halo cells
  integer :: nr               = 5         ! Number of cells for lateral relaxation (irlx, irly)
  integer :: ibcx             = 1         ! X/Y boundary conditions (0=cyclic,1<= different relaxation options)
  integer :: ibcy             = 1         ! X/Y boundary conditions (0=cyclic,1<= different relaxation options)
  integer :: ibc_top          = 7825      ! Lateral Top boundary conditions (0=cyclic,1<=4 different relaxation options,5=symmetric)
  integer :: ibc_bot          = 7825      ! Lateral Bottom boundary conditions (0=cyclic,1<=4 different relaxation options,5=symmetric)
  integer :: ibc_left         = 7825      ! Lateral Left boundary conditions (0=cyclic,1<=4 different relaxation options,5=symmetric)
  integer :: ibc_right        = 7825      ! Lateral Right boundary conditions (0=cyclic,1<=4 different relaxation options,5=symmetric)
  integer :: ibctopopen       = 0         ! Open outlet, feeslip option at top
  integer :: iab              = 1         ! Options for top absorption initialization
  real(wp):: zab              = 400._wp   ! Bottom of top damper (m)
  real(wp):: tow              = 10._wp    ! Magnitude parameter for top relaxation
  real(wp):: zStabBot         = 400._wp   ! Height of bottom stability layer (m) (common values are 400. and 800.)
  real(wp):: zStabMid         = 450._wp   ! Height of middle stability layer (m) (common values are 450. and 1000.)
  real(wp):: thetaStabMid     = 0._wp     ! Theta increase in the middle stability layer (K/m) (common values are 8./50. and 8.0/200.0)
  real(wp):: thetaStabTop     = 0._wp     ! Theta increase in the top stability layer (K/m) (common values are 6./1000. and 3.0/1000.0)
  integer :: iturb            = 2         ! Turbulence resolution options
  integer :: iemissions       = 0         ! Option for tracking particle emissions
  integer :: iheatsource      = 0         ! Option for adding artificial heat source
  integer :: ifbrand          = 0         ! Option for tracking firebrands
  integer :: icorio           = 0         ! Option fors coriolis effect
  integer :: ilspgf           = 0         ! Large scale pressure gradient forcing options
  integer :: izlspgf          = 1         ! Reference for max flux is whole domain or a slice at height zu (1)
  integer :: frqlspgf         = 10000     ! Large scale pressure gradient forcing frequency
  integer :: islip            = 0         ! Option for no-slip boundary condition forcing
  real(wp):: Tambient         = 300.      ! Ambient temperature (K)
  real(wp):: pressground      = 101325.   ! Reference pressure for ground (Pa)
  real(wp):: zgroundref       = 0.        ! Reference elevation for ground (m)
  integer :: iperturb         = 0         ! Perturbation option for cyclic runs
  real(wp):: u0               = 3._wp     ! Initial ambient x direction wind magnitude (m/s)
  integer :: uswitch          = 2         ! Options for x direction ambient wind profiles
  real(wp):: v0               = 0._wp     ! Initial ambient y direction wind magnitude (m/s)
  integer :: vswitch          = 0         ! Options for y direction ambient wind profiles
  real(wp):: zu               = 25._wp    ! Height of referenced ambient winds
  integer :: ius              = 0         ! Lower x indice for reference zone in lai computation (when 0 whole domain is used) for uswitch=2
  integer :: iue              = 0         ! Upper x indice for reference zone in lai computation (when 0 whole domain is used) for uswitch=2
  integer :: jus              = 0         ! Lower y indice for reference zone in lai computation (when 0 whole domain is used) for vswitch=2
  integer :: jue              = 0         ! Upper y indice for reference zone in lai computation (when 0 whole domain is used) for vswitch=2
  real(wp):: slopeangle       = 0._wp     ! Archaic slope usage (use topo now) for rotating gravity
  real(wp):: slopeazimuth     = 0._wp     ! Archaic slope usage (use topo now) for rotating gravity
  integer :: iord             = 2         ! Advection order of resolution
  integer :: nonos            = 1         ! Advection for non-oscillatory flow option
  integer :: idiv             = 1         ! Divergent flows option
  integer :: nfct             = 1
  integer :: nonosold         = 1         ! Non-oscillatory option in mpdataold
  integer :: frqoutput        = 100       ! Frequency (large time step) of output
  integer :: frqfilstr        = 0         ! Frequency (large time step) of file filtering
  integer :: icfmeflag        = 0         ! Option to read-in slices of the grid
  integer :: ihdf             = 1         ! Use HDF5 for I/O (uneffective if you compile without -DHDF)
  integer :: windfieldout     = 0         ! Option for outputting a windfield
  integer :: windfieldin      = 0         ! Option for outputting a windfield
  integer :: ixevariation     = 0         ! Option for outputting a windfield
  integer :: isensor          = 0         ! Switch for using sensors

  character(len=257) :: outname     = 'comp.out' ! Name of output files
  character(len=257) :: restartfile = ''         ! Name of restart file
  character(len=257) :: topofile    = ''         ! Name of topo file

end module gridlist_variables  

!-----------------------------------------------------------------------
! MSGA arrays and associated variables
!-----------------------------------------------------------------------
module msga_variables
  use mpi
  Implicit None
  
#ifdef DBL_PREC
  integer, parameter :: mpi_wp = mpi_double_precision
#else
  integer, parameter :: mpi_wp = mpi_real
#endif

  integer :: mpi_rank,npos,mpos,numprocs
  integer :: blocktype,blocktype0
  integer :: rightedge=0
  integer :: leftedge=0
  integer :: botedge=0
  integer :: topedge=0
  integer :: peleft,peright,peabove,pebelow
  integer :: perightabove,perightbelow,peleftbelow,peleftabove
  integer :: ierror=0

end module msga_variables

!-----------------------------------------------------------------------
! Forcing arrays from physics and boundary forcings
!-----------------------------------------------------------------------
module forcings
  use gridlist_variables, only : wp

  Implicit None
      
  real(wp),allocatable :: forceSE_xv(:,:,:,:),forceLE_xv(:,:,:,:)
  real(wp),allocatable :: forceSI_xv(:,:,:,:)
  real(wp),allocatable :: forceSE_xvfuel(:,:,:,:,:)
  real(wp),allocatable :: forceLE_xvfuel(:,:,:,:,:)
  real(wp),allocatable :: forceSI_xvfuel(:,:,:,:,:)
  real(wp),allocatable :: xvLim(:,:),xvfuelLim(:,:) ! limits on forcings
  real(wp),allocatable :: tkewght(:)
  real(wp),allocatable :: uavg(:,:,:),vavg(:,:,:),wavg(:,:,:)
  real(wp),allocatable :: oavg(:,:,:)
  real(wp),allocatable :: wt(:),wtf(:)
  real(wp),allocatable :: u1(:,:,:),u2(:,:,:),u3(:,:,:) ! face contravariant velocities*dt/di/gi

end module forcings

!-----------------------------------------------------------------------
! Xv arrays and associated variables/pointers
!-----------------------------------------------------------------------
module xvall
  use gridlist_variables, only : wp

  Implicit None
      
  real(wp),allocatable :: xe(:,:,:,:)           ! environmental field
  real(wp),allocatable :: xv(:,:,:,:)           ! transported field
  real(wp),allocatable :: xvrho(:,:,:,:)        ! xv/rho field
  real(wp),allocatable :: xvfuel(:,:,:,:,:)     ! non-advected field
  real(wp),allocatable :: xvtmp(:,:,:,:)        ! MOA xv field
  real(wp),allocatable :: xvrhotmp(:,:,:,:)     ! xv/rho field
  real(wp),allocatable :: xvfueltmp(:,:,:,:,:)  ! MOA xvfuel field
  real(wp),allocatable :: relaxxv(:,:,:,:)      ! relaxation coeff for bc
  character(len=12),allocatable :: xv_list(:)
  character(len=12),allocatable :: xvfuel_list(:)

  ! Following are the pointers to the xvarray making it more dynamic
  ! xv array
  integer :: iuvel=1,ivvel=1,iwvel=1
  integer :: itheta=1,ika=1,ikb=1,iO2=1,iH2O=1,irho=1
  integer :: iM0=1,iM1=1
  integer :: imixfrac=1
  integer :: iEmitStart=1,iEmitStop=0
  integer :: nv

  ! xvfuel array
  integer :: irhof=1,irhow=1,isies=1,ipsiw=1
  integer :: nvfuel

end module xvall

!-----------------------------------------------------------------------
! Thermodynamic arrays and associated variables
!-----------------------------------------------------------------------
module thermo_variables
  use gridlist_variables, only : wp

  Implicit None
  
  real(wp) :: cpwater=4186.0_wp ! J/kg*K
  real(wp),allocatable :: pr(:,:,:),pre(:,:,:)
  real(wp),allocatable :: cp_gas(:,:,:),cv_gas(:,:,:)
  real(wp),allocatable :: mw_gas(:,:,:)
  real(wp),allocatable :: tempg(:,:,:)

end module thermo_variables

!-----------------------------------------------------------------------
! Metric arrays and associated variables
!-----------------------------------------------------------------------
module metric_variables
  use gridlist_variables, only : wp

  Implicit None
 
  real(wp) :: aa2,aa3 
  real(wp),allocatable :: x(:),y(:),z(:),zedge(:)
  real(wp),allocatable :: zs(:,:)
  real(wp),allocatable :: gmul(:),c13(:,:),c23(:,:)
  real(wp),allocatable :: h(:,:,:),gi(:,:,:)

  real(wp) :: zb

end module metric_variables

!-----------------------------------------------------------------------
! Large pressure gradient arrays and associated variables
!-----------------------------------------------------------------------
module lspgf_variables
  use gridlist_variables, only : wp

  Implicit None
  
  real(wp), allocatable:: massFluxTemp(:) ! mass flux in the direction of geostrophic wind (ug,vg), temporary on local domain
  real(wp), allocatable:: massFlux(:) ! mass flux in the direction of geostrophic wind (ug,vg), after mpi_allreduce
  real(wp) :: totMassFlux, targMassFlux

  real(wp) :: cosg,sing ! geostrophic wind cos and sin
  integer :: iwindx 
  ! if iwindx==1 mass flux is computed for each j (as if wind aligned with x axis) : n*vg<m*ug
  ! if iwindx==0 mass flux is computed for each i (as if wind aligned with y axis) : n*vg>m*ug
  integer :: nMassFlux ! size of massFluxArray (either n or m depending of iwindx)
  real(wp) :: intsintheta ! vertical integral of sintheta
  real(wp),allocatable:: sintheta(:,:,:) ! evaluation of v/sqrt(u2+v2)
  real(wp),allocatable:: sinthetaf(:,:,:) ! evaluation of f*v/sqrt(u2+v2)
  real(wp),allocatable:: flspgf(:,:) ! evaluation of f
  ! based on ekman theory to impose a pressure gradient (ilspgf)

end module lspgf_variables

!-----------------------------------------------------------------------
! Gridsetup variables include parallelization variables and 
! general temporal and spatial grid variables.
!-----------------------------------------------------------------------
module gridsetup
  use gridlist_variables, only : wp

  Implicit None

  integer :: timestep=0 ! Temporal counter in iterations
  real(wp):: time=0.0_wp ! Temporal counter in seconds
  real(wp):: dts
  integer :: np,mp
  integer :: itrestart=0 ! Iteration of restart
  integer :: ittot=0
  real(wp):: dxi,dyi,dzi,dti

end module gridsetup

!-----------------------------------------------------------------------
! Constant variables do not change throughout the simulation.
!-----------------------------------------------------------------------
module constants
  use gridlist_variables, only : wp

  Implicit None

  real(wp) :: pi    = acos(-1.0_wp)     ! pi, yuuummmmm!!!
  real(wp) :: g     = 9.81_wp          ! gravity acceleration (m/s2)
  real(wp) :: Rgas  = 8314.5_wp        ! ideal gas constant (J/kmol K)
  real(wp) :: pref  = 101325.0_wp
  real(wp) :: tolerance = 1.e-8_wp  ! Tolerance of error
  real(wp) :: Na    = 6.0221409e26_wp  ! Avogadro's number (#/kmol)
  real(wp) :: kB    = 1.38064852e-23_wp ! Boltzmann's constant (J/K)
  real(wp) :: sigma = 5.670374419e-8_wp ! Stefan-Boltzmann constant (W/(m^2K^4))
  real(wp) :: fcor2 = 1.4584e-4_wp*cos(acos(-1.0_wp)+55.0_wp/180.0_wp) ! I don't know what this is
  real(wp) :: fcor3 = 1.4584e-4_wp*sin(acos(-1.0_wp)+55.0_wp/180.0_wp) ! I don't know what this is

end module constants

!-----------------------------------------------------------------------
! Ignition arrays and associated variables
!-----------------------------------------------------------------------
module ign_variables
  use gridlist_variables, only : wp

  Implicit None

  integer :: nIgn  
  real(wp) :: targetTemp,rampRate
  integer,allocatable :: ignLocation(:,:)
  real(wp),allocatable :: ignTime(:)

end module ign_variables

!-----------------------------------------------------------------------
! namefile subroutine derives names for specific files from macro
! data and iteration number
!-----------------------------------------------------------------------
subroutine namefile(itn,fnamein,fnameout)
  Implicit None

  ! Local Variables
  integer,intent(in) :: itn
  character(*),intent(in) :: fnamein
  character(*),intent(inout) :: fnameout
  
  integer :: indx
  character(len=32) :: fend  
       
  ! Executable Code 
  write(fend,"(i0)") itn
  indx=len_trim(fnamein)
  fnameout=fnamein(1:indx)//'.'//fend

end subroutine namefile
