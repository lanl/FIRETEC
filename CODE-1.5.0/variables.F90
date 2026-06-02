!-----------------------------------------------------------------------
! Gridlist variables set to defaults here and changed only in the 
! gridlist itself
!-----------------------------------------------------------------------
module gridlist_variables
  
  Implicit None

#ifdef DBL_PREC
  integer,parameter :: prec = 8
#else
  integer,parameter :: prec = 4
#endif
  
  integer :: irst             = 0       ! Restart flag
  integer :: nt               = 10000   ! Total number of large timesteps
  integer :: nts              = 10      ! Number of small timesteps (higrad) in a large timestep
  real(prec) :: dts           = 0.001   ! Temporal size of a small timestep (higrad)
  integer :: ntp              = 10      ! Number of small timesteps (firetec) in a large timestep
  integer :: n                = 200     ! Number of cells in the x dimension
  integer :: m                = 200     ! Number of cells in the y dimension
  integer :: l                = 40      ! Number of cells in the z dimension
  real(prec) :: dx            = 2.      ! Size of cell in x dimension
  real(prec) :: dy            = 2.      ! Size of cell in y dimension
  real(prec) :: dz            = 15.     ! Size of cell in z dimension before stretching
  real(prec) :: aa1           = 0.1     ! Stretching parameter for z direction
  integer :: nprocx           = 8       ! Number of processors in the x dimension
  integer :: nprocy           = 8       ! Number of processors in the y dimension
  integer :: ih               = 3       ! Number of halo cells
  integer :: nr               = 5       ! Number of cells for lateral relaxation (irlx, irly)
  integer :: ibcx             = 1       ! X/Y boundary conditions (0=cyclic,1<= different relaxation options)
  integer :: ibcy             = 1       ! X/Y boundary conditions (0=cyclic,1<= different relaxation options)
  integer :: ibclatopen       = 0       ! Top boundary relaxation options
  integer :: ibctopopen       = 0       ! Open outlet, feeslip option at top
  integer :: iab              = 1       ! Options for top absorption initialization
  real(prec) :: zab           = 400.    ! Bottom of top damper (m)
  real(prec) :: zabt          = 400.    ! Bottom of top damper for iab=3 option
  real(prec) :: tow           = 10.     ! Magnitude parameter for top relaxation
  integer :: itheta           = 0       ! Stable layer at domain top
  integer :: ifire            = 1       ! Firetec capability
  integer :: iturb            = 2       ! Turbulence resolution options
  integer :: irad             = 2       ! Radiation options
  integer :: icallrad         = 5       ! Frequency of radiation call
  real(prec) :: crad          = 50.     ! Soot radiation constant
  integer :: irandseed        = 0       ! Random function seed
  integer :: iseed            = 0       ! Seed for randomization function
  integer :: isootmodel       = 0       ! Options for radiating soot model
  integer :: iradeastflux     = 0       ! Specific radiative flux output through
  integer :: inonlocal        = 0       ! Option for non-local burning
  integer :: iemissions       = 0       ! Option for tracking particle emissions
  integer :: iheatsource      = 0       ! Option for adding artificial heat source
  integer :: idiffsies        = 0       ! Option for diffusing ignitions for creeping fires
  real(prec) :: relativeHumidity = 0.   ! Relative humidity if tracking water vapor
  integer :: ifbrand          = 0       ! Option for tracking firebrands
  integer :: icorio           = 0       ! Option fors coriolis effect
  integer :: ilspgf           = 0       ! Large scale pressure gradient forcing options
  integer :: izlspgf          = 1       ! Reference for max flux is whole domain or a slice at height zu (1)
  integer :: frqlspgf         = 10000   ! Large scale pressure gradient forcing frequency
  integer :: islip            = 0       ! Option for no-slip boundary condition forcing
  real(prec) :: tambient      = 300.    ! Ambient temperature (K)
  real(prec) :: pressground   = 101325. ! Reference pressure for ground (Pa)
  real(prec) :: zgroundref    = 0.      ! Reference elevation for ground (m)
  integer :: iperturb         = 0       ! Perturbation option for cyclic runs
  real(prec) :: u0            = 3.      ! Initial ambient x direction wind magnitude (m/s)
  integer :: uswitch          = 2       ! Options for x direction ambient wind profiles
  real(prec) :: v0            = 0.      ! Initial ambient y direction wind magnitude (m/s)
  integer :: vswitch          = 0       ! Options for y direction ambient wind profiles
  real(prec) :: zu            = 25.     ! Height of referenced ambient winds
  integer :: ius              = 0       ! Lower x indice for reference zone in lai computation (when 0 whole domain is used) for uswitch=2
  integer :: iue              = 0       ! Upper x indice for reference zone in lai computation (when 0 whole domain is used) for uswitch=2
  integer :: jus              = 0       ! Lower y indice for reference zone in lai computation (when 0 whole domain is used) for vswitch=2
  integer :: jue              = 0       ! Upper y indice for reference zone in lai computation (when 0 whole domain is used) for vswitch=2
  real(prec) :: slopeangle    = 0.      ! Archaic slope usage (use topo now) for rotating gravity
  real(prec) :: slopeazimuth  = 0.      ! Archaic slope usage (use topo now) for rotating gravity
  integer :: iord             = 2       ! Advection order of resolution
  integer :: nonos            = 1       ! Advection for non-oscillatory flow option
  integer :: idiv             = 1       ! Divergent flows option
  integer :: nfct             = 1
  integer :: nonosold         = 1       ! Non-oscillatory option in mpdataold
  integer :: ifuel            = 1       ! Option for burning fuel pdf choice
  integer :: nfuel            = 1       ! Number of fuel types for multiple-fuels
  integer :: ivegread         = 0       ! Option for fuel bed read-in
  real(prec) :: rhoMicro      = 500.    ! Fuel's true density (kg/m3)
  real(prec) :: cpwood        = 2500.   ! Fuel's heat capacity (J/kg K)
  integer :: frqoutput        = 100     ! Frequency (large time step) of output
  integer :: frqfilstr        = 0       ! Frequency (large time step) of file filtering
  integer :: icfmeflag        = 0       ! Option to read-in slices of the grid
  integer :: ihdf             = 1       ! Use HDF5 for I/O (uneffective if you compile without -DHDF)
  integer :: windfieldout     = 0       ! Option for outputting a windfield
  integer :: windfieldin      = 0       ! Option for outputting a windfield
  integer :: ixevariation     = 0       ! Option for outputting a windfield
  integer :: isensor          = 0       ! Switch for using sensors

  character(len=60)  :: ignfile     = 'ignite.dat'    ! Name of ignition file
  character(len=257) :: outname     = 'comp.out'      ! Name of output files
  character(len=257) :: restartfile = 'comp.out.h5'   ! Name of restart file
  character(len=60)  :: topofile    = ''              ! Name of topo file

end module gridlist_variables  

!-----------------------------------------------------------------------
! MSGA arrays and associated variables
!-----------------------------------------------------------------------
module msga_variables

  use mpi
  Implicit None
  
  integer :: mpi_rank,npos,mpos,numprocs
  integer :: blocktype,blocktype0
  integer :: rightedge,leftedge,botedge,topedge
  integer :: rightdedge=0
  integer :: leftdedge=0
  integer :: botdedge=0
  integer :: topdedge=0
  integer :: peleft,peright,peabove,pebelow
  integer :: perightabove,perightbelow,peleftbelow,peleftabove
  integer :: ierror=0

end module msga_variables

!-----------------------------------------------------------------------
! Forcing arrays from physics and boundary forcings
!-----------------------------------------------------------------------
module forcings
  use gridlist_variables, only : prec

  Implicit None
      
  real(prec),allocatable :: force(:,:,:,:)
  real(prec),allocatable :: forceSE_xv(:,:,:,:),forceLE_xv(:,:,:,:)
  real(prec),allocatable :: forceSI_xv(:,:,:,:)
  real(prec),allocatable :: forceSE_xvfuel(:,:,:,:,:)
  real(prec),allocatable :: forceLE_xvfuel(:,:,:,:,:)
  real(prec),allocatable :: forceSI_xvfuel(:,:,:,:,:)
  real(prec),allocatable :: xvLim(:,:),xvfuelLim(:,:) ! limits on forcings
  real(prec),allocatable :: tkewght(:)
  real(prec),allocatable :: uavg(:,:,:),vavg(:,:,:),oavg(:,:,:)
  real(prec),allocatable :: wt(:),wtf(:)
  real(prec),allocatable :: u1(:,:,:),u2(:,:,:),u3(:,:,:) ! face contravariant velocities*dt/di/gi

end module forcings

!-----------------------------------------------------------------------
! Xv arrays and associated variables/pointers
!-----------------------------------------------------------------------
module xvall
  use gridlist_variables, only : prec

  Implicit None
      
  real(prec),allocatable :: xv(:,:,:,:)       ! transported field
  real(prec),allocatable :: xvfuel(:,:,:,:,:) ! non-advected field
  real(prec),allocatable :: xvfuelLim(:,:)    ! limits on non-advected field
  real(prec),allocatable :: xvrho(:,:,:,:)    ! xv/rho field
  real(prec),allocatable :: xe(:,:,:,:)       ! environmental field
  real(prec),allocatable :: xvtmp(:,:,:,:)        ! MOA xv field
  real(prec),allocatable :: xvfueltmp(:,:,:,:,:)  ! MOA xvfuel field
  real(prec),allocatable :: relaxxv(:,:,:,:)  ! relaxation coeff for bc
  real(prec),allocatable :: uprofile(:,:,:) ! uprofile from LAI
  character(len=12),allocatable :: xv_list(:)
  character(len=12),allocatable :: xvfuel_list(:)

  ! Following are the pointers to the xvarray making it more dynamic
  ! xv array
  integer :: iuvel=1,ivvel=1,iwvel=1
  integer :: itemp=1,ika=1,ikb=1,iO2=1,iH2O=1,irho=1
  integer :: iM0=1,iM1=1
  integer :: imixfrac=1
  integer :: iEmitStart=1,iEmitStop=1
  integer :: nv

  ! xvfuel array
  integer :: irhof=1,irhow=1,isies=3,ipsiw=4
  integer :: nvfuel

end module xvall

!-----------------------------------------------------------------------
! Thermodynamic arrays and associated variables
!-----------------------------------------------------------------------
module thermo_variables
  use gridlist_variables, only : prec

  Implicit None
  
  real(prec) :: cpwater=4186. ! J/kg*K
  real(prec) :: specifichumidity
  real(prec),allocatable :: pr(:,:,:),pre(:,:,:)
  real(prec),allocatable :: cp_gas(:,:,:),cv_gas(:,:,:)
  real(prec),allocatable :: mw_gas(:,:,:)
  real(prec),allocatable :: tempg(:,:,:)

end module thermo_variables

!-----------------------------------------------------------------------
! Fuel arrays and associated variables
!-----------------------------------------------------------------------
module fuel_variables
  use gridlist_variables, only : prec

  Implicit None

  integer :: lfuel 
  real(prec) :: min_rhoFuel=1.e-6 
  real(prec) :: rnfuel=0.4552 ! stoichiometric fuel constant (Drysdale pp. 179)
  real(prec) :: rno=0.5448    ! stoichiometric oxygen constant (Drysdale pp. 179)
  real(prec) :: tcrit=600.    ! temperature of pyrolysis
  real(prec) :: twvap=373.    ! temperature of water vaporization
  real(prec) :: tfstep=310.   ! tail of the psi ramp curve
  real(prec) :: tstep
  real(prec) :: hf=8913.48e3   ! heat of reaction for simple wood  (J/Kg of products)
  real(prec) :: hwevap=2.257e6 ! energy needed for water evap at 373 K (Moran)(J/Kg)
  real(prec),allocatable :: sizeScale(:,:,:,:),actualFuelDepth(:,:,:)
  real(prec),allocatable :: temps(:,:,:,:)
  real(prec),allocatable :: convht(:,:,:,:)
  ! real(prec),allocatable :: fcorr(:,:)
  real(prec),allocatable :: rhoFuelInitial(:,:,:,:) 

end module fuel_variables

!-----------------------------------------------------------------------
! Metric arrays and associated variables
!-----------------------------------------------------------------------
module metric_variables
  use gridlist_variables, only : prec

  Implicit None
 
  real(prec) :: aa2,aa3 
  real(prec),allocatable :: x(:),y(:),z(:),zedge(:)
  real(prec),allocatable :: zs(:,:)
  real(prec),allocatable :: gmul(:),c13(:,:),c23(:,:),gi(:,:,:)

  real(prec) :: zb

end module metric_variables

!-----------------------------------------------------------------------
! Large pressure gradient arrays and associated variables
!-----------------------------------------------------------------------
module lspgf_variables
  use gridlist_variables, only : prec

  Implicit None
  
  real(prec), allocatable:: massFluxTemp(:) ! mass flux in the direction of geostrophic wind (ug,vg), temporary on local domain
  real(prec), allocatable:: massFlux(:) ! mass flux in the direction of geostrophic wind (ug,vg), after mpi_allreduce
  real(prec) :: totMassFlux, targMassFlux

  real(prec) :: cosg,sing ! geostrophic wind cos and sin
  integer :: iwindx 
  ! if iwindx==1 mass flux is computed for each j (as if wind aligned with x axis) : n*vg<m*ug
  ! if iwindx==0 mass flux is computed for each i (as if wind aligned with y axis) : n*vg>m*ug
  integer :: nMassFlux ! size of massFluxArray (either n or m depending of iwindx)
  real(prec) :: intsintheta ! vertical integral of sintheta
  real(prec),allocatable:: sintheta(:,:,:) ! evaluation of v/sqrt(u2+v2)
  real(prec),allocatable:: sinthetaf(:,:,:) ! evaluation of f*v/sqrt(u2+v2)
  real(prec),allocatable:: flspgf(:,:) ! evaluation of f
  ! based on ekman theory to impose a pressure gradient (ilspgf)

end module lspgf_variables

!-----------------------------------------------------------------------
! Gridsetup variables include parallelization variables and 
! general temporal and spatial grid variables.
!-----------------------------------------------------------------------
module gridsetup
  use gridlist_variables, only : prec

  Implicit None

  integer :: timestep=0 ! Temporal counter in iterations
  real(prec) :: time=0 ! Temporal counter in seconds
  real(prec) :: dt,dtp
  integer :: np,mp
  integer :: itrestart=0 ! Iteration of restart
  integer :: ittot=0
  real(prec) :: dxi,dyi,dzi,dti

end module gridsetup

!-----------------------------------------------------------------------
! Constant variables do not change throughout the simulation.
!-----------------------------------------------------------------------
module constants
  use gridlist_variables, only : prec

  Implicit None

  real(prec) :: pi    = acos(-1.)     ! pi, yuuummmmm!!!
  real(prec) :: g     = 9.81          ! gravity acceleration (m/s2)
  real(prec) :: Rgas  = 8314.5        ! ideal gas constant (J/kmol K)
  real(prec) :: pref  = 101325.
  real(prec) :: tolerance = 0.000001  ! Tolerance of error
#if DBL_PREC
  real(prec) :: Na    = 6.0221409d26  ! Avogadro's number (#/kmol)
  real(prec) :: kB    = 1.38064852d-23! Boltzmann's constant (J/K)
  real(prec) :: sigma = 5.678d-8      ! Stefan-Boltzmann constant (W/(m^2K^4))
  real(prec) :: fcor2 = 1.4584d-4*cos(acos(-1.)+55./180.) ! I don't know what this is
  real(prec) :: fcor3 = 1.4584d-4*sin(acos(-1.)+55./180.) ! I don't know what this is
#else
  real(prec) :: Na    = 6.0221409e26  ! Avogadro's number (#/kmol)
  real(prec) :: kB    = 1.38064852e-23! Boltzmann's constant (J/K)
  real(prec) :: sigma = 5.678e-8      ! Stefan-Boltzmann constant (W/(m^2K^4))
  real(prec) :: fcor2 = 1.4584e-4*cos(acos(-1.)+55./180.) ! I don't know what this is
  real(prec) :: fcor3 = 1.4584e-4*sin(acos(-1.)+55./180.) ! I don't know what this is
#endif

end module constants

!-----------------------------------------------------------------------
! Radiation arrays and associated variables
!-----------------------------------------------------------------------
module radiation_variables
  use gridlist_variables, only : prec

  Implicit None
  
  integer :: lmc
  integer :: nphotbatch=1e6
  real(prec) :: Emin=1000.    ! Minimal energy for a photon
  real(prec) :: Etot
  real(prec),allocatable :: sourcesol(:,:,:,:),sourcegas(:,:,:)
  real(prec),allocatable :: papvift(:,:,:,:),papvgas(:,:,:)
  real(prec),allocatable :: papvtotAZ(:,:,:),volumeAZ(:,:,:)
  real(prec),allocatable :: zCartEdgeAZGlobal(:,:,:,:)
  real(prec),allocatable :: EAZ(:,:,:)
  real(prec),allocatable :: eastFlux(:,:,:)
  integer,allocatable :: kfindex(:,:,:)
  integer,allocatable :: nphotemisAZsplit(:,:,:)
  integer,allocatable :: nphotonAZgather(:,:,:),nphotonAZ(:,:,:)
  integer,allocatable :: nphotonEastFlux(:,:,:)

  real(prec),allocatable :: firad (:,:,:),fsiesrad(:,:,:,:)

end module radiation_variables

!-----------------------------------------------------------------------
! Ignition arrays and associated variables
!-----------------------------------------------------------------------
module ign_variables
  use gridlist_variables, only : prec

  Implicit None

  integer :: nIgn  
  real(prec) :: targetTemp,rampRate
  integer,allocatable :: ignLocation(:,:)
  real(prec),allocatable :: ignTime(:)

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
  character(len=10) :: fend  
       
  ! Executable Code 
  write(fend,"(i0)") itn
  indx=index(fnamein,' ')
  fnameout=fnamein(1:indx-1)//'.'//fend

end subroutine namefile

!2345678***************************************************
! module workavg contains variables that are computed in MOA (inner loop)
module workavg
  use gridlist_variables, only : prec

  Implicit None
  real(prec),allocatable::  wavg(:,:,:)
!                   ! uavg, vavg, wavg and oavg are cell centered average
!                     velocities computed in the inner loop   (m/s)
  real(prec),allocatable:: f1avg(:,:,:),f2avg(:,:,:), &
                     f3avg(:,:,:)
!                  ! f1avg, f2avg, and f3avg are cell centered forcing 
!                   funcions that are used in inner loop and in turbulence
end module workavg

!2345678***************************************************
! module forcinner contains forcing functions for the inner loop
module forcinner
  use gridlist_variables, only : prec
  
  Implicit None
  real(prec), allocatable::f1(:,:,:),f2(:,:,:),f3(:,:,:)
  real(prec), allocatable::rg_over_prrcp(:,:,:) ! for higrad
  real(prec), allocatable::cp_over_cv(:,:,:) ! for higrad
end module forcinner

!2345678***************************************************
! module metric_variables_old contains arrays for the coordinate transformations
module metric_variables_old
  use gridlist_variables, only : prec
  
  Implicit None
  real(prec),allocatable:: h(:, :,:)
  real(prec),allocatable:: zsio(:,:)

end module metric_variables_old
