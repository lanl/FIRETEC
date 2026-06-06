!----------------------------------------------------------------
! Gridlist variables set to defaults here and changed only in the 
! gridlist itself
!----------------------------------------------------------------
module gridlist_variables
  
  Implicit None
  
  integer :: irst             = 0       ! Restart flag
  integer :: nt               = 10000   ! Total number of large timesteps
  integer :: nts              = 10      ! Number of small timesteps (higrad) in a large timestep
  real    :: dts              = 0.001   ! Temporal size of a small timestep (higrad)
  integer :: ntp              = 10      ! Number of small timesteps (firetec) in a large timestep
  integer :: n                = 200     ! Number of cells in the x dimension
  integer :: m                = 200     ! Number of cells in the y dimension
  integer :: l                = 40      ! Number of cells in the z dimension
  real    :: dx               = 2.      ! Size of cell in x dimension
  real    :: dy               = 2.      ! Size of cell in y dimension
  real    :: dz               = 15.     ! Size of cell in z dimension before stretching
  real    :: aa1              = 0.1     ! Stretching parameter for z direction
  integer :: nprocx           = 8       ! Number of processors in the x dimension
  integer :: nprocy           = 8       ! Number of processors in the y dimension
  integer :: ih               = 3       ! Number of halo cells
  integer :: nr               = 5       ! Number of cells for lateral relaxation (irlx, irly)
  integer :: ibcx             = 1       ! X/Y boundary conditions (0=cyclic,1<= different relaxation options)
  integer :: ibcy             = 1       ! X/Y boundary conditions (0=cyclic,1<= different relaxation options)
  integer :: ibclatopen       = 0       ! Top boundary relaxation options
  integer :: ibctopopen       = 0       ! Open outlet, feeslip option at top
  integer :: iab              = 1       ! Options for top absorption initialization
  real    :: zab              = 400.    ! Bottom of top damper (m)
  real    :: zabt             = 400.    ! Bottom of top damper for iab=3 option
  real    :: tow              = 10.     ! Magnitude parameter for top relaxation
  integer :: itheta           = 0       ! Stable layer at domain top
  integer :: ifire            = 1       ! Firetec capability
  integer :: iturb            = 2       ! Turbulence resolution options
  integer :: isa              = 2       ! Calculation of sa
  real    :: rturbprandtl     = 2.      ! Inverse of turb prandtl number for scalar diffusion
  integer :: irad             = 2       ! Radiation options
  integer :: icallrad         = 5       ! Frequency of radiation call
  real    :: crad             = 50.     ! Soot radiation constant
  integer :: irandseed        = 0       ! Random function seed
  integer :: iseed            = 0       ! Seed for randomization function
  integer :: isootmodel       = 0       ! Options for radiating soot model
  integer :: iradeastflux     = 0       ! Specific radiative flux output through
  integer :: inonlocal        = 0       ! Option for non-local burning
  integer :: iemissions       = 0       ! Option for tracking particle emissions
  integer :: idiffsies        = 0       ! Option for diffusing ignitions for creeping fires
  integer :: irhovapor        = 0       ! Option for tracking water vapor
  integer :: ifbrand          = 0       ! Option for tracking firebrands
  integer :: icorio           = 0       ! Option fors coriolis effect
  integer :: ilspgf           = 0       ! Large scale pressure gradient forcing options
  integer :: izlspgf          = 1       ! Reference for max flux is whole domain or a slice at height zu (1)
  integer :: frqlspgf         = 10000   ! Large scale pressure gradient forcing frequency
  integer :: islip            = 0       ! Option for no-slip boundary condition forcing
  real    :: tambient         = 300.    ! Ambient temperature (K)
  real    :: relativehumidity = 0.      ! Relative humidity if tracking water vapor
  real    :: pressground      = 1.e5    ! Reference pressure for ground (Pa)
  real    :: zgroundref       = 0.      ! Reference elevation for ground (m)
  integer :: iperturb         = 0       ! Perturbation option for cyclic runs
  real    :: u0               = 3.      ! Initial ambient x direction wind magnitude (m/s)
  real    :: uramp            = 10.     ! Target of ambient x direction wind magnitude (m/s)
  real    :: uramptime        = 30.     ! Time to reach target ambient x direction wind (s)
  integer :: uswitch          = 2       ! Options for x direction ambient wind profiles
  real    :: v0               = 0.      ! Initial ambient y direction wind magnitude (m/s)
  real    :: vramp            = 0.      ! Target of ambient y direction wind magnitude (m/s)
  real    :: vramptime        = 0.      ! Time to reach target ambient y direction wind (s)
  integer :: vswitch          = 0       ! Options for y direction ambient wind profiles
  real    :: zu               = 25.     ! Height of referenced ambient winds
  integer :: ius              = 0       ! Lower x indice for reference zone in lai computation (when 0 whole domain is used) for uswitch=2
  integer :: iue              = 0       ! Upper x indice for reference zone in lai computation (when 0 whole domain is used) for uswitch=2
  integer :: jus              = 0       ! Lower y indice for reference zone in lai computation (when 0 whole domain is used) for vswitch=2
  integer :: jue              = 0       ! Upper y indice for reference zone in lai computation (when 0 whole domain is used) for vswitch=2
  real    :: slopeangle       = 0.      ! Archaic slope usage (use topo now) for rotating gravity
  real    :: slopeazimuth     = 0.      ! Archaic slope usage (use topo now) for rotating gravity
  integer :: iord             = 2       ! Advection order of resolution
  integer :: nonos            = 1       ! Advection for non-oscillatory flow option
  integer :: idiv             = 1       ! Divergent flows option
  integer :: nfct             = 1
  integer :: nonosold         = 1       ! Non-oscillatory option in mpdataold
  integer :: ifuel            = 1       ! Option for burning fuel pdf choice
  integer :: nfuel            = 1       ! Number of fuel types for multiple-fuels
  integer :: ivegread         = 0       ! Option for fuel bed read-in
  real    :: rhoMicro         = 500.    ! Fuel's true density (kg/m3)
  real    :: cpwood           = 2500.   ! Fuel's heat capacity (J/kg K)
  real    :: Water2WoodRatio  = 0.56    ! Water product to wood mass ratio
  integer :: frqoutput        = 100     ! Frequency (large time step) of output
  integer :: frqfilstr        = 0       ! Frequency (large time step) of file filtering
  integer :: ipotflow         = 0       ! Potflow testing scenario
  integer :: ignVertExtent    = 0       ! Top cell of vertical ignition
  integer :: icfmeflag        = 0       ! Option to read-in slices of the grid
  integer :: iwindfieldout    = 0       ! Option for outputting a windfield
  integer :: iwindfieldin     = 0       ! Option for reading-in a windfield
  integer :: windspeedupfactor= 1       ! Factor of speed up between wind run and fire run
  integer :: itwindfield      = 80000   ! Timestep to start saving wind data
  integer :: itinterp         = 10      ! # of timesteps between windfield interpolation
  integer :: ibcells          = 5       ! # of x cells saved on boundary
  integer :: jbcells          = 5       ! # of y cells saved on boundary
  integer :: is               = 0       ! Lower x indice for production xvdata file on subdomain (0 = whole domain is used)
  integer :: ie               = 0       ! Upper x indice for production xvdata file on subdomain (0 = whole domain is used)
  integer :: js               = 0       ! Lower y indice for production xvdata file on subdomain (0 = whole domain is used)
  integer :: je               = 0       ! Upper y indice for production xvdata file on subdomain (0 = whole domain is used)

  character(len=60)  :: ignfile           = 'ignite.dat'    ! Name of ignition file
  character(len=257) :: outname           = 'comp.out'      ! Name of output files
  character(len=257) :: restartfile       = 'comp.out.100'  ! Name of restart file
  character(len=60)  :: topofile          = ''              ! Name of topo file
  character(len=257) :: windfieldstartfile= 'windfieldstart'! Name of windfield start file
  character(len=257) :: xvdataname        = 'xvdata'        ! Name of xvdata files

end module gridlist_variables  

!----------------------------------------------------------------
! MSGA arrays and associated variables
!----------------------------------------------------------------
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

!----------------------------------------------------------------
! Forcing arrays from physics and boundary forcings
!----------------------------------------------------------------
module forcings

  Implicit None
      
  real,allocatable :: force(:,:,:,:)
  real,allocatable :: tkewght(:)
  real,allocatable :: uavg(:,:,:),vavg(:,:,:),oavg(:,:,:)
  real,allocatable :: wt(:),wtf(:)
  real,allocatable :: u1(:,:,:),u2(:,:,:),u3(:,:,:) ! face contravariant velocities*dt/di/gi

end module forcings

!----------------------------------------------------------------
! Xv arrays and associated variables/pointers
!----------------------------------------------------------------
module xvall

  Implicit None
      
  real,allocatable :: xv(:,:,:,:)       ! transported field
  real,allocatable :: xvrho(:,:,:,:)    ! xv/rho field
  real,allocatable :: xe(:,:,:,:)       ! environmental field
  real,allocatable :: xvtmp(:,:,:,:)    ! MOA xv field
  real,allocatable :: relaxxv(:,:,:,:)  ! relaxation coeff for bc
  real,allocatable :: xvdataold(:,:,:,:),xvdatanew(:,:,:,:) ! Data arrays from windrun
  real,allocatable :: uprofile(:,:,:) ! uprofile from LAI
  integer :: ixevariation = 0 ! flag that environmental field will vary with time
  character(len=12),allocatable :: xv_list(:)

  ! Following are the pointers to the xvarray making it more dynamic
  integer :: iuvel=1,ivvel=1,iwvel=1
  integer :: itemp=1,ika=1,ikb=1,iO2=1,irho=1
  integer :: ivapor=1
  integer :: iM0=1,iM1=1
  integer :: imixfrac=1
  integer :: nv

end module xvall

!----------------------------------------------------------------
! Thermodynamic arrays and associated variables
!----------------------------------------------------------------
module thermo_variables

  Implicit None
  
  real :: cpwater=4186. ! J/kg*K
  real :: cp_gas,cp_over_cv_gas,rg_over_cp_gas,rg_over_prrcp_gas
  real :: specifichumidity
  real,allocatable :: pr(:,:,:),pre(:,:,:)
  real,allocatable :: tempg(:,:,:)

end module thermo_variables

!----------------------------------------------------------------
! Fuel arrays and associated variables
!----------------------------------------------------------------
module fuel_variables

  Implicit None

  integer :: lfuel 
  real :: min_rhoFuel=1.e-6 
  real :: rnfuel=0.4552 ! stoichiometric fuel constant (Drysdale pp. 179)
  real :: rno=0.5448    ! stoichiometric oxygen constant (Drysdale pp. 179)
  real :: tcrit=600.    ! temperature of pyrolysis
  real :: twvap=373.    ! temperature of water vaporization
  real :: tfstep=310.   ! tail of the psi ramp curve
  real :: tstep
  real,allocatable :: rhoFuel(:,:,:,:),rhowater(:,:,:,:)
  real,allocatable :: sizeScale(:,:,:,:),actualFuelDepth(:,:,:)
  real,allocatable :: sies(:,:,:,:),psiwmax(:,:,:,:)
  real,allocatable :: temps(:,:,:,:)
  real,allocatable :: convht(:,:,:,:)
  real,allocatable :: fcorr(:,:)
  real,allocatable :: rhoFuelInitial(:,:,:,:)
  real,allocatable :: siesDiff(:,:,:,:)

end module fuel_variables

!----------------------------------------------------------------
! Metric arrays and associated variables
!----------------------------------------------------------------
module metric_variables

  Implicit None
 
  real :: aa2,aa3 
  real,allocatable :: x(:),y(:),z(:),zedge(:)
  real,allocatable :: zs(:,:)
  real,allocatable :: gmul(:),c13(:,:),c23(:,:),gi(:,:,:)

  real :: zb

end module metric_variables

!----------------------------------------------------------------
! Large pressure gradient arrays and associated variables
!----------------------------------------------------------------
module lspgf_variables

  Implicit None
  
  real, allocatable:: massFluxTemp(:) ! mass flux in the direction of geostrophic wind (ug,vg), temporary on local domain
  real, allocatable:: massFlux(:) ! mass flux in the direction of geostrophic wind (ug,vg), after mpi_allreduce
  real :: totMassFlux, targMassFlux

  real :: cosg,sing ! geostrophic wind cos and sin
  integer :: iwindx 
  ! if iwindx==1 mass flux is computed for each j (as if wind aligned with x axis) : n*vg<m*ug
  ! if iwindx==0 mass flux is computed for each i (as if wind aligned with y axis) : n*vg>m*ug
  integer :: nMassFlux ! size of massFluxArray (either n or m depending of iwindx)
  real :: intsintheta ! vertical integral of sintheta
  real,allocatable:: sintheta(:,:,:) ! evaluation of v/sqrt(u2+v2)
  real,allocatable:: sinthetaf(:,:,:) ! evaluation of f*v/sqrt(u2+v2)
  real,allocatable:: flspgf(:,:) ! evaluation of f
  ! based on ekman theory to impose a pressure gradient (ilspgf)

end module lspgf_variables

!----------------------------------------------------------------
! Gridsetup variables include parallelization variables and 
! general temporal and spatial grid variables.
!----------------------------------------------------------------
module gridsetup

  Implicit None

  integer :: timestep=0 ! Temporal counter in iterations
  real :: time=0 ! Temporal counter in seconds
  real :: dt,dtp
  integer :: np,mp
  integer :: itrestart=0 ! Iteration of restart
  integer :: ittot=0
  real :: dxi,dyi,dzi,dti

end module gridsetup

!----------------------------------------------------------------
! Constant variables do not change throughout the simulation.
!----------------------------------------------------------------
module constants

  Implicit None

  real :: pi    = acos(-1.) ! pi, yuuummmmm!!!
  real :: g     = 9.81      ! gravity acceleration (m/s2) 
  real :: Rgas  = 8314.5    ! ideal gas constant (J/kmol K)
  real :: sigma = 5.678e-8  ! Stefan-Boltzmann constant               (W/(m^2K^4))
  real :: fcor2 = 1.4584e-4*cos(acos(-1.)+55./180.) ! I don't know what this is
  real :: fcor3 = 1.4584e-4*sin(acos(-1.)+55./180.) ! I don't know what this is

end module constants

!----------------------------------------------------------------
! Radiation arrays and associated variables
!----------------------------------------------------------------
module radiation_variables

  Implicit None
  
  integer :: lmc
  integer :: nphotbatch=1e6
  real :: Emin=1000.    ! Minimal energy for a photon
  real :: Etot
  real,allocatable :: sourcesol(:,:,:,:),sourcegas(:,:,:)
  real,allocatable :: papvift(:,:,:,:),papvgas(:,:,:)
  real,allocatable :: papvtotAZ(:,:,:),volumeAZ(:,:,:)
  real,allocatable :: zCartEdgeAZGlobal(:,:,:,:)
  real,allocatable :: EAZ(:,:,:)
  real,allocatable :: eastFlux(:,:,:)
  integer,allocatable :: kfindex(:,:,:)
  integer,allocatable :: nphotemisAZsplit(:,:,:)
  integer,allocatable :: nphotonAZgather(:,:,:),nphotonAZ(:,:,:)
  integer,allocatable :: nphotonEastFlux(:,:,:)

  real,allocatable :: firad (:,:,:),fsiesrad(:,:,:,:)

end module radiation_variables

!----------------------------------------------------------------
! Turbulence arrays and associated variables
!----------------------------------------------------------------
module turb_variables

  Implicit None

  real :: cd=1.0
  real :: sc=0.1
  real,allocatable :: sa(:,:,:),sb(:,:,:)
  real,allocatable :: saxy(:,:,:),saz(:,:,:)

  real,allocatable :: zonehts(:),zonedzs(:),zonerhos(:)
  real,allocatable :: iftwght(:),fuelinds(:),zoneus(:),zonevs(:)
  real,allocatable :: deltaus(:),deltavs(:)

  real,allocatable :: sqrtG_Kxy(:,:,:),sqrtG_Kz(:,:,:)
  real,allocatable :: sqrtG_KG33(:,:,:)
  real,allocatable :: K_axy(:,:,:),K_az(:,:,:),K_b(:,:,:)
  real,allocatable :: rtke_abc(:,:,:)

end module turb_variables

!----------------------------------------------------------------
! Ignition arrays and associated variables
!----------------------------------------------------------------
module ign_variables

  Implicit None

  integer :: nIgn  
  real :: targetTemp,rampRate
  integer,allocatable :: ignLocation(:,:)
  real,allocatable :: ignTime(:)

end module ign_variables

!----------------------------------------------------------------
! namefile subroutine derives names for specific files from macro
! data and iteration number
!----------------------------------------------------------------
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
  Implicit None
  real,allocatable::  wavg(:,:,:)
!                   ! uavg, vavg, wavg and oavg are cell centered average
!                     velocities computed in the inner loop   (m/s)
  real,allocatable:: f1avg(:,:,:),f2avg(:,:,:), &
                     f3avg(:,:,:)
!                  ! f1avg, f2avg, and f3avg are cell centered forcing 
!                   funcions that are used in inner loop and in turbulence
end module workavg

!2345678***************************************************
! module forcinner contains forcing functions for the inner loop
module forcinner
  Implicit None
  real, allocatable::f1(:,:,:),f2(:,:,:),f3(:,:,:)
  real, allocatable::rg_over_prrcp(:,:,:) ! for higrad
  real, allocatable::cp_over_cv(:,:,:) ! for higrad
end module forcinner

!2345678***************************************************
! module metric_variables_old contains arrays for the coordinate transformations
module metric_variables_old
  Implicit None
  real,allocatable:: h(:, :,:)
  real,allocatable:: zsio(:,:)

end module metric_variables_old

!2345678***************************************************
module turba
  Implicit None
  !stresrij arrays:
  real,allocatable::AB11c(:,:,:),AB12c(:,:,:),AB13c(:,:,:), &
    AB21c(:,:,:),AB22c(:,:,:),AB23c(:,:,:), &
    AB31c(:,:,:),AB32c(:,:,:),AB33c(:,:,:)
  real, allocatable::Ci(:,:,:),Cj(:,:,:),Ck(:,:,:)

  ! source terms
  real,allocatable:: rkc(:, :,:) 
  integer::kfuelmax ! max height of the fuel
  integer:: isoturb
  integer:: ilapdo
  real::diffcst ! diffusion cst (=0.09)
end module turba

!2345678***************************************************
module xvbin    !JMC  10/21/6  This module is for reading windfields
  Implicit None
  integer :: itabsold,itabsnew
  real,allocatable:: xvbdataold(:,:,:,:)
  real,allocatable:: xvbdatanew(:,:,:,:)
  character(len=257)::fxvbdataname=''        !JMC
   
  contains
  !---------------------------------------------------
  real function lin_interp(xi,xl,xh,yl,yh)
    Implicit None
    real :: xi,xl,xh,yl,yh

    lin_interp=(xi-xl)/real(xh-xl)*(yh-yl)+yl

  end function lin_interp
  !---------------------------------------------------
end module xvbin
