!-----------------------------------------------------------------------
! Gridlist variables related to heatsource set to defaults here 
! and changed only in the gridlist itself
!-----------------------------------------------------------------------
module heatsource_gridlist_variables
  use gridlist_variables, only : prec
   
  Implicit None
  
  ! Executable Code
  integer :: isource_function  !0 const; 1 step funct; 2 sin(pi*omega*t)
  integer :: nsource_hs
  real(prec),allocatable :: av_hs(:)      ! m2/m3 surf area to source volume 
  real(prec),allocatable :: flux_hs(:)    ! W/m2
  real(prec),allocatable :: depth_hs(:)   ! m
  real(prec),allocatable :: freq_hs(:)    ! 1/s

  ! Circle heat source
  real(prec),allocatable :: radius_hs(:)  ! m
  real(prec),allocatable :: xcen_hs(:)    ! m
  real(prec),allocatable :: ycen_hs(:)    ! m

  ! Rectangluar heat source
  real(prec),allocatable :: xl_hs(:)    ! m
  real(prec),allocatable :: xu_hs(:)    ! m
  real(prec),allocatable :: yl_hs(:)    ! m
  real(prec),allocatable :: yu_hs(:)    ! m

  ! mass source
  integer :: imass_source         ! Source gas-phase mass
  integer :: MLRHRRinputs = 0        ! If 1, provide tab-delimited file titled "MLRHRRinputs.dat.#" for each heat source
  real(prec),allocatable :: fm_hs(:)    ! kg/m3 s
  real(prec),allocatable :: TME(:,:)   ! timestamp for inputs
  real(prec),allocatable :: MLR(:,:)   ! mass loss rate (must be at least 2 values but can be the same)
  real(prec),allocatable :: HRR(:,:)   ! heat release rate (must be at least 2 values but can be the same)
  character(len=257) :: inputName = 'MLRHRRinputs.dat'  ! Name of sensor data files

end module heatsource_gridlist_variables

!-----------------------------------------------------------------------
! General emission variables used throughout the emissions module
!-----------------------------------------------------------------------
module heatsource_general_variables
  
  Implicit None

  ! Executable Code

end module heatsource_general_variables
