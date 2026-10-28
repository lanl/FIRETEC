!-----------------------------------------------------------------------
! Gridlist variables related to heatsource set to defaults here 
! and changed only in the gridlist itself
!-----------------------------------------------------------------------
module heatsource_gridlist_variables
  use gridlist_variables, only : wp
   
  Implicit None
  
  ! Executable Code
  integer :: nsource_hs
  integer,allocatable :: isource_function_hs(:)  !0 const; 1 step funct; 2 sin(pi*omega*t)
  real(wp),allocatable :: starttime_hs(:)      ! Time to begin each heat source
  real(wp),allocatable :: endtime_hs(:)        ! Time to end each heat source
  real(wp),allocatable :: av_hs(:)      ! m2/m3 surf area to source volume 
  real(wp),allocatable :: flux_hs(:)    ! W/m2
  real(wp),allocatable :: depth_hs(:)   ! m
  real(wp),allocatable :: freq_hs(:)    ! 1/s

  ! Circle heat source
  real(wp),allocatable :: radius_hs(:)  ! m
  real(wp),allocatable :: xcen_hs(:)    ! m
  real(wp),allocatable :: ycen_hs(:)    ! m

  ! Rectangluar heat source
  real(wp),allocatable :: xl_hs(:)    ! m
  real(wp),allocatable :: xu_hs(:)    ! m
  real(wp),allocatable :: yl_hs(:)    ! m
  real(wp),allocatable :: yu_hs(:)    ! m

  ! mass source
  integer,allocatable :: imass_source_hs(:)    ! Source gas-phase mass option per heat source
  integer,allocatable :: MLRHRRinputs_hs(:)    ! If 1, provide MLRHRRinputs.dat.# for that heat source
  real(wp),allocatable :: fm_hs(:)    ! kg/m3 s
  real(wp),allocatable :: TME(:,:)   ! timestamp for inputs
  real(wp),allocatable :: MLR(:,:)   ! mass loss rate (must be at least 2 values but can be the same)
  real(wp),allocatable :: HRR(:,:)   ! heat release rate (must be at least 2 values but can be the same)
  logical,allocatable :: MLRHRR_eof_hs(:) ! true after final MLR/HRR record is reached
  character(len=257) :: inputName = 'MLRHRRinputs.dat'  ! Name of sensor data files

end module heatsource_gridlist_variables

!-----------------------------------------------------------------------
! General emission variables used throughout the emissions module
!-----------------------------------------------------------------------
module heatsource_general_variables
  
  Implicit None

  ! Executable Code

end module heatsource_general_variables
