!-----------------------------------------------------------------------
! Gridlist variables related to emissions set to defaults here 
! and changed only in the gridlist itself
!-----------------------------------------------------------------------
module emission_gridlist_variables
  use gridlist_variables, only : wp 
  Implicit None
  
  ! Executable Code
  integer :: nEmit = 0 ! Number of resolved emitted gas species 
  integer :: nAero = 0 ! Number of resolved emitted aerosol species
  integer :: nMAero= 2 ! Number of resolved moments for aerosol size distributions
  real(wp):: relativeHumidity = 0._wp     ! Relative humidity if tracking water vapor

end module emission_gridlist_variables

!-----------------------------------------------------------------------
! General emission variables used throughout the emissions module
!-----------------------------------------------------------------------
module emission_general_variables
  use gridlist_variables, only : wp
  
  Implicit None

  ! Executable Code
  character(len=6),allocatable :: spEmit(:)
  real(wp),allocatable :: cpEmit(:)
  real(wp),allocatable :: mwEmit(:)

  character(len=12),allocatable :: spAero(:)
  real(wp),allocatable :: rhoAero(:)

  real(wp) :: specificHumidity

  contains

  !---------------------------------------------------------------------
  ! Lagrangian Interpolation function
  !---------------------------------------------------------------------
  real(wp) function LagrangeInterp(n,p,array)
    Implicit None
    integer,intent(in) :: n
    real(wp),intent(in) :: p
    real(wp),intent(in) :: array(n)
    
    integer :: i,j
    real(wp) :: iw

    LagrangeInterp=0._wp
    if(any(array.le.0._wp)) return

    do i=1,n
      iw=1
      do j=1,n
        if(i.ne.j) iw=iw*(p-j)/(i-j)
      enddo
      LagrangeInterp=LagrangeInterp+log10(array(i))*iw
    enddo
    LagrangeInterp=10._wp**LagrangeInterp
  end function LagrangeInterp

end module emission_general_variables

!-----------------------------------------------------------------------
! Emission variables used for the point source emissions model
!-----------------------------------------------------------------------
module emission_pointSource_variables
  use gridlist_variables, only : wp

  Implicit None

  ! Executable Code
  integer :: nEmitPoints
  real(wp),allocatable :: emitPointLocation(:,:)
  real(wp),allocatable :: emitPointRate(:)
  character(len=6),allocatable :: emitPointSpecies(:)
  
  integer :: nAeroPoints
  real(wp),allocatable :: aeroPointLocation(:,:)
  real(wp),allocatable :: aeroPointRate(:,:)
  character(len=12),allocatable :: aeroPointSpecies(:)

end module emission_pointSource_variables

!-----------------------------------------------------------------------
! Emission variables used for the emission factor emissions model
!-----------------------------------------------------------------------
module emission_factor_variables
  use gridlist_variables, only : wp

  Implicit None

  ! Executable Code
  real(wp),allocatable :: efEmit(:)
  real(wp),allocatable :: efAero(:)

end module emission_factor_variables
