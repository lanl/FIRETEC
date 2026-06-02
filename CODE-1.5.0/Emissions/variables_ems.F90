!-----------------------------------------------------------------------
! Gridlist variables related to emissions set to defaults here 
! and changed only in the gridlist itself
!-----------------------------------------------------------------------
module emission_gridlist_variables
   
  Implicit None
  
  ! Executable Code
  integer :: nEmit = 0 ! Number of resolved emitted gas species 
  integer :: nAero = 0 ! Number of resolved emitted aerosol species
  integer :: nMAero= 2 ! Number of resolved moments for aerosol size distributions

end module emission_gridlist_variables

!-----------------------------------------------------------------------
! General emission variables used throughout the emissions module
!-----------------------------------------------------------------------
module emission_general_variables
  use gridlist_variables, only : prec
  
  Implicit None

  ! Executable Code
  character(len=6),allocatable :: spEmit(:)
  real(prec),allocatable :: cpEmit(:)
  real(prec),allocatable :: mwEmit(:)

  character(len=12),allocatable :: spAero(:)
  real(prec),allocatable :: rhoAero(:)

  contains

  ! Lagrangian Interpolation function
  real(prec) function LagrangeInterp(n,p,array)
    Implicit None
    integer,intent(in) :: n
    real(prec),intent(in) :: p
    real(prec),intent(in) :: array(n)
    
    integer :: i,j
    real(prec) :: iw

    LagrangeInterp=0.
    do i=1,n
      iw=1
      do j=1,n
        if(i.ne.j) iw=iw*(p-j)/(i-j)
      enddo
      LagrangeInterp=log10(array(i))*iw
    enddo
    LagrangeInterp=10**LagrangeInterp
  end function LagrangeInterp

end module emission_general_variables

!-----------------------------------------------------------------------
! Emission variables used for the point source emissions model
!-----------------------------------------------------------------------
module emission_pointSource_variables
  use gridlist_variables, only : prec

  Implicit None

  ! Executable Code
  integer :: nEmitPoints
  real(prec),allocatable :: emitPointLocation(:,:)
  real(prec),allocatable :: emitPointRate(:)
  character(len=6),allocatable :: emitPointSpecies(:)
  
  integer :: nAeroPoints
  real(prec),allocatable :: aeroPointLocation(:,:)
  real(prec),allocatable :: aeroPointRate(:,:)
  character(len=12),allocatable :: aeroPointSpecies(:)

end module emission_pointSource_variables

!-----------------------------------------------------------------------
! Emission variables used for the emission factor emissions model
!-----------------------------------------------------------------------
module emission_factor_variables
  use gridlist_variables, only : prec

  Implicit None

  ! Executable Code
  real(prec),allocatable :: efEmit(:)
  real(prec),allocatable :: efAero(:)

end module emission_factor_variables
