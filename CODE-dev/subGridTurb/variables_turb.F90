!-----------------------------------------------------------------------
! Turbulence gridlist variables
!-----------------------------------------------------------------------
module turb_gridlist_variables
  use gridlist_variables, only : wp

  Implicit None

  integer   :: isa              = 2       ! Calculation of sa
  real(wp)  :: rturbprandtl     = 2._wp   ! Inverse of turb prandtl number for scalar diffusion

end module turb_gridlist_variables

!-----------------------------------------------------------------------
! Turbulence arrays and associated variables
!-----------------------------------------------------------------------
module linn_turb_variables
  use gridlist_variables, only : wp

  Implicit None

  real(wp) :: cd=1.0_wp
  real(wp) :: sc=0.1_wp
  real(wp),allocatable :: sa(:,:,:),sb(:,:,:)
  real(wp),allocatable :: saxy(:,:,:),saz(:,:,:)

  real(wp),allocatable :: zonehts(:),zonedzs(:),zonerhos(:)
  real(wp),allocatable :: iftwght(:),fuelinds(:),zoneus(:),zonevs(:)
  real(wp),allocatable :: deltaus(:),deltavs(:)

  real(wp),allocatable :: sqrtG_Kxy(:,:,:),sqrtG_Kz(:,:,:)
  real(wp),allocatable :: sqrtG_KG33(:,:,:)
  real(wp),allocatable :: K_axy(:,:,:),K_az(:,:,:),K_b(:,:,:)
  real(wp),allocatable :: rtke_abc(:,:,:)

end module linn_turb_variables
