!-----------------------------------------------------------------------
! Turbulence gridlist variables
!-----------------------------------------------------------------------
module turb_gridlist_variables
  use gridlist_variables, only : prec

  Implicit None

  integer     :: isa              = 2       ! Calculation of sa
  real(prec)  :: rturbprandtl     = 2.      ! Inverse of turb prandtl number for scalar diffusion

end module turb_gridlist_variables

!-----------------------------------------------------------------------
! Turbulence arrays and associated variables
!-----------------------------------------------------------------------
module linn_turb_variables
  use gridlist_variables, only : prec

  Implicit None

  real(prec) :: cd=1.0
  real(prec) :: sc=0.1
  real(prec),allocatable :: sa(:,:,:),sb(:,:,:)
  real(prec),allocatable :: saxy(:,:,:),saz(:,:,:)

  real(prec),allocatable :: zonehts(:),zonedzs(:),zonerhos(:)
  real(prec),allocatable :: iftwght(:),fuelinds(:),zoneus(:),zonevs(:)
  real(prec),allocatable :: deltaus(:),deltavs(:)

  real(prec),allocatable :: sqrtG_Kxy(:,:,:),sqrtG_Kz(:,:,:)
  real(prec),allocatable :: sqrtG_KG33(:,:,:)
  real(prec),allocatable :: K_axy(:,:,:),K_az(:,:,:),K_b(:,:,:)
  real(prec),allocatable :: rtke_abc(:,:,:)

end module linn_turb_variables
