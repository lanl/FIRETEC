!-----------------------------------------------------------------------
! Gridlist variables related to emissions set to defaults here 
! and changed only in the gridlist itself
!-----------------------------------------------------------------------
module firebrand_gridlist_variables
  use gridlist_variables, only : wp
  Implicit None

  integer :: shapeFB = 1  ! Flag to determine characteristic shape of firebrands (1=disk, 2=cylinder, 3=sphere)
  integer :: dsizeFB = 1  ! Flag to determine how size changes as burning commences (1=radius, 2=height)
  real(wp) :: radius = 0.001_wp  ! initial radius of firebrands (m)
  real(wp) :: height = 0.001_wp  ! initial height of firebrands (m)
  real(wp) :: rhoFB  = 100.0_wp    ! initial density of firebrands (kg/m3)
  real(wp) :: TempHot= 500.0_wp   ! temperature threshold for launching firebrands (K)
  character(len=257) :: outname_fb     = 'brand.out'   ! Name of firebrand output files
  character(len=257) :: restartfile_fb = '' ! Name of restart file

end module firebrand_gridlist_variables

!-----------------------------------------------------------------------
! General emission variables used throughout the firebrand module
!-----------------------------------------------------------------------
module firebrand_general_variables
  use gridlist_variables, only : wp
  Implicit None

  integer :: nvbrand = 12
    ! Deterministic characteristics of firebrands
    ! 1,2,3 = x,y,z Cartisian location
    ! 4,5,6 = u,v,w firebrand velocities
    ! 7,8   = characteristic sizes (radius, height)
    ! 9,10  = temperature and density of firebrand
    ! 11,12 = angle of orientation (radians from x, radians from z)
  integer :: ix,iy,iz,iuvel,ivvel,iwvel
  integer :: irad,ihght,itemp,irho,itheta,iphi
  
  integer :: numBrands = 0
  real(wp),allocatable :: xvbrand(:,:)
  character(len=12),allocatable :: xvbrand_list(:)

  real(wp) :: Cd_dsk = 1.1_wp
  real(wp) :: Cd_cyl = 1.2_wp
  real(wp) :: Cd_sph = 0.4_wp
  
end module firebrand_general_variables
