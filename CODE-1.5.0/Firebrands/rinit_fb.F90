!-----------------------------------------------------------------------
! rinitFirebrands initializes all arrays and variables related to 
! firebrand evolution throughout the simulation
!-----------------------------------------------------------------------
subroutine rinitFirebrands(irst)
  use firebrand_general_variables, only : ix,iy,iz,iuvel,ivvel,iwvel, &
    irad,ihght,itemp,irho,itheta,iphi,nvbrand,xvbrand_list,xvbrand, &
    numBrands
  Implicit None

  ! Local Variables
  integer,intent(in) :: irst

  ! Executable Code
  ! Set of pointers for xvbrand
  ix=1
  iy=2
  iz=3
  iuvel=4
  ivvel=5
  iwvel=6
  irad=7
  ihght=8
  itemp=9
  irho=10
  itheta=11
  iphi=12

  allocate(xvbrand_list(nvbrand))
  xvbrand_list(ix) = 'x_fb'
  xvbrand_list(iy) = 'y_fb'
  xvbrand_list(iz) = 'z_fb'
  xvbrand_list(iuvel) = 'u_fb'
  xvbrand_list(ivvel) = 'v_fb'
  xvbrand_list(iwvel) = 'w_fb'
  xvbrand_list(irad)  = 'radius_fb'
  xvbrand_list(ihght) = 'height_fb'
  xvbrand_list(itemp) = 'temp_fb'
  xvbrand_list(irho)  = 'rho_fb'
  xvbrand_list(itheta)= 'theta_fb'
  xvbrand_list(iphi)  = 'phi_fb'

  ! TODO Restart
  if(irst.eq.1) then
    call irstreadio_fb
  else
    allocate(xvbrand(numBrands,nvbrand))
  endif

end subroutine rinitFirebrands
