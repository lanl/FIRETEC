!-----------------------------------------------------------------------
! Commonly used functions in the firebrand module
!-----------------------------------------------------------------------
subroutine GetWind(x_p,y_p,z_p,i,j,k,u_c,v_c,w_c,il,iu,jl,ju,klh,lls, &
    u_p,v_p,w_p)
  use gridlist_variables, only : prec,dx,dy
  use zcart_function, only : zcart
  use msga_variables, only : npos,mpos
  use metric_variables, only : zs,z,zedge
  use gridsetup, only : np,mp
  Implicit None

  ! Local Variables
  integer,intent(in) :: i,j,k                       ! Cell indices for point of interest
  integer,intent(in) :: il,iu,jl,ju,klh,lls
  real(prec),intent(in) :: x_p,y_p,z_p              ! Cartesian coordinates for point of interest
  real(prec),intent(in) :: u_c(il:iu,jl:ju,klh:lls) ! Cell velocities for before, self, and after cell of interest
  real(prec),intent(in) :: v_c(il:iu,jl:ju,klh:lls) ! Cell velocities for before, self, and after cell of interest
  real(prec),intent(in) :: w_c(il:iu,jl:ju,klh:lls) ! Cell velocities for before, self, and after cell of interest
  real(prec),intent(inout) :: u_p,v_p,w_p           ! Output point velocities 

  integer :: im,jm,km,ip,jp,kp
  real(prec) :: x_c,y_c,z_c  ! Cartesian coordinates of cell center
  real(prec) :: xd,yd,zd
  real(prec) :: u00,u01,u10,u11,u0,u1
  real(prec) :: v00,v01,v10,v11,v0,v1
  real(prec) :: w00,w01,w10,w11,w0,w1

  ! Executable Code
  x_c=((npos-1)*np+i-0.5)*dx
  y_c=((mpos-1)*mp+j-0.5)*dy
  z_c=zcart(z(k),i,j)

  ! Find quadrant of cell containing point of interest
  if(x_p.gt.x_c)then
    ip=i+1
    im=i
    xd=(x_p-x_c)/dx
  else
    ip=i
    im=i-1
    xd=(x_c-x_p)/dx
  endif
  if(y_p.gt.y_c)then
    jp=j+1
    jm=j
    yd=(y_p-y_c)/dy
  else
    jp=j
    jm=j-1
    yd=(y_c-y_p)/dy
  endif
  if(z_p.gt.z_c)then
    km=k
    if(k.ne.lls)then
      kp=k+1
      zd=(z_p-z_c)/(zcart(z(kp),i,j)-z_c)
    else
      kp=lls
      zd=(z_p-z_c)/(zcart(zedge(lls),i,j)-z_c)
    endif
  else
    kp=k
    km=k-1
    if(k.ne.1)then
      zd=(z_c-z_p)/(z_c-zcart(z(km),i,j))
    else
      zd=(z_c-z_p)/(z_c-zs(i,j))
    endif
  endif
  
  ! Triple interpolation (en.wikipedia.org/wiki/Trilinear_interpolation)
  u00=u_c(im,jm,km)*(1-xd)+u_c(ip,jm,km)
  u01=u_c(im,jm,kp)*(1-xd)+u_c(ip,jm,kp)
  u10=u_c(im,jp,km)*(1-xd)+u_c(ip,jp,km)
  u11=u_c(im,jp,kp)*(1-xd)+u_c(ip,jp,kp)
  u0 =u00*(1-yd)+u10*yd
  u1 =u01*(1-yd)+u11*yd
  u_p=u0*(1-zd)+u1*zd

  v00=v_c(im,jm,km)*(1-xd)+v_c(ip,jm,km)
  v01=v_c(im,jm,kp)*(1-xd)+v_c(ip,jm,kp)
  v10=v_c(im,jp,km)*(1-xd)+v_c(ip,jp,km)
  v11=v_c(im,jp,kp)*(1-xd)+v_c(ip,jp,kp)
  v0 =v00*(1-yd)+v10*yd
  v1 =v01*(1-yd)+v11*yd
  v_p=v0*(1-zd)+v1*zd
  
  w00=w_c(im,jm,km)*(1-xd)+w_c(ip,jm,km)
  w01=w_c(im,jm,kp)*(1-xd)+w_c(ip,jm,kp)
  w10=w_c(im,jp,km)*(1-xd)+w_c(ip,jp,km)
  w11=w_c(im,jp,kp)*(1-xd)+w_c(ip,jp,kp)
  w0 =w00*(1-yd)+w10*yd
  w1 =w01*(1-yd)+w11*yd
  w_p=w0*(1-zd)+w1*zd

end subroutine GetWind
