!-----------------------------------------------------------------------
! Commonly used functions in the firebrand module
!-----------------------------------------------------------------------
subroutine GetWind(x_p,y_p,z_p,i,j,k,u_c,v_c,w_c,il,iu,jl,ju,klh,lls, &
    u_p,v_p,w_p)
  use gridlist_variables, only : wp,dx,dy
  use msga_variables, only : npos,mpos
  use metric_variables, only : zs,z,zedge
  use gridsetup, only : np,mp
  Implicit None

  ! Local Variables
  integer,intent(in) :: i,j,k                       ! Cell indices for point of interest
  integer,intent(in) :: il,iu,jl,ju,klh,lls
  real(wp),intent(in) :: x_p,y_p,z_p              ! Cartesian coordinates for point of interest
  real(wp),intent(in) :: u_c(il:iu,jl:ju,klh:lls) ! Cell velocities for before, self, and after cell of interest
  real(wp),intent(in) :: v_c(il:iu,jl:ju,klh:lls) ! Cell velocities for before, self, and after cell of interest
  real(wp),intent(in) :: w_c(il:iu,jl:ju,klh:lls) ! Cell velocities for before, self, and after cell of interest
  real(wp),intent(inout) :: u_p,v_p,w_p           ! Output point velocities 

  integer :: im,jm,km,ip,jp,kp
  real(wp) :: x_c,y_c,z_c  ! Cartesian coordinates of cell center
  real(wp) :: wxm,wxp,wym,wyp,wzm,wzp
  real(wp) :: u00,u01,u10,u11,u0,u1
  real(wp) :: v00,v01,v10,v11,v0,v1
  real(wp) :: w00,w01,w10,w11,w0,w1
  real(wp),external :: zcart

  ! Executable Code
  x_c=((npos-1)*np+i-0.5_wp)*dx
  y_c=((mpos-1)*mp+j-0.5_wp)*dy
  z_c=zcart(z(k),i,j)

  ! Find octant of cell containing point of interest
  if(x_p.gt.x_c)then
    ip=i+1
    im=i
    wxp=(x_p-x_c)/dx
    wxm=1._wp-wxp
  else
    ip=i
    im=i-1
    wxm=(x_c-x_p)/dx
    wxp=1._wp-wxm
  endif
  if(y_p.gt.y_c)then
    jp=j+1
    jm=j
    wyp=(y_p-y_c)/dy
    wym=1._wp-wyp
  else
    jp=j
    jm=j-1
    wym=(y_c-y_p)/dy
    wyp=1._wp-wym
  endif
  if(z_p.gt.z_c)then
    km=k
    if(k.ne.lls)then
      kp=k+1
      wzp=(z_p-z_c)/(zcart(z(kp),i,j)-z_c)
    else
      kp=lls
      wzp=(z_p-z_c)/(zcart(zedge(lls),i,j)-z_c)
    endif
    wzm=1._wp-wzp
  else
    kp=k
    km=k-1
    if(k.ne.1)then
      wzm=(z_c-z_p)/(z_c-zcart(z(km),i,j))
    else
      wzm=(z_c-z_p)/(z_c-zs(i,j))
    endif
    wzp=1._wp-wzm
  endif
  
  ! Triple interpolation (en.wikipedia.org/wiki/Trilinear_interpolation)
  u00=u_c(im,jm,km)*wxm+u_c(ip,jm,km)*wxp
  u01=u_c(im,jm,kp)*wxm+u_c(ip,jm,kp)*wxp
  u10=u_c(im,jp,km)*wxm+u_c(ip,jp,km)*wxp
  u11=u_c(im,jp,kp)*wxm+u_c(ip,jp,kp)*wxp
  u0 =u00*wym+u10*wyp
  u1 =u01*wym+u11*wyp
  u_p=u0*wzm+u1*wzp

  v00=v_c(im,jm,km)*wxm+v_c(ip,jm,km)*wxp
  v01=v_c(im,jm,kp)*wxm+v_c(ip,jm,kp)*wxp
  v10=v_c(im,jp,km)*wxm+v_c(ip,jp,km)*wxp
  v11=v_c(im,jp,kp)*wxm+v_c(ip,jp,kp)*wxp
  v0 =v00*wym+v10*wyp
  v1 =v01*wym+v11*wyp
  v_p=v0*wzm+v1*wzp
  
  w00=w_c(im,jm,km)*wxm+w_c(ip,jm,km)*wxp
  w01=w_c(im,jm,kp)*wxm+w_c(ip,jm,kp)*wxp
  w10=w_c(im,jp,km)*wxm+w_c(ip,jp,km)*wxp
  w11=w_c(im,jp,kp)*wxm+w_c(ip,jp,kp)*wxp
  w0 =w00*wym+w10*wyp
  w1 =w01*wym+w11*wyp
  w_p=w0*wzm+w1*wzp

end subroutine GetWind
