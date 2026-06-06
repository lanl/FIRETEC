!----------------------------------------------------------------
! rijgradu routine calculates changes and evolution of turbulence 
! parameters due to shear production
!----------------------------------------------------------------
subroutine rijgradu(fk,K_xy,K_z,rk,il,iu,jl,ju,lls)
  use gridlist_variables, only : l,dx,dy,dz
  use gridsetup, only : np,mp,dt,dxi,dyi,dzi
  use metric_variables, only : gmul,c13,c23,gi
  use xvall, only : xvrho,iuvel,ivvel,iwvel
  Implicit None    
  
  ! Local Variables
  integer,intent(in) :: il,iu,jl,ju,lls
  real,intent(inout) :: fk(il:iu,jl:ju,lls)
  real,intent(in) :: K_xy(il:iu,jl:ju,lls)  ! K_axy or K_b
  real,intent(in) :: K_z(il:iu,jl:ju,lls)   ! K_az or K_b
  real,intent(in) :: rk(il:iu,jl:ju,lls)    ! rho*ka or rho*kb
 
  integer :: i,j,k,kp1,kv
  real :: rijgrad
  real :: other
  real :: tmp
  real :: hdxi,hdyi,hdzi,r11,r22,r33,r12,r13,r23,r31,r21,r32
  real :: J31,J32 ! jacobians
  real,dimension(3) :: gradx,grady,gradz

  ! Executable Code
  hdxi=0.5*dxi
  hdyi=0.5*dyi
  hdzi=0.5*dzi

  do k=1,l
    kp1=min(k+1,l)  
    do j=1,mp
      do i=1,np
        J31=c13(i,j)*gmul(k)
        J32=c23(i,j)*gmul(k)
        ! cartesian derivatives
        ! cell centered cartesian-derivative, computed at cell center and divergence term 
        do kv=iuvel,iwvel
          tmp=(xvrho(i,j,kp1,kv)-xvrho(i,j,k-1,kv))*hdzi
          gradx(kv)=(xvrho(i+1,j,k,kv)-xvrho(i-1,j,k,kv))*hdxi+J31*tmp
          grady(kv)=(xvrho(i,j+1,k,kv)-xvrho(i,j-1,k,kv))*hdyi+J32*tmp
          gradz(kv)=gi(i,j,k)*tmp
        enddo
        other=2./3.*(K_xy(i,j,k)*(gradx(1)+grady(2))+K_z(i,j,k)*gradz(3)+rk(i,j,k))
        r11=-2*K_xy(i,j,k)*gradx(1)+other
        r21=-K_xy(i,j,k)*(grady(1)+gradx(2))
        r31=-K_xy(i,j,k)*(gradz(1)+gradx(3))
        r12=-K_xy(i,j,k)*(grady(1)+gradx(2))  !=r21...
        r22=-2*K_xy(i,j,k)*grady(2)+other
        r32=-K_xy(i,j,k)*(grady(3)+gradz(2))
        r13=-K_z(i,j,k)*(gradz(1)+gradx(3))
        r23=-K_z(i,j,k)*(grady(3)+gradz(2))
        r33=-2*K_z(i,j,k)*gradz(3)+other
        rijgrad=-(r11*gradx(1)+r12*grady(1)+r13*gradz(1)+r21*gradx(2)+r22*grady(2) &
          +r23*gradz(2)+r31*gradx(3)+r32*grady(3)+r33*gradz(3))
        fk(i,j,k)=fk(i,j,k)+rijgrad*2*dt
      enddo
    enddo
  enddo
  
end subroutine rijgradu
