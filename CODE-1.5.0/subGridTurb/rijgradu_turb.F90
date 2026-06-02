!-----------------------------------------------------------------------
! rijgradu routine calculates changes and evolution of turbulence 
! parameters due to shear production
!-----------------------------------------------------------------------
subroutine rijgradu(fk,K_xy,K_z,rk,il,iu,jl,ju,lls)
  use gridlist_variables, only : prec
  use gridsetup, only : np,mp,dxi,dyi,dzi
  use metric_variables, only : gmul,c13,c23,gi
  use xvall, only : xvrho,iuvel,ivvel,iwvel
  Implicit None    
  
  ! Local Variables
  integer,intent(in) :: il,iu,jl,ju,lls
  real(prec),intent(inout) :: fk(il:iu,jl:ju,lls)
  real(prec),intent(in) :: K_xy(il:iu,jl:ju,lls)  ! K_axy or K_b
  real(prec),intent(in) :: K_z(il:iu,jl:ju,lls)   ! K_az or K_b
  real(prec),intent(in) :: rk(il:iu,jl:ju,lls)    ! rho*ka or rho*kb
 
  integer :: i,j,k,kp1,kv
  real(prec) :: rijgrad
  real(prec) :: other
  real(prec) :: tmp
  real(prec) :: hdxi,hdyi,hdzi,r11,r22,r33,r12,r13,r23,r31,r21,r32
  real(prec) :: J31,J32 ! jacobians
  real(prec),dimension(iuvel:iwvel) :: gradx,grady,gradz

  ! Executable Code
  hdxi=0.5*dxi
  hdyi=0.5*dyi
  hdzi=0.5*dzi

  do k=1,lls
    kp1=min(k+1,lls)  
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
        other=2./3.*(K_xy(i,j,k)*(gradx(iuvel)+grady(ivvel)) &
          +K_z(i,j,k)*gradz(iwvel)+rk(i,j,k))
        r11=-2*K_xy(i,j,k)*gradx(iuvel)+other
        r21=-K_xy(i,j,k)*(grady(iuvel)+gradx(ivvel))
        r31=-K_xy(i,j,k)*(gradz(iuvel)+gradx(iwvel))
        r12=-K_xy(i,j,k)*(grady(iuvel)+gradx(ivvel))  !=r21...
        r22=-2*K_xy(i,j,k)*grady(ivvel)+other
        r32=-K_xy(i,j,k)*(grady(iwvel)+gradz(ivvel))
        r13=-K_z(i,j,k)*(gradz(iuvel)+gradx(iwvel))
        r23=-K_z(i,j,k)*(grady(iwvel)+gradz(ivvel))
        r33=-2*K_z(i,j,k)*gradz(iwvel)+other
        rijgrad=-(r11*gradx(iuvel)+r12*grady(iuvel)+r13*gradz(iuvel) &
          +r21*gradx(ivvel)+r22*grady(ivvel) &
          +r23*gradz(ivvel)+r31*gradx(iwvel) &
          +r32*grady(iwvel)+r33*gradz(iwvel))
        fk(i,j,k)=fk(i,j,k)+rijgrad
      enddo
    enddo
  enddo
  
end subroutine rijgradu
