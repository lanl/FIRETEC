      subroutine rijgradu(fk,K_xy,K_z,rk,il,iu,jl,ju,lls)
      use gridsetup
      use metryic
      use bc
      use turba
      use updatedFields
      use msga

      Implicit None
      integer,intent(in) :: il,iu,jl,ju,lls
      real fk(il:iu,jl:ju,lls),
     .     K_xy(il:iu,jl:ju,lls),  ! K_axy or K_b
     .     K_z(il:iu,jl:ju,lls),  !K_az or K_b
     .     rk(il:iu,jl:ju,lls)   ! rho*ka or rho*kb
    
      real :: rijgrad
      real::ux,uy,uz,vx,vy,vz,wx,wy,wz ! derivative on cartesian coordinates
      real::other
      !cell centered cartesian-derivative, computed at cell center and divergence term 
      integer :: i,j,k,kp1
      real :: tmp
      real :: hdxi,hdyi,hdzi,r11,r22,r33,r12,r13,r23,r31,r21,r32
      real :: J31,J32 ! jacobians
      hdxi=0.5*dxi
      hdyi=0.5*dyi
      hdzi=0.5*dzi

      do k=1,l
      kp1=k+1
      if (k==l) kp1=l
      do j=1,mp
      do i=1,np
        J31 = c13(i,j)*gmul(k)
        J32 = c23(i,j)*gmul(k)
        ! cartesian derivatives
        tmp=(u(i,j,kp1)-u(i,j,k-1))*hdzi
        ux=(u(i+1,j,k)-u(i-1,j,k))*hdxi+J31*tmp
        uy=(u(i,j+1,k)-u(i,j-1,k))*hdyi+J32*tmp
        uz=gi(i,j,k)*tmp
        tmp=(v(i,j,kp1)-v(i,j,k-1))*hdzi
        vx=(v(i+1,j,k)-v(i-1,j,k))*hdxi+J31*tmp
        vy=(v(i,j+1,k)-v(i,j-1,k))*hdyi+J32*tmp
        vz=gi(i,j,k)*tmp
        tmp=(w(i,j,kp1)-w(i,j,k-1))*hdzi
        wx=(w(i+1,j,k)-w(i-1,j,k))*hdxi+J31*tmp
        wy=(w(i,j+1,k)-w(i,j-1,k))*hdyi+J32*tmp
        wz=gi(i,j,k)*tmp
        other=2./3.*(K_xy(i,j,k)*(ux+vy)+K_z(i,j,k)*wz+rk(i,j,k))
        r11=-2*K_xy(i,j,k)*ux+other
        r21=-K_xy(i,j,k)*(uy+vx)
        r31=-K_xy(i,j,k)*(uz+wx)
        r12=-K_xy(i,j,k)*(uy+vx)  !=r21...
        r22=-2*K_xy(i,j,k)*vy+other
        r32=-K_xy(i,j,k)*(wy+vz)
        r13=-K_z(i,j,k)*(uz+wx)
        r23=-K_z(i,j,k)*(wy+vz)
        r33=-2*K_z(i,j,k)*wz+other
        rijgrad=-(r11*ux+r12*uy+r13*uz+r21*vx+r22*vy+r23*vz
     +   +r31*wx+r32*wy+r33*wz)
        fk(i,j,k)=fk(i,j,k)+rijgrad*2*dt
      enddo
      enddo
      enddo
      !call updated(fk,fk,np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      return
      end
