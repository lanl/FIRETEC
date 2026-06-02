!-----------------------------------------------------------------------
! Diffuse a variable phi that should be defined on [0,np+1]*[0,mp+1]*[0,l],
! This requires that the variable have been updated with corners and bottom bc such as done in (see fieldUpdate)
! Method: we compute directly the contravariant terms (fluxes normal to the surfaces of a cell)
! we specify a matrix [K] where the diagonal elements are the simple diffusivities in this direction 
! and the off diagonal elements would potentially be off-diagonal coef. if fully nonisotropic diffusion 
! with cross diffusion terms.  For the current stage of the code off-diagonal terms are zero
! we compute a matrix KG that combines the elements of the diffusion coeficient matrix Kij with metric tensor G using 
! the train of thought that [G]=[J][J]transpose therefore [KG]=[J][K][J]transpose
! d(Kij(dphi/dxi))/dxj (x and y are physical space) =1/sqrtG*d(sqrtG*KGij*dphi/dxj)/dxi (x and y are model space)
! hic=sumoverj(sqrtG*KGij*dphi/dxj) is the contravariant flux in the i direction times sqrt G  at the i-1/2 face therefore:

! the sqrtG_GKij are the following (based on arrays computed in fieldUpdate.f 
! sqrtG_GK11 = sqrtG_GK22 =sqrtG_Kxy(i,j,k)
! sqrtG_GK21 = sqrtG_GK12 = 0
! sqrtG_GK31 = sqrtG_GK13 =sqrtG_Kxy(i,j,k) * J31
! sqrtG_GK32 = sqrtG_GK23 =sqrtG_Kxy(i,j,k) * J32
! sqrtG_GK33 is already an array (computed in fieldUpdate.f) = sqrtG_Kxy * (J31**2+J32**2)+sqrtG_Kz * gi**
! as a result: with px, py and pz the derivative of phi in model space:
! hxc=sqrtG_Kxy * px + (sqrtG_Kxy * J31) * pz
! hyc=sqrtG_Kxy * py + (sqrtG_Kxy * J32) * pz
! hzc= (sqrtG_Kxy * J31) * px + (sqrtG_Kxy * J31) * py + sqrtG_KG33 * pz
! when K=k*Identity is isotropic KGij=k*Gij...

! The following arrays are defined on [0,np+1]*[0,mp+1]*[1,l] (without corners):
! rtkte_abc, sqrtG_Kxy,sqrtG_Kz,sqrtG_GK33
!-----------------------------------------------------------------------

subroutine diffuse(fphi,phi,il,iu,jl,ju,lls,cdiff)
  use gridlist_variables, only : prec
  use gridsetup, only : np,mp,dxi,dyi,dzi
  use metric_variables, only : c13,c23,gmul,gi
  use linn_turb_variables, only : sqrtG_Kxy,sqrtG_KG33
  Implicit None

  ! Local Variables
  integer,intent(in) :: il,iu,jl,ju,lls
  real(prec),intent(inout) :: fphi(il:iu,jl:ju,lls) ! rhs of variable to diffuse
  real(prec),intent(in) :: phi(il:iu,jl:ju,0:lls) ! variable to diffuse (defined on 0:l)
  real(prec),intent(in) :: cdiff  ! ratio between vt of variable iv and simple vt (ie rturbprandtl for theta and o2)
  
  real(prec),allocatable::    &
           hxc(:, :,:), & ! contravariant coordinate of x phi-flux  component*sqrt(G) at i-1/2
           hyc(:, :,:), & ! contravariant coordinate of y phi-flux  component*sqrt(G) at j-1/2
           hzc(:, :,:)    ! contravariant coordinate of z phi-flux  component*sqrt(G)  at k-1/2
  integer :: i,j,k,kp1,kk,km1
  real(prec) :: hdxi,hdyi,hdzi
  real(prec)::  sqrtG_KG11a,sqrtG_KG13a,sqrtG_KG22a,sqrtG_KG23a,sqrtG_KG31a &
   ,sqrtG_KG32a,sqrtG_KG33a !  edge values of sqrtG_KGij
  real(prec):: px,py,pz ! phi derivative on model grid
  real(prec) ::r
 
  !noise filter
  real(prec) :: phi_diff1, phi_diff2, phi_diff3, phi_diff4

  ! Executable Code 
  hdxi=0.5*dxi
  hdyi=0.5*dyi
  hdzi=0.5*dzi
  allocate(hxc(np+1,mp,lls)); hxc=0.0
  allocate(hyc(np,mp+1,lls)); hyc=0.0
  allocate(hzc(np,mp,lls+1)); hzc=0.0

  ! compute contravariant component x-flux at (i-1/2,j,k) hxc (sqrt(g) burried  in...) top bc of phi is hardcoded kp1=l when k=l
  do k=1,lls
    kp1=min(lls,k+1)
    do j=1,mp
      do i=1,np+1
        !at i-1/2
        sqrtG_KG11a=0.5*(sqrtG_Kxy(i,j,k)+sqrtG_Kxy(i-1,j,k))
        ! sqrtG_KG12a=0
        sqrtG_KG13a=0.5*(sqrtG_Kxy(i,j,k)*c13(i,j)*gmul(k)+sqrtG_Kxy(i-1,j,k)*c13(i-1,j)*gmul(k))
        
        !noise filter for px
        phi_diff1 = phi(i,j,k)-phi(i-1,j,k)
        if(abs(phi_diff1*1d+6).LT.phi(i,j,k)) phi_diff1=0.0
        px = phi_diff1 * dxi        
        
        pz = 0.5*(phi(i,j,kp1)-phi(i,j,k-1)+phi(i-1,j,kp1)-phi(i-1,j,k-1))*hdzi !at i-1/2
        ! compute contravariant flux*sqrt(g):
        hxc(i,j,k) = sqrtG_KG11a * px + sqrtG_KG13a * pz  
      enddo
    enddo
  enddo

  ! compute contravariant component y-flux at (i,j-1/2,k) hyc (sqrt(g) burried in...) top bc of phi is hardcoded kp1=l when k=l
  do k=1,lls
    kp1=min(lls,k+1)
    do i=1,np
      do j=1,mp+1                 
        !at j-1/2
        sqrtG_KG22a=0.5*(sqrtG_Kxy(i,j,k)+sqrtG_Kxy(i,j-1,k))
        sqrtG_KG23a=0.5*(sqrtG_Kxy(i,j,k)*c23(i,j)*gmul(k)+sqrtG_Kxy(i,j-1,k)*c23(i,j-1)*gmul(k))
        
        !noise filter for py
        phi_diff1 = phi(i,j,k)-phi(i,j-1,k)
        if(abs(phi_diff1*1d+6).LT.phi(i,j,k)) phi_diff1=0.0
        py= phi_diff1 * dyi
        
        pz = 0.5*(phi(i,j,kp1)-phi(i,j,k-1)+phi(i,j-1,kp1)-phi(i,j-1,k-1))*hdzi
        ! compute contravariant flux*sqrt(g):
        hyc(i,j,k) = sqrtG_KG22a * py + sqrtG_KG23a * pz 
      enddo
    enddo
  enddo
      
  ! compute contravariant component z-flux at i,j,k-1/2 hzc (sqrtg burried in)   top boundary condition on phi, rvtxy, rvtya is burried here (kk=l when k=l+1)
  do k=1,lls+1
    km1=max(1,k-1)
    kk=min(k,lls)
    do j=1,mp
      do i=1,np
        !at k-1/2
        sqrtG_KG31a=0.5*(sqrtG_Kxy(i,j,kk)*c13(i,j)*gmul(kk)+sqrtG_Kxy(i,j,km1)*c13(i,j)*gmul(km1))
        sqrtG_KG32a=0.5*(sqrtG_Kxy(i,j,kk)*c23(i,j)*gmul(kk)+sqrtG_Kxy(i,j,km1)*c23(i,j)*gmul(km1))
        sqrtG_KG33a=0.5*(sqrtG_KG33(i,j,kk)+sqrtG_KG33(i,j,km1))
        
        !noise filter for px and py
        phi_diff1 = phi(i+1,j,kk)-phi(i-1,j,kk)
        if(abs(phi_diff1*1d+6).LT.phi(i,j,kk)) phi_diff1=0.0
        phi_diff2 = phi(i+1,j,km1)-phi(i-1,j,km1)
        if(abs(phi_diff2*1d+6).LT.phi(i,j,kk)) phi_diff2=0.0
        px = 0.5 * (phi_diff1 + phi_diff2) * hdxi

        phi_diff3 = phi(i,j+1,kk)-phi(i,j-1,kk)
        if(abs(phi_diff3*1d+6).LT.phi(i,j,kk)) phi_diff3=0.0
        phi_diff4 = phi(i,j+1,km1)-phi(i,j-1,km1)
        if(abs(phi_diff4*1d+6).LT.phi(i,j,kk)) phi_diff4=0.0
        py = 0.5 * (phi_diff3 + phi_diff4) * hdyi
        
        pz = (phi(i,j,kk)-phi(i,j,k-1))*dzi
        ! compute contravariant flux *sqrt(g):
        hzc(i,j,k)=sqrtG_KG31a * px +sqrtG_KG32a * py + sqrtG_KG33a * pz
      enddo
    enddo
  enddo
 
  ! compute Laplacian term
  do k=1,lls
    do j=1,mp
      do i=1,np
        r=cdiff*gi(i,j,k)*(dxi*(hxc(i+1,j,k)-hxc(i,j,k))+dyi*(hyc(i,j+1,k)-hyc(i,j,k)) &
          +dzi*(hzc(i,j,k+1)-hzc(i,j,k)))
        fphi(i,j,k)=fphi(i,j,k)+r
      enddo
    enddo
  enddo

  deallocate (hxc,hyc,hzc)

end subroutine diffuse
