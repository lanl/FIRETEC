!-----------------------------------------------------------------------
! compute cell centered quantities including bc that are used for turb and firetec
!
! u,v,w are defined on [0,np+1]*[0,mp+1]*[0,l], updated with corners
! value for k=0 is hardcoded here, assuming wcontravariant=0 at ground for w, and 
! gradient normal to the ground is zero for u,v, and any scalars
! NB : for k=0, quantities are not defined in the last halo cell (i=1-ih, i=np+ih, j=1-ih, j=l+ih) !!!!!

! For rijgradu, the following arrays are required (no update require...)
! - K_axy,K_az,K_b : diffusivity coefficient at scale a (in horizontal and vertical direction) and at scale b
! For "diffuse" and "stressrij", the following arrays are defined on [0,np+1]*[0,mp+1]*[1,l] 
! for the following arrays corners are required for computation of bottom bc  :
! - rtke_abc : total modeled tke times rho   ! TODO : check but can probably be removed
! - tkea, tkeb : xv(5)/xv(nv), xv(6)/xv(nv)
! For the following arrays corners are not required
! - K_axy,K_az,K_b : diffusivity coefficient at scale a (in horizontal and vertical direction) and at scale b
! - sqrtG_Kxy = sqrtG*(K_axy+(1+sc*sqrt(kbcratio))*K_b)
! - sqrtG_Kz = sqrtG*(K_az+(1+sc*sqrt(kbcratio))*K_b)
! - sqrtG_GK33 = sqrtG*(Kxy * J31**2 + Kxy * J32**2 +Kz * J33**2) ! generalised metric tensor term
!-----------------------------------------------------------------------
subroutine turbFieldUpdate()
  use gridlist_variables, only : l,ih,iturb,prec
  use xvall, only : xv,irho,ika,ikb
  use metric_variables, only : gmul,gi,c13,c23
  use gridsetup, only : np,mp
  use linn_turb_variables, only : K_axy,K_az,K_b,sqrtG_Kxy,sqrtG_Kz, &
    sqrtG_KG33,saxy,saz,sb,sc,rtke_abc
  use SubGround_function
  Implicit none

  ! Local Variables
  integer::i,j,k
  real(prec)::vtbc=0
  real(prec) :: kbcratio=0.2
  real(prec) :: diffcst=0.09
  real(prec):: J31, J32 ! jacobians      
  real(prec):: g13, g23, g33 ! metric tensor 

  ! Executable Code
  if (iturb.eq.2) vtbc = 1+sc*sqrt(kbcratio)
  do k=1,l
    do j=1,mp
      do i=1,np
        K_axy(i,j,k) = diffcst * xv(i,j,k,irho) * &
          saxy(i,j,k) * sqrt(xv(i,j,k,ika)/xv(i,j,k,irho))
        K_az(i,j,k) = diffcst * xv(i,j,k,irho) * &
          saz(i,j,k) * sqrt(xv(i,j,k,ika)/xv(i,j,k,irho))
        K_b(i,j,k) = diffcst * xv(i,j,k,irho) * &
          sb(i,j,k) * sqrt(xv(i,j,k,ikb)/xv(i,j,k,irho))
        ! diffusion arrays    
        rtke_abc(i,j,k)=xv(i,j,k,ika) + (1+kbcratio) * xv(i,j,k,ikb)
        sqrtG_Kxy(i,j,k) = 1/gi(i,j,k) * (K_axy(i,j,k)+vtbc*K_b(i,j,k))
        sqrtG_Kz(i,j,k) = 1/gi(i,j,k) * (K_az(i,j,k)+vtbc*K_b(i,j,k))
        J31 = c13(i,j)*gmul(k)
        J32 = c23(i,j)*gmul(k)
        sqrtG_KG33(i,j,k) = sqrtG_Kxy (i,j,k) * (J31**2+J32**2) &
          +sqrtG_Kz(i,j,k) * gi(i,j,k)**2
      enddo
    enddo
  enddo         

  ! update including corners
  call update(rtke_abc(:,:,1:l),rtke_abc(:,:,1:l),np,mp,l, &
    1-ih,np+ih,1-ih,mp+ih,1,0)

  ! update without corners
  call update(sqrtG_Kxy,sqrtG_Kxy,np,mp,l,1-ih,np+ih,1-ih,mp+ih,0,0)
  call update(sqrtG_Kz,sqrtG_Kz,np,mp,l,1-ih,np+ih,1-ih,mp+ih,0,0)
  call update(sqrtG_KG33,sqrtG_KG33,np,mp,l,1-ih,np+ih,1-ih,mp+ih,0,0)
  
  k=1 !bottombc
  do j=1-ih+1,mp+ih-1
    do i=1-ih+1,np+ih-1
      g13 = c13(i,j)*gmul(1)
      g23 = c23(i,j)*gmul(1)
      g33 = g13**2+g23**2+gi(i,j,1)**2
      rtke_abc(i,j,0) = getSubGroundValue(rtke_abc(i,j,k), &
       rtke_abc(i-1,j,k),rtke_abc(i+1,j,k), &
       rtke_abc(i,j-1,k),rtke_abc(i,j+1,k),g13,g23,g33)
    enddo
  enddo
      
end subroutine turbFieldUpdate
