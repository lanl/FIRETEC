!----------------------------------------------------------------
! stressij routine computes the contribution to momentum equations from the derivative of 
! the the reynold stress tensor Rij
! Method : The main interest of this method is to compute with the most of accuracy
! the terms that are used because the average used to compute a given term at a given place
! are done globally (avoiding products of averages...)
! Rij=-Aij-Bij+deltaij*C with (cartesian coordinates):
!    Aij= Kj*dui/dxj
!    Bij= Kj*duj/dxi
!    C=2/3*(rtkte_abc+k1*du/dx+k2*dv/dy+k3*dw/dz)  (deltaij=1 if i=j, 0 otherwhise)
! because fu=-div(R1j), fv=-div(R2j) and fw=-div(R3j) and linearity of divergence, we can write:
!    fu=div(A1j)+div(B1j)-dC/dx 
!    fv=div(A2j)+div(B2j)-dC/dy 
!    fw=div(A3j)+div(B3j)-dC/dz 
! div(A1j)=1/sqrtG*d(A1jc)/dxj with A1jc=sumoverl(sqrtG*KGjl*du/dxl) (xl are model space)
! div(A2j)=1/sqrtG*d(A2jc)/dxj with A2jc=sumoverl(sqrtG*KGjl*dv/dxl) (xl are model space)
! div(A3j)=1/sqrtG*d(A3jc)/dxj with A3jc=sumoverl(sqrtG*KGjl*dw/dxl) (xl are model space)
! for each i, div(Bij)=1/sqrtG*d(Bijc)/dxj with
!   Bijc=sumoverlm(sqrtG*Kl*Jjl*Jmi*dul/dxm) (xm are model space)
!       later tlm=sqrtG*Kl*Jjl*Jmi
!       because Bijc = sqrtG*(Jjl*Bil) and Bil=Kl*dul/dxm*Jmi
! dC/dxi=J1i*Cx+J2i*Cy+J3i*Cz, where Cx,Cy and Cz are derivatives of C on model grid
!
! using linearity again :
!  fu=1/sqrtG*d(AB1jc)/dxj-Cx-J31*Cz, with AB1jc=A1jc+B1jc
!  fv=1/sqrtG*d(AB2jc)/dxj-Cy-J32*Cz, with AB2jc=A2jc+B2jc
!  fw=1/sqrtG*d(AB3jc)/dxj-J33*Cz, with AB3jc=A3jc+B3jc
!
! the sqrtG_KGij are the following (based on arrays computed in fieldUpdate.f 
! sqrtG_KG11 = sqrtG_KG22 =sqrtG_Kxy(i,j,k)
! sqrtG_KG21 = sqrtG_KG12 = 0
! sqrtG_KG31 = sqrtG_KG13 =sqrtG_Kxy(i,j,k) * J31
! sqrtG_KG32 = sqrtG_KG23 =sqrtG_Kxy(i,j,k) * J32
! sqrtG_KG33 is already an array (computed in fieldUpdate.f) = sqrtG_Kxy * (J31**2+J32**2)+sqrtG_Kz * gi**2	  
! these sqrtG_GKij are computed at cell face (and labeled sqrtG_GKija)
! for computation of C terms, we also need to compute cell face values Kxy_J31a,Kxy_J32a,Kz_J33a
! for computation of B terms, we need to compute cell face values sqrtG_Kxya, sqrtG_Kxy_J31a, 
! sqrt_Kxy_J32a,sqrtG_Kxy_J31_J31a, sqrtG_Kxy_J31_J32a,sqrtG_Kxy_J32_J32a, sqrtG_Kz_J33a,
! sqrtG_Kz_J31_J33a, sqrtG_Kz_J32_J33a, sqrtG_Kxy_J33a,sqrtG_Kz_J33_J33a  
! we use the fact that sqrtG=1/J33 and values of sqrtG_GK13, sqrtG_GK31, sqrtG_GK32 and sqrtGK23 
! already computed for Aijc terms to reduce the number of average to compute (and simplify notations...)
!
! Here we assume that:
! u,v,w are defined on [0,np+1]*[0,mp+1]*[0,l], meaning that they have been updated 
! and that bottom bc is implemented (top bc is hard coded here)
! for u, v, w : NB corners are required!
! The following arrays are defined on [0,np+1]*[0,mp+1]*[1,l] (without corners):
! rtkte_abc, sqrtG_Kxy,sqrtG_Kz,sqrtG_GK33
!----------------------------------------------------------------
subroutine stressrij(fx,fy,fz,il,iu,jl,ju,lls)
  use gridlist_variables, only : l,dx,dy,dz
  use metric_variables, only : gmul,c13,c23,gi
  use turb_variables, only : sqrtG_Kxy,sqrtG_Kz,sqrtG_KG33,rtke_abc
  use gridsetup, only : np,mp,dt,dxi,dyi,dzi
  use forcings, only : force
  use xvall, only : xvrho,xv,iuvel,ivvel,iwvel,ika,ikb
  use turba
  use metric_variables_old
  Implicit None

  ! Local Variables
  integer,intent(in) :: il,iu,jl,ju,lls
  real,intent(inout) :: fx(il:iu,jl:ju,lls)
  real,intent(inout) :: fy(il:iu,jl:ju,lls)
  real,intent(inout) :: fz(il:iu,jl:ju,lls)

  integer :: i,j,k,kp1,kk,km1
  real hdxi,hdyi,hdzi
  real:: J31,J32 ! jacobians
  ! cell face quantities for Aijc
  real:: sqrtG_KG11a,sqrtG_KG22a, sqrtG_KG13a,sqrtG_KG31a, &
   sqrtG_KG23a,sqrtG_KG32a, sqrtG_KG33a 
  ! cell face for Bijc
  real :: sqrtG_Kxya &
     ,sqrtG_Kxy_J31_J31a, sqrtG_Kxy_J31_J32a,&
     sqrtG_Kxy_J32_J32a, Kza,Kz_J31a, Kz_J32a, Kxya,Kz_J33a   
  ! cell face model derivative
  real::ux,uy,uz,vx,vy,vz,wx,wy,wz
  ! contravariant coordinates of Aij= Kj*dui/dxj
  real::A11c,A12c,A13c,A21c,A22c,A23c,A31c,A32c,A33c
  ! contravariant coordinates of Bij= Kj*duj/dxi
  real::B11c,B12c,B13c,B21c,B22c,B23c,B31c,B32c,B33c
  real::t11a,t12a,t13a,t21a,t22a,t23a,t31a,t32a,t33a 
  !  C=2/3*(rtkte_abc+k1*du/dx+k2*dv/dy+k3*dw/dz)
  ! cell face quantity for C
  real:: rtke_abca, Kxy_uxa, Kxy_vya, Kz_wza,Kxy_J31a,Kxy_J32a
  ! C arrays at i-1/2 (Ci), j-1/2 (Cj), k-1/2 (Ck)
  !real, allocatable::Ci(:,:,:),Cj(:,:,:),Ck(:,:,:)
  real::Cx,Cy,Cz ! C model grid derivatives
  real::fu,fv,fw ! contribution of reynold stres to momentum u, v,w

  ! Executable Code
  ! compute some local constants_old
  hdxi=0.5*dxi
  hdyi=0.5*dyi
  hdzi=0.5*dzi

  ! computation of Ai1c, Bi1c and C in i-1/2,j,k     
  ! top boundary condition on u,v,w, for uz, vz and wz is burried here (kp1=l when k=l)
  do k=1,l
    kp1=min(l,k+1)
    do j=1,mp
      do i=1,np+1
        ! 0/ computation of required i-1/2 face averages:
        sqrtG_KG11a=0.5*(sqrtG_Kxy(i,j,k)+sqrtG_Kxy(i-1,j,k))
        sqrtG_KG13a=0.5*(sqrtG_Kxy(i,j,k)*c13(i,j)*gmul(k) &
          +sqrtG_Kxy(i-1,j,k)*c13(i-1,j)*gmul(k))
        sqrtG_Kxya =0.5*(sqrtG_Kxy(i,j,k)+sqrtG_Kxy(i-1,j,k))
        sqrtG_KG23a=0.5*(sqrtG_Kxy(i,j,k)*c23(i,j)*gmul(k) &
          +sqrtG_Kxy(i-1,j,k)*c23(i-1,j)*gmul(k))
        Kxya    =0.5*(sqrtG_Kxy(i,j,k) * gi(i,j,k) &
          +sqrtG_Kxy(i-1,j,k) * gi(i-1,j,k))     
        Kxy_J31a=0.5*(sqrtG_Kxy(i,j,k)*gi(i,j,k)*c13(i,j)*gmul(k) &
          +sqrtG_Kxy(i-1,j,k)*gi(i-1,j,k)*c13(i-1,j)*gmul(k))
        Kxy_J32a=0.5*(sqrtG_Kxy(i,j,k)*gi(i,j,k)*c23(i,j)*gmul(k) &
          +sqrtG_Kxy(i-1,j,k)*gi(i-1,j,k)*c23(i-1,j)*gmul(k))
        Kz_J33a =0.5*(sqrtG_Kz(i,j,k)*gi(i,j,k)**2 &
          +sqrtG_Kz(i-1,j,k)*gi(i-1,j,k)**2)     
        rtke_abca=0.5*(rtke_abc(i,j,k) + rtke_abc(i-1,j,k))
        ux = (xvrho(i,j,k,iuvel)-xvrho(i-1,j,k,iuvel))*dxi
        vx = (xvrho(i,j,k,ivvel)-xvrho(i-1,j,k,ivvel))*dxi
        wx = (xvrho(i,j,k,iwvel)-xvrho(i-1,j,k,iwvel))*dxi
        uy = 0.5*(xvrho(i,j+1,k,iuvel)-xvrho(i,j-1,k,iuvel)+xvrho(i-1,j+1,k,iuvel)-xvrho(i-1,j-1,k,iuvel))*hdyi
        vy = 0.5*(xvrho(i,j+1,k,ivvel)-xvrho(i,j-1,k,ivvel)+xvrho(i-1,j+1,k,ivvel)-xvrho(i-1,j-1,k,ivvel))*hdyi
        uz = 0.5*(xvrho(i,j,kp1,iuvel)-xvrho(i,j,k-1,iuvel)+xvrho(i-1,j,kp1,iuvel)-xvrho(i-1,j,k-1,iuvel))*hdzi
        vz = 0.5*(xvrho(i,j,kp1,ivvel)-xvrho(i,j,k-1,ivvel)+xvrho(i-1,j,kp1,ivvel)-xvrho(i-1,j,k-1,ivvel))*hdzi
        wz = 0.5*(xvrho(i,j,kp1,iwvel)-xvrho(i,j,k-1,iwvel)+xvrho(i-1,j,kp1,iwvel)-xvrho(i-1,j,k-1,iwvel))*hdzi
      
        ! 1/ computation of Ai1c :
        A11c = sqrtG_KG11a * ux + sqrtG_KG13a * uz
        A21c = sqrtG_KG11a * vx + sqrtG_KG13a * vz
        A31c = sqrtG_KG11a * wx + sqrtG_KG13a * wz
      
        ! 2/ computation of Bi1c : NB J1l=1 only for l=1 (0 otherwhise), Bi1c=Bi1 !
        ! B11c : i=1,j=1; tlm=sqrtG*Kl*J1l*Jm1 
        t11a = sqrtG_Kxya            ! t11=sqrtG*Kxy*J11*J11			
        !t12a = 0                    ! t12=sqrtG*Kxy*J11*J21			
        t13a = sqrtG_KG13a           ! t13=sqrtG*Kxy*J11*J31
        B11c = t11a * ux + t13a * uz ! the others are 0 cause J12=J13=0
        !B21c : i=2,j=1; tlm=sqrtG*Kl*J1l*Jm2 
        !t11a = 0                     ! t11=sqrtG*Kxy*J11*J12			
        t12a = sqrtG_Kxya            ! t12=sqrtG*Kxy*J11*J22			
        t13a = sqrtG_KG23a           ! t13=sqrtG*Kxy*J11*J32
        B21c = t12a * uy + t13a * uz ! the others are 0 cause J12=J13=0
        ! B31c : i=3,j=1; tlm=sqrtG*J1l*Jm3 
        !t11a = 0                    ! t11=sqrtG*kxy*J11*J13=0			
        !t12a = 0                    ! t12=sqrtG*kxy*J11*J23=0			
        t13a = Kxya                    ! t13=sqrtG*kxy*J11*J33=kxy          
        B31c = t13a * uz ! the others are 0 cause J12=J13=0 

        ! 3/ ABi1c=Ai1c+Bi1c
        AB11c(i,j,k) = A11c+B11c
        AB21c(i,j,k) = A21c+B21c
        AB31c(i,j,k) = A31c+B31c
 
        ! 4/ Ci(i,j,k) is i-1/2 face value of C:
        ! C=2/3*(rtkte_abc+k1*du/dx+k2*dv/dy+k3*dw/dz)  
        Kxy_uxa = Kxya * ux + Kxy_J31a * uz ! (cartesian derivatives)
        Kxy_vya = Kxya * vy + Kxy_J32a * vz ! (cartesian derivatives) 
        Kz_wza = Kz_J33a * wz ! (cartesian derivatives)         
        Ci(i,j,k) = 2./3.*(rtke_abca + Kxy_uxa + Kxy_vya + Kz_wza)
      enddo
    enddo
  enddo

  ! computation of Ai2c, Bi2c and C in i,j-1/2,k     
  ! top boundary condition on u,v,w, for uz, vz and wz is burried here (kp1=l when k=l)
  do k=1,l
    kp1=min(l,k+1)
    do j=1,mp+1
      do i=1,np
        ! 0/ computation of required j-1/2 face averages:
        sqrtG_KG22a=0.5*(sqrtG_Kxy(i,j,k)+sqrtG_Kxy(i,j-1,k))
        sqrtG_KG23a=0.5*(sqrtG_Kxy(i,j,k)*c23(i,j)*gmul(k) &
                     +sqrtG_Kxy(i,j-1,k)*c23(i,j-1)*gmul(k))
        sqrtG_Kxya = 0.5*(sqrtG_Kxy(i,j,k)+sqrtG_Kxy(i,j-1,k))
        sqrtG_KG13a=0.5*(sqrtG_Kxy(i,j,k)*c13(i,j)*gmul(k) &
                    +sqrtG_Kxy(i,j-1,k)*c13(i,j-1)*gmul(k))
        Kxya = 0.5*(sqrtG_Kxy(i,j,k) * gi(i,j,k) &
                 + sqrtG_Kxy(i,j-1,k) * gi(i,j-1,k))     
        Kxy_J31a=0.5*(sqrtG_Kxy(i,j,k)*gi(i,j,k)*c13(i,j)*gmul(k) &
               +sqrtG_Kxy(i,j-1,k)*gi(i,j-1,k)*c13(i,j-1)*gmul(k))
        Kxy_J32a=0.5*(sqrtG_Kxy(i,j,k)*gi(i,j,k)*c23(i,j)*gmul(k) &
                +sqrtG_Kxy(i,j-1,k)*gi(i,j-1,k)*c23(i,j-1)*gmul(k))
        Kz_J33a = 0.5*(sqrtG_Kz(i,j,k) * gi(i,j,k)**2 + &
                       sqrtG_Kz(i,j-1,k) * gi(i,j-1,k)**2)     
        rtke_abca = 0.5 * (rtke_abc(i,j,k) + rtke_abc(i,j-1,k))
    
        ux = 0.5*(xvrho(i+1,j,k,iuvel)-xvrho(i-1,j,k,iuvel)+xvrho(i+1,j-1,k,iuvel)-xvrho(i-1,j-1,k,iuvel))*hdxi
        vx = 0.5*(xvrho(i+1,j,k,ivvel)-xvrho(i-1,j,k,ivvel)+xvrho(i+1,j-1,k,ivvel)-xvrho(i-1,j-1,k,ivvel))*hdxi
        uy = (xvrho(i,j,k,iuvel)-xvrho(i,j-1,k,iuvel))*dyi
        vy = (xvrho(i,j,k,ivvel)-xvrho(i,j-1,k,ivvel))*dyi
        wy = (xvrho(i,j,k,iwvel)-xvrho(i,j-1,k,iwvel))*dyi
        uz = 0.5*(xvrho(i,j,kp1,iuvel)-xvrho(i,j,k-1,iuvel)+xvrho(i,j-1,kp1,iuvel)-xvrho(i,j-1,k-1,iuvel))*hdzi
        vz = 0.5*(xvrho(i,j,kp1,ivvel)-xvrho(i,j,k-1,ivvel)+xvrho(i,j-1,kp1,ivvel)-xvrho(i,j-1,k-1,ivvel))*hdzi
        wz = 0.5*(xvrho(i,j,kp1,iwvel)-xvrho(i,j,k-1,iwvel)+xvrho(i,j-1,kp1,iwvel)-xvrho(i,j-1,k-1,iwvel))*hdzi
 
        ! 1/ computation of Ai2c :
        A12c = sqrtG_KG22a * uy + sqrtG_KG23a * uz
        A22c = sqrtG_KG22a * vy + sqrtG_KG23a * vz
        A32c = sqrtG_KG22a * wy + sqrtG_KG23a * wz
      
        ! 2/ computation of Bi2c :NB J1l=1 only for l=1 (0 otherwhise), Bi2c=Bi2 ! 
        ! B12c : i=1,j=2; tlm=sqrtG*Kl*J2l*Jm1 
        t21a = sqrtG_Kxya            ! t21=sqrtG*kxy*J22*J11			
        !t22a = 0                    ! t22=sqrtG*kxy*J22*J21			
        t23a = sqrtG_KG13a           ! t23=sqrtG*kxy*J22*J31
        B12c = t21a * vx + t23a * vz ! the others are 0 cause J21=J23=0
        ! B22c : i=2,j=2; tlm=sqrtG*Kl*J2l*Jm2 
        !t21a = 0                    ! t21=sqrtG*Kxy*J22*J12			
        t22a = sqrtG_Kxya            ! t22=sqrtG*Kxy*J22*J22			
        t23a = sqrtG_KG23a           ! t23=sqrtG*Kxy*J22*J32
        B22c = t22a * vy + t23a * vz ! the others are 0 cause J21=J23=0
        ! B32c : i=3,j=2; tlm=sqrtG*J2l*Jm3 
        !t21a = 0                    ! t21=sqrtG*kxy*J22*J13=0			
        !t22a = 0                    ! t22=sqrtG*kxy*J22*J23=0			
        t23a = Kxya                    ! t23=sqrtG*kxy*J22*J33          
        B23c = t23a * vz    ! the others are 0 cause J21=J23=0 

        ! 3/ ABi2c=Ai2c+Bi2c
        AB12c(i,j,k) = A12c+B12c
        AB22c(i,j,k) = A22c+B22c
        AB32c(i,j,k) = A32c+B32c
 
        ! 4/ Cj(i,j,k) is j-1/2 face value of C:
        ! C=2/3*(rtkte_abc+k1*du/dx+k2*dv/dy+k3*dw/dz)  
        Kxy_uxa = Kxya * ux + Kxy_J31a * uz ! (cartesian derivatives)
        Kxy_vya = Kxya * vy + Kxy_J32a * vz ! (cartesian derivatives) 
        Kz_wza = Kz_J33a * wz ! (cartesian derivatives)         
        Cj(i,j,k) = 2/3 * (rtke_abca + Kxy_uxa + Kxy_vya + Kz_wza)
      enddo
    enddo
  enddo


  ! computation of Ai3c, Bi3c and C in i,j,k-1/2     
  ! top boundary condition on u,v,w, for uz, vz and wz is burried here (kp1=l when k=l)
  do k=1,l+1
    km1=max(k-1,1)
    kk =min(k,l)    
    do j=1,mp
      do i=1,np
        ! 0/ computation of required k-1/2 face averages:
        sqrtG_KG31a=0.5*(sqrtG_Kxy(i,j,kk)*c13(i,j)*gmul(kk) &
                        +sqrtG_Kxy(i,j,km1)*c13(i,j)*gmul(km1))
        sqrtG_KG32a=0.5*(sqrtG_Kxy(i,j,kk)*c23(i,j)*gmul(kk) &
                        +sqrtG_Kxy(i,j,km1)*c23(i,j)*gmul(km1))
        sqrtG_KG33a=0.5*(sqrtG_KG33(i,j,kk)+sqrtG_KG33(i,j,km1))
        Kxya = 0.5*(sqrtG_Kxy(i,j,kk) * gi(i,j,kk) +  sqrtG_Kxy(i,j,km1) * gi(i,j,km1))     
        Kza = 0.5*(sqrtG_Kz(i,j,kk) * gi(i,j,kk) +  sqrtG_Kz(i,j,km1) * gi(i,j,km1))     
        Kz_J31a = 0.5*(sqrtG_Kz(i,j,kk) * gi(i,j,kk) * c13(i,j) * gmul(kk) &
            + sqrtG_Kz(i,j,km1) * gi(i,j,km1)* c13(i,j)* gmul(km1))   
        Kz_J32a = 0.5*(sqrtG_Kz(i,j,kk) * gi(i,j,kk) * c23(i,j) * gmul(kk) &
            + sqrtG_Kz(i,j,km1) * gi(i,j,km1)* c23(i,j)* gmul(km1))   
        Kz_J33a = 0.5*(sqrtG_Kz(i,j,kk) * gi(i,j,kk)**2 + sqrtG_Kz(i,j,km1) * gi(i,j,km1)**2)  
        sqrtG_Kxy_J31_J31a=0.5*(sqrtG_Kxy(i,j,kk)*(c13(i,j)*gmul(kk))**2+ &
           sqrtG_Kxy(i,j,km1)*(c13(i,j)*gmul(km1))**2) 
        sqrtG_Kxy_J31_J32a=0.5*(sqrtG_Kxy(i,j,kk)*c13(i,j)*c23(i,j)*(gmul(kk))**2+ &
           sqrtG_Kxy(i,j,km1)*c13(i,j)*c23(i,j)*(gmul(km1))**2)
        sqrtG_Kxy_J32_J32a=0.5*(sqrtG_Kxy(i,j,kk)*(c23(i,j)*gmul(kk))**2+ &
           sqrtG_Kxy(i,j,km1)*(c23(i,j)*gmul(km1))**2)
        Kxy_J31a=0.5*(sqrtG_Kxy(i,j,kk)*gi(i,j,kk)*c13(i,j)*gmul(kk) &
             +sqrtG_Kxy(i,j,km1)*gi(i,j,km1)*c13(i,j)*gmul(km1))
        Kxy_J32a=0.5*(sqrtG_Kxy(i,j,kk)*gi(i,j,kk)*c23(i,j)*gmul(kk) &
             +sqrtG_Kxy(i,j,km1)*gi(i,j,km1)*c23(i,j)*gmul(km1))
        rtke_abca = 0.5 * (rtke_abc(i,j,kk) + rtke_abc(i,j,km1))
     
        ux = 0.5*(xvrho(i+1,j,kk,iuvel)-xvrho(i-1,j,kk,iuvel)+xvrho(i+1,j,k-1,iuvel)-xvrho(i-1,j,k-1,iuvel))*hdxi
        vx = 0.5*(xvrho(i+1,j,kk,ivvel)-xvrho(i-1,j,kk,ivvel)+xvrho(i+1,j,k-1,ivvel)-xvrho(i-1,j,k-1,ivvel))*hdxi
        wx = 0.5*(xvrho(i+1,j,kk,iwvel)-xvrho(i-1,j,kk,iwvel)+xvrho(i+1,j,k-1,iwvel)-xvrho(i-1,j,k-1,iwvel))*hdxi
        uy = 0.5*(xvrho(i,j+1,kk,iuvel)-xvrho(i,j-1,kk,iuvel)+xvrho(i,j+1,k-1,iuvel)+xvrho(i,j-1,k-1,iuvel))*hdyi
        vy = 0.5*(xvrho(i,j+1,kk,ivvel)-xvrho(i,j-1,kk,ivvel)+xvrho(i,j+1,k-1,ivvel)+xvrho(i,j-1,k-1,ivvel))*hdyi
        wy = 0.5*(xvrho(i,j+1,kk,iwvel)-xvrho(i,j-1,kk,iwvel)+xvrho(i,j+1,k-1,iwvel)+xvrho(i,j-1,k-1,iwvel))*hdyi
        uz = (xvrho(i,j,kk,iuvel)-xvrho(i,j,k-1,iuvel))*hdzi
        vz = (xvrho(i,j,kk,ivvel)-xvrho(i,j,k-1,ivvel))*hdzi
        wz = (xvrho(i,j,kk,iwvel)-xvrho(i,j,k-1,iwvel))*hdzi
  
        ! 1/ computation of Ai3c :
        A13c=sqrtG_KG31a * ux +sqrtG_KG32a * uy + sqrtG_KG33a * uz
        A23c=sqrtG_KG31a * vx +sqrtG_KG32a * vy + sqrtG_KG33a * vz
        A33c=sqrtG_KG31a * wx +sqrtG_KG32a * wy + sqrtG_KG33a * wz

        ! 2/ computation of Bi3c =Bi1*J31+Bi2*J32+Bi3*J33 
        ! B13c : i=1,j=3; tlm=sqrtG*Kl*J3l*Jm1 
        t11a = sqrtG_KG31a           ! t11=sqrtG*Kxy*J31*J11			
        !t12a = 0                    ! t12=sqrtG*Kxy*J31*J21			
        t13a = sqrtG_Kxy_J31_J31a       ! t13=sqrtG*Kxy*J31*J31
        t21a = sqrtG_KG32a           ! t21=sqrtG*Kxy*J32*J11			
        !t22a = 0                    ! t22=sqrtG*Kxy*J32*J21			
        t23a = sqrtG_Kxy_J31_J32a       ! t23=sqrtG*Kxy*J32*J31
        t31a = Kza                    ! t31=sqrtG*kz*J33*J11=Kz			
        !t32a = 0                    ! t32=sqrtG*kz*J33*J21=0			
        t33a = Kz_J31a                 ! t33=sqrtG*kz*J33*J31=Kz*J31
        B13c = t11a * ux + t13a * uz + &
          t21a * vx + t23a * vz + t31a * wx + t33a * wz 
        ! B23c : i=2,j=3; tlm=sqrtG*Kl*J3l*Jm2 
        !t11a = 0                    ! t11=sqrtG*Kxy*J31*J12			
        t12a = sqrtG_KG31a            ! t12=sqrtG*Kxy*J31*J22			
        t13a = sqrtG_Kxy_J31_J32a       ! t13=sqrtG*Kxy*J31*J32
        !t21a = 0                    ! t21=sqrtG*Kxy*J32*J12			
        t22a = sqrtG_KG32a           ! t22=sqrtG*Kxy*J32*J22			
        t23a = sqrtG_Kxy_J32_J32a       ! t23=sqrtG*Kxy*J32*J32
        !t31a = 0                    ! t31=sqrtG*Kz*J33*J12			
        t32a = Kza                   ! t32=sqrtG*Kz*J33*J22=Kz			
        t33a = Kz_J32a                 ! t33=sqrtG*Kz*J33*J32=Kz*J32
        B23c = t12a * uy + t13a * uz + &
          t22a * vy + t23a * vz + t32a * wy + t33a * wz 
        ! B33c : i=3,j=3; tlm=sqrtG*Kl*J3l*Jm3 
        !t11a = 0                    ! t11=sqrtG*Kxy*J31*J13	=0	
        !t12a = 0                     ! t12=sqrtG*Kxy*J31*J23=0			
        t13a = Kxy_J31a              ! t13=sqrtG*Kxy*J31*J33=kxy*J31
        !t21a = 0                    ! t21=sqrtG*Kxy*J32*J13=0			
        !t22a = 0                    ! t22=sqrtG*Kxy*J32*J23	=0		
        t23a = Kxy_J32a             ! t23=sqrtG*Kxy*J32*J33=Kxy*J32
        !t31a = 0                    ! t31=sqrtG*Kz*J33*J13=0			
        !t32a = 0                   ! t32=sqrtG*Kz*J33*J23=0			
        t33a = Kz_J33a                 ! t33=sqrtG*Kz*J33*J33=Kxy*J33
        B33c = t13a * uz + t23a * vz + t33a * wz 

        ! 3/ ABi3c=Ai3c+Bi3c
        AB13c(i,j,k) = A13c+B13c
        AB23c(i,j,k) = A23c+B23c
        AB33c(i,j,k) = A33c+B33c
 
        ! 4/ Ck(i,j,k) is k-1/2 face value of C:
        ! C=2/3*(rtkte_abc+k1*du/dx+k2*dv/dy+k3*dw/dz)  
        Kxy_uxa = Kxya * ux + Kxy_J31a * uz ! (cartesian derivatives)
        Kxy_vya = Kxya * vy + Kxy_J32a * vz ! (cartesian derivatives) 
        Kz_wza = Kz_J33a * wz ! (cartesian derivatives)         
        Ck(i,j,k) = 2./3. * (rtke_abca + Kxy_uxa + Kxy_vya + Kz_wza)
      enddo
    enddo
  enddo

  !  computation of all derivatives of contravariant components ABijc
  ! and model grid Cx,Cy,Cz of C
  do k=1,l
    do j=1,mp
      do i=1,np  
        ! model derivatives of C :
        Cx = (Ci(i+1,j,k) - Ci(i,j,k)) * dxi
        Cy = (Cj(i,j+1,k) - Cj(i,j,k)) * dyi 
        Cz = (Ck(i,j,k+1) - Ck(i,j,k)) * dzi 
        J31 = c13(i,j) * gmul(k)
        J32 = c23(i,j) * gmul(k)
        fu = gi(i,j,k)*((AB11c(i+1,j,k)-AB11c(i,j,k))*dxi + & 
           (AB12c(i,j+1,k)-AB12c(i,j,k))*dyi + (AB13c(i,j,k+1)-AB13c(i,j,k))*dzi) &
           - Cx - J31 * Cz      
        fv = gi(i,j,k)*((AB21c(i+1,j,k)-AB21c(i,j,k))*dxi + &
           (AB22c(i,j+1,k)-AB22c(i,j,k))*dyi + (AB23c(i,j,k+1)-AB23c(i,j,k))*dzi) &
           - Cy - J32 * Cz      
        fw = gi(i,j,k)*((AB31c(i+1,j,k)-AB31c(i,j,k))*dxi + &
           (AB32c(i,j+1,k)-AB32c(i,j,k))*dyi + (AB33c(i,j,k+1)-AB33c(i,j,k))*dzi) &
           - gi(i,j,k) * Cz          
        fx(i,j,k)=fx(i,j,k)+2.*fu*dt
        fy(i,j,k)=fy(i,j,k)+2.*fv*dt
        fz(i,j,k)=fz(i,j,k)+2.*fw*dt
      enddo
    enddo
  enddo

end subroutine stressrij
