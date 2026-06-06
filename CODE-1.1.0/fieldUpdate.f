!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
! - tkea, tkeb : xvb(5)/xvb(nv), xvb(6)/xvb(nv)
! - theta, ox, vap
! For the following arrays corners are not required
! - K_axy,K_az,K_b : diffusivity coefficient at scale a (in horizontal and vertical direction) and at scale b
! - sqrtG_Kxy = sqrtG*(K_axy+(1+sc*sqrt(kbcratio))*K_b)
! - sqrtG_Kz = sqrtG*(K_az+(1+sc*sqrt(kbcratio))*K_b)
! - sqrtG_GK33 = sqrtG*(Kxy * J31**2 + Kxy * J32**2 +Kz * J33**2) ! generalised metric tensor term

      subroutine fieldUpdate()
      use gridsetup
      use updatedFields
      use turba
      use metryic
      use msga
      use xvo
      integer::i,j,k 
      real::vtbc=0
      real:: J31, J32 ! jacobians      
      real:: g13, g23, g33 ! metric tensor 
      real::uground,vground,wground ! cartesian at ground level     
      real, external::getSubGroundValue ! value subground assuming freeslip

      if (iturb.eq.2) vtbc = 1+sc*sqrt(kbcratio)
      do k=1,l
      do j=1,mp
      do i=1,np
        u(i,j,k) = xvb(i,j,k,1)/xvb(i,j,k,nv)
        v(i,j,k) = xvb(i,j,k,2)/xvb(i,j,k,nv)
        w(i,j,k) = xvb(i,j,k,3)/xvb(i,j,k,nv)

        K_axy(i,j,k) = diffcst * xvb(i,j,k,nv) * 
     +        saxy(i,j,k) * sqrt(xvb(i,j,k,5)/xvb(i,j,k,nv))
        K_az(i,j,k) = diffcst * xvb(i,j,k,nv) * 
     +        saz(i,j,k) * sqrt(xvb(i,j,k,5)/xvb(i,j,k,nv))
        K_b(i,j,k) = diffcst * xvb(i,j,k,nv) * 
     +         sb(i,j,k) * sqrt(xvb(i,j,k,6)/xvb(i,j,k,nv))
        ! diffusion arrays    
        theta(i,j,k)=xvb(i,j,k,4)/xvb(i,j,k,nv)
        tke_a(i,j,k)=xvb(i,j,k,5)/xvb(i,j,k,nv)
        tke_b(i,j,k)=xvb(i,j,k,6)/xvb(i,j,k,nv)
        ox(i,j,k)=xvb(i,j,k,7)/xvb(i,j,k,nv)
        rtke_abc(i,j,k)=xvb(i,j,k,5) + (1+kbcratio) * xvb(i,j,k,6)
        sqrtG_Kxy(i,j,k) = 1/gi(i,j,k) * (K_axy(i,j,k)+vtbc*K_b(i,j,k))
        sqrtG_Kz(i,j,k) = 1/gi(i,j,k) * (K_az(i,j,k)+vtbc*K_b(i,j,k))
        J31 = c13(i,j)*gmul(k)
        J32 = c23(i,j)*gmul(k)
        sqrtG_KG33(i,j,k) = sqrtG_Kxy (i,j,k) * (J31**2+J32**2)
     +           +sqrtG_Kz(i,j,k) * gi(i,j,k)**2
        ! the other components of sqrtGKij can easily be computed on the fly:
        ! sqrtG_GK11 = sqrtG_GK22 =sqrtG_Kxy(i,j,k)
        ! sqrtG_GK21 = sqrtG_GK12 = 0
        ! sqrtG_GK31 = sqrtG_GK13 =sqrtG_Kxy(i,j,k) * J31
        ! sqrtG_GK32 = sqrtG_GK23 =sqrtG_Kxy(i,j,k) * J32
        ! when K=k*Identity is isotropic KGij=k*Gij...

        !firetec arrays
        sqrtk(i,j,k) = sqrt(rtke_abc(i,j,k)/xvb(i,j,k,nv))

      enddo
      enddo
      enddo   

      ! update including corners
      call updated(u(:,:,1:l),u(:,:,1:l),np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      call updated(v(:,:,1:l),v(:,:,1:l),np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      call updated(w(:,:,1:l),w(:,:,1:l),np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      call updated(theta(:,:,1:l),theta(:,:,1:l),np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      call updated(tke_a(:,:,1:l),tke_a(:,:,1:l),np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      call updated(tke_b(:,:,1:l),tke_b(:,:,1:l),np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      call updated(ox(:,:,1:l),ox(:,:,1:l),np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
      call updated(rtke_abc(:,:,1:l),rtke_abc(:,:,1:l),np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)

      ! update without corners
      !call updated(K_axy,K_axy,np,mp,l,1-ih,np+ih,1-ih,mp+ih,0)
      !call updated(K_az,K_az,np,mp,l,1-ih,np+ih,1-ih,mp+ih,0)
      !call updated(K_b,K_b,np,mp,l,1-ih,np+ih,1-ih,mp+ih,0)
      !call updated(sqrtk,sqrtk,np,mp,l,1-ih,np+ih,1-ih,mp+ih,0)
      call updated(sqrtG_Kxy,sqrtG_Kxy,np,mp,l,1-ih,np+ih,1-ih,mp+ih,0)
      call updated(sqrtG_Kz,sqrtG_Kz,np,mp,l,1-ih,np+ih,1-ih,mp+ih,0)
      call updated(sqrtG_KG33,sqrtG_KG33,np,mp,l,1-ih,np+ih,1-ih,mp+ih,0)

      ! bottom boundary of uz,vz,wz... (NB : should be updated first!) 
      if (ih.lt.3) then
         write(6,*) 
     +    "last halo cell will not be set correctly (see fieldUpdate.f)"
         stop
      endif
      
      k=1 !bottombc
      do j=1-ih+1,mp+ih-1
      do i=1-ih+1,np+ih-1
        g13 = c13(i,j)*gmul(1)
        g23 = c23(i,j)*gmul(1)
        g33 = g13**2+g23**2+gi(i,j,1)**2
        !uz : grad(u) 's normal to ground component is 0
        u(i,j,0) = getSubGroundValue(u(i,j,k),u(i-1,j,k),u(i+1,j,k),u(i,j-1,k),u(i,j+1,k)
     +              ,g13,g23,g33)
        !vz : grad(v) 's normal to ground component is 0
        v(i,j,0) = getSubGroundValue(v(i,j,k),v(i-1,j,k),v(i+1,j,k),v(i,j-1,k),v(i,j+1,k)
     +                ,g13,g23,g33)
        !wz : contravariant w is 0 at ground
        uground=0.5*(u(i,j,1)+u(i,j,0))
        vground=0.5*(v(i,j,1)+v(i,j,0))
        wground = - (g13 * uground + g23 * vground)/gi(i,j,1)
        w(i,j,0) = 2*wground-w(i,j,1)
        theta(i,j,0) = getSubGroundValue(theta(i,j,k),theta(i-1,j,k),theta(i+1,j,k)
     +        ,theta(i,j-1,k),theta(i,j+1,k),g13,g23,g33)
        tke_a(i,j,0) = getSubGroundValue(tke_a(i,j,k),tke_a(i-1,j,k),tke_a(i+1,j,k)
     +        ,tke_a(i,j-1,k),tke_a(i,j+1,k),g13,g23,g33)
        tke_b(i,j,0) = getSubGroundValue(tke_b(i,j,k),tke_b(i-1,j,k),tke_b(i+1,j,k)
     +        ,tke_b(i,j-1,k),tke_b(i,j+1,k),g13,g23,g33)
        ox(i,j,0) = getSubGroundValue(ox(i,j,k),ox(i-1,j,k),ox(i+1,j,k),ox(i,j-1,k)
     +        ,ox(i,j+1,k),g13,g23,g33)
        rtke_abc(i,j,0) = getSubGroundValue(rtke_abc(i,j,k),rtke_abc(i-1,j,k),rtke_abc(i+1,j,k)
     +        ,rtke_abc(i,j-1,k),rtke_abc(i,j+1,k),g13,g23,g33)
      enddo
      enddo
      
      ! below same stuf for irhovapor
      if (irhovapor.eq.1) then
         do k=1,l
         do j=1,mp
         do i=1,np
            vap(i,j,k)=xvb(i,j,k,8)/xvb(i,j,k,nv) 
         enddo
         enddo
         enddo
         call updated(vap(:,:,1:l),vap(:,:,1:l),np,mp,l,1-ih,np+ih,1-ih,mp+ih,1)
         k=1 !bottombc
         do j=1-ih+1,mp+ih-1
         do i=1-ih+1,np+ih-1
           g13 = c13(i,j)*gmul(1)
           g23 = c23(i,j)*gmul(1)
           g33 = g13**2+g23**2+gi(i,j,1)**2
           vap(i,j,0) = getSubGroundValue(vap(i,j,k),vap(i-1,j,k),
     +      vap(i+1,j,k),vap(i,j-1,k),vap(i,j+1,k),g13,g23,g33)
         enddo
         enddo 
      endif   ! irhovapor
      
          
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      !free slip bc
     
      real function getSubGroundValue(u,uim1,uip1,ujm1,ujp1,g13,g23,g33)
        use gridsetup
        Implicit None
        real u,uim1,uip1,ujm1,ujp1,g13,g23,g33
        real ux1,uy1 ! 2dx - u derivatives in cell 1 on the model grid
        real uzground ! gradient at the wall on the model grid
        ux1=0.5*(uip1-uim1)*dxi
        uy1=0.5*(ujp1-ujm1)*dyi
        uzground = -(g13 * ux1 + g23 * uy1) / g33 !uz : grad(u) 's normal to ground component is 0
        getSubGroundValue=u - dz*uzground
        return
      end function getSubGroundValue
      
