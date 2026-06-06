subroutine turb()
  use gridlist_variables, only : iturb,ih,l,rturbprandtl,ifire
  use gridsetup, only : np,mp,dti
  use forcings, only : force
  use xvall, only : itemp,nv,ika,ikb,xvrho,xv,iuvel,ivvel,iwvel
  use turb_variables, only : K_axy,K_az,K_b,sa,sb
  use turba
  use workavg
  Implicit None
 
  ! Local Variables
  integer :: i,j,k,kv

  ! Executable Code
  force=0.

! 1/ compute forces on momentum
  call stressrij(force(:,:,:,1),force(:,:,:,2),force(:,:,:,3),1-ih,np+ih,1-ih,mp+ih,l)
  ! compute drag contribution
  call dragm(xv,1-ih,np+ih,1-ih,mp+ih,l,nv)
     
  ! compute forces for ka and kb
  ! shear production
  call rijgradu(force(:,:,:,ika),K_axy,K_az,xv(:,:,:,ika),1-ih,np+ih,1-ih,mp+ih,l)
  ! turbulence dissipation+vegetation effect
  call dragtk(force(:,:,:,ika),xv(:,:,:,ika),sa,0.*xv(:,:,:,ika), &
    1-ih,np+ih,1-ih,mp+ih,l)
  if (iturb.eq.2) then
    ! shear production
    call rijgradu(force(:,:,:,ikb),K_b,K_b,xv(:,:,:,ikb),1-ih,np+ih,1-ih,mp+ih,l)
    ! turbulence dissipation+vegetation effect
    call dragtk(force(:,:,:,ikb),xv(:,:,:,ikb),sb,xv(:,:,:,ika), &
      1-ih,np+ih,1-ih,mp+ih,l)
  endif
     
  ! diffusion of all non-wind fields
  do kv=itemp,nv-1
    if(kv.eq.ika.or.kv.eq.ikb)then
      call diffuse(force(:,:,:,kv),xvrho(:,:,:,kv),1-ih,np+ih,1-ih,mp+ih,l,1.)
    else
      call diffuse(force(:,:,:,kv),xvrho(:,:,:,kv),1-ih,np+ih,1-ih,mp+ih,l,rturbprandtl)
    endif
  enddo

  ! prepare contrubution of these terms to rhs
  if(ifire.eq.1) then
    do k=1,l
      do j=1,mp
        do i=1,np
          do kv=1,nv
            force(i,j,k,kv)=force(i,j,k,kv)*dti
          enddo
        enddo
      enddo
    enddo
  else  ! ifire.eq.0
    do k=1,l
      do j=1,mp
        do i=1,np
          xv(i,j,k,1)=xv(i,j,k,1)+0.5*force(i,j,k,iuvel)
          xv(i,j,k,2)=xv(i,j,k,2)+0.5*force(i,j,k,ivvel)
          xv(i,j,k,3)=xv(i,j,k,3)+0.5*force(i,j,k,iwvel)
          xv(i,j,k,5)=xv(i,j,k,5)+0.5*force(i,j,k,ika)
          xv(i,j,k,6)=xv(i,j,k,6)+0.5*force(i,j,k,ikb)
         enddo
       enddo
     enddo
  endif 

end subroutine turb
