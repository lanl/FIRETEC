!-----------------------------------------------------------------------
! small_explicit_forcings calculates all the forcing terms that will be 
! applied in the implicit stepping method but will only be calculated 
! once on the large time step
!-----------------------------------------------------------------------
subroutine calcForceTurb(xv,force_xv,il,iu,jl,ju,lls,nv)
  use gridlist_variables, only : iturb,wp
  use turb_gridlist_variables, only : rturbprandtl
  use xvall, only : iuvel,ivvel,iwvel,itheta,ika,ikb,xvrho
  use linn_turb_variables, only : K_axy,K_az,K_b,sa,sb
  Implicit None

  ! Local Variables
  integer,intent(in) :: il,iu,jl,ju,lls,nv
  real(wp),intent(in) :: xv(il:iu,jl:ju,lls,nv)
  real(wp),intent(inout) :: force_xv(il:iu,jl:ju,lls,nv)

  integer :: kv

  ! Executable Code
  if(iturb.eq.1)then ! Diffusive Radiation Scheme
    ! 'Smagorinski Not implemented yet'
  elseif (iturb.eq.2) then ! Linn Subgrid Turbulence Model
    
    ! compute forces on momentum
    call stressrij(force_xv(:,:,:,iuvel),force_xv(:,:,:,ivvel), &
      force_xv(:,:,:,iwvel),il,iu,jl,ju,lls)
    
    ! shear production
    call rijgradu(force_xv(:,:,:,ika),K_axy,K_az,xv(:,:,:,ika), &
      il,iu,jl,ju,lls)
    call rijgradu(force_xv(:,:,:,ikb),K_b,K_b,xv(:,:,:,ikb), &
      il,iu,jl,ju,lls)
    
    ! turbulence dissipation+vegetation effect
    call dragtk(force_xv(:,:,:,ika),xv(:,:,:,ika),sa, &
      0._wp*xv(:,:,:,ika),il,iu,jl,ju,lls)
    call dragtk(force_xv(:,:,:,ikb),xv(:,:,:,ikb),sb, &
      xv(:,:,:,ika),il,iu,jl,ju,lls)
     
    ! turbulent diffusion of all non-wind fields
    do kv=itheta,nv-1
      if(kv.eq.ika.or.kv.eq.ikb)then
        call diffuse(force_xv(:,:,:,kv),xvrho(:,:,:,kv), &
          il,iu,jl,ju,lls,1._wp)
      else
        call diffuse(force_xv(:,:,:,kv),xvrho(:,:,:,kv), &
          il,iu,jl,ju,lls,rturbprandtl)
      endif
    enddo

  endif
     
end subroutine calcForceTurb
