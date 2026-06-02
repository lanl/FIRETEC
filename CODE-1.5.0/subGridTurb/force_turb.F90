!-----------------------------------------------------------------------
! small_explicit_forcings calculates all the forcing terms that will be 
! applied in the implicit stepping method but will only be calculated 
! once on the large time step
!-----------------------------------------------------------------------
subroutine calcForceTurb(xv,force_xv,il,iu,jl,ju,lls,nv)
  use gridlist_variables, only : iturb,ih,l,prec
  use turb_gridlist_variables, only : rturbprandtl
  use gridsetup, only : np,mp
  use xvall, only : iuvel,ivvel,iwvel,itemp,ika,ikb,xvrho
  use linn_turb_variables, only : K_axy,K_az,K_b,sa,sb
  use msga_variables, only : mpi_rank,ierror
  Implicit None

  ! Local Variables
  integer,intent(in) :: il,iu,jl,ju,lls,nv
  real(prec),intent(in) :: xv(il:iu,jl:ju,lls,nv)
  real(prec),intent(inout) :: force_xv(il:iu,jl:ju,lls,nv)

  integer :: kv
  real(prec) :: one=1.0

  ! Executable Code
  if(iturb.eq.1)then ! Diffusive Radiation Scheme
    if(mpi_rank.eq.0) print*,'Smagorinski Not implemented'
    call mpi_finalize(ierror)
    STOP
  elseif (iturb.eq.2) then ! Linn Subgrid Turbulence Model
    
    ! compute forces on momentum
    call stressrij(force_xv(:,:,:,iuvel),force_xv(:,:,:,ivvel), &
      force_xv(:,:,:,iwvel),1-ih,np+ih,1-ih,mp+ih,l)
    
    ! shear production
    call rijgradu(force_xv(:,:,:,ika),K_axy,K_az,xv(:,:,:,ika), &
      1-ih,np+ih,1-ih,mp+ih,l)
    call rijgradu(force_xv(:,:,:,ikb),K_b,K_b,xv(:,:,:,ikb), &
      1-ih,np+ih,1-ih,mp+ih,l)
    
    ! turbulence dissipation+vegetation effect
    call dragtk(force_xv(:,:,:,ika),xv(:,:,:,ika),sa, &
      0.*xv(:,:,:,ika),1-ih,np+ih,1-ih,mp+ih,l)
    call dragtk(force_xv(:,:,:,ikb),xv(:,:,:,ikb),sb, &
      xv(:,:,:,ika),1-ih,np+ih,1-ih,mp+ih,l)
     
    ! turbulent diffusion of all non-wind fields
    do kv=itemp,nv-1
      if(kv.eq.ika.or.kv.eq.ikb)then
        call diffuse(force_xv(:,:,:,kv),xvrho(:,:,:,kv), &
          1-ih,np+ih,1-ih,mp+ih,l,one)
      else
        call diffuse(force_xv(:,:,:,kv),xvrho(:,:,:,kv), &
          1-ih,np+ih,1-ih,mp+ih,l,rturbprandtl)
      endif
    enddo

  endif
     
end subroutine calcForceTurb
