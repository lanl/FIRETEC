!-----------------------------------------------------------------------
! the subroutine convection computes convective heat exchange forcings
! between solid and gas phases
!-----------------------------------------------------------------------
subroutine convection(force_theta,force_sies,xvrho,xvfuel, &
    ilh,iuh,jlh,juh,il,iu,jl,ju,kl,ku,lfuel,nvp,nvf,nfuel)
  use gridlist_variables, only : wp
  use vegetation_gridlist_variables, only : rhoMicro
  use xvall, only : irhof,irhow,ika,ikb,iuvel,ivvel,iwvel
  use thermo_variables, only : tempg
  use fuel_variables, only : sizeScale,temps,actualFuelDepth
  use constants, only : pi,Pref
  use thermo_variables, only : pr,cp_gas,cv_gas
  Implicit None

  ! Local Variables
  integer,intent(in) :: ilh,iuh,jlh,juh,il,iu,jl,ju,kl,ku,lfuel, &
    nvp,nvf,nfuel
  real(wp),intent(in) :: xvrho(ilh:iuh,jlh:juh,kl:ku,nvp-1), &
    xvfuel(nfuel,il:iu,jl:ju,lfuel,nvf)
  real(wp),intent(inout) :: force_theta(ilh:iuh,jlh:juh,ku), &
    force_sies(nfuel,il:iu,jl:ju,lfuel)
  
  integer :: i,j,k,ift
  real(wp) :: rktemp
  real(wp) :: sp,re,h,av
  real(wp) :: thermCondAir=0.0338_wp
  real(wp) :: kinViscAir=0.000015_wp
  real(wp) :: convht
  real(wp) :: u_corr(nfuel),v_corr(nfuel),w_corr(nfuel)

  ! Executable Code
  do k=1,lfuel
    do j=jl,ju
      do i=il,iu
        call fuelDragCorrection(i,j,nfuel,actualFuelDepth(:,i,j), &
          xvfuel(:,i,j,k,irhof),xvfuel(:,i,j,k,irhow), &
          xvrho(i,j,k,iuvel),xvrho(i,j,k,ivvel),xvrho(i,j,k,iwvel), &
          u_corr,v_corr,w_corr,k)  
          
        do ift=1,nfuel
          ! computation of sp (velocity in fuel + turbulence) 
          rktemp=xvrho(i,j,k,ika)+1.2_wp*xvrho(i,j,k,ikb)
          sp=sqrt(rktemp)+ &
            sqrt(u_corr(ift)**2+v_corr(ift)**2+w_corr(ift)**2)

          ! computation of h and av TODO account for free convection
          re=sizescale(ift,i,j,k)*sp/kinViscAir
          h=0.683_wp/pi*re**0.466_wp*thermCondAir/sizeScale(ift,i,j,k) ! AJJ optical average over full range of rotation is 1/pi
          av=2._wp*(xvfuel(ift,i,j,k,irhof)/rhoMicro) &
            /sizeScale(ift,i,j,k)

          ! computation of convht forcings
          convht=h*av*(tempg(i,j,k)-temps(ift,i,j,k))    ! wss
          
          force_sies(ift,i,j,k)=force_sies(ift,i,j,k) &
            +convht/(xvfuel(ift,i,j,k,irhof)+xvfuel(ift,i,j,k,irhow))
          force_theta(i,j,k)=force_theta(i,j,k) &
            -convht/cp_gas(i,j,k)* &
            (Pref/pr(i,j,k))**(1._wp-cv_gas(i,j,k)/cp_gas(i,j,k))
        enddo
      enddo
    enddo
  enddo

end subroutine convection
