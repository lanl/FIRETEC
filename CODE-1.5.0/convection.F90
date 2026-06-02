!-----------------------------------------------------------------------
! the subroutine convection computes convective heat exchange forcings
! between solid and gas phases
!-----------------------------------------------------------------------
subroutine convection(force_theta,force_sies,xvrho,xvfuel, &
    il,iu,jl,ju,lls,nvp,nvf,nfuel)
  use gridlist_variables, only : rhoMicro,prec
  use xvall, only : irhof,irhow,ika,ikb,iuvel,ivvel,iwvel
  use thermo_variables, only : tempg
  use fuel_variables, only : sizeScale,temps,actualFuelDepth
  use constants, only : pi,Pref
  use thermo_variables, only : pr,cp_gas,cv_gas
  Implicit None

  ! Local Variables
  integer,intent(in) :: il,iu,jl,ju,lls,nvp,nvf,nfuel
  real(prec),intent(in) :: xvrho(il:iu,jl:ju,lls,nvp), &
    xvfuel(nfuel,il:iu,jl:ju,lls,nvf)
  real(prec),intent(inout) :: force_theta(il:iu,jl:ju,lls), &
    force_sies(nfuel,il:iu,jl:ju,lls)
  
  integer :: i,j,k,ift
  real(prec) :: rktemp
  real(prec) :: sp,re,h,av
  real(prec) :: thermCondAir=33.8e-3
  real(prec) :: convht
  real(prec) :: u_corr(nfuel),v_corr(nfuel),w_corr(nfuel)

  ! Executable Code
  do k=1,lls
    do j=jl,ju
      do i=il,iu
        call fuelDragCorrection(i,j,nfuel,actualFuelDepth(:,i,j), &
          xvfuel(:,i,j,k,irhof),xvfuel(:,i,j,k,irhow), &
          xvrho(i,j,k:k+2,iuvel),xvrho(i,j,k:k+2,ivvel), &
          xvrho(i,j,k:k+2,iwvel),u_corr,v_corr,w_corr,k)  
          
        do ift=1,nfuel
          ! computation of sp (velocity in fuel + turbulence) 
          rktemp=(xvrho(i,j,k,ika)+1.2*xvrho(i,j,k,ikb)) 
          sp=sqrt(rktemp)+ &
            sqrt(u_corr(ift)**2+v_corr(ift)**2+w_corr(ift)**2)

          ! computation of h and av
          re=sizescale(ift,i,j,k)*sp/2.e-05
          h=0.683/pi*re**0.466*thermCondAir/sizeScale(ift,i,j,k) ! AJJ optical average over full range of rotation is 2/pi
          av=2.*(xvfuel(ift,i,j,k,irhof)/rhoMicro)/sizeScale(ift,i,j,k)

          ! computation of convht forcings
          convht=h*av*(tempg(i,j,k)-temps(ift,i,j,k))    !wss
          
          force_sies(ift,i,j,k)=force_sies(ift,i,j,k) &
            +convht/(xvfuel(ift,i,j,k,irhof)+xvfuel(ift,i,j,k,irhow))
          force_theta(i,j,k)=force_theta(i,j,k) &
            -convht/cp_gas(i,j,k)* &
            (Pref/pr(i,j,k))**(1.-cv_gas(i,j,k)/cp_gas(i,j,k))
        enddo
      enddo
    enddo
  enddo

end subroutine convection
