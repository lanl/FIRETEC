!----------------------------------------------------------------
! the subroutine convection() computes convht
!----------------------------------------------------------------
subroutine convection(ift,i,j,k)
  use gridlist_variables, only : rhoMicro
  use xvall, only : xv,ika,ikb,iuvel,ivvel,iwvel,irho,itemp,nv
  use thermo_variables, only : pr,rg_over_prrcp_gas,cp_over_cv_gas, &
    rg_over_cp_gas,tempg
  use fuel_variables, only : fcorr,sizeScale,rhoFuel,convht,temps
  use constants, only : pi
  use turba
  Implicit None

  ! Local Variables
  integer,intent(in) :: i,j,k,ift
  real :: rktemp
  real :: sp,re,h,av
  real :: thermCondAir=33.8e-3

  ! Executable Code
  ! computation of sp (velocity in fuel + turbulence)
  rktemp=(xv(i,j,k,ika)+1.2*xv(i,j,k,ikb))*fcorr(ift,2)/xv(i,j,k,irho)
  sp=sqrt(rktemp)+sqrt(xv(i,j,k,iuvel)**2.+xv(i,j,k,ivvel)**2.+xv(i,j,k,iwvel)**2.)*fcorr(ift,1)/xv(i,j,k,irho)

  ! computation of h and av
  re=sizescale(ift,i,j,k)*sp/2.e-05
  h=0.683/pi*re**0.466*thermCondAir/sizescale(ift,i,j,k) ! AJJ optical average over full range of rotation is 2/pi
  av=2.*(rhoFuel(ift,i,j,k)/rhoMicro)/sizescale(ift,i,j,k)
  ! computation of tempg
  pr(i,j,k)=(xv(i,j,k,itemp)*rg_over_prrcp_gas)**cp_over_cv_gas
  tempg(i,j,k)=xv(i,j,k,itemp)/xv(i,j,k,nv)*(pr(i,j,k)*1.e-5)**rg_over_cp_gas

  ! computation of convht
  convht(ift,i,j,k)=h*av*(tempg(i,j,k)-temps(ift,i,j,k))    !wss

end subroutine convection
