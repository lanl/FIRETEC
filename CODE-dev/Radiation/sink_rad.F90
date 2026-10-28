!-----------------------------------------------------------------------
! radiationSink estimates a heat loss due to radiation based on 
! comparing local temperatures against environmental temperatures
!-----------------------------------------------------------------------
subroutine radiationSink(force,force_fuel,rhoFuel,rhoWater,tempg,temps,&
    ilh,iuh,jlh,juh,il,iu,jl,ju,lls,lfuel,nfuel)
  use gridlist_variables, only : dx,dy,wp
  use vegetation_gridlist_variables, only : rhoMicro
  use fuel_variables, only : sizeScale
  use thermo_variables, only : pre,pr,cp_gas,cv_gas
  use constants, only : sigma,pref
  use xvall, only : xe,itheta,irho
  Implicit none
  
  ! Local Variables  
  integer,intent(in) :: ilh,iuh,jlh,juh,il,iu,jl,ju,lls,lfuel,nfuel
  real(wp),intent(in) :: tempg(ilh:iuh,jlh:juh,lls)
  real(wp),intent(in) :: temps(nfuel,ilh:iuh,jlh:juh,lfuel)
  real(wp),intent(in) :: rhoFuel(nfuel,il:iu,jl:ju,lfuel)
  real(wp),intent(in) :: rhoWater(nfuel,il:iu,jl:ju,lfuel)
  real(wp),intent(inout) :: force(ilh:iuh,jlh:juh,lls)
  real(wp),intent(inout) :: force_fuel(nfuel,il:iu,jl:ju,lfuel)

  integer :: i,j,k,it 
  real(wp) :: tenv,av,firad_gas,firad_solid

  ! Executable Code
  do k=1,lls
    do j=jl,ju
      do i=il,iu
        tenv=xe(i,j,k,itheta)/xe(i,j,k,irho) &
          *(pre(i,j,k)/pref)**(1.0_wp-cv_gas(i,j,k)/cp_gas(i,j,k))
        firad_gas=-2.0_wp*sigma/sqrt(dx*dy)*(tempg(i,j,k)**4-tenv**4)
        force(i,j,k)=force(i,j,k)+firad_gas/cp_gas(i,j,k)* &
          (pref/pr(i,j,k))**(1.0_wp-cv_gas(i,j,k)/cp_gas(i,j,k))
      enddo
    enddo
  enddo
  do k=1,lfuel
    do j=jl,ju
      do i=il,iu
        tenv=xe(i,j,k,itheta)/xe(i,j,k,irho) &
          *(pre(i,j,k)/pref)**(1.0_wp-cv_gas(i,j,k)/cp_gas(i,j,k))
        do it=1,nfuel
          av=2.0_wp*(rhoFuel(it,i,j,k)/rhoMicro)/sizeScale(it,i,j,k)
          firad_solid=-0.95_wp*sigma*av*(temps(it,i,j,k)**4-tenv**4)
          force_fuel(it,i,j,k)=force_fuel(it,i,j,k)+ &
            firad_solid/(rhoFuel(it,i,j,k)+rhoWater(it,i,j,k))
        enddo
      enddo
    enddo
  enddo
        
end subroutine radiationSink
        
