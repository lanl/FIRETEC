!-----------------------------------------------------------------------
! dragm routine calculates changes and evolution of wind due to 
! drag effects
!-----------------------------------------------------------------------
subroutine dragm(force,xv,il,iu,jl,ju,lls,nvp)
  use gridlist_variables, only : wp
  use vegetation_gridlist_variables, only : rhoMicro,nfuel
  use xvall, only: xvrho,iuvel,ivvel,iwvel,xvfuel,irhof
  use fuel_variables, only : lfuel,sizeScale,actualFuelDepth
  use metric_variables, only : zs,z
  use gridsetup, only : np,mp
  use constants, only : pi,tolerance
  use linn_turb_variables, only : cd
  Implicit None

  ! Local Variables
  integer,intent(in) :: il,iu,jl,ju,lls,nvp
  real(wp),intent(inout) :: force(il:iu,jl:ju,lls,nvp)
  real(wp),intent(in) :: xv(il:iu,jl:ju,lls,nvp)

  integer :: ift,i,j,k
  real(wp) :: sp,av
  real(wp) :: zla2,zla3,rinterp
  real(wp) :: dragxt,dragyt,dragzt
  real(wp) :: uFuel,vFuel,wFuel
  real(wp) :: bendh,bendx,bendy,bendz
  real(wp),external :: zcart

  ! Executable Code
  do k=1,lfuel
    do j=1,mp
      do i=1,np
        dragxt=0.0_wp
        dragyt=0.0_wp
        dragzt=0.0_wp
        if(k.eq.1) then
          do ift=1,nfuel
            ! Fuel velocities, interprets velocity at fuel top based on 2nd and 3rd cell velocities
            zla3=zcart(z(3),i,j)-zs(i,j)
            zla2=zcart(z(2),i,j)-zs(i,j)
            rinterp=(actualFuelDepth(ift,i,j)-zla2)/(zla3-zla2)
            uFuel=(xvrho(i,j,2,iuvel)+rinterp* &
              (xvrho(i,j,3,iuvel)-xvrho(i,j,2,iuvel)))/2._wp
            vFuel=(xvrho(i,j,2,ivvel)+rinterp* &
              (xvrho(i,j,3,ivvel)-xvrho(i,j,2,ivvel)))/2._wp
            wFuel=(xvrho(i,j,2,iwvel)+rinterp* &
              (xvrho(i,j,3,iwvel)-xvrho(i,j,2,iwvel)))/2._wp
            sp=sqrt(uFuel**2+vFuel**2+wFuel**2)
            
            ! Bending grass correlations
            bendh=cos(atan((uFuel**2+vFuel**2)/(3._wp+wFuel**2)))
            bendx=1._wp-(1._wp-bendh)*(xvrho(i,j,k,iuvel)**2/ &
              (xvrho(i,j,k,iuvel)**2+xvrho(i,j,k,ivvel)**2+tolerance))
            bendy=1._wp-(1._wp-bendh)*(xvrho(i,j,k,ivvel)**2/ &
              (xvrho(i,j,k,iuvel)**2+xvrho(i,j,k,ivvel)**2+tolerance))
            bendz=sin(atan((xvrho(i,j,k,iuvel)**2+ &
              xvrho(i,j,k,ivvel)**2)/(3._wp+xvrho(i,j,k,iwvel)**2)))
            
            av=2._wp/pi/sizeScale(ift,i,j,k)* &
              xvfuel(ift,i,j,k,irhof)/rhoMicro !correct for cylinder if ss is radius
            dragxt=dragxt+0.5_wp*cd*av*sp*bendx/sqrt(2._wp)
            dragyt=dragyt+0.5_wp*cd*av*sp*bendy/sqrt(2._wp)
            dragzt=dragzt+0.5_wp*cd*av*sp*(1._wp/sqrt(2._wp)+ &
              (1._wp-1._wp/sqrt(2._wp))*bendz)
          enddo
        else
          sp=sqrt(xvrho(i,j,k,iuvel)**2+xvrho(i,j,k,ivvel)**2+ &
            xvrho(i,j,k,iwvel)**2)
          do ift=1,nfuel
            av=2._wp/pi/sizeScale(ift,i,j,k)* &
              xvfuel(ift,i,j,k,irhof)/rhoMicro !correct for cylinder if ss is radius
            
            dragxt=dragxt+0.5_wp*cd*av*sp
            dragyt=dragyt+0.5_wp*cd*av*sp
            dragzt=dragzt+0.5_wp*cd*av*sp
          enddo
        endif
        force(i,j,k,iuvel)=force(i,j,k,iuvel)- &
          dragxt*xv(i,j,k,iuvel)
        force(i,j,k,ivvel)=force(i,j,k,ivvel)- &
          dragyt*xv(i,j,k,ivvel)
        force(i,j,k,iwvel)=force(i,j,k,iwvel)- &
          dragzt*xv(i,j,k,iwvel)
      enddo !i
    enddo !j
  enddo !k
 
end subroutine dragm

!-----------------------------------------------------------------------
! fuelDragCorrection calculates a correction term for wind velocities
! and turbulent kinetic energy in the fuel layers
!-----------------------------------------------------------------------
subroutine fuelDragCorrection(i,j,nfuel,afd,rhoFuel,rhoWater,u,v,w, &
  u_corr,v_corr,w_corr,iflag) 
  use gridlist_variables, only : wp 
  use metric_variables, only : zs,zedge 
  use fuel_variables, only : kfindex,kfrank
  Implicit None

  ! Local Variables
  integer,intent(in) :: i,j,nfuel,iflag   
  real(wp),intent(in) :: u,v,w 
  real(wp),intent(in) :: afd(nfuel),rhoFuel(nfuel),rhoWater(nfuel)
  real(wp),intent(inout) :: u_corr(nfuel),v_corr(nfuel),w_corr(nfuel) 

  integer :: it,it2
  real(wp) :: afd0
  real(wp) :: layer_ms(nfuel)  
  real(wp) :: sum_layer_ms(nfuel)  
  real(wp) :: cell_ht,layer_ht(nfuel) 
  real(wp) :: uFuelBot,vFuelBot,wFuelBot,uFuelTop,vFuelTop,wFuelTop 
  real(wp) :: u_fuel_layer(nfuel),v_fuel_layer(nfuel)
  real(wp) :: w_fuel_layer(nfuel) 
  real(wp),external :: zcart

  ! Executable Code
  if(iflag.ne.1) then ! not bot cell; no need for correction 
    u_corr = u
    v_corr = v
    w_corr = w
    return
  endif

  ! actual mass in each fuel layer
  afd0=0
  do it=1,nfuel
    layer_ht(it)=afd(kfindex(i,j,it))-afd0
    layer_ms(it)=0._wp
    do it2=it,nfuel
      layer_ms(it)=layer_ms(it)+layer_ht(it)/afd(kfindex(i,j,it2))* &
        (rhoWater(kfindex(i,j,it2))+rhoFuel(kfindex(i,j,it2)))
    enddo
    sum_layer_ms(it)=sum(layer_ms(1:it))
    afd0=afd(kfindex(i,j,it))
  enddo

  ! Calculate velocity at top of fuel bed  
  ! The average of the entire extent of the cell must be u 
  cell_ht=zcart(zedge(2),i,j)-zcart(zedge(1),i,j)
  uFuelTop=u*cell_ht/(cell_ht-layer_ht(nfuel)/2._wp)
  vFuelTop=v*cell_ht/(cell_ht-layer_ht(nfuel)/2._wp)
  wFuelTop=w*cell_ht/(cell_ht-layer_ht(nfuel)/2._wp)

  ! Iterate down through the layers
  do it=nfuel,1,-1
    uFuelBot=uFuelTop*(1-layer_ms(it)/sum_layer_ms(it))
    u_fuel_layer(it)=(uFuelTop+uFuelBot)/2._wp
    uFuelTop=uFuelBot
    vFuelBot=vFuelTop*(1-layer_ms(it)/sum_layer_ms(it))
    v_fuel_layer(it)=(vFuelTop+vFuelBot)/2._wp
    vFuelTop=vFuelBot
    wFuelBot=wFuelTop*(1-layer_ms(it)/sum_layer_ms(it))
    w_fuel_layer(it)=(wFuelTop+wFuelBot)/2._wp
    wFuelTop=wFuelBot
  enddo

  ! Parse layers back into the fuel types
  do it=1,nfuel
    u_corr(it)=sum(u_fuel_layer(1:kfrank(i,j,it)) &
      *layer_ht(1:kfrank(i,j,it)))/afd(it)
    v_corr(it)=sum(v_fuel_layer(1:kfrank(i,j,it)) &
      *layer_ht(1:kfrank(i,j,it)))/afd(it)
    w_corr(it)=sum(w_fuel_layer(1:kfrank(i,j,it)) &
      *layer_ht(1:kfrank(i,j,it)))/afd(it)
  enddo

end subroutine fuelDragCorrection 
