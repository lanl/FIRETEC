!-----------------------------------------------------------------------
! dragm routine calculates changes and evolution of wind due to 
! drag effects
!-----------------------------------------------------------------------
subroutine dragm(force,xv,il,iu,jl,ju,lls,nvp)
  use gridlist_variables, only : nfuel,rhoMicro,prec
  use xvall, only: xvrho,iuvel,ivvel,iwvel,xvfuel,irhof 
  use fuel_variables, only : lfuel,sizeScale,actualFuelDepth
  use metric_variables, only : zs,z
  use gridsetup, only : np,mp
  use constants, only : pi
  use linn_turb_variables, only : cd
  use zcart_function
  Implicit None

  ! Local Variables
  integer,intent(in) :: il,iu,jl,ju,lls,nvp
  real(prec),intent(inout) :: force(il:iu,jl:ju,lls,nvp)
  real(prec),intent(in) :: xv(il:iu,jl:ju,lls,nvp)

  integer :: ift,i,j,k
  real(prec) :: sp,av
  real(prec) :: zla2,zla3,rinterp
  real(prec) :: dragxt,dragyt,dragzt
  real(prec) :: uFuel,vFuel,wFuel
  real(prec) :: bendh,bendx,bendy,bendz

  ! Executable Code
  do k=1,lfuel
    do j=1,mp
      do i=1,np
        dragxt=0.
        dragyt=0.
        dragzt=0.
        if(k.eq.1) then
          do ift=1,nfuel
            ! Fuel velocities, interprets velocity at fuel top based on 2nd and 3rd cell velocities
            zla3=zcart(z(3),i,j)-zs(i,j)
            zla2=zcart(z(2),i,j)-zs(i,j)
            rinterp=(actualFuelDepth(ift,i,j)-zla2)/(zla3-zla2)
            uFuel=(xvrho(i,j,2,iuvel)+rinterp* &
              (xvrho(i,j,3,iuvel)-xvrho(i,j,2,iuvel)))/2.
            vFuel=(xvrho(i,j,2,ivvel)+rinterp* &
              (xvrho(i,j,3,ivvel)-xvrho(i,j,2,ivvel)))/2.
            wFuel=(xvrho(i,j,2,iwvel)+rinterp* &
              (xvrho(i,j,3,iwvel)-xvrho(i,j,2,iwvel)))/2.
            sp=sqrt(uFuel**2.+vFuel**2.+wFuel**2.)
            
            ! Bending grass correlations
            bendh=cos(atan((uFuel**2.+vFuel**2.)/(3.+wFuel**2.)))
#if DBL_EC
            bendx=1.-(1.-bendh)*(xvrho(i,j,k,iuvel)**2./ &
              (xvrho(i,j,k,iuvel)**2.+xvrho(i,j,k,ivvel)**2.+1d-6))
            bendy=1.-(1.-bendh)*(xvrho(i,j,k,ivvel)**2./ &
              (xvrho(i,j,k,iuvel)**2+xvrho(i,j,k,ivvel)**2+1d-6))
#else
            bendx=1.-(1.-bendh)*(xvrho(i,j,k,iuvel)**2./ &
              (xvrho(i,j,k,iuvel)**2.+xvrho(i,j,k,ivvel)**2.+1e-6))
            bendy=1.-(1.-bendh)*(xvrho(i,j,k,ivvel)**2./ &
              (xvrho(i,j,k,iuvel)**2+xvrho(i,j,k,ivvel)**2+1e-6))
#endif
            bendz=sin(atan((xvrho(i,j,k,iuvel)**2.+ &
              xvrho(i,j,k,ivvel)**2.)/(3.+xvrho(i,j,k,iwvel)**2.)))
            
            av=2./pi/sizeScale(ift,i,j,k)* &
              xvfuel(ift,i,j,k,irhof)/rhoMicro !correct for cylinder if ss is radius
            dragxt=dragxt+.5*cd*av*sp*bendx/sqrt(2.)
            dragyt=dragyt+.5*cd*av*sp*bendy/sqrt(2.)
            dragzt=dragzt+.5*cd*av*sp*(1./sqrt(2.)+ &
              (1.-1./sqrt(2.))*bendz)
          enddo
        else
          sp=sqrt(xvrho(i,j,k,iuvel)**2.+xvrho(i,j,k,ivvel)**2.+ &
            xvrho(i,j,k,iwvel)**2.)
          do ift=1,nfuel
            av=2./pi/sizeScale(ift,i,j,k)* &
              xvfuel(ift,i,j,k,irhof)/rhoMicro !correct for cylinder if ss is radius
            
            dragxt=dragxt+.5*cd*av*sp
            dragyt=dragyt+.5*cd*av*sp
            dragzt=dragzt+.5*cd*av*sp
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
  use gridlist_variables, only : prec 
  use metric_variables, only : zs,z 
  use zcart_function  
  Implicit None

  ! Local Variables
  integer,intent(in) :: i,j,nfuel,iflag   
  real(prec),intent(in) :: u(3),v(3),w(3) 
  real(prec),intent(in) :: afd(nfuel),rhoFuel(nfuel),rhoWater(nfuel)
  real(prec),intent(inout) :: u_corr(nfuel),v_corr(nfuel),w_corr(nfuel) 

  integer :: it
  integer :: itfindex(nfuel) 
  real(prec) :: zla2,zla3,rinterp 
  real(prec) :: afdtmp(nfuel)  
  real(prec) :: layer_ms(nfuel)  
  real(prec) :: sum_layer_ms(nfuel)  
  real(prec) :: layer_ht(nfuel) 
  real(prec) :: uFuelBot,vFuelBot,wFuelBot,uFuelTop,vFuelTop,wFuelTop 
  real(prec) :: u_fuel_layer(nfuel),v_fuel_layer(nfuel)
  real(prec) :: w_fuel_layer(nfuel) 

  ! Executable Code
  if(iflag.ne.1) then ! not bot cell; no need for correction 
    u_corr = u(1)
    v_corr = v(1)
    w_corr = w(1)
    return
  endif

  ! Sort individual fuel zones by afd
  afdtmp=afd(:)
  do it=1,nfuel
    itfindex(minloc(afdtmp,dim=1))=it 
    afdtmp(minloc(afdtmp,dim=1))=maxval(afd)+1
  enddo

  ! actual mass in each fuel layer   
  layer_ht(1) = afd(itfindex(1))
  layer_ms(1) = sum((rhoWater+rhoFuel)*layer_ht(1)/afd)
  sum_layer_ms(1) = layer_ms(1)
  if(nfuel.gt.1)then 
    do it=2,nfuel
      layer_ht(it)=afd(itfindex(it))-afd(itfindex(it-1))
      layer_ms(it)= &
        sum(rhoWater(itfindex(it:nfuel))+rhoFuel(itfindex(it:nfuel))) &
        *layer_ht(it)
      sum_layer_ms(it) = sum_layer_ms(it-1) + layer_ms(it) 
    enddo
  endif

  ! Calculate velocity at top of fuel bed  
  ! Fuel velocities, interprets velocity at fuel top based on 2nd and 3rd cell velocities 
  zla3=zcart(z(3),i,j)-zs(i,j) 
  zla2=zcart(z(2),i,j)-zs(i,j) 
  rinterp=(afd(itfindex(nfuel))-zla2)/(zla3-zla2) 
  uFuelTop=(u(2)+rinterp*(u(3)-u(2)))/2. 
  vFuelTop=(v(2)+rinterp*(v(3)-v(2)))/2.
  wFuelTop=(w(2)+rinterp*(w(3)-w(2)))/2.

  ! Iterate down through the layers
  do it=nfuel,1,-1
    uFuelBot=uFuelTop*(1-layer_ms(it)/sum_layer_ms(it))
    u_fuel_layer(it)=(uFuelTop+uFuelBot)/2.
    uFuelTop=uFuelBot
    vFuelBot=vFuelTop*(1-layer_ms(it)/sum_layer_ms(it))
    v_fuel_layer(it)=(vFuelTop+vFuelBot)/2.
    vFuelTop=vFuelBot
    wFuelBot=wFuelTop*(1-layer_ms(it)/sum_layer_ms(it))
    w_fuel_layer(it)=(wFuelTop+wFuelBot)/2.
    wFuelTop=wFuelBot
  enddo

  ! Parse layers back into the fuel types
  do it=1,nfuel
    u_corr(it)= &
      sum(u_fuel_layer(1:itfindex(it))*layer_ht(1:itfindex(it)))/afd(it)
    v_corr(it)= &
      sum(v_fuel_layer(1:itfindex(it))*layer_ht(1:itfindex(it)))/afd(it)
    w_corr(it)= &
      sum(w_fuel_layer(1:itfindex(it))*layer_ht(1:itfindex(it)))/afd(it)
  enddo

end subroutine fuelDragCorrection 
