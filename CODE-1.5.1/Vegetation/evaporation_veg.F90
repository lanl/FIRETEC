!-----------------------------------------------------------------------
! evapPDF computes evaporation rates due to rapid heating on the small 
! timestep
!-----------------------------------------------------------------------
subroutine evapPDF(force_xv,force_xvfuel,xvfuel,ilh,iuh,jlh,juh,ku, &
    il,iu,jl,ju,lls,nvp,nvf,nfuel)
  use gridlist_variables, only : tambient,wp
  use gridsetup, only : dts
  use xvall, only : irhof,irhow,isies,ipsiw,irho,itheta,iH2O
  use thermo_variables, only : cpwater,cp_gas,cv_gas,pr
  use fuel_variables, only : Temps,twvap,tcrit,tfstep,hwevap
  use constants, only : Pref
  Implicit None

  ! Local Variables  
  integer,intent(in) :: ilh,iuh,jlh,juh,ku,il,iu,jl,ju,lls,nvp,nvf,nfuel
  real(wp),intent(in) :: xvfuel(nfuel,il:iu,jl:ju,lls,nvf)
  real(wp),intent(inout) :: force_xv(ilh:iuh,jlh:juh,ku,nvp), &
    force_xvfuel(nfuel,il:iu,jl:ju,lls,nvf)

  integer :: i,j,k,it
  real(wp) :: tapsw,th
  real(wp) :: reactWaterSolid,waterMassLoss
  real(wp) :: energyGas,energySolid
  real(wp) :: fw
  real(wp) :: psiw

  ! Executable Code
  tapsw=tambient+(twvap-tambient)*(tfstep-tambient) &
    /(tcrit-tambient) ! for 2 m res. step used in water ramp functions (K)
  th=2.0_wp*twvap-tapsw
  do k=1,lls
    do j=jl,ju
      do i=il,iu
        do it=1,nfuel
          ! Water Evaporation 
          if(Temps(it,i,j,k).le.tapsw)then
            psiw=0._wp
          elseif(Temps(it,i,j,k).ge.th)then
            psiw=1._wp
          else
            psiw=(Temps(it,i,j,k)-tapsw)/(th-tapsw)
          endif
          psiw=max(xvfuel(it,i,j,k,ipsiw),psiw)
          force_xvfuel(it,i,j,k,ipsiw)=force_xvfuel(it,i,j,k,ipsiw) &
            +(psiw-xvfuel(it,i,j,k,ipsiw))/dts
          if(xvfuel(it,i,j,k,ipsiw).lt.1)then
            fw=min(xvfuel(it,i,j,k,irhow)/dts,xvfuel(it,i,j,k,irhow) &
              *(psiw-xvfuel(it,i,j,k,ipsiw)) &
              /(1.0_wp-xvfuel(it,i,j,k,ipsiw))/dts)
          else
            fw=0._wp
          endif
          
          ! Energy Changes
          reactWaterSolid=-hwevap*fw
          waterMassLoss=twvap*cpwater*fw

          energySolid=reactWaterSolid-waterMassLoss
          energyGas=waterMassLoss

          ! Forcings for fuel physical terms
          force_xvfuel(it,i,j,k,irhow)=force_xvfuel(it,i,j,k,irhow)-fw
          force_xvfuel(it,i,j,k,isies)=force_xvfuel(it,i,j,k,isies) &
            +energySolid &
            /(xvfuel(it,i,j,k,irhof)+xvfuel(it,i,j,k,irhow))

          ! Forcings for gas physical terms
          force_xv(i,j,k,itheta)=force_xv(i,j,k,itheta)+ &
            energyGas/cp_gas(i,j,k)* &
            (Pref/pr(i,j,k))**(1.0_wp-cv_gas(i,j,k)/cp_gas(i,j,k))
          force_xv(i,j,k,irho)=force_xv(i,j,k,irho)+fw
          if(iH2O.ne.1) force_xv(i,j,k,iH2O)=force_xv(i,j,k,iH2O)+fw
        enddo
      enddo
    enddo
  enddo
 
end subroutine evapPDF

!-----------------------------------------------------------------------
! evapBulk computes evaporation rates due to slow heating on the small 
! timestep
!-----------------------------------------------------------------------
subroutine evapBulk(force_xv,force_xvfuel,xv,xvrho,xvfuel, &
    ilh,iuh,jlh,juh,kl,ku,il,iu,jl,ju,lls,nvp,nvf,nfuel)
  use gridlist_variables, only : wp
  use vegetation_gridlist_variables, only : rhoMicro
  use xvall, only : irhof,irhow,isies,irho,iH2O,itheta,ika,ikb, &
    iuvel,ivvel,iwvel
  use thermo_variables, only : cpwater,cp_gas,cv_gas,pr
  use fuel_variables, only : Temps,twvap,hwevap,sizeScale, &
    actualFuelDepth
  use constants, only : Pref,pi
  Implicit None

  ! Local Variables  
  integer,intent(in) :: ilh,iuh,jlh,juh,kl,ku,nvp
  integer,intent(in) :: il,iu,jl,ju,lls,nvf
  integer,intent(in) :: nfuel
  real(wp),intent(in) :: xv(ilh:iuh,jlh:juh,ku,nvp), &
    xvrho(ilh:iuh,jlh:juh,kl:ku,nvp-1),xvfuel(nfuel,il:iu,jl:ju,lls,nvf)
  real(wp),intent(inout) :: force_xv(ilh:iuh,jlh:juh,ku,nvp), &
    force_xvfuel(nfuel,il:iu,jl:ju,lls,nvf)

  integer :: i,j,k,it
  real(wp) :: gasVapPres,surfVapPres,satVapPres ! Vapor pressure (Pa)
  real(wp) :: humid_regime,RHeqv,RHa,RHb,RHc,TempFah,fuelMC
  real(wp) :: rktemp,sp
  real(wp) :: re,le,h,c,m ! Reynolds & Lewis numbers
  real(wp) :: evapR,airR ! Evaporation resistance (m/s)
  real(wp) :: av,fw
  real(wp) :: reactWaterSolid,waterMassLoss
  real(wp) :: energyGas,energySolid
  real(wp) :: kinViscAir=0.00005_wp ! Kinematic viscosity of air (m2/s)
  real(wp) :: waterDiffAir=0.0000235_wp  ! Water diffusivity in Air (m2/s)
  real(wp) :: thermDiffAir=0.000024_wp ! thermal diffusivity of air (m2/s)
  real(wp) :: thermCondAir=0.0338_wp
  real(wp) :: u_corr(nfuel),v_corr(nfuel),w_corr(nfuel)

  ! Executable Code
  le=thermDiffAir/waterDiffAir ! Lewis number
  do k=1,lls
    do j=jl,ju
      do i=il,iu
        ! computation of sp (velocity in fuel + turbulence)
        call fuelDragCorrection(i,j,nfuel, &
          actualFuelDepth(:,i,j), &
          xvfuel(:,i,j,k,irhof),xvfuel(:,i,j,k,irhow), &
          xvrho(i,j,k,iuvel),xvrho(i,j,k,ivvel), &
          xvrho(i,j,k,iwvel),u_corr,v_corr,w_corr,k)  
        
        ! Evaporation based on vapor pressure gradients
        if(iH2O.ne.1)then
          gasVapPres=pr(i,j,k)*xv(i,j,k,iH2O)/xv(i,j,k,irho)
        else
          gasVapPres=0.0_wp
        endif
        do it=1,nfuel
          if (Temps(it,i,j,k).gt.335)then
            cycle
          else
            ! Calculate vapor pressure at the surface of the fuels by computing the humidity percent at which
            ! the fuel would be equilibriated (Fosberg, Rothermel, and Andrews, 1981) and multiplying that 
            ! fraction by the saturated vapor pressure at the fuel surface (Bolton 1980)
            satVapPres=611.2_wp*exp(17.67_wp* &
              (Temps(it,i,j,k)-273.15_wp)/(Temps(it,i,j,k)-29.65_wp)) ! Saturated vapor pressure, Bolton, 1980
            humid_regime=gasVapPres/satVapPres
            fuelMC=xvfuel(it,i,j,k,irhow)/xvfuel(it,i,j,k,irhof)*100
            TempFah=(Temps(it,i,j,k)-273.15_wp)*9._wp/5._wp+32._wp
            if(humid_regime.lt.0.1_wp)then
              RHeqv=(fuelMC-0.03229_wp)/ &
                (0.281073_wp-0.0005578_wp*TempFah)/100
            elseif(humid_regime.lt.0.5_wp)then
              RHeqv=(fuelMC-2.22749_wp+0.014784_wp*TempFah)/ &
                0.160107_wp/100
            else
              RHa=0.005565_wp
              RHb=-0.00063_wp*TempFah-0.483199_wp
              RHc=21.0606_wp-fuelMC
              RHeqv=(-RHb+sqrt(RHb**2-4._wp*RHa*RHc))/(2._wp*RHa)/100
              if(RHeqv.gt.1._wp) &
                RHeqv=(-RHb-sqrt(RHb**2-4._wp*RHa*RHc))/(2._wp*RHa)/100
            endif
            surfVapPres=satVapPres*max(0._wp,RHeqv)
            
            rktemp=xvrho(i,j,k,ika)+1.2_wp*xvrho(i,j,k,ikb) 
            sp=sqrt(rktemp)+ &
              sqrt(u_corr(it)**2+v_corr(it)**2+w_corr(it)**2)
        
            ! computation of h and av
            re=sizescale(it,i,j,k)*sp/kinViscAir
            if(re.lt.40._wp)then
              c=0.911_wp
              m=0.385_wp
            else
              c=0.683_wp
              m=0.466_wp
            endif
            h=c/pi*re**m*thermCondAir/sizeScale(it,i,j,k) ! AJJ optical average over full range of rotation is 2/pi
            av=2._wp*(xvfuel(it,i,j,k,irhof)/rhoMicro) &
              /sizeScale(it,i,j,k)
            
            airR=xv(i,j,k,irho)*cp_gas(i,j,k)*le**(2._wp/3._wp)/h
            !fuelR =
            evapR=1._wp/airR!+fuelR) ! Evaporation resistance considers only air at this point
            
            fw=0.622_wp*av*xv(i,j,k,irho)*evapR &
              *(surfVapPres-gasVapPres)/pr(i,j,k) 

            ! Energy Changes
            reactWaterSolid=-hwevap*fw
            waterMassLoss=twvap*cpwater*fw

            energySolid=reactWaterSolid-waterMassLoss
            energyGas=waterMassLoss

            ! Forcings for fuel physical terms
            force_xvfuel(it,i,j,k,irhow)=force_xvfuel(it,i,j,k,irhow)-fw
            force_xvfuel(it,i,j,k,isies)=force_xvfuel(it,i,j,k,isies)+ &
             energySolid/(xvfuel(it,i,j,k,irhof)+xvfuel(it,i,j,k,irhow))

            ! Forcings for gas physical terms
            force_xv(i,j,k,itheta)=force_xv(i,j,k,itheta)+ &
              energyGas/cp_gas(i,j,k)* &
              (Pref/pr(i,j,k))**(1._wp-cv_gas(i,j,k)/cp_gas(i,j,k))
            force_xv(i,j,k,irho)=force_xv(i,j,k,irho)+fw
            if(iH2O.ne.1) force_xv(i,j,k,iH2O)=force_xv(i,j,k,iH2O)+fw
          endif
        enddo
      enddo
    enddo
  enddo

end subroutine evapBulk
