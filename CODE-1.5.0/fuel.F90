!-----------------------------------------------------------------------
! burn computes reaction terms of fuel on the small timestep
!-----------------------------------------------------------------------
subroutine burnFuel(force_xv,force_xvfuel,xv,xvfuel, &
    il,iu,jl,ju,lls,nvp,nvf,nfuel)
  use gridlist_variables, only : cpwood,tambient,dts,prec
  use xvall, only : ikb,iO2,irho,irhof,irhow,isies,ipsiw,itemp,xvrho 
  use thermo_variables, only : cpwater,cp_gas,cv_gas,pr
  use fuel_variables, only : temps,rhoFuelInitial,rno,rnfuel,twvap, &
    tcrit,tfstep,sizeScale,hf,hwevap
  use radiation_variables, only : fsiesrad
  use linn_turb_variables, only : sc
  use constants, only : Pref
  use computePsi_function
  Implicit None

  ! Local Variables  
  integer,intent(in) :: il,iu,jl,ju,lls,nvp,nvf,nfuel
  real(prec),intent(in) :: xv(il:iu,jl:ju,lls,nvp), &
    xvfuel(nfuel,il:iu,jl:ju,lls,nvf)
  real(prec),intent(inout) :: force_xv(il:iu,jl:ju,lls,nvp), &
    force_xvfuel(nfuel,il:iu,jl:ju,lls,nvf)

  integer :: i,j,k,ift
  real(prec) :: sigmac,th 
  real(prec) :: percHydroRemaining,hydroFactor,tapsw
  real(prec) :: slambdaof
  real(prec) :: reactFuelGas,reactFuelSolid,fuelMassLoss
  real(prec) :: reactWaterSolid,waterMassLoss
  real(prec) :: energyGas,energySolid
  real(prec) :: ff,fw
  real(prec) :: rkctemp
  real(prec) :: psif,psiw,thetaSolid
  real(prec) :: cf 
  real(prec) :: hydroThresh=0.4
  real(prec) :: thetag=0.75 ! fraction of energy that is deposited in the gas phase
                      ! for the computation of thetasolid
  real(prec) :: cfhydro=.30*3. ! Coef. of mass reaction rates.  (found through numerical experiments with
                         ! grass) cfhydro of .8 seems a little low with code of 5/24/01, flat 
                         ! ground, no wind, rhof=1 kg/m^3 and 5% moisture
  real(prec) :: cfchar=.015*6. ! Coef. of mass reaction rates.  (found through 

  ! Executable Code
  do k=1,lls
    do j=jl,ju
      do i=il,iu
        do ift=1,nfuel
          ! Turbulent mixing
          rkctemp=.2*xvrho(i,j,k,ikb)
          ! the .2 factor is included as an extrapolation for turbulence at the c
          ! scale.  The .2 factor should be replaced by a function of sc and sb.
          sigmac=sc*0.5*sqrt(rkctemp)

          ! Flame heat, reaction extent, and length correlations
          psif=computePsi(temps(ift,i,j,k))
          percHydroRemaining=max(0.,(xvfuel(ift,i,j,k,irhof) &
            -hydroThresh*rhoFuelinitial(ift,i,j,k)) &
            /(rhoFuelinitial(ift,i,j,k)*(1.-hydroThresh)))
          cf=cfhydro*percHydroRemaining+cfchar*(1-percHydroRemaining)
          hydroFactor=exp(-rno/rnfuel*psif* &
            xvfuel(ift,i,j,k,irhof)/xv(i,j,k,iO2))
          thetaSolid=percHydroRemaining*(1.-thetag)*hydroFactor &
            +(1.-percHydroRemaining)*thetag

          ! Water Evaporation 
          tapsw=tambient+(twvap-tambient)*(tfstep-tambient) &
            /(tcrit-tambient) ! for 2 m res. step used in water ramp functions (K)
          th=2.*twvap-tapsw
          if(temps(ift,i,j,k).le.tapsw)then
            psiw=0.
          elseif(temps(ift,i,j,k).ge.th)then
            psiw=1.
          else
            psiw=(temps(ift,i,j,k)-tapsw)/(th-tapsw)
          endif
          psiw=max(xvfuel(ift,i,j,k,ipsiw),psiw)
          force_xvfuel(ift,i,j,k,ipsiw)=force_xvfuel(ift,i,j,k,ipsiw) &
            +(psiw-xvfuel(ift,i,j,k,ipsiw))/dts
          if(xvfuel(ift,i,j,k,ipsiw).lt.1)then
            fw=min(xvfuel(ift,i,j,k,irhow)/dts,xvfuel(ift,i,j,k,irhow) &
              *(psiw-xvfuel(ift,i,j,k,ipsiw)) &
              /(1.-xvfuel(ift,i,j,k,ipsiw))/dts)
          else
            fw=0.
          endif

          ! Fuel Consumption
          slambdaof=sum(xvfuel(:,i,j,k,irhof))*xv(i,j,k,iO2) &
            /(sum(xvfuel(:,i,j,k,irhof))/rnfuel+xv(i,j,k,iO2)/rno)**2
          ff=cf*xvfuel(ift,i,j,k,irhof)*xv(i,j,k,iO2)*sigmac*psif* &
            slambdaof/(100.*sizescale(ift,i,j,k)**2.)

          ! Energy Changes
          reactFuelGas=(1.-thetaSolid)*hf*ff
          reactFuelSolid=thetaSolid*hf*ff
          fuelMassLoss=tcrit*rnfuel*cpwood*ff
          reactWaterSolid=-hwevap*fw
          waterMassLoss=twvap*cpwater*fw

          energySolid=reactFuelSolid-fuelMassLoss &
            +reactWaterSolid-waterMassLoss
          energyGas=reactFuelGas+fuelMassLoss+waterMassLoss

          ! Forcings for fuel physical terms
          force_xvfuel(ift,i,j,k,irhof)=force_xvfuel(ift,i,j,k,irhof) &
            -ff*rnfuel
          force_xvfuel(ift,i,j,k,irhow)=force_xvfuel(ift,i,j,k,irhow)-fw
          force_xvfuel(ift,i,j,k,isies)=force_xvfuel(ift,i,j,k,isies)+ &
            (energySolid+fsiesrad(ift,i,j,k))/ &
            (xvfuel(ift,i,j,k,irhof)+xvfuel(ift,i,j,k,irhow))

          ! Forcings for gas physical terms
          force_xv(i,j,k,itemp)=force_xv(i,j,k,itemp)+ &
            energyGas/cp_gas(i,j,k)* &
            (Pref/pr(i,j,k))**(1.-cv_gas(i,j,k)/cp_gas(i,j,k))
          force_xv(i,j,k,iO2)=force_xv(i,j,k,iO2)-ff*rno
          force_xv(i,j,k,irho)=force_xv(i,j,k,irho)+ff*rnfuel+fw
        enddo
      enddo
    enddo
  enddo
 
end subroutine burnFuel
