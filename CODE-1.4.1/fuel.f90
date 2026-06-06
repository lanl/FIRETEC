!----------------------------------------------------------------
! burn computes reaction terms of fuel on the small timestep
!----------------------------------------------------------------
subroutine fuel(ift,i,j,k,rneteng,ff_sum,fw_sum,rhoFuelold)
  use gridlist_variables, only : iemissions,idiffsies,nfuel, &
    cpwood,tambient
  use xvall, only : xv,xvrho,ikb,iO2,irho
  use thermo_variables, only : cpwater
  use fuel_variables, only : temps,rhoFuel,rhoWater,rhoFuelInitial, &
    rno,rnfuel,twvap,tcrit,tfstep,psiwmax,sizeScale,sies, &
    convht,siesdiff,fcorr
  use radiation_variables, only : fsiesrad
  use turb_variables, only : sc
  use gridsetup, only : dtp
  Implicit None

  ! Local Variables  
  integer,intent(in) :: i,j,k,ift
  real,intent(inout) :: rneteng ! defined for gas phase
  real,intent(inout) :: ff_sum 
  real,intent(inout) :: fw_sum
  real,intent(in) :: rhoFuelold
  
  real :: sigmac, th 
  real :: percHydroRemaining,hydroFactor,tapsw
  real :: slambdaof,psiwmaxold
  real :: rhoWaterold,rhos
  real :: reactht,rmassloss,waterevpht,rwaterlossht
  real :: ff,fw
  real :: rkctemp
  real :: frhoFuel,frhosies,frhoWater,psif,psiw,thetaSolid
  real :: cf 
  real :: hydroThresh=0.4
  real :: thetag=0.75 ! fraction of energy that is deposited in the gas phase
                      ! for the computation of thetasolid
  real :: cfhydro=.30*3. ! Coef. of mass reaction rates.  (found through numerical experiments with
                         ! grass) cfhydro of .8 seems a little low with code of 5/24/01, flat 
                         ! ground, no wind, rhof=1 kg/m^3 and 5% moisture
  real :: cfchar=.015*6. ! Coef. of mass reaction rates.  (found through 
  real :: hf=8913.48e3   ! heat of reaction for simple wood  (J/Kg of products)
  real :: hwevap=2.257e6 ! energy needed for water evap at 373 K (Moran)(J/Kg)
  real,external :: computePsi      

  ! Executable Code
  ! Turbulent mixing
  rkctemp=.2*xv(i,j,k,ikb)*fcorr(ift,2)/xv(i,j,k,irho)
  ! the .2 factor is included as an extrapolation for turbulence at the c
  ! scale.  The .2 factor should be replaced by a function of sc and sb.
  sigmac=sc*0.5*sqrt(rkctemp)

  ! Flame heat, reaction extent, and length correlations
  psif=computePsi(temps(ift,i,j,k))
  percHydroRemaining=max(0.,(rhoFuel(ift,i,j,k)-hydroThresh*rhoFuelinitial(ift,i,j,k)) &
    /(rhoFuelinitial(ift,i,j,k)*(1.-hydroThresh)))
  cf=cfhydro*percHydroRemaining+cfchar*(1-percHydroRemaining)
  hydroFactor=exp(-rno/rnfuel*psif*rhoFuelold/xv(i,j,k,iO2))
  thetaSolid=percHydroRemaining*(1.-thetag)*hydroFactor+(1.-percHydroRemaining)*thetag

  ! Water Evaporation 
  tapsw=tambient+(twvap-tambient)*(tfstep-tambient)/(tcrit-tambient) ! for 2 m res. step used in water ramp functions (K)
  th=2.*twvap-tapsw
  if(temps(ift,i,j,k).le.tapsw)then
    psiw=0.
  elseif(temps(ift,i,j,k).ge.th)then
    psiw=1.
  else
    psiw=(temps(ift,i,j,k)-tapsw)/(th-tapsw)
  endif
  psiwmaxold=psiwmax(ift,i,j,k)
  psiwmax(ift,i,j,k)=max(psiwmax(ift,i,j,k),psiw)
  if(psiwmax(ift,i,j,k).lt.1.and.rhoWater(ift,i,j,k).gt.0.)then
    fw=min(rhoWater(ift,i,j,k)*(psiwmax(ift,i,j,k)-psiwmaxold) &
      /(1.-psiwmaxold)/dtp,rhoWater(ift,i,j,k)/dtp)
    fw_sum=fw+fw_sum
  else
    fw=0.
  endif
  frhoWater=-fw*dtp

  ! Fuel Consumption
  slambdaof=rhoFuelold*xv(i,j,k,iO2)/(rhoFuelold/rnfuel+xv(i,j,k,iO2)/rno)**2
  ff=cf*rhoFuel(ift,i,j,k)*xv(i,j,k,iO2)*sigmac*psif*slambdaof &
    /(100.*sizescale(ift,i,j,k)**2.)
  ff_sum=ff_sum+ff
  frhoFuel=-ff*dtp*rnfuel

  ! Change in sies
  reactht=thetaSolid*hf*ff
  rmassloss=-tcrit*rnfuel*cpwood*ff
  waterevpht=-fw*hwevap
  rwaterlossht=-fw*cpwater*twvap
  frhosies=(convht(ift,i,j,k)+reactht+rmassloss+waterevpht+rwaterlossht)*dtp
  
  ! Update fuel physical terms
  rhos=rhoFuel(ift,i,j,k)+rhoWater(ift,i,j,k)
  rhoFuel(ift,i,j,k)=rhoFuel(ift,i,j,k)+frhoFuel
  if(idiffsies.eq.1) rhoFuel(ift,i,j,k)=rhoFuel(ift,i,j,k)/(1+dtp*siesdiff(ift,i,j,k)/hf) ! addition of siesdiff mass loss for energy balance
  rhoWater(ift,i,j,k)=max(0.0,rhoWater(ift,i,j,k)+frhoWater)
  sies(ift,i,j,k)=(sies(ift,i,j,k)*rhos+frhosies+fsiesrad(ift,i,j,k)*dtp)/(rhoFuel(ift,i,j,k)+rhoWater(ift,i,j,k))
  if(idiffsies.eq.1) sies(ift,i,j,k)=sies(ift,i,j,k)+dtp*siesdiff(ift,i,j,k)
  call updateSolidThermProps(rhoFuel(ift,i,j,k),rhoWater(ift,i,j,k),sies(ift,i,j,k),temps(ift,i,j,k))

  ! Variable for gas phase
  rneteng=rneteng+ff*hf*(1.-thetaSolid)
  if (iemissions.eq.1) call soot(ift,i,j,k,ff)
 
end subroutine fuel

!----------------------------------------------------------------
! computePsi(temps) computes the pdf with several methods
!----------------------------------------------------------------
real function computePsi(temps)
  use gridlist_variables, only : ifuel
  use fuel_variables, only : tcrit,tfstep
  Implicit None

  ! Local Variables
  real,intent(in) :: temps

  real :: rlcrit,rlramp,rljoint1,scalefactor,rpsi
  real :: rltop,rljoint2,rpsi2,rl
  real :: xerf,terf,er,sumf
  real :: tramp=500. ! the tail of the ramp straight line (if extrapolated to psi-0)
  real :: psijoint=0.30 ! height of the ramp joint
  real :: c1=0.5
  real :: c2=0.0079
  real :: c3=1.
  real :: p=0.47047
  real :: a1=0.3480242
  real :: a2=-0.0958798
  real :: a3=0.7478556

  ! Executable Code  
  if(ifuel.eq.1)then  !Use RRL's method of calculating computePsi!
    rlcrit=(tcrit-tfstep)
    rlramp=(tramp-tfstep)
    rljoint1=psijoint*2.*(rlcrit-rlramp)+rlramp
    scalefactor=sqrt(rlramp**2-(rljoint1-rlramp)**2)/psijoint
    rpsi=2.0/scalefactor*(rlcrit-rlramp)*rljoint1+scalefactor*psijoint
    rltop=2.*rlcrit
    rljoint2=rltop-rljoint1
    rpsi2=scalefactor-rpsi
       
    ! temps is used here
    rl=temps-tfstep
    if (rl.le.0.0 ) then
      computePsi=0.0
    elseif (rl.le.rljoint1) then
      computePsi=1./scalefactor*(rpsi-sqrt(rpsi**2-rl**2))
    elseif (rl.le.rljoint2) then
      computePsi=.5/(rlcrit-rlramp)*(rl-rlramp)
    elseif (rl.le.rltop) then
      computePsi=1./scalefactor*(rpsi2+sqrt(rpsi**2-(rl-rltop)**2))
    elseif (rl.gt.rltop) then
      computePsi=1.
    endif

  elseif(ifuel.eq.2)then  !Use Michael Clark's method of calculating computePsi!
    xerf = c2*(temps-tcrit)
    terf = 1./(1.+p*abs(xerf))
    sumf = a1*terf+a2*terf*terf+a3*terf*terf*terf
    er   = 1-sumf*exp(-xerf*xerf)
    computePsi = c1*(c3+er)
    if (xerf.lt.0.0) computePsi=-computePsi+c3
      
  elseif(ifuel.eq.3)then  !Use JAS's method of calculating computePsi!
    if(temps.le.tfstep)then 
      computePsi=0.0
    elseif(temps.le.(2*(tcrit-tfstep)+tfstep))then
      er = erf(c2*(temps-tcrit))                         
      computePsi = c1*(c3+er)
    else 
      computePsi = 1.0
    endif  
  endif    !end if(ifuel==1,2,or 3)        
        
end function computePsi
