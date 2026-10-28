!-----------------------------------------------------------------------
! burn computes reaction terms of fuel on the small timestep
!-----------------------------------------------------------------------
subroutine burnFuel(force_xv,force_xvfuel,xv,xvfuel, &
    ilh,iuh,jlh,juh,ku,il,iu,jl,ju,lls,nvp,nvf,nfuel)
  use gridlist_variables, only : wp
  use vegetation_gridlist_variables, only : cpwood
  use xvall, only : ikb,iO2,irho,irhof,irhow,isies,itheta
  use thermo_variables, only : cp_gas,cv_gas,pr
  use fuel_variables, only : temps,rhoFuelInitial,rno,rnfuel, &
    sizeScale,hf,tcrit
  use linn_turb_variables, only : sc
  use constants, only : Pref
  Implicit None

  ! Local Variables  
  integer,intent(in) :: ilh,iuh,jlh,juh,ku,il,iu,jl,ju,lls,nvp,nvf,nfuel
  real(wp),intent(in) :: xv(ilh:iuh,jlh:juh,ku,nvp), &
    xvfuel(nfuel,il:iu,jl:ju,lls,nvf)
  real(wp),intent(inout) :: force_xv(ilh:iuh,jlh:juh,ku,nvp), &
    force_xvfuel(nfuel,il:iu,jl:ju,lls,nvf)

  integer :: i,j,k,ift
  real(wp) :: sigmac 
  real(wp) :: percHydroRemaining,hydroFactor
  real(wp) :: slambdaof
  real(wp) :: reactFuelGas,reactFuelSolid,fuelMassLoss
  real(wp) :: energyGas,energySolid
  real(wp) :: ff
  real(wp) :: rkctemp
  real(wp) :: psif,thetaSolid
  real(wp) :: cf 
  real(wp) :: hydroThresh=0.4_wp
  real(wp) :: thetag=0.75_wp ! fraction of energy that is deposited in the gas phase
                      ! for the computation of thetasolid
  real(wp) :: cfhydro=0.30_wp*3.0_wp ! Coef. of mass reaction rates.  (found through numerical experiments with
                         ! grass) cfhydro of .8 seems a little low with code of 5/24/01, flat 
                         ! ground, no wind, rhof=1 kg/m^3 and 5% moisture
  real(wp) :: cfchar=0.015_wp*6.0_wp ! Coef. of mass reaction rates.  (found through 

  ! Executable Code
  do k=1,lls
    do j=jl,ju
      do i=il,iu
        do ift=1,nfuel
          ! Turbulent mixing
          rkctemp=0.2_wp*xv(i,j,k,ikb)/xv(i,j,k,irho)
          ! the .2 factor is included as an extrapolation for turbulence at the c
          ! scale.  The .2 factor should be replaced by a function of sc and sb.
          sigmac=sc*0.5_wp*sqrt(rkctemp)

          ! Fuel Consumption
          call computePsif(psif,temps(ift,i,j,k))
          percHydroRemaining=max(0.0_wp,(xvfuel(ift,i,j,k,irhof) &
            -hydroThresh*rhoFuelinitial(ift,i,j,k)) &
            /(rhoFuelinitial(ift,i,j,k)*(1.0_wp-hydroThresh)))
          cf=cfhydro*percHydroRemaining+cfchar*(1-percHydroRemaining)
          hydroFactor=exp(-rno/rnfuel*psif* &
            xvfuel(ift,i,j,k,irhof)/xv(i,j,k,iO2))
          thetaSolid=percHydroRemaining*(1.0_wp-thetag)*hydroFactor &
            +(1.0_wp-percHydroRemaining)*thetag
          slambdaof=sum(xvfuel(:,i,j,k,irhof))*xv(i,j,k,iO2) &
            /(sum(xvfuel(:,i,j,k,irhof))/rnfuel+xv(i,j,k,iO2)/rno)**2
          ff=cf*xvfuel(ift,i,j,k,irhof)*xv(i,j,k,iO2)*sigmac*psif* &
            slambdaof/(100.0_wp*sizescale(ift,i,j,k)**2.0_wp)

          ! Energy Changes
          reactFuelGas=(1.0_wp-thetaSolid)*hf*ff
          reactFuelSolid=thetaSolid*hf*ff
          fuelMassLoss=tcrit*rnfuel*cpwood*ff

          energySolid=reactFuelSolid-fuelMassLoss
          energyGas=reactFuelGas+fuelMassLoss

          ! Forcings for fuel physical terms
          force_xvfuel(ift,i,j,k,irhof)=force_xvfuel(ift,i,j,k,irhof) &
            -ff*rnfuel
          force_xvfuel(ift,i,j,k,isies)=force_xvfuel(ift,i,j,k,isies)+ &
            energySolid/ &
            (xvfuel(ift,i,j,k,irhof)+xvfuel(ift,i,j,k,irhow))

          ! Forcings for gas physical terms
          force_xv(i,j,k,itheta)=force_xv(i,j,k,itheta)+ &
            energyGas/cp_gas(i,j,k)* &
            (Pref/pr(i,j,k))**(1._wp-cv_gas(i,j,k)/cp_gas(i,j,k))
          force_xv(i,j,k,irho)=force_xv(i,j,k,irho)+ff*rnfuel
        enddo
      enddo
    enddo
  enddo
 
end subroutine burnFuel

!-----------------------------------------------------------------------
! computePsi(temps) computes the pdf with several methods
!-----------------------------------------------------------------------
subroutine computePsif(psif,temps)
  use gridlist_variables, only : wp
  use vegetation_gridlist_variables, only : ifuel
  use fuel_variables, only : tcrit,tfstep
  Implicit None
  
  ! Local Variables
  real(wp),intent(inout) :: psif
  real(wp),intent(in) :: temps
  
  real(wp) :: rlcrit,rlramp,rljoint1,scalefactor,rpsi
  real(wp) :: rltop,rljoint2,rpsi2,rl
  real(wp) :: xerf,terf,er,sumf
  real(wp) :: tramp=500.0_wp ! the tail of the ramp straight line (if extrapolated to psi-0)
  real(wp) :: psijoint=0.30_wp ! height of the ramp joint
  real(wp) :: c1=0.5_wp,c2=0.0079_wp,c3=1.0_wp
  real(wp) :: p=0.47047_wp
  real(wp) :: a1=0.3480242_wp,a2=-0.0958798_wp,a3=0.7478556_wp
  
  ! Executable Code  
  if(ifuel.eq.1)then  !Use RRL's method of calculating psif!
    rlcrit=(tcrit-tfstep)
    rlramp=(tramp-tfstep)
    rljoint1=psijoint*2.0_wp*(rlcrit-rlramp)+rlramp
    scalefactor=sqrt(rlramp**2-(rljoint1-rlramp)**2)/psijoint
    rpsi=2.0_wp/scalefactor*(rlcrit-rlramp)*rljoint1+scalefactor*psijoint
    rltop=2.0_wp*rlcrit
    rljoint2=rltop-rljoint1
    rpsi2=scalefactor-rpsi

    ! temps is used here
    rl=temps-tfstep
    if (rl.le.0.0_wp ) then
      psif=0.0_wp
    elseif (rl.le.rljoint1) then
      psif=1.0_wp/scalefactor*(rpsi-sqrt(rpsi**2-rl**2))
    elseif (rl.le.rljoint2) then
      psif=0.5_wp/(rlcrit-rlramp)*(rl-rlramp)
    elseif (rl.le.rltop) then
      psif=1.0_wp/scalefactor*(rpsi2+sqrt(rpsi**2-(rl-rltop)**2))
    elseif (rl.gt.rltop) then
      psif=1.0_wp
    endif
  
  elseif(ifuel.eq.2)then  !Use Michael Clark's method of calculating psif!
    xerf = c2*(temps-tcrit)
    terf = 1.0_wp/(1.0_wp+p*abs(xerf))
    sumf = a1*terf+a2*terf*terf+a3*terf*terf*terf
    er   = 1-sumf*exp(-xerf*xerf)
    psif = c1*(c3+er)
    if (xerf.lt.0.0_wp) psif=-psif+c3

  elseif(ifuel.eq.3)then  !Use JAS's method of calculating psif!
    if(temps.le.tfstep)then 
      psif=0.0_wp
    elseif(temps.le.(2*(tcrit-tfstep)+tfstep))then
      er = erf(c2*(temps-tcrit))                         
      psif = c1*(c3+er)
    else 
      psif = 1.0_wp
    endif  
  endif    !end if(ifuel==1,2,or 3)        

end subroutine computePsif
