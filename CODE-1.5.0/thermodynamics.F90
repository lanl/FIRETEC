!-----------------------------------------------------------------------
! updates gas thermal dynamic properties due to changes in gas 
! composition, such as the presence of water vapor
!-----------------------------------------------------------------------
subroutine updateGasThermProps(yEmit,cpEmit,mwEmit,nEmit,theta,rho, &
    cp_gas,cv_gas,mw_gas,pr,tempg)
  use gridlist_variables, only : prec
  use constants, only : Rgas,Pref
  Implicit None
  
  ! Local Variables
  real(prec),intent(inout) :: cp_gas,cv_gas,mw_gas,pr,tempg
  real(prec),intent(in) :: yEmit(nEmit),cpEmit(nEmit),mwEmit(nEmit)
  real(prec),intent(in) :: theta,rho
  integer,intent(in) :: nEmit
  
  real(prec) :: cp_N2=1040.  ! J/kg*K
  real(prec) :: MW_N2=28.014 ! kg/kmol
  real(prec) :: rg_over_prrcp_gas

  ! Executable Code
  cp_gas=sum(cpEmit*yEmit)+cp_N2*(1.-sum(yEmit))
  mw_gas=sum(mwEmit*yEmit)+MW_N2*(1.-sum(yEmit))
  cv_gas=cp_gas-Rgas/mw_gas
  rg_over_prrcp_gas=Rgas/mw_gas*Pref**(-Rgas/cp_gas/mw_gas)
  pr=(theta*rho*rg_over_prrcp_gas)**(cp_gas/cv_gas)
  tempg=theta*(pr/Pref)**(1.-cv_gas/cp_gas)

end subroutine updateGasThermProps

!-----------------------------------------------------------------------
! updates solid thermal dynamic properties due to changes in 
! solid composition, such as from burning
!-----------------------------------------------------------------------
subroutine updateSolidThermProps(rhoFuel,rhoWater,sies,temps)
  use gridlist_variables, only : cpwood,prec
  use thermo_variables, only : cpwater
  Implicit None
  
  ! Local Variables
  real(prec),intent(in) :: rhoFuel,rhoWater,sies
  real(prec),intent(inout) :: temps
  real(prec) :: cpsolid

  ! Executable Code
  cpsolid=(rhoFuel*cpwood+rhoWater*cpwater)/(rhoFuel+rhoWater)
  temps=sies/cpsolid

end subroutine updateSolidThermProps
