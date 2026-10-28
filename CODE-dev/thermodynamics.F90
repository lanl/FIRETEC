!-----------------------------------------------------------------------
! updates gas thermal dynamic properties due to changes in gas 
! composition, such as the presence of water vapor
!-----------------------------------------------------------------------
subroutine updateGasThermProps(yEmit,cpEmit,mwEmit,nEmit,theta,rho, &
    cp_gas,cv_gas,mw_gas,pr,tempg)
  use gridlist_variables, only : wp
  use constants, only : Rgas,Pref
  Implicit None
  
  ! Local Variables
  integer,intent(in) :: nEmit
  real(wp),intent(inout) :: cp_gas,cv_gas,mw_gas,pr,tempg
  real(wp),intent(in) :: yEmit(nEmit),cpEmit(nEmit),mwEmit(nEmit)
  real(wp),intent(in) :: theta,rho
  
  real(wp),parameter :: cp_N2=1040.0_wp  ! J/kg*K
  real(wp),parameter :: MW_N2=28.014_wp ! kg/kmol
  real(wp) :: rg_over_prrcp_gas

  ! Executable Code
  cp_gas=sum(cpEmit*yEmit)+cp_N2*(1._wp-sum(yEmit))
  mw_gas=1._wp/(sum(yEmit/mwEmit)+(1._wp-sum(yEmit))/MW_N2)
  cv_gas=cp_gas-Rgas/mw_gas
  rg_over_prrcp_gas=Rgas/mw_gas*Pref**(-Rgas/cp_gas/mw_gas)
  pr=(theta*rho*rg_over_prrcp_gas)**(cp_gas/cv_gas)
  tempg=theta*(pr/Pref)**(Rgas/mw_gas/cp_gas)

end subroutine updateGasThermProps

!-----------------------------------------------------------------------
! updates solid thermal dynamic properties due to changes in 
! solid composition, such as from burning
!-----------------------------------------------------------------------
subroutine updateSolidThermProps(rhoFuel,rhoWater,sies,temps)
  use gridlist_variables, only : wp
  use vegetation_gridlist_variables, only : cpwood
  use thermo_variables, only : cpwater
  use constants, only : tolerance
  Implicit None
  
  ! Local Variables
  real(wp),intent(in) :: rhoFuel,rhoWater,sies
  real(wp),intent(inout) :: temps
  real(wp) :: cpsolid

  ! Executable Code
  cpsolid=(rhoFuel*cpwood+rhoWater*cpwater)/(rhoFuel+rhoWater+tolerance)
  temps=sies/cpsolid

end subroutine updateSolidThermProps
