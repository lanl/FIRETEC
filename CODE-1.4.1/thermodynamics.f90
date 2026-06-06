!----------------------------------------------------------------
! updates gas thermal dynamic properties due to changes in gas 
! composition, such as the presence of water vapor
!----------------------------------------------------------------
subroutine updateGasThermProps(yO2,yH2O)
  use thermo_variables, only : cp_gas,cp_over_cv_gas, &
    rg_over_cp_gas,rg_over_prrcp_gas
  use constants, only : Rgas 
  Implicit None
  
  ! Local Variables
  real,intent(in) :: yO2,yH2O
  real :: cp_O2=918.   ! J/kg*K
  real :: cp_H2O=1996. ! J/kg*K
  real :: cp_N2=1040.  ! J/kg*K
  real :: MW_O2=31.998 ! kg/kmol
  real :: MW_N2=28.014 ! kg/kmol
  real :: MW_H2O=18.015! kg/kmol
  real :: MW_gas,cv_gas

  ! Executable Code
  cp_gas=cp_O2*yO2+cp_H2O*yH2O+cp_N2*(1.-yO2-yH2O)
  MW_gas=MW_O2*yO2+MW_H2O*yH2O+MW_N2*(1.-yO2-yH2O)
  cv_gas=cp_gas-Rgas/MW_gas
  cp_over_cv_gas=cp_gas/cv_gas
  rg_over_cp_gas=1.-1./cp_over_cv_gas
  rg_over_prrcp_gas=(cp_gas-cv_gas)/(1.e5**rg_over_cp_gas)

end subroutine updateGasThermProps

!----------------------------------------------------------------
! updates solid thermal dynamic properties due to changes in 
! solid composition, such as from burning
!----------------------------------------------------------------
subroutine updateSolidThermProps(rhoFuel,rhoWater,sies,temps)
  use gridlist_variables, only : cpwood
  use thermo_variables, only : cpwater
  Implicit None
  
  ! Local Variables
  real,intent(in) :: rhoFuel,rhoWater,sies
  real,intent(inout) :: temps
  real :: cpsolid

  ! Executable Code
  cpsolid=(rhoFuel*cpwood+rhoWater*cpwater)/(rhoFuel+rhoWater)
  temps=sies/cpsolid

end subroutine updateSolidThermProps
