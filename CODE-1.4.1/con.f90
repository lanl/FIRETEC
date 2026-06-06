!234567**********************************************
subroutine con
! this subroutine sets the values of the constants_old that
! are going to be used throughout HIGRAD and FIRETEC
  use gridlist_variables, only : iemissions,inonlocal,tambient, &
    Water2WoodRatio
  use thermo_variables, only : cpwater,cp_gas,cp_over_cv_gas, &
    rg_over_cp_gas,rg_over_prrcp_gas,specifichumidity
  use fuel_variables, only : tfstep,tcrit,tstep
  use constants, only : Rgas
  use gridsetup_old
  use constants_old
  use soot_constants_old
  use nonlocal
!  use turba, only :rturbprandtl
  Implicit None

  return
end
