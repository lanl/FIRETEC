!-----------------------------------------------------------------------
! Reads basic user input for turbulent subgrid modeling in 
! HIGRAD/FIRETEC
!-----------------------------------------------------------------------
subroutine input_subGridTurb(iturb)
  use turb_gridlist_variables, only : isa,rturbprandtl
  Implicit None
  
  ! Local Variables
  integer,intent(in) :: iturb

  ! Executable Code
  if(iturb.eq.2)then ! Rod Linn subgrid Turbulence Model
    call QueryGridlist_integer('isa',isa,48)
    call QueryGridlist_real   ('rturbprandtl',rturbprandtl,48)
  endif

end subroutine input_subGridTurb
