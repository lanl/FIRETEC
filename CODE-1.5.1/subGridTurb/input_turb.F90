!-----------------------------------------------------------------------
! Reads basic user input for turbulent subgrid modeling in 
! HIGRAD/FIRETEC
!-----------------------------------------------------------------------
subroutine input_subGridTurb(iturb)
  use gridlist_variables, only : wp
  use turb_gridlist_variables, only : isa,rturbprandtl
  use QueryGridlist_mod, only : QueryGridlist_variable
  Implicit None
  
  ! Local Variables
  integer,intent(in) :: iturb

  ! Executable Code
  if(iturb.eq.2)then ! Rod Linn subgrid Turbulence Model
    call QueryGridlist_variable('isa',isa,48)
    call QueryGridlist_variable('rturbprandtl',rturbprandtl,48)
    if(rturbprandtl.lt.0.0_wp)then
      print*,'ERROR: rturbprandtl must be non-negative for scalar turbulence diffusion'
      STOP
    endif
    if(isa.lt.1.or.isa.gt.5)then
      print*,'Unsupported turbulence length-scale option isa=',isa
      print*,'Supported values are isa=1,2,3,4,5'
      STOP
    endif
  endif

end subroutine input_subGridTurb
